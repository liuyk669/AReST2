"""AReST GSE182127 reproducibility workflow.

The public notebook calls this module to recompute AReST scores from the
released AnnData file, apply max-statistic permutation thresholding, extract a
DBSCAN region, and evaluate the result against the reference lesion label.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from sklearn.cluster import DBSCAN
from sklearn.decomposition import PCA
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.neighbors import NearestNeighbors


HEALTHY_SAMPLES = (
    "GSM5519054_control",
    "GSM5519060_sham",
    "GSM5519055_hemeHpx_1000",
)
ABNORMAL_SAMPLES = (
    "GSM5519056_heme_0030",
    "GSM5519057_heme_0125",
    "GSM5519058_heme_0500",
    "GSM5519059_heme_1000",
)
TARGET_SAMPLE = "GSM5519059_heme_1000"


@dataclass(frozen=True)
class MethodConfig:
    n_top_hvg: int = 2000
    n_pc: int = 30
    k_expr: int = 80
    sigma_quantile: float = 0.60
    spatial_k: int = 6
    spatial_iters: int = 1
    n_perm: int = 1000
    alpha: float = 0.05
    seed: int = 0


@dataclass(frozen=True)
class RegionConfig:
    min_samples: int = 7
    min_cluster_size: int = 20


def _to_dense_nonnegative(matrix: Any) -> np.ndarray:
    if sparse.issparse(matrix):
        matrix = matrix.toarray()
    arr = np.asarray(matrix, dtype=float)
    arr[~np.isfinite(arr)] = 0.0
    arr[arr < 0.0] = 0.0
    return arr


def _normalize_log1p(adata: ad.AnnData) -> ad.AnnData:
    out = adata.copy()
    sc.pp.normalize_total(out, target_sum=1e4)
    sc.pp.log1p(out)
    return out


def _prepare_pca(adata: ad.AnnData, config: MethodConfig) -> np.ndarray:
    raw = adata.copy()
    top_n = min(int(config.n_top_hvg), raw.n_vars)
    sc.pp.highly_variable_genes(raw, flavor="seurat_v3", n_top_genes=top_n)
    if int(raw.var["highly_variable"].sum()) < 50:
        raise ValueError("Too few HVGs to build PCA.")
    hvg_mask = raw.var["highly_variable"].to_numpy()

    ad_proc = _normalize_log1p(adata)
    matrix = ad_proc[:, hvg_mask].X
    if sparse.issparse(matrix):
        matrix = matrix.toarray()
    matrix = np.asarray(matrix, dtype=float)
    matrix[~np.isfinite(matrix)] = np.nan
    if np.isnan(matrix).any():
        col_means = np.nanmean(matrix, axis=0)
        col_means = np.where(np.isfinite(col_means), col_means, 0.0)
        ii, jj = np.where(np.isnan(matrix))
        matrix[ii, jj] = col_means[jj]
    n_comp = int(min(config.n_pc, matrix.shape[1], matrix.shape[0] - 1))
    return PCA(n_components=n_comp, random_state=0).fit_transform(matrix)


def _build_expr_knn(x_pca: np.ndarray, k_expr: int) -> np.ndarray:
    n = x_pca.shape[0]
    k = int(min(k_expr, n - 1))
    idx = NearestNeighbors(n_neighbors=k + 1, algorithm="auto").fit(x_pca).kneighbors(return_distance=False)
    out = np.zeros((n, k), dtype=int)
    for i, row in enumerate(idx):
        row = row[row != i]
        if row.size < k:
            row = np.pad(row, (0, k - row.size), mode="edge")
        out[i] = row[:k]
    return out


def _score_from_local_target_fraction(p_hat: np.ndarray, rho: float) -> np.ndarray:
    rho = float(np.clip(rho, 1e-9, 1.0 - 1e-9))
    f1 = p_hat / rho
    f0 = (1.0 - p_hat) / (1.0 - rho)
    return np.clip((f1 - f0) / (f1 + f0), -1.0, 1.0)


def _moran_i(values: np.ndarray, coords: np.ndarray, k: int = 6) -> float:
    values = np.asarray(values, dtype=float).reshape(-1)
    coords = np.asarray(coords, dtype=float)
    ok = np.isfinite(values) & np.all(np.isfinite(coords), axis=1)
    values = values[ok]
    coords = coords[ok]
    n = values.size
    if n <= max(k + 1, 3):
        return float("nan")
    centered = values - float(np.mean(values))
    denom = float(np.sum(centered**2))
    if denom <= 1e-12:
        return float("nan")
    idx = NearestNeighbors(n_neighbors=min(k + 1, n), algorithm="kd_tree").fit(coords).kneighbors(return_distance=False)
    w = np.zeros((n, n), dtype=np.uint8)
    for i in range(n):
        w[i, idx[i, 1:]] = 1
    w = np.maximum(w, w.T)
    sum_w = float(w.sum())
    if sum_w <= 0:
        return float("nan")
    return float((n / sum_w) * (np.sum(w * np.outer(centered, centered)) / denom))


def _spatial_average(values: np.ndarray, coords: np.ndarray, config: MethodConfig) -> np.ndarray:
    # This is part of the region-calling statistic used in the paper example.
    coords = np.asarray(coords, dtype=float)
    out = np.asarray(values, dtype=float).copy()
    k_use = int(min(config.spatial_k + 1, coords.shape[0]))
    _, idx = NearestNeighbors(n_neighbors=k_use, algorithm="kd_tree").fit(coords).kneighbors(coords)
    for _ in range(max(1, int(config.spatial_iters))):
        out = np.mean(out[idx], axis=1)
    return out


def load_gse182127_h5ad(h5ad_path: str | Path) -> ad.AnnData:
    keep = set(HEALTHY_SAMPLES + ABNORMAL_SAMPLES)
    adata = ad.read_h5ad(h5ad_path)
    adata.obs_names_make_unique()
    adata.obs["sample"] = adata.obs["sample"].astype(str)
    adata = adata[adata.obs["sample"].isin(keep).to_numpy()].copy()
    if "spatial_aligned" not in adata.obsm:
        raise ValueError("Expected coordinates in adata.obsm['spatial_aligned'].")
    return adata


def compute_arest_scores(adata: ad.AnnData, config: MethodConfig) -> dict[str, Any]:
    coords = np.asarray(adata.obsm["spatial_aligned"], dtype=float)
    sample_vec = adata.obs["sample"].astype(str).to_numpy()
    is_target_condition = np.isin(sample_vec, ABNORMAL_SAMPLES)
    y01 = is_target_condition.astype(int)
    rho = float(np.mean(y01))

    x_pca = _prepare_pca(adata, config)
    knn_idx = _build_expr_knn(x_pca, config.k_expr)
    nb_coords = coords[knn_idx]
    d2 = np.sum((nb_coords - coords[:, None, :]) ** 2, axis=2).astype(float)
    d = np.sqrt(np.maximum(d2[np.isfinite(d2)], 0.0))
    sigma = float(np.quantile(d, config.sigma_quantile))
    if not np.isfinite(sigma) or sigma <= 0:
        sigma = 1.0

    weights = np.exp(-d2 / (2.0 * sigma**2 + 1e-12))
    weight_sum = np.sum(weights, axis=1) + 1e-12
    p_hat = np.sum(weights * y01[knn_idx].astype(float), axis=1) / weight_sum
    score_raw = _score_from_local_target_fraction(p_hat, rho)

    tables = []
    score_spatial = np.zeros_like(score_raw)
    for sample in HEALTHY_SAMPLES + ABNORMAL_SAMPLES:
        mask = sample_vec == sample
        sample_score = _spatial_average(score_raw[mask], coords[mask], config)
        score_spatial[mask] = sample_score
        tables.append(
            pd.DataFrame(
                {
                    "obs_name": adata.obs_names[mask].astype(str),
                    "sample": sample,
                    "group": "target" if sample in ABNORMAL_SAMPLES else "reference",
                    "x": coords[mask, 0],
                    "y": coords[mask, 1],
                    "score_raw": score_raw[mask],
                    "score": sample_score,
                }
            )
        )

    return {
        "coords": coords,
        "sample_vec": sample_vec,
        "y01": y01,
        "knn_idx": knn_idx,
        "weights": weights,
        "weight_sum": weight_sum,
        "rho": rho,
        "score_raw": score_raw,
        "score": score_spatial,
        "score_table": pd.concat(tables, ignore_index=True),
        "sigma": sigma,
    }


def permutation_max_threshold(score_cache: dict[str, Any], config: MethodConfig, target_sample: str = TARGET_SAMPLE) -> dict[str, Any]:
    sample_vec = np.asarray(score_cache["sample_vec"], dtype=str)
    target_mask = sample_vec == target_sample
    target_coords = np.asarray(score_cache["coords"], dtype=float)[target_mask]
    target_raw = np.asarray(score_cache["score_raw"], dtype=float)[target_mask]
    observed = _spatial_average(target_raw, target_coords, config)

    rng = np.random.default_rng(int(config.seed))
    y01 = np.asarray(score_cache["y01"], dtype=int)
    knn_idx = np.asarray(score_cache["knn_idx"], dtype=int)
    weights = np.asarray(score_cache["weights"], dtype=float)
    weight_sum = np.asarray(score_cache["weight_sum"], dtype=float)
    max_stats = np.zeros(int(config.n_perm), dtype=float)
    for b in range(int(config.n_perm)):
        y_perm = rng.permutation(y01)
        p_hat = np.sum(weights * y_perm[knn_idx].astype(float), axis=1) / weight_sum
        raw_perm = _score_from_local_target_fraction(p_hat, float(score_cache["rho"]))
        max_stats[b] = float(np.nanmax(_spatial_average(raw_perm[target_mask], target_coords, config)))
    tau = float(np.quantile(max_stats[np.isfinite(max_stats)], 1.0 - float(config.alpha)))
    p_adj = np.mean(max_stats[:, None] >= observed[None, :], axis=0)
    return {
        "observed": observed,
        "max_stats": max_stats,
        "tau": tau,
        "p_adj": p_adj,
        "sig_mask": observed > tau,
    }


def dbscan_region(coords: np.ndarray, sig_mask: np.ndarray, config: RegionConfig) -> dict[str, Any]:
    coords = np.asarray(coords, dtype=float)
    sig_mask = np.asarray(sig_mask, dtype=bool)
    sig_idx = np.where(sig_mask)[0]
    sig_coords = coords[sig_idx]
    if sig_coords.shape[0] < 7:
        raise ValueError("Need at least 7 significant spots to compute the sixth-neighbor radius.")
    dists, _ = NearestNeighbors(n_neighbors=7, algorithm="kd_tree").fit(sig_coords).kneighbors(sig_coords)
    eps = float(np.median(dists[:, 6]))
    labels_sig = DBSCAN(eps=eps, min_samples=int(config.min_samples)).fit_predict(sig_coords)
    unique, counts = np.unique(labels_sig[labels_sig >= 0], return_counts=True)
    size_map = {int(label): int(size) for label, size in zip(unique, counts)}
    kept = [label for label, size in size_map.items() if size >= int(config.min_cluster_size)]
    kept = sorted(kept, key=lambda label: (-size_map[label], label))
    labels = np.full(coords.shape[0], -1, dtype=int)
    region = np.zeros(coords.shape[0], dtype=bool)
    for new_label, old_label in enumerate(kept):
        global_idx = sig_idx[labels_sig == old_label]
        labels[global_idx] = new_label
        region[global_idx] = True
    return {
        "region_mask": region,
        "labels": labels,
        "eps": eps,
        "n_sig": int(sig_mask.sum()),
        "n_clusters_total": int(len(unique)),
        "n_clusters_retained": int(len(kept)),
        "retained_cluster_sizes": [size_map[label] for label in kept],
        "n_noise_sig": int(np.sum(labels_sig < 0)),
    }


def evaluate_predictions(score: np.ndarray, region: np.ndarray, truth: np.ndarray, coords: np.ndarray) -> pd.DataFrame:
    score = np.asarray(score, dtype=float)
    region = np.asarray(region, dtype=bool)
    truth = np.asarray(truth, dtype=bool)
    tp = int(np.sum(region & truth))
    fp = int(np.sum(region & ~truth))
    fn = int(np.sum(~region & truth))
    precision = tp / (tp + fp) if tp + fp else 0.0
    recall = tp / (tp + fn) if tp + fn else 0.0
    jaccard = tp / (tp + fp + fn) if tp + fp + fn else 0.0
    dice = 2 * tp / (2 * tp + fp + fn) if 2 * tp + fp + fn else 0.0
    return pd.DataFrame(
        [
            {"metric": "AUROC", "value": float(roc_auc_score(truth.astype(int), score))},
            {"metric": "AUPRC", "value": float(average_precision_score(truth.astype(int), score))},
            {"metric": "Moran_I", "value": float(_moran_i(score, coords, k=6))},
            {"metric": "Region_precision", "value": float(precision)},
            {"metric": "Region_recall", "value": float(recall)},
            {"metric": "Region_Jaccard", "value": float(jaccard)},
            {"metric": "Region_Dice", "value": float(dice)},
            {"metric": "Region_size", "value": int(region.sum())},
            {"metric": "Reference_region_size", "value": int(truth.sum())},
            {"metric": "Region_overlap", "value": int(tp)},
        ]
    )


def _plot_score(df: pd.DataFrame, out: Path) -> None:
    fig = plt.figure(figsize=(5.0, 5.2), dpi=220)
    ax = fig.add_axes([0.08, 0.18, 0.84, 0.76])
    cb_ax = fig.add_axes([0.08, 0.07, 0.84, 0.035])
    vmin, vmax = np.nanpercentile(df["score"], [5, 95])
    span = max(float(vmax - vmin), 1e-9)
    vmin = float(vmin - 0.08 * span)
    vmax = float(vmax + 0.08 * span)
    sca = ax.scatter(df["x"], df["y"], c=df["score"], cmap="coolwarm", vmin=vmin, vmax=vmax, s=12, linewidths=0)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.invert_yaxis()
    for spine in ax.spines.values():
        spine.set_linewidth(1.0)
        spine.set_color("black")
    cb = fig.colorbar(sca, cax=cb_ax, orientation="horizontal")
    cb.set_ticks(np.linspace(vmin, vmax, 5))
    cb.ax.tick_params(labelsize=8, length=3, pad=2)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out.with_suffix(".png"), dpi=300)
    fig.savefig(out.with_suffix(".pdf"))
    plt.close(fig)


def _plot_region(df: pd.DataFrame, mask_col: str, out: Path) -> None:
    palette = {
        "cortex": "#e9d8a6",
        "striatum": "#b8d4e8",
        "diencephalon": "#c7e3bf",
        "globus_pallidus": "#dec9e9",
        "corpus_callosum": "#f1c0a8",
        "plexus": "#d9d9d9",
        "NA": "#eceff2",
    }
    fig, ax = plt.subplots(figsize=(5.0, 4.8), dpi=220)
    labels = df["anatomy_paper6"].fillna("NA").astype(str)
    for label in sorted(labels.unique()):
        sub = df[labels == label]
        ax.scatter(sub["x"], sub["y"], s=10, c=palette.get(label, "#eceff2"), alpha=0.65, linewidths=0, label=label)
    mask = df[mask_col].astype(bool).to_numpy()
    ax.scatter(df.loc[mask, "x"], df.loc[mask, "y"], s=14, c="#ff7a1a", alpha=0.97, linewidths=0)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.invert_yaxis()
    for spine in ax.spines.values():
        spine.set_linewidth(1.0)
        spine.set_color("black")
    ax.legend(frameon=False, fontsize=6, loc="lower left", bbox_to_anchor=(0.0, -0.02), ncol=2, handletextpad=0.2)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def _paper6_anatomy(anatomy: pd.Series) -> pd.Series:
    mapping = {
        "cortex": "cortex",
        "caudate_putamen": "striatum",
        "globus_pallidus": "globus_pallidus",
        "thalamus": "diencephalon",
        "hypothalamus": "diencephalon",
        "corpus_callosum": "corpus_callosum",
        "plexus": "plexus",
    }
    return anatomy.astype("object").map(mapping).fillna("NA")


def _prepare_author_annotation(annotation_csv: str | Path) -> pd.DataFrame:
    ann = pd.read_csv(annotation_csv)
    if "barcode" not in ann.columns or "cluster" not in ann.columns:
        raise ValueError("Annotation CSV must include at least 'barcode' and 'cluster' columns.")
    out = ann.copy()
    out["barcode"] = out["barcode"].astype(str)
    out["author_cluster"] = out["cluster"].astype(int)
    out["author10_12"] = out["author_cluster"].isin([10, 12])
    if "anatomy" not in out.columns:
        out["anatomy"] = "NA"
    out["anatomy_paper6"] = _paper6_anatomy(out["anatomy"])
    return out[["barcode", "author_cluster", "author10_12", "anatomy", "anatomy_paper6"]]


def run_gse182127_reproduction(
    h5ad_path: str | Path,
    annotation_csv: str | Path,
    output_dir: str | Path,
    method_config: MethodConfig | None = None,
    region_config: RegionConfig | None = None,
) -> dict[str, Any]:
    method_config = method_config or MethodConfig()
    region_config = region_config or RegionConfig()
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "figures").mkdir(exist_ok=True)
    (output_dir / "tables").mkdir(exist_ok=True)

    adata = load_gse182127_h5ad(h5ad_path)
    score_cache = compute_arest_scores(adata, method_config)
    score_cache["score_table"].to_csv(output_dir / "tables" / "all_sample_scores.csv", index=False)

    target_scores = score_cache["score_table"].loc[score_cache["score_table"]["sample"] == TARGET_SAMPLE].reset_index(drop=True)
    perm = permutation_max_threshold(score_cache, method_config, target_sample=TARGET_SAMPLE)
    target_scores["p_adj_maxT"] = perm["p_adj"]
    target_scores["sig_positive"] = perm["sig_mask"].astype(int)

    ann = _prepare_author_annotation(annotation_csv)
    target_scores["barcode"] = target_scores["obs_name"].astype(str).str.split("__", n=1).str[0]
    target = target_scores.merge(ann, on="barcode", how="left")
    if target["author10_12"].isna().any():
        raise ValueError("Annotation CSV did not cover all target spots.")
    target["author10_12"] = target["author10_12"].astype(bool)

    # Region calling uses the saved target table, matching the released analysis
    # trace and making the notebook outputs byte-for-byte reproducible.
    region_input = output_dir / "tables" / "target_scores_for_region_calling.csv"
    target.to_csv(region_input, index=False)
    target = pd.read_csv(region_input)
    target["author10_12"] = target["author10_12"].astype(bool)

    region = dbscan_region(target[["x", "y"]].to_numpy(float), target["sig_positive"].astype(bool).to_numpy(), region_config)
    target["arest_region"] = region["region_mask"].astype(int)
    target["dbscan_label"] = region["labels"]
    target.to_csv(output_dir / "tables" / "target_scores_region.csv", index=False)

    metrics = evaluate_predictions(
        target["score"].to_numpy(float),
        target["arest_region"].astype(bool).to_numpy(),
        target["author10_12"].astype(bool).to_numpy(),
        target[["x", "y"]].to_numpy(float),
    )
    metrics.to_csv(output_dir / "tables" / "metrics.csv", index=False)

    _plot_score(target, output_dir / "figures" / "arest_score_map")
    _plot_region(target, "author10_12", output_dir / "figures" / "reference_region_map")
    _plot_region(target, "arest_region", output_dir / "figures" / "arest_region_map")

    summary = {
        "method_config": asdict(method_config),
        "region_config": asdict(region_config),
        "n_total_spots": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "target_sample": TARGET_SAMPLE,
        "sigma": float(score_cache["sigma"]),
        "tau_high": float(perm["tau"]),
        "n_significant": int(perm["sig_mask"].sum()),
        **{row["metric"]: row["value"] for _, row in metrics.iterrows()},
        **region,
    }
    # Convert arrays before writing JSON.
    summary.pop("region_mask", None)
    summary.pop("labels", None)
    pd.Series(summary).to_json(output_dir / "summary.json", indent=2)
    return {"target": target, "metrics": metrics, "summary": summary, "output_dir": output_dir}
