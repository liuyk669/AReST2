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
    n_perm: int = 1000
    alpha: float = 0.05
    seed: int = 0


@dataclass(frozen=True)
class RegionConfig:
    min_samples: int = 7
    min_cluster_size: int = 20


def _prepare_pca(adata: ad.AnnData, config: MethodConfig) -> np.ndarray:
    raw = adata.copy()
    top_n = min(int(config.n_top_hvg), raw.n_vars)
    sc.pp.highly_variable_genes(raw, flavor="seurat_v3", n_top_genes=top_n)
    hvg_mask = raw.var["highly_variable"].to_numpy()

    ad_proc = adata.copy()
    sc.pp.normalize_total(ad_proc, target_sum=1e4)
    sc.pp.log1p(ad_proc)
    matrix = ad_proc[:, hvg_mask].X
    if sparse.issparse(matrix):
        matrix = matrix.toarray()
    matrix = np.asarray(matrix, dtype=float)
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
    n = values.size
    centered = values - float(np.mean(values))
    denom = float(np.sum(centered**2))
    idx = NearestNeighbors(n_neighbors=k + 1, algorithm="kd_tree").fit(coords).kneighbors(return_distance=False)
    w = np.zeros((n, n), dtype=np.uint8)
    for i in range(n):
        w[i, idx[i, 1:]] = 1
    sum_w = float(w.sum())
    return float((n / sum_w) * (np.sum(w * np.outer(centered, centered)) / denom))


def load_gse182127_h5ad(h5ad_path: str | Path) -> ad.AnnData:
    keep = set(HEALTHY_SAMPLES + ABNORMAL_SAMPLES)
    adata = ad.read_h5ad(h5ad_path)
    adata.obs_names_make_unique()
    adata.obs["sample"] = adata.obs["sample"].astype(str)
    adata = adata[adata.obs["sample"].isin(keep).to_numpy()].copy()
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
    d = np.sqrt(d2.ravel())
    sigma = float(np.quantile(d, config.sigma_quantile))

    weights = np.exp(-d2 / (2.0 * sigma**2 + 1e-12))
    weight_sum = np.sum(weights, axis=1) + 1e-12
    p_hat = np.sum(weights * y01[knn_idx].astype(float), axis=1) / weight_sum
    score_raw = _score_from_local_target_fraction(p_hat, rho)

    tables = []
    score_spatial = np.zeros_like(score_raw)
    for sample in HEALTHY_SAMPLES + ABNORMAL_SAMPLES:
        mask = sample_vec == sample
        sample_coords = coords[mask]
        _, smooth_idx = NearestNeighbors(n_neighbors=7, algorithm="kd_tree").fit(sample_coords).kneighbors(sample_coords)
        sample_score = np.mean(score_raw[mask][smooth_idx], axis=1)
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
    _, smooth_idx = NearestNeighbors(n_neighbors=7, algorithm="kd_tree").fit(target_coords).kneighbors(target_coords)
    observed = np.mean(target_raw[smooth_idx], axis=1)

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
        max_stats[b] = float(np.max(np.mean(raw_perm[target_mask][smooth_idx], axis=1)))
    tau = float(np.quantile(max_stats, 1.0 - float(config.alpha)))
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
    precision = tp / (tp + fp)
    recall = tp / (tp + fn)
    jaccard = tp / (tp + fp + fn)
    dice = 2 * tp / (2 * tp + fp + fn)
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
    fig = plt.figure(figsize=(4.2, 4.45), dpi=220)
    ax = fig.add_axes([0.08, 0.18, 0.84, 0.76])
    cb_ax = fig.add_axes([0.08, 0.075, 0.84, 0.035])
    vmin, vmax = -0.43, 1.00
    sca = ax.scatter(df["x"], df["y"], c=df["score"], cmap="coolwarm", vmin=vmin, vmax=vmax, s=10, linewidths=0)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.invert_yaxis()
    for spine in ax.spines.values():
        spine.set_linewidth(1.0)
        spine.set_color("black")
    cb = fig.colorbar(sca, cax=cb_ax, orientation="horizontal")
    cb.set_ticks([-0.43, -0.07, 0.28, 0.64, 1.00])
    cb.set_ticklabels(["-0.43", "-0.07", "0.28", "0.64", "1.00"])
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
    fig = plt.figure(figsize=(4.2, 4.55), dpi=220)
    ax = fig.add_axes([0.08, 0.25, 0.84, 0.70])
    labels = df["anatomy_paper6"].fillna("NA").astype(str)
    for label in sorted(labels.unique()):
        sub = df[labels == label]
        ax.scatter(sub["x"], sub["y"], s=8, c=palette.get(label, "#eceff2"), alpha=0.65, linewidths=0, label=label)
    mask = df[mask_col].astype(bool).to_numpy()
    ax.scatter(df.loc[mask, "x"], df.loc[mask, "y"], s=11, c="#ff7a1a", alpha=0.97, linewidths=0)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.invert_yaxis()
    for spine in ax.spines.values():
        spine.set_linewidth(1.0)
        spine.set_color("black")
    ax.legend(
        frameon=False,
        fontsize=6,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.06),
        ncol=3,
        handletextpad=0.2,
        columnspacing=0.8,
    )
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out.with_suffix(".png"), dpi=300, bbox_inches="tight", pad_inches=0.03)
    fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight", pad_inches=0.03)
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
    out = ann.copy()
    out["barcode"] = out["barcode"].astype(str)
    out["author_cluster"] = out["cluster"].astype(int)
    out["author10_12"] = out["author_cluster"].isin([10, 12])
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
    target["author10_12"] = target["author10_12"].astype(bool)

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
    summary.pop("region_mask", None)
    summary.pop("labels", None)
    pd.Series(summary).to_json(output_dir / "summary.json", indent=2)
    return {"target": target, "metrics": metrics, "summary": summary, "output_dir": output_dir}
