from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
import seaborn as sns

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from compute_pairwise_donor_overlap import (
    build_variant_index,
    contains_with_hamming_leq_one,
    load_significant_donor_clonotypes,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Plot pairwise cluster-dispersion summaries and a detailed donor-pair "
            "view showing how matched clonotypes are distributed across clusters."
        )
    )
    parser.add_argument(
        "--pairwise",
        type=Path,
        default=Path("results") / "yfv_redcea" / "pairwise_overlap" / "yfv_significant_overlap_pairwise.tsv",
        help="Pairwise overlap TSV from compute_pairwise_donor_overlap.py",
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=Path("results") / "yfv_redcea" / "raw_dbscan",
        help="Directory with YFV summary and cluster-members tables.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results") / "yfv_redcea" / "pairwise_overlap",
        help="Directory for output figures.",
    )
    parser.add_argument("--donor-x", default="Q1", help="Source donor for detailed pair plot.")
    parser.add_argument("--donor-y", default="Q2", help="Target donor for detailed pair plot.")
    parser.add_argument("--prefix", default="yfv_", help="Filename prefix filter.")
    parser.add_argument("--method", default="dbscan", help="Method filter inferred from filenames.")
    parser.add_argument("--fdr-threshold", type=float, default=0.05)
    parser.add_argument("--fold-threshold", type=float, default=0.0)
    return parser.parse_args()


def plot_heatmap(matrix: pd.DataFrame, path: Path, title: str, cmap: str, center: float | None = None) -> None:
    fig, ax = plt.subplots(figsize=(8.0, 6.2))
    vmax = float(np.nanmax(matrix.to_numpy()))
    vmin = float(np.nanmin(matrix.to_numpy()))
    sns.heatmap(
        matrix,
        cmap=cmap,
        annot=True,
        fmt=".2f",
        linewidths=0.5,
        linecolor="white",
        cbar_kws={"label": title},
        vmin=vmin,
        vmax=vmax,
        center=center,
        ax=ax,
    )
    ax.set_title(title)
    ax.set_xlabel("donor Y")
    ax.set_ylabel("donor X")
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(path.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)


def build_source_cluster_counts(
    donor_to_clusters: dict[str, dict[str, set[str]]],
    donor_x: str,
    donor_y: str,
) -> pd.DataFrame:
    y_index = build_variant_index(
        {
            clonotype
            for cluster_clonotypes in donor_to_clusters[donor_y].values()
            for clonotype in cluster_clonotypes
        }
    )
    rows: list[dict[str, object]] = []
    for cluster_id, clonotypes in donor_to_clusters[donor_x].items():
        matched = sum(1 for clonotype in clonotypes if contains_with_hamming_leq_one(clonotype, y_index))
        if matched == 0:
            continue
        rows.append(
            {
                "cluster_id": cluster_id,
                "cluster_size_x": len(clonotypes),
                "matched_clonotypes_to_y": matched,
                "matched_fraction_in_cluster_x": matched / len(clonotypes),
            }
        )
    return pd.DataFrame(rows).sort_values(
        ["matched_clonotypes_to_y", "matched_fraction_in_cluster_x", "cluster_size_x"],
        ascending=[False, False, False],
    )


def build_target_cluster_counts(
    donor_to_clusters: dict[str, dict[str, set[str]]],
    donor_x: str,
    donor_y: str,
) -> pd.DataFrame:
    x_index = build_variant_index(
        {
            clonotype
            for cluster_clonotypes in donor_to_clusters[donor_x].values()
            for clonotype in cluster_clonotypes
        }
    )
    rows: list[dict[str, object]] = []
    for cluster_id, clonotypes in donor_to_clusters[donor_y].items():
        matched = sum(1 for clonotype in clonotypes if contains_with_hamming_leq_one(clonotype, x_index))
        if matched == 0:
            continue
        rows.append(
            {
                "cluster_id": cluster_id,
                "cluster_size_y": len(clonotypes),
                "matched_clonotypes_from_x": matched,
                "matched_fraction_in_cluster_y": matched / len(clonotypes),
            }
        )
    return pd.DataFrame(rows).sort_values(
        ["matched_clonotypes_from_x", "matched_fraction_in_cluster_y", "cluster_size_y"],
        ascending=[False, False, False],
    )


def plot_pair_detail(
    source_df: pd.DataFrame,
    target_df: pd.DataFrame,
    donor_x: str,
    donor_y: str,
    path: Path,
) -> None:
    top_n = 30
    fig, axes = plt.subplots(1, 2, figsize=(15.5, 6.2))

    src = source_df.head(top_n).copy().reset_index(drop=True)
    src["rank"] = np.arange(1, len(src) + 1)
    sns.barplot(data=src, x="rank", y="matched_clonotypes_to_y", color="#2a6f97", ax=axes[0])
    axes[0].set_title(f"{donor_x} source clusters hitting {donor_y}")
    axes[0].set_xlabel("ranked source cluster")
    axes[0].set_ylabel("matched clonotypes from source cluster")
    axes[0].text(
        0.99,
        0.98,
        f"hit clusters: {len(source_df)}\nmedian matched/cluster: {source_df['matched_clonotypes_to_y'].median():.1f}",
        transform=axes[0].transAxes,
        ha="right",
        va="top",
        bbox={"boxstyle": "round,pad=0.3", "facecolor": "white", "edgecolor": "#999999"},
    )

    tgt = target_df.head(top_n).copy().reset_index(drop=True)
    tgt["rank"] = np.arange(1, len(tgt) + 1)
    sns.barplot(data=tgt, x="rank", y="matched_clonotypes_from_x", color="#c85c39", ax=axes[1])
    axes[1].set_title(f"{donor_y} target clusters receiving {donor_x}")
    axes[1].set_xlabel("ranked target cluster")
    axes[1].set_ylabel("matched clonotypes into target cluster")
    axes[1].text(
        0.99,
        0.98,
        f"hit clusters: {len(target_df)}\nmedian matched/cluster: {target_df['matched_clonotypes_from_x'].median():.1f}",
        transform=axes[1].transAxes,
        ha="right",
        va="top",
        bbox={"boxstyle": "round,pad=0.3", "facecolor": "white", "edgecolor": "#999999"},
    )

    fig.suptitle(
        f"Cluster dispersion for {donor_x} -> {donor_y}: few matched clonotypes can still spread across many clusters",
        y=1.02,
    )
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(path.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    if hasattr(sns, "set_theme"):
        sns.set_theme(style="white", font_scale=1.0)
    else:
        sns.set(style="white", font_scale=1.0)

    pairwise = pd.read_csv(args.pairwise, sep="\t")
    pairwise["cluster_coverage_over_clonotype_coverage"] = (
        pairwise["share_clusters_x_in_y"] / pairwise["share_x_in_y"]
    )
    pairwise["matched_clusters_per_100_clonotypes"] = (
        100.0 * pairwise["n_matched_clusters_x_in_y"] / pairwise["n_matched_x_in_y"]
    )

    ratio_matrix = pairwise.pivot(
        index="donor_x", columns="donor_y", values="cluster_coverage_over_clonotype_coverage"
    ).sort_index(axis=0).sort_index(axis=1)
    spread_matrix = pairwise.pivot(
        index="donor_x", columns="donor_y", values="matched_clusters_per_100_clonotypes"
    ).sort_index(axis=0).sort_index(axis=1)

    plot_heatmap(
        ratio_matrix,
        args.output_dir / "yfv_cluster_vs_clonotype_coverage_ratio",
        "cluster coverage / clonotype coverage",
        cmap="viridis",
        center=1.0,
    )
    plot_heatmap(
        spread_matrix,
        args.output_dir / "yfv_matched_clusters_per_100_clonotypes",
        "matched source clusters per 100 matched clonotypes",
        cmap="mako",
    )

    donor_to_clonotypes, donor_to_clusters, _ = load_significant_donor_clonotypes(
        args.input_dir,
        prefix=args.prefix,
        method_filter=args.method,
        donors_filter={args.donor_x, args.donor_y},
        fdr_threshold=args.fdr_threshold,
        fold_threshold=args.fold_threshold,
    )
    _ = donor_to_clonotypes
    source_df = build_source_cluster_counts(donor_to_clusters, args.donor_x, args.donor_y)
    target_df = build_target_cluster_counts(donor_to_clusters, args.donor_x, args.donor_y)
    plot_pair_detail(
        source_df,
        target_df,
        args.donor_x,
        args.donor_y,
        args.output_dir / f"yfv_{args.donor_x}_to_{args.donor_y}_cluster_dispersion",
    )

    print("Saved cluster-coverage ratio heatmap")
    print("Saved matched-clusters-per-100-clonotypes heatmap")
    print(f"Saved detailed cluster-dispersion plot for {args.donor_x} -> {args.donor_y}")


if __name__ == "__main__":
    main()
