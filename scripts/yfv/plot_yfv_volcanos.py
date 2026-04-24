from __future__ import annotations

import argparse
import math
import re
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
import seaborn as sns

matplotlib.use("Agg")
import matplotlib.pyplot as plt


SUMMARY_SUFFIXES = ("_summary_tcrempnet.tsv", "_summary_tcrempnet.tsv.gz")
CLONOTYPES_SUFFIXES = (
    "_tcremp_clusters.tsv",
    "_tcremp_clusters.tsv.gz",
    "_enriched_clonotypes_tcremp.tsv",
    "_enriched_clonotypes_tcremp.tsv.gz",
)
KNOWN_METHODS = ("leiden", "vdbscan")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Reproduce YFV volcano plots from RedCEA summary tables and optionally "
            "highlight clusters with exact YFV VDJdb clonotype matches."
        )
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=Path("results") / "yfv_redcea" / "raw",
        help="Directory with copied RedCEA YFV result tables.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results") / "yfv_redcea" / "volcanos",
        help="Directory for volcano figures and annotated summary tables.",
    )
    parser.add_argument(
        "--vdjdb",
        type=Path,
        default=Path("data") / "vdjdb" / "vdjdb.slim.txt",
        help="VDJdb table used to annotate YFV-matching clusters.",
    )
    parser.add_argument(
        "--pval-threshold",
        type=float,
        default=0.05,
        help="FDR threshold used to define significant clusters.",
    )
    parser.add_argument(
        "--fold-threshold",
        type=float,
        default=0.0,
        help="Minimum log fold-change threshold used to define enriched clusters.",
    )
    parser.add_argument(
        "--prefix",
        default="yfv_",
        help="Only process summary files whose basename starts with this prefix.",
    )
    return parser.parse_args()


def load_yfv_cdr3s(vdjdb_path: Path) -> set[str]:
    vdjdb = pd.read_csv(vdjdb_path, sep="\t")
    required = {"cdr3", "antigen.species"}
    missing = required.difference(vdjdb.columns)
    if missing:
        raise ValueError(f"VDJdb file is missing required columns: {sorted(missing)}")

    mask = vdjdb["antigen.species"].fillna("").str.upper().eq("YFV")
    if "gene" in vdjdb.columns:
        mask &= vdjdb["gene"].fillna("").str.upper().eq("TRB")

    cdr3s = vdjdb.loc[mask, "cdr3"].dropna().astype(str).str.strip().str.upper()
    return set(cdr3s)


def infer_sample_label(summary_path: Path, prefix: str) -> str:
    name = summary_path.name
    summary_suffix = next((suffix for suffix in SUMMARY_SUFFIXES if name.endswith(suffix)), None)
    if not name.startswith(prefix) or summary_suffix is None:
        return summary_path.stem
    core = name[len(prefix) : -len(summary_suffix)]
    core = re.sub(r"_(leiden|vdbscan)$", "", core)
    return core


def infer_method(summary_path: Path) -> str:
    stem = summary_path.name
    for suffix in SUMMARY_SUFFIXES:
        if stem.endswith(suffix):
            stem = stem[: -len(suffix)]
            break
    for method in KNOWN_METHODS:
        if stem.endswith(f"_{method}"):
            return method
    return "unknown"


def find_clonotypes_file(summary_path: Path) -> Path:
    for summary_suffix in SUMMARY_SUFFIXES:
        if not summary_path.name.endswith(summary_suffix):
            continue
        prefix = summary_path.name[: -len(summary_suffix)]
        for clonotypes_suffix in CLONOTYPES_SUFFIXES:
            expected = summary_path.with_name(f"{prefix}{clonotypes_suffix}")
            if expected.exists():
                return expected
    return summary_path.with_name(summary_path.stem)


def resolve_cdr3_column(df: pd.DataFrame) -> str | None:
    candidates = (
        "cdr3aa_beta",
        "junction_aa",
        "cdr3",
        "cdr3aa",
        "cdr3_b_aa",
        "cdr3_beta",
    )
    for candidate in candidates:
        if candidate in df.columns:
            return candidate
    return None


def annotate_with_vdjdb(summary_df: pd.DataFrame, clonotypes_path: Path, yfv_cdr3s: set[str]) -> pd.DataFrame:
    annotated = summary_df.copy()
    annotated["vdjdb"] = False

    if not clonotypes_path.exists():
        return annotated

    clonotypes = pd.read_csv(clonotypes_path, sep="\t")
    if "cluster_id" not in clonotypes.columns:
        return annotated

    cdr3_col = resolve_cdr3_column(clonotypes)
    if cdr3_col is None:
        return annotated

    matches = (
        clonotypes[cdr3_col]
        .fillna("")
        .astype(str)
        .str.strip()
        .str.upper()
        .isin(yfv_cdr3s)
    )
    matched_clusters = set(clonotypes.loc[matches, "cluster_id"].tolist())
    annotated["vdjdb"] = annotated["cluster_id"].isin(matched_clusters)
    return annotated


def prepare_volcano_frame(summary_df: pd.DataFrame, pval_threshold: float, fold_threshold: float) -> pd.DataFrame:
    df = summary_df.copy()
    required = {"cluster_id", "enrichment_pvalue_zbinom", "enrichment_fdr_zbinom", "log_fold_change"}
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(f"Summary table is missing required columns: {sorted(missing)}")

    eps = 1e-10
    df["log10_pval"] = -np.log10(df["enrichment_pvalue_zbinom"].clip(lower=0) + eps)
    df["significant"] = (df["enrichment_fdr_zbinom"] < pval_threshold) & (df["log_fold_change"] > fold_threshold)

    def assign_label(row: pd.Series) -> str:
        if row["significant"] and row["vdjdb"]:
            return "sign_vdjdb"
        if row["significant"]:
            return "sign"
        if row["vdjdb"]:
            return "vdjdb"
        return "no"

    df["label"] = df.apply(assign_label, axis=1)
    return df


def plot_volcano(df: pd.DataFrame, sample_name: str, pval_threshold: float, fold_threshold: float, ax) -> None:
    colors = {
        "sign_vdjdb": "#084594",
        "sign": "#4292C6",
        "vdjdb": "#FDBB84",
        "no": "#FFFFB2",
    }
    markers = {True: "s", False: "o"}

    sns.scatterplot(
        data=df,
        x="log_fold_change",
        y="log10_pval",
        hue="label",
        style="vdjdb",
        palette=colors,
        markers=markers,
        edgecolor="black",
        linewidth=0.3,
        s=70,
        alpha=0.85,
        ax=ax,
    )

    sig_mask = df["enrichment_fdr_zbinom"] < pval_threshold
    adj_pval = float(df.loc[sig_mask, "enrichment_pvalue_zbinom"].max()) if sig_mask.any() else pval_threshold

    ax.axvline(fold_threshold, ls="--", color="black")
    ax.axhline(-np.log10(pval_threshold), ls="--", color="black")
    ax.axhline(-np.log10(adj_pval + 1e-10), ls=":", color="blue")
    ax.set_xlabel("log10(Fold Enrichment)")
    ax.set_ylabel("-log10(p-value)")
    ax.set_title(f"{sample_name}, FDR<{pval_threshold}: {int(df['significant'].sum())}")
    ax.legend(title="", loc="best")


def process_sample(
    summary_path: Path,
    output_dir: Path,
    yfv_cdr3s: set[str],
    pval_threshold: float,
    fold_threshold: float,
    prefix: str,
) -> tuple[str, pd.DataFrame]:
    sample_label = infer_sample_label(summary_path, prefix=prefix)
    method = infer_method(summary_path)
    summary_df = pd.read_csv(summary_path, sep="\t")
    clonotypes_path = find_clonotypes_file(summary_path)
    annotated = annotate_with_vdjdb(summary_df, clonotypes_path, yfv_cdr3s)
    annotated["method"] = method
    volcano_df = prepare_volcano_frame(annotated, pval_threshold=pval_threshold, fold_threshold=fold_threshold)

    fig, ax = plt.subplots(figsize=(8, 6))
    title = sample_label if method == "unknown" else f"{sample_label} [{method}]"
    plot_volcano(volcano_df, sample_name=title, pval_threshold=pval_threshold, fold_threshold=fold_threshold, ax=ax)
    fig.tight_layout()
    base_name = sample_label if method == "unknown" else f"{sample_label}_{method}"
    fig.savefig(output_dir / f"{base_name}_volcano.png", dpi=200, bbox_inches="tight")
    fig.savefig(output_dir / f"{base_name}_volcano.svg", bbox_inches="tight")
    plt.close(fig)

    volcano_df.to_csv(output_dir / f"{base_name}_summary_annotated.tsv", sep="\t", index=False)
    return base_name, volcano_df


def plot_panel(
    sample_frames: list[tuple[str, pd.DataFrame]],
    output_dir: Path,
    pval_threshold: float,
    fold_threshold: float,
    panel_name: str,
) -> None:
    if not sample_frames:
        return

    n = len(sample_frames)
    ncols = min(3, n)
    nrows = math.ceil(n / ncols)
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(6.5 * ncols, 5.5 * nrows), squeeze=False)

    for ax, (sample_name, df) in zip(axes.flat, sample_frames):
        plot_volcano(df, sample_name=sample_name, pval_threshold=pval_threshold, fold_threshold=fold_threshold, ax=ax)

    for ax in axes.flat[n:]:
        ax.axis("off")

    fig.tight_layout()
    fig.savefig(output_dir / f"{panel_name}.png", dpi=200, bbox_inches="tight")
    fig.savefig(output_dir / f"{panel_name}.svg", bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    args = parse_args()
    input_dir = args.input_dir.resolve()
    output_dir = args.output_dir.resolve()
    vdjdb_path = args.vdjdb.resolve()

    if not input_dir.exists():
        raise FileNotFoundError(f"Input directory does not exist: {input_dir}")
    if not vdjdb_path.exists():
        raise FileNotFoundError(f"VDJdb file does not exist: {vdjdb_path}")

    summary_paths = sorted(
        path
        for suffix in SUMMARY_SUFFIXES
        for path in input_dir.glob(f"{args.prefix}*{suffix}")
        if path.is_file()
    )
    if not summary_paths:
        expected = ", ".join(f"{args.prefix}*{suffix}" for suffix in SUMMARY_SUFFIXES)
        print(f"No summary files matching [{expected}] were found in {input_dir}")
        return 1

    output_dir.mkdir(parents=True, exist_ok=True)
    if hasattr(sns, "set_theme"):
        sns.set_theme(style="whitegrid")
    else:
        sns.set_style("whitegrid")
    yfv_cdr3s = load_yfv_cdr3s(vdjdb_path)

    sample_frames: list[tuple[str, pd.DataFrame]] = []
    for summary_path in summary_paths:
        sample_label, volcano_df = process_sample(
            summary_path=summary_path,
            output_dir=output_dir,
            yfv_cdr3s=yfv_cdr3s,
            pval_threshold=args.pval_threshold,
            fold_threshold=args.fold_threshold,
            prefix=args.prefix,
        )
        sample_frames.append((sample_label, volcano_df))
        print(f"Rendered volcano for {sample_label}")

    by_method: dict[str, list[tuple[str, pd.DataFrame]]] = {}
    for sample_name, df in sample_frames:
        method = str(df["method"].iloc[0]) if "method" in df.columns and not df.empty else "unknown"
        by_method.setdefault(method, []).append((sample_name, df))

    plot_panel(
        sample_frames=sample_frames,
        output_dir=output_dir,
        pval_threshold=args.pval_threshold,
        fold_threshold=args.fold_threshold,
        panel_name="yfv_volcano_panel_all",
    )
    for method, frames in sorted(by_method.items()):
        plot_panel(
            sample_frames=frames,
            output_dir=output_dir,
            pval_threshold=args.pval_threshold,
            fold_threshold=args.fold_threshold,
            panel_name=f"yfv_volcano_panel_{method}",
        )
    print(f"Saved outputs to {output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
