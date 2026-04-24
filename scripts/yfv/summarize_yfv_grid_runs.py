from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd


SUMMARY_SUFFIX = "_summary_tcrempnet.tsv"
CLUSTERS_SUFFIX = "_tcremp_clusters.tsv"
ENRICHED_SUFFIX = "_enriched_clonotypes_tcremp.tsv"
LEIDEN_RUN_RE = re.compile(
    r"^yfv_(?P<donor>[A-Z]\d+)_(?P<replica>F\d+)_leiden_"
    r"k(?P<k_neighbors>\d+)_res(?P<resolution>[A-Za-z0-9p.-]+)_min(?P<min_samples>\d+)$"
)
DBSCAN_RUN_RE = re.compile(
    r"^yfv_(?P<donor>[A-Z]\d+)_(?P<replica>F\d+)_dbscan_"
    r"k(?P<k_neighbors>\d+)_ek(?P<eps_k_neighbors>\d+)_min(?P<min_samples>\d+)$"
)
VDBSCAN_RUN_RE = re.compile(
    r"^yfv_(?P<donor>[A-Z]\d+)_(?P<replica>F\d+)_vdbscan_"
    r"k(?P<k_neighbors>\d+)_ek(?P<eps_k_neighbors>\d+)_min(?P<min_samples>\d+)"
    r"_eps(?P<eps_estimation_based_on>[A-Za-z0-9_-]+)_sym(?P<vdbscan_sym_rule>[A-Za-z0-9_-]+)$"
)
CDR3_CANDIDATES = ("cdr3aa_beta", "junction_aa", "cdr3", "cdr3aa", "cdr3_b_aa", "cdr3_beta")


def find_repo_root() -> Path:
    current = Path.cwd().resolve()
    for candidate in [current] + list(current.parents):
        if (candidate / "scripts").exists() and (candidate / "data").exists():
            return candidate
    raise FileNotFoundError("Could not locate repository root from current working directory.")


def parse_args() -> argparse.Namespace:
    root = find_repo_root()
    parser = argparse.ArgumentParser(
        description=(
            "Summarize YFV run directories against VDJdb LLW clonotypes. "
            "Supports Leiden, DBSCAN, and vDBSCAN run naming schemes."
        )
    )
    parser.add_argument(
        "--runs-root",
        type=Path,
        default=root / "results" / "yfv_redcea" / "runs",
        help="Directory containing YFV run subdirectories.",
    )
    parser.add_argument(
        "--vdjdb",
        type=Path,
        default=root / "data" / "vdjdb" / "vdjdb.slim.txt",
        help="VDJdb table.",
    )
    parser.add_argument(
        "--epitope",
        default="LLWNGPMAV",
        help="VDJdb epitope to track.",
    )
    parser.add_argument(
        "--species",
        default="YFV",
        help="VDJdb antigen.species value to track.",
    )
    parser.add_argument(
        "--fdr-threshold",
        type=float,
        default=0.05,
        help="FDR threshold used to define enriched vs not enriched.",
    )
    parser.add_argument(
        "--donors",
        default="",
        help="Optional comma-separated donor filter, e.g. S1,Q1",
    )
    parser.add_argument(
        "--algos",
        default="",
        help="Optional comma-separated algo filter, e.g. dbscan or leiden,dbscan,vdbscan",
    )
    parser.add_argument(
        "--output-prefix",
        type=Path,
        default=root / "results" / "yfv_redcea" / "yfv_grid_llw_summary",
        help="Prefix for output TSV files.",
    )
    return parser.parse_args()


def normalize_segment(series: pd.Series) -> pd.Series:
    return (
        series.fillna("")
        .astype(str)
        .str.split(",")
        .str[0]
        .str.strip()
        .str.split("*")
        .str[0]
        .str.upper()
    )


def resolve_cdr3_column(df: pd.DataFrame) -> str:
    for column in CDR3_CANDIDATES:
        if column in df.columns:
            return column
    raise ValueError(f"Could not find a CDR3 column among {CDR3_CANDIDATES}")


def parse_resolution(tag: str) -> float:
    return float(tag.replace("p", "."))


def parse_run_name(run_name: str) -> dict[str, object] | None:
    match = LEIDEN_RUN_RE.match(run_name)
    if match is not None:
        return {
            "cluster_algo": "leiden",
            "donor": match.group("donor"),
            "replica": match.group("replica"),
            "k_neighbors": int(match.group("k_neighbors")),
            "cluster_min_samples": int(match.group("min_samples")),
            "leiden_resolution": parse_resolution(match.group("resolution")),
            "resolution_tag": match.group("resolution"),
            "eps_k_neighbors": None,
            "eps_estimation_based_on": None,
            "vdbscan_sym_rule": None,
            "param_id": (
                f"leiden|k{match.group('k_neighbors')}|res{match.group('resolution')}|"
                f"min{match.group('min_samples')}"
            ),
        }

    match = DBSCAN_RUN_RE.match(run_name)
    if match is not None:
        return {
            "cluster_algo": "dbscan",
            "donor": match.group("donor"),
            "replica": match.group("replica"),
            "k_neighbors": int(match.group("k_neighbors")),
            "cluster_min_samples": int(match.group("min_samples")),
            "leiden_resolution": None,
            "resolution_tag": None,
            "eps_k_neighbors": int(match.group("eps_k_neighbors")),
            "eps_estimation_based_on": None,
            "vdbscan_sym_rule": None,
            "param_id": (
                f"dbscan|k{match.group('k_neighbors')}|ek{match.group('eps_k_neighbors')}|"
                f"min{match.group('min_samples')}"
            ),
        }

    match = VDBSCAN_RUN_RE.match(run_name)
    if match is not None:
        return {
            "cluster_algo": "vdbscan",
            "donor": match.group("donor"),
            "replica": match.group("replica"),
            "k_neighbors": int(match.group("k_neighbors")),
            "cluster_min_samples": int(match.group("min_samples")),
            "leiden_resolution": None,
            "resolution_tag": None,
            "eps_k_neighbors": int(match.group("eps_k_neighbors")),
            "eps_estimation_based_on": match.group("eps_estimation_based_on"),
            "vdbscan_sym_rule": match.group("vdbscan_sym_rule"),
            "param_id": (
                f"vdbscan|k{match.group('k_neighbors')}|ek{match.group('eps_k_neighbors')}|"
                f"min{match.group('min_samples')}|eps{match.group('eps_estimation_based_on')}|"
                f"sym{match.group('vdbscan_sym_rule')}"
            ),
        }
    return None


def load_llw_vdjdb(vdjdb_path: Path, *, epitope: str, species: str) -> pd.DataFrame:
    vdjdb = pd.read_csv(vdjdb_path, sep="\t")
    mask = (
        vdjdb["antigen.epitope"].fillna("").astype(str).str.upper().eq(epitope.upper())
        & vdjdb["antigen.species"].fillna("").astype(str).str.upper().eq(species.upper())
        & vdjdb["gene"].fillna("").astype(str).str.upper().eq("TRB")
    )
    llw = vdjdb.loc[mask].copy()
    if llw.empty:
        raise ValueError(f"No VDJdb TRB rows found for antigen.epitope={epitope!r}, antigen.species={species!r}")

    llw["cdr3_norm"] = llw["cdr3"].fillna("").astype(str).str.strip().str.upper()
    llw["v_norm"] = normalize_segment(llw["v.segm"])
    llw["j_norm"] = normalize_segment(llw["j.segm"])
    return llw


def discover_run_dirs(
    runs_root: Path,
    *,
    donors_filter: set[str] | None,
    algos_filter: set[str] | None,
) -> list[tuple[Path, dict[str, object]]]:
    if not runs_root.exists():
        raise FileNotFoundError(f"Runs root does not exist: {runs_root}")

    run_dirs: list[tuple[Path, dict[str, object]]] = []
    for path in sorted(runs_root.iterdir()):
        if not path.is_dir():
            continue
        parsed = parse_run_name(path.name)
        if parsed is None:
            continue
        if donors_filter is not None and str(parsed["donor"]) not in donors_filter:
            continue
        if algos_filter is not None and str(parsed["cluster_algo"]) not in algos_filter:
            continue
        if any(path.glob(f"*{SUMMARY_SUFFIX}")):
            run_dirs.append((path, parsed))

    if not run_dirs:
        raise FileNotFoundError(f"No matching YFV run directories found under {runs_root}")
    return run_dirs


def classify_zone(fdr: float, logfc: float, *, fdr_threshold: float) -> str:
    enriched = pd.notna(fdr) and float(fdr) < fdr_threshold
    positive = pd.notna(logfc) and float(logfc) >= 0.0
    if enriched and positive:
        return "enriched_positive"
    if enriched and not positive:
        return "enriched_negative"
    if not enriched and positive:
        return "not_enriched_positive"
    return "not_enriched_negative"


def summarize_run(
    run_dir: Path,
    run_meta: dict[str, object],
    *,
    llw_vdjdb: pd.DataFrame,
    fdr_threshold: float,
) -> dict[str, object]:
    summary_path = next(run_dir.glob(f"*{SUMMARY_SUFFIX}"))
    prefix = summary_path.name[: -len(SUMMARY_SUFFIX)]
    clusters_path = run_dir / f"{prefix}{CLUSTERS_SUFFIX}"
    enriched_path = run_dir / f"{prefix}{ENRICHED_SUFFIX}"

    if not clusters_path.exists():
        raise FileNotFoundError(f"Clusters file not found for {run_dir.name}: {clusters_path}")
    if not enriched_path.exists():
        raise FileNotFoundError(f"Enriched clonotypes file not found for {run_dir.name}: {enriched_path}")

    summary_df = pd.read_csv(summary_path, sep="\t")
    clusters_df = pd.read_csv(clusters_path, sep="\t")
    enriched_df = pd.read_csv(enriched_path, sep="\t")

    cdr3_clusters = resolve_cdr3_column(clusters_df)
    cdr3_enriched = resolve_cdr3_column(enriched_df)

    summary_df["cluster_id"] = summary_df["cluster_id"].astype(str)
    clusters_df["cluster_id"] = clusters_df["cluster_id"].astype(str)
    enriched_df["cluster_id"] = enriched_df["cluster_id"].astype(str)
    clusters_df["cdr3_norm"] = clusters_df[cdr3_clusters].fillna("").astype(str).str.strip().str.upper()
    enriched_df["cdr3_norm"] = enriched_df[cdr3_enriched].fillna("").astype(str).str.strip().str.upper()

    llw_cdr3s = set(llw_vdjdb["cdr3_norm"])
    matched_cluster_rows = clusters_df[clusters_df["cdr3_norm"].isin(llw_cdr3s)].copy()
    matched_clusters = set(matched_cluster_rows["cluster_id"])
    matched_enriched_rows = enriched_df[enriched_df["cdr3_norm"].isin(llw_cdr3s)].copy()

    vdjdb_cluster_summary = summary_df[summary_df["cluster_id"].isin(matched_clusters)].copy()
    vdjdb_cluster_summary["zone"] = [
        classify_zone(fdr, logfc, fdr_threshold=fdr_threshold)
        for fdr, logfc in zip(
            vdjdb_cluster_summary["enrichment_fdr_zbinom"],
            vdjdb_cluster_summary["log_fold_change"],
        )
    ]
    zone_counts = vdjdb_cluster_summary["zone"].value_counts().to_dict()

    result = {
        "run_dir": str(run_dir),
        "run_name": run_dir.name,
        "prefix": prefix,
        "n_clusters_total": int(summary_df["cluster_id"].nunique()),
        "mean_cluster_size": float(summary_df["cluster_size"].mean()),
        "median_cluster_size": float(summary_df["cluster_size"].median()),
        "n_enriched_positive_clusters_total": int(
            ((summary_df["enrichment_fdr_zbinom"] < fdr_threshold) & (summary_df["log_fold_change"] > 0)).sum()
        ),
        "n_llw_vdjdb_trb_clonotypes_total": int(len(llw_cdr3s)),
        "n_llw_vdjdb_trb_clonotypes_found_in_enriched": int(matched_enriched_rows["cdr3_norm"].nunique()),
        "n_llw_vdjdb_trb_enriched_rows": int(len(matched_enriched_rows)),
        "n_llw_vdjdb_clusters": int(len(matched_clusters)),
        "mean_log_fold_change_all_clusters": float(summary_df["log_fold_change"].mean()),
        "mean_log_fold_change_llw_vdjdb_clusters": float(vdjdb_cluster_summary["log_fold_change"].mean())
        if not vdjdb_cluster_summary.empty
        else float("nan"),
        "mean_cluster_size_llw_vdjdb_clusters": float(vdjdb_cluster_summary["cluster_size"].mean())
        if not vdjdb_cluster_summary.empty
        else float("nan"),
        "llw_vdjdb_clusters_enriched_positive": int(zone_counts.get("enriched_positive", 0)),
        "llw_vdjdb_clusters_enriched_negative": int(zone_counts.get("enriched_negative", 0)),
        "llw_vdjdb_clusters_not_enriched_positive": int(zone_counts.get("not_enriched_positive", 0)),
        "llw_vdjdb_clusters_not_enriched_negative": int(zone_counts.get("not_enriched_negative", 0)),
    }
    result.update(run_meta)
    return result


def build_param_summary(detailed: pd.DataFrame) -> pd.DataFrame:
    detailed = detailed.copy()
    for column, fill_value in [
        ("leiden_resolution", "NA"),
        ("resolution_tag", "NA"),
        ("eps_k_neighbors", "NA"),
        ("eps_estimation_based_on", "NA"),
        ("vdbscan_sym_rule", "NA"),
    ]:
        detailed[column] = detailed[column].where(pd.notna(detailed[column]), fill_value)

    group_cols = [
        "cluster_algo",
        "k_neighbors",
        "cluster_min_samples",
        "leiden_resolution",
        "resolution_tag",
        "eps_k_neighbors",
        "eps_estimation_based_on",
        "vdbscan_sym_rule",
        "param_id",
    ]
    summary = (
        detailed.groupby(group_cols)
        .agg(
            runs=("run_name", "count"),
            donors=("donor", lambda x: ",".join(sorted(set(x)))),
            mean_n_clusters_total=("n_clusters_total", "mean"),
            mean_cluster_size=("mean_cluster_size", "mean"),
            mean_n_enriched_positive_clusters_total=("n_enriched_positive_clusters_total", "mean"),
            mean_llw_found_in_enriched=("n_llw_vdjdb_trb_clonotypes_found_in_enriched", "mean"),
            max_llw_found_in_enriched=("n_llw_vdjdb_trb_clonotypes_found_in_enriched", "max"),
            mean_llw_vdjdb_clusters=("n_llw_vdjdb_clusters", "mean"),
            mean_log_fold_change_all_clusters=("mean_log_fold_change_all_clusters", "mean"),
            mean_log_fold_change_llw_vdjdb_clusters=("mean_log_fold_change_llw_vdjdb_clusters", "mean"),
            mean_llw_clusters_enriched_positive=("llw_vdjdb_clusters_enriched_positive", "mean"),
            mean_llw_clusters_enriched_negative=("llw_vdjdb_clusters_enriched_negative", "mean"),
            mean_llw_clusters_not_enriched_positive=("llw_vdjdb_clusters_not_enriched_positive", "mean"),
            mean_llw_clusters_not_enriched_negative=("llw_vdjdb_clusters_not_enriched_negative", "mean"),
        )
        .reset_index()
        .sort_values(
            [
                "mean_llw_found_in_enriched",
                "mean_llw_clusters_enriched_positive",
                "mean_log_fold_change_llw_vdjdb_clusters",
            ],
            ascending=[False, False, False],
        )
        .reset_index(drop=True)
    )
    summary.insert(0, "rank", range(1, len(summary) + 1))
    return summary


def parse_filter_list(raw: str) -> set[str] | None:
    values = [item.strip() for item in raw.split(",") if item.strip()]
    return set(values) if values else None


def main() -> None:
    args = parse_args()
    llw_vdjdb = load_llw_vdjdb(args.vdjdb, epitope=args.epitope, species=args.species)
    donors_filter = parse_filter_list(args.donors)
    algos_filter = parse_filter_list(args.algos)
    run_dirs = discover_run_dirs(args.runs_root, donors_filter=donors_filter, algos_filter=algos_filter)

    rows = [
        summarize_run(run_dir, run_meta, llw_vdjdb=llw_vdjdb, fdr_threshold=args.fdr_threshold)
        for run_dir, run_meta in run_dirs
    ]
    detailed = pd.DataFrame(rows).sort_values(
        ["cluster_algo", "donor", "k_neighbors", "leiden_resolution", "eps_k_neighbors", "cluster_min_samples"],
        na_position="last",
    ).reset_index(drop=True)
    summary = build_param_summary(detailed)

    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    detailed_path = args.output_prefix.with_name(args.output_prefix.name + "_per_run.tsv")
    summary_path = args.output_prefix.with_name(args.output_prefix.name + "_by_params.tsv")
    detailed.to_csv(detailed_path, sep="\t", index=False)
    summary.to_csv(summary_path, sep="\t", index=False)

    print(f"Analyzed {len(detailed)} runs under {args.runs_root}")
    print(f"Tracked {llw_vdjdb['cdr3_norm'].nunique()} unique VDJdb TRB clonotypes for {args.epitope}/{args.species}")
    print(f"Saved {detailed_path}")
    print(f"Saved {summary_path}")
    print("\nTop parameter sets:")
    print(summary.head(10).to_string(index=False))


if __name__ == "__main__":
    main()
