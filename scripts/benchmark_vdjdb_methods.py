#!/usr/bin/env python
from __future__ import annotations

import argparse
import csv
import math
import subprocess
from collections import defaultdict
from pathlib import Path
from io import StringIO

import pandas as pd


EPITOPES = [
    ("GLCTLVAML", "GLC"),
    ("YLQPRTFLL", "YLQ"),
]
PADJ_THRESHOLD = 1e-5
ALL_DATASET_NAME = "ALL"


def normalize_segment(series: pd.Series, *, drop_allele: bool = True) -> pd.Series:
    normalized = series.fillna("").astype(str).str.split(",").str[0].str.strip()
    if drop_allele:
        normalized = normalized.str.split("*").str[0]
    return normalized


def load_vdjdb(path: Path, *, species: str, gene: str, epitopes: list[str] | None = None) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t").copy()
    mask = df["species"].eq(species) & df["gene"].eq(gene) & df["v.segm"].fillna("").ne("") & df["j.segm"].fillna("").ne("")
    if epitopes is not None:
        mask &= df["antigen.epitope"].isin(epitopes)
    df = df[mask].copy()
    df["cdr3aa"] = df["cdr3"].astype(str)
    # Keep allele annotations for method inputs such as tcrdist3, which expects
    # IMGT-style gene names like TRBV29-1*01 in its reference database.
    df["v.segm"] = normalize_segment(df["v.segm"], drop_allele=False)
    df["j.segm"] = normalize_segment(df["j.segm"], drop_allele=False)
    return (
        df[
            [
                "gene",
                "species",
                "antigen.epitope",
                "antigen.gene",
                "antigen.species",
                "mhc.a",
                "mhc.b",
                "mhc.class",
                "cdr3aa",
                "v.segm",
                "j.segm",
                "v.end",
                "j.start",
            ]
        ]
        .drop_duplicates()
        .reset_index(drop=True)
    )


def load_truth(path: Path, *, epitopes: list[str]) -> pd.DataFrame:
    truth = pd.read_csv(path).drop(columns=["Unnamed: 0"], errors="ignore")
    truth = truth[truth["epitope_aa"].isin(epitopes)].dropna(subset=["padj"]).copy()
    truth["valid"] = truth["padj"] < PADJ_THRESHOLD
    truth["TRBV"] = normalize_segment(truth["TRBV"])
    truth["TRBJ"] = normalize_segment(truth["TRBJ"])
    return truth


def write_giana_input(df: pd.DataFrame, output_path: Path) -> None:
    rows = []
    for idx, row in df.reset_index(drop=True).iterrows():
        rows.append([row["cdr3aa"], row["v.segm"], row["j.segm"], f"row_{idx}"])
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerows(rows)


def load_giana_output(path: Path) -> pd.DataFrame:
    valid_lines: list[str] = []
    with path.open() as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped:
                continue
            parts = stripped.split("\t")
            if len(parts) < 4:
                continue
            valid_lines.append(stripped)
    if not valid_lines:
        return pd.DataFrame(columns=["cdr3aa", "cluster", "v.segm", "j.segm"])
    raw = pd.read_csv(StringIO("\n".join(valid_lines)), sep="\t", header=None)
    if raw.shape[1] < 4:
        raise ValueError(f"Unexpected GIANA output shape after filtering: {raw.shape}")
    return raw


def mock_back_translate(seq: str) -> str:
    codons = {
        "A": "GCT",
        "C": "TGT",
        "D": "GAT",
        "E": "GAA",
        "F": "TTT",
        "G": "GGT",
        "H": "CAT",
        "I": "ATT",
        "K": "AAA",
        "L": "TTA",
        "M": "ATG",
        "N": "AAT",
        "P": "CCT",
        "Q": "CAA",
        "R": "CGT",
        "S": "TCT",
        "T": "ACT",
        "V": "GTT",
        "W": "TGG",
        "Y": "TAT",
    }
    return "".join(codons[aa] for aa in seq)


def write_vdjtools_input(df: pd.DataFrame, output_path: Path) -> None:
    tmp = df.copy()
    tmp["count"] = 1
    tmp["freq"] = 1.0 / max(len(tmp), 1)
    tmp["cdr3nt"] = tmp["cdr3aa"].map(mock_back_translate)
    tmp["v"] = tmp["v.segm"]
    tmp["d"] = ""
    tmp["j"] = tmp["j.segm"]
    tmp["vend"] = tmp["v.end"].fillna(-1).astype(int)
    tmp["dstart"] = -1
    tmp["dend"] = -1
    tmp["jstart"] = tmp["j.start"].fillna(-1).astype(int)
    cols = ["count", "freq", "cdr3nt", "cdr3aa", "v", "d", "j", "vend", "dstart", "dend", "jstart"]
    output_path.parent.mkdir(parents=True, exist_ok=True)
    tmp[cols].to_csv(output_path, sep="\t", index=False)


def prepare_inputs(args: argparse.Namespace) -> None:
    epitopes = [ep for ep, _ in EPITOPES]
    df = load_vdjdb(args.vdjdb, species=args.species, gene=args.gene)
    truth = load_truth(args.truth, epitopes=epitopes)

    input_root = args.work_dir / "inputs"
    input_root.mkdir(parents=True, exist_ok=True)
    truth.to_csv(args.work_dir / "truth_glc_ylq.csv", index=False)

    generic_path = input_root / "generic" / f"{ALL_DATASET_NAME}.tsv"
    giana_path = input_root / "giana" / f"{ALL_DATASET_NAME}.tsv"

    generic_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(generic_path, sep="\t", index=False)
    write_giana_input(df, giana_path)

    print(f"Prepared inputs in {input_root}")
    print(f"Saved truth table to {args.work_dir / 'truth_glc_ylq.csv'}")


def hamming_distance(seq1: str, seq2: str) -> int | None:
    if len(seq1) != len(seq2):
        return None
    return sum(ch1 != ch2 for ch1, ch2 in zip(seq1, seq2))


def connected_components(n_nodes: int, edges: list[tuple[int, int]]) -> list[list[int]]:
    graph: dict[int, list[int]] = {idx: [] for idx in range(n_nodes)}
    for left, right in edges:
        graph[left].append(right)
        graph[right].append(left)

    seen: set[int] = set()
    comps: list[list[int]] = []
    for start in range(n_nodes):
        if start in seen:
            continue
        stack = [start]
        seen.add(start)
        comp: list[int] = []
        while stack:
            node = stack.pop()
            comp.append(node)
            for nxt in graph[node]:
                if nxt not in seen:
                    seen.add(nxt)
                    stack.append(nxt)
        comps.append(sorted(comp))
    return comps


def save_cluster_members(
    input_df: pd.DataFrame,
    cluster_map: dict[int, str],
    output_path: Path,
    *,
    dataset_name: str,
    method: str,
) -> None:
    columns = ["gene", "cdr3aa", "v.segm", "j.segm", "cid", "antigen.epitope", "method", "dataset_name"]
    rows = []
    for idx, cid in cluster_map.items():
        row = input_df.iloc[idx]
        rows.append(
            {
                "gene": row["gene"],
                "cdr3aa": row["cdr3aa"],
                "v.segm": row["v.segm"],
                "j.segm": row["j.segm"],
                "cid": cid,
                "antigen.epitope": row["antigen.epitope"],
                "method": method,
                "dataset_name": dataset_name,
            }
        )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows, columns=columns).to_csv(output_path, sep="\t", index=False)


def run_tcrdist3(args: argparse.Namespace) -> None:
    valid_genes: set[str] | None = None
    from tcrdist.repertoire import TCRrep
    try:
        from tcrdist.all_genes import all_genes

        valid_genes = set(all_genes["human"].keys())
    except Exception:
        valid_genes = None

    input_df = pd.read_csv(args.input, sep="\t").copy()
    input_df = input_df.dropna(subset=["cdr3aa", "v.segm", "j.segm"]).copy()
    input_df["cdr3aa"] = input_df["cdr3aa"].astype(str).str.strip()
    input_df["v.segm"] = input_df["v.segm"].astype(str).str.strip()
    input_df["j.segm"] = input_df["j.segm"].astype(str).str.strip()

    # tcrdist3 cannot handle rows whose V/J genes are not present in its reference
    # database; they later surface as None values during sparse distance building.
    # Filter aggressively to canonical IMGT-looking beta-chain names with alleles.
    valid_gene_pattern = r"^TRB[VDJ][A-Z0-9/-]*\*\d+$"
    mask = (
        input_df["cdr3aa"].ne("")
        & input_df["v.segm"].str.match(valid_gene_pattern, na=False)
        & input_df["j.segm"].str.match(valid_gene_pattern, na=False)
    )
    if valid_genes is not None:
        mask &= input_df["v.segm"].isin(valid_genes) & input_df["j.segm"].isin(valid_genes)
    filtered_out = int((~mask).sum())
    if filtered_out:
        print(f"Filtered out {filtered_out} rows before tcrdist3 due to missing/non-canonical V/J annotations")
    input_df = input_df[mask].reset_index(drop=True)
    if input_df.empty:
        save_cluster_members(input_df, {}, args.output, dataset_name=args.dataset_name, method="tcrdist3")
        print(f"Saved {args.output}")
        return

    tcrdist_df = (
        input_df.rename(
            columns={
                "cdr3aa": "cdr3_b_aa",
                "v.segm": "v_b_gene",
                "j.segm": "j_b_gene",
            }
        )[["cdr3_b_aa", "v_b_gene", "j_b_gene"]]
        .copy()
    )
    tcrdist_df["count"] = 1

    tr = TCRrep(
        cell_df=tcrdist_df,
        organism="human",
        chains=["beta"],
        db_file=args.db_file,
        deduplicate=False,
        compute_distances=False,
    )
    tr.cpus = args.cpus
    tr.compute_sparse_rect_distances(radius=args.radius, chunk_size=args.chunk_size)
    sparse = tr.rw_beta.tocoo()

    edges: list[tuple[int, int]] = []
    for i, j, value in zip(sparse.row, sparse.col, sparse.data):
        if i >= j:
            continue
        # In tcrdist3 sparse matrices, exact zero distances are encoded as -1.
        if value == -1 or value <= args.radius:
            edges.append((int(i), int(j)))

    cluster_map: dict[int, str] = {}
    cid_counter = 1
    for comp in connected_components(len(input_df), edges):
        if len(comp) < args.min_cluster_size:
            continue
        cid = f"H.B.{cid_counter}"
        cid_counter += 1
        for idx in comp:
            cluster_map[idx] = cid

    save_cluster_members(input_df, cluster_map, args.output, dataset_name=args.dataset_name, method="tcrdist3")
    print(f"Saved {args.output}")


def convert_giana(args: argparse.Namespace) -> None:
    input_df = pd.read_csv(args.input, sep="\t").copy()
    raw = load_giana_output(args.giana_output)

    raw = raw.rename(columns={0: "cdr3aa", 1: "cluster", 2: "v.segm", 3: "j.segm"})
    raw["cluster"] = raw["cluster"].astype(str)
    valid = raw[raw["cluster"].str.strip().ne("") & ~raw["cluster"].isin(["-1", "nan", "NA"])]
    counts = valid["cluster"].value_counts()
    keep = set(counts[counts >= args.min_cluster_size].index)

    cluster_map: dict[int, str] = {}
    for idx, row in valid.iterrows():
        cluster_id = row["cluster"]
        if cluster_id not in keep:
            continue
        cluster_map[idx] = f"H.B.{cluster_id}"

    save_cluster_members(input_df, cluster_map, args.output, dataset_name=args.dataset_name, method="GIANA")
    print(f"Saved {args.output}")


def convert_tcrnet(args: argparse.Namespace) -> None:
    input_df = pd.read_csv(args.input, sep="\t").copy()
    annotated = pd.read_csv(args.annotated, sep="\t").copy()
    if "p.value.g" not in annotated.columns or "degree.s" not in annotated.columns:
        raise ValueError("TCRNet annotated file must contain 'p.value.g' and 'degree.s' columns")

    annotated["enriched"] = (annotated["degree.s"] >= args.min_degree) & (annotated["p.value.g"] < args.pvalue_threshold)
    seqs = annotated["cdr3aa"].astype(str).tolist()
    enriched_idx = [idx for idx, value in enumerate(annotated["enriched"].tolist()) if value]

    edge_set: set[tuple[int, int]] = set()
    for i in enriched_idx:
        for j in range(len(seqs)):
            if i == j:
                continue
            dist = hamming_distance(seqs[i], seqs[j])
            if dist == 1:
                edge_set.add(tuple(sorted((i, j))))

    expanded_nodes = sorted({node for edge in edge_set for node in edge})
    index_remap = {old: new for new, old in enumerate(expanded_nodes)}
    comps = connected_components(len(expanded_nodes), [(index_remap[a], index_remap[b]) for a, b in edge_set])

    cluster_map: dict[int, str] = {}
    cid_counter = 1
    for comp in comps:
        original = [expanded_nodes[idx] for idx in comp]
        if len(original) < args.min_cluster_size:
            continue
        cid = f"H.B.{cid_counter}"
        cid_counter += 1
        for idx in original:
            cluster_map[idx] = cid

    save_cluster_members(input_df, cluster_map, args.output, dataset_name=args.dataset_name, method="TCRnet")
    print(f"Saved {args.output}")


def compute_metrics(y_true: pd.Series, y_pred: pd.Series) -> dict[str, float | int]:
    y_true = y_true.astype(bool)
    y_pred = y_pred.astype(bool)
    tp = int((y_true & y_pred).sum())
    fp = int((~y_true & y_pred).sum())
    fn = int((y_true & ~y_pred).sum())
    precision = tp / (tp + fp) if tp + fp else 0.0
    recall = tp / (tp + fn) if tp + fn else 0.0
    f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
    return {
        "precision": precision,
        "recall": recall,
        "f1": f1,
        "tp": tp,
        "fp": fp,
        "fn": fn,
        "predicted_positive": int(y_pred.sum()),
        "true_positive_total": int(y_true.sum()),
    }


def evaluate(args: argparse.Namespace) -> None:
    truth = load_truth(args.truth, epitopes=[ep for ep, _ in EPITOPES])
    rows: list[dict[str, object]] = []

    for method_dir in sorted(args.results_root.iterdir()):
        if not method_dir.is_dir():
            continue
        members_path = method_dir / "cluster_members_TRB.txt"
        if not members_path.exists():
            continue

        members = pd.read_csv(members_path, sep="\t").copy()
        members["v.segm"] = normalize_segment(members["v.segm"])
        members["j.segm"] = normalize_segment(members["j.segm"])

        for epitope, short_name in EPITOPES:
            truth_ep = truth[truth["epitope_aa"].eq(epitope)].copy()
            cluster_ids = set(members.loc[members["antigen.epitope"].eq(epitope), "cid"].dropna().astype(str))
            cluster_scope = members[members["cid"].astype(str).isin(cluster_ids)].copy()

            cdr3_signatures = set(cluster_scope["cdr3aa"].astype(str))
            vj_signatures = set(
                cluster_scope["cdr3aa"].astype(str)
                + "|"
                + cluster_scope["v.segm"].astype(str)
                + "|"
                + cluster_scope["j.segm"].astype(str)
            )

            beta_key = truth_ep["cdr3_beta_aa"].fillna("").astype(str)
            beta_key_vj = beta_key + "|" + truth_ep["TRBV"].fillna("").astype(str) + "|" + truth_ep["TRBJ"].fillna("").astype(str)

            for with_vj, y_pred in [
                (False, beta_key.isin(cdr3_signatures)),
                (True, beta_key_vj.isin(vj_signatures)),
            ]:
                rows.append(
                    {
                        "method": method_dir.name,
                        "epitope": epitope,
                        "epitope_short": short_name,
                        "n_clusters": int(len(cluster_ids)),
                        "cluster_member_total": int(len(cluster_scope)),
                        "cluster_member_target_epitope": int(cluster_scope["antigen.epitope"].eq(epitope).sum()),
                        "cluster_member_other_epitope": int(cluster_scope["antigen.epitope"].ne(epitope).sum()),
                        "match_mode": "cdr3+vj" if with_vj else "cdr3",
                        "total": int(len(truth_ep)),
                        "positives": int(truth_ep["valid"].sum()),
                        **compute_metrics(truth_ep["valid"], y_pred),
                    }
                )

    detailed = pd.DataFrame(rows).sort_values(["method", "match_mode", "epitope_short"]).reset_index(drop=True)
    summary = (
        detailed.groupby(["method", "match_mode"], as_index=False)
        .agg(
            mean_precision=("precision", "mean"),
            mean_recall=("recall", "mean"),
            mean_f1=("f1", "mean"),
            min_f1=("f1", "min"),
        )
        .sort_values(["mean_f1", "mean_precision", "mean_recall"], ascending=False)
        .reset_index(drop=True)
    )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    detailed_path = args.output.with_name(args.output.stem + "_detailed.csv")
    summary_path = args.output.with_name(args.output.stem + "_summary.csv")
    detailed.to_csv(detailed_path, index=False)
    summary.to_csv(summary_path, index=False)
    print(f"Saved {detailed_path}")
    print(f"Saved {summary_path}")
    if not summary.empty:
        print(summary.to_string(index=False))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Benchmark VDJdb clustering methods on GLC/YLQ TRB data.")
    subparsers = parser.add_subparsers(dest="command", required=True)

    prepare_parser = subparsers.add_parser("prepare")
    prepare_parser.add_argument("--vdjdb", type=Path, required=True)
    prepare_parser.add_argument("--truth", type=Path, required=True)
    prepare_parser.add_argument("--work-dir", type=Path, required=True)
    prepare_parser.add_argument("--species", default="HomoSapiens")
    prepare_parser.add_argument("--gene", default="TRB")

    tcrdist_parser = subparsers.add_parser("run-tcrdist3")
    tcrdist_parser.add_argument("--input", type=Path, required=True)
    tcrdist_parser.add_argument("--output", type=Path, required=True)
    tcrdist_parser.add_argument("--dataset-name", default=ALL_DATASET_NAME)
    tcrdist_parser.add_argument("--radius", type=int, default=24)
    tcrdist_parser.add_argument("--min-cluster-size", type=int, default=3)
    tcrdist_parser.add_argument("--db-file", default="alphabeta_gammadelta_db.tsv")
    tcrdist_parser.add_argument("--cpus", type=int, default=1)
    tcrdist_parser.add_argument("--chunk-size", type=int, default=100)

    giana_parser = subparsers.add_parser("convert-giana")
    giana_parser.add_argument("--input", type=Path, required=True)
    giana_parser.add_argument("--giana-output", type=Path, required=True)
    giana_parser.add_argument("--output", type=Path, required=True)
    giana_parser.add_argument("--dataset-name", default=ALL_DATASET_NAME)
    giana_parser.add_argument("--min-cluster-size", type=int, default=3)

    tcrnet_parser = subparsers.add_parser("convert-tcrnet")
    tcrnet_parser.add_argument("--input", type=Path, required=True)
    tcrnet_parser.add_argument("--annotated", type=Path, required=True)
    tcrnet_parser.add_argument("--output", type=Path, required=True)
    tcrnet_parser.add_argument("--dataset-name", default=ALL_DATASET_NAME)
    tcrnet_parser.add_argument("--min-cluster-size", type=int, default=5)
    tcrnet_parser.add_argument("--min-degree", type=int, default=2)
    tcrnet_parser.add_argument("--pvalue-threshold", type=float, default=0.05)

    eval_parser = subparsers.add_parser("evaluate")
    eval_parser.add_argument("--truth", type=Path, required=True)
    eval_parser.add_argument("--results-root", type=Path, required=True)
    eval_parser.add_argument("--output", type=Path, required=True)

    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.command == "prepare":
        prepare_inputs(args)
    elif args.command == "run-tcrdist3":
        run_tcrdist3(args)
    elif args.command == "convert-giana":
        convert_giana(args)
    elif args.command == "convert-tcrnet":
        convert_tcrnet(args)
    elif args.command == "evaluate":
        evaluate(args)
    else:
        raise ValueError(f"Unsupported command: {args.command}")


if __name__ == "__main__":
    main()
