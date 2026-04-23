#!/usr/bin/env python
from __future__ import annotations

import argparse
import csv
import math
import subprocess
from collections import defaultdict
from pathlib import Path

import pandas as pd


EPITOPES = [
    ("GLCTLVAML", "GLC"),
    ("YLQPRTFLL", "YLQ"),
]
PADJ_THRESHOLD = 1e-5


def normalize_segment(series: pd.Series, *, drop_allele: bool = True) -> pd.Series:
    normalized = series.fillna("").astype(str).str.split(",").str[0].str.strip()
    if drop_allele:
        normalized = normalized.str.split("*").str[0]
    return normalized


def load_vdjdb(path: Path, *, species: str, gene: str, epitopes: list[str]) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t").copy()
    df = df[
        df["species"].eq(species)
        & df["gene"].eq(gene)
        & df["antigen.epitope"].isin(epitopes)
        & df["v.segm"].fillna("").ne("")
        & df["j.segm"].fillna("").ne("")
    ].copy()
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
    df = load_vdjdb(args.vdjdb, species=args.species, gene=args.gene, epitopes=epitopes)
    truth = load_truth(args.truth, epitopes=epitopes)

    input_root = args.work_dir / "inputs"
    input_root.mkdir(parents=True, exist_ok=True)
    truth.to_csv(args.work_dir / "truth_glc_ylq.csv", index=False)

    for epitope, short_name in EPITOPES:
        ep_df = df[df["antigen.epitope"].eq(epitope)].copy().reset_index(drop=True)
        generic_path = input_root / "generic" / f"{short_name}.tsv"
        giana_path = input_root / "giana" / f"{short_name}.tsv"

        generic_path.parent.mkdir(parents=True, exist_ok=True)
        ep_df.to_csv(generic_path, sep="\t", index=False)
        write_giana_input(ep_df, giana_path)

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
    short_name: str,
    method: str,
) -> None:
    columns = ["gene", "cdr3aa", "v.segm", "j.segm", "cid", "antigen.epitope", "method", "epitope_short"]
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
                "epitope_short": short_name,
            }
        )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows, columns=columns).to_csv(output_path, sep="\t", index=False)


def run_tcrdist3(args: argparse.Namespace) -> None:
    from tcrdist.repertoire import TCRrep

    input_df = pd.read_csv(args.input, sep="\t").copy()
    short_name = args.short_name
    epitope = input_df["antigen.epitope"].iat[0]

    tr = TCRrep(
        cell_df=input_df.rename(
            columns={
                "cdr3aa": "cdr3_b_aa",
                "v.segm": "v_b_gene",
                "j.segm": "j_b_gene",
            }
        )[["cdr3_b_aa", "v_b_gene", "j_b_gene"]],
        organism="human",
        chains=["beta"],
        db_file=args.db_file,
    )
    matrix = pd.DataFrame(tr.pw_beta)

    edges: list[tuple[int, int]] = []
    for i in range(len(matrix)):
        for j in range(i + 1, len(matrix)):
            value = matrix.iat[i, j]
            if value <= args.radius:
                edges.append((i, j))

    cluster_map: dict[int, str] = {}
    cid_counter = 1
    for comp in connected_components(len(input_df), edges):
        if len(comp) < args.min_cluster_size:
            continue
        cid = f"H.B.{epitope}.{cid_counter}"
        cid_counter += 1
        for idx in comp:
            cluster_map[idx] = cid

    save_cluster_members(input_df, cluster_map, args.output, short_name=short_name, method="tcrdist3")
    print(f"Saved {args.output}")


def convert_giana(args: argparse.Namespace) -> None:
    input_df = pd.read_csv(args.input, sep="\t").copy()
    raw = pd.read_csv(args.giana_output, sep="\t", header=None)
    if raw.shape[1] < 4:
        raise ValueError(f"Unexpected GIANA output shape: {raw.shape}")

    raw = raw.rename(columns={0: "cdr3aa", 1: "cluster", 2: "v.segm", 3: "j.segm"})
    raw["cluster"] = raw["cluster"].astype(str)
    valid = raw[raw["cluster"].str.strip().ne("") & ~raw["cluster"].isin(["-1", "nan", "NA"])]
    counts = valid["cluster"].value_counts()
    keep = set(counts[counts >= args.min_cluster_size].index)

    cluster_map: dict[int, str] = {}
    epitope = input_df["antigen.epitope"].iat[0]
    for idx, row in valid.iterrows():
        cluster_id = row["cluster"]
        if cluster_id not in keep:
            continue
        cluster_map[idx] = f"H.B.{epitope}.{cluster_id}"

    save_cluster_members(input_df, cluster_map, args.output, short_name=args.short_name, method="GIANA")
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
    epitope = input_df["antigen.epitope"].iat[0]
    cid_counter = 1
    for comp in comps:
        original = [expanded_nodes[idx] for idx in comp]
        if len(original) < args.min_cluster_size:
            continue
        cid = f"H.B.{epitope}.{cid_counter}"
        cid_counter += 1
        for idx in original:
            cluster_map[idx] = cid

    save_cluster_members(input_df, cluster_map, args.output, short_name=args.short_name, method="TCRnet")
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
            members_ep = members[members["antigen.epitope"].eq(epitope)].copy()

            cdr3_signatures = set(members_ep["cdr3aa"].astype(str))
            vj_signatures = set(
                members_ep["cdr3aa"].astype(str)
                + "|"
                + members_ep["v.segm"].astype(str)
                + "|"
                + members_ep["j.segm"].astype(str)
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
    tcrdist_parser.add_argument("--short-name", required=True)
    tcrdist_parser.add_argument("--radius", type=int, default=24)
    tcrdist_parser.add_argument("--min-cluster-size", type=int, default=3)
    tcrdist_parser.add_argument("--db-file", default="combo_xcr_2024-03-05.tsv")

    giana_parser = subparsers.add_parser("convert-giana")
    giana_parser.add_argument("--input", type=Path, required=True)
    giana_parser.add_argument("--giana-output", type=Path, required=True)
    giana_parser.add_argument("--output", type=Path, required=True)
    giana_parser.add_argument("--short-name", required=True)
    giana_parser.add_argument("--min-cluster-size", type=int, default=3)

    tcrnet_parser = subparsers.add_parser("convert-tcrnet")
    tcrnet_parser.add_argument("--input", type=Path, required=True)
    tcrnet_parser.add_argument("--annotated", type=Path, required=True)
    tcrnet_parser.add_argument("--output", type=Path, required=True)
    tcrnet_parser.add_argument("--short-name", required=True)
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
