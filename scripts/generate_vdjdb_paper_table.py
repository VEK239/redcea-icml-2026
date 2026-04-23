#!/usr/bin/env python
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


EPITOPE_ORDER = ["GLC", "YLQ"]
METHOD_ORDER = ["GIANA", "TCRdist3", "TCRnet", "RedCEA", "RedCEA+VJ"]


def find_repo_root() -> Path:
    current = Path.cwd().resolve()
    for candidate in [current] + list(current.parents):
        if (candidate / "scripts").exists() and (candidate / "results").exists():
            return candidate
    raise FileNotFoundError("Could not locate repository root from current working directory.")


def build_table(benchmark_detailed: Path, redcea_tuned_detailed: Path) -> pd.DataFrame:
    bench = pd.read_csv(benchmark_detailed)
    red = pd.read_csv(redcea_tuned_detailed)

    baseline = (
        bench[
            bench["method"].isin(["giana", "tcrdist3", "tcrnet"])
            & bench["match_mode"].eq("cdr3")
        ][["method", "epitope_short", "precision", "recall", "f1"]]
        .copy()
    )
    baseline["Method"] = baseline["method"].map(
        {
            "giana": "GIANA",
            "tcrdist3": "TCRdist3",
            "tcrnet": "TCRnet",
        }
    )

    tuned_rows = []
    for epitope in EPITOPE_ORDER:
        for with_vj, label in [(False, "RedCEA"), (True, "RedCEA+VJ")]:
            sub = red[(red["epitope_short"].eq(epitope)) & (red["with_vj"].eq(with_vj))].copy()
            sub = sub.sort_values(["f1", "precision", "recall"], ascending=[False, False, False]).head(1)
            if sub.empty:
                continue
            row = sub.iloc[0]
            tuned_rows.append(
                {
                    "Method": label,
                    "epitope_short": epitope,
                    "precision": row["precision"],
                    "recall": row["recall"],
                    "f1": row["f1"],
                }
            )
    tuned = pd.DataFrame(tuned_rows)

    out = pd.concat(
        [
            baseline[["Method", "epitope_short", "precision", "recall", "f1"]],
            tuned[["Method", "epitope_short", "precision", "recall", "f1"]],
        ],
        ignore_index=True,
    )
    out["epitope_order"] = out["epitope_short"].map({ep: i for i, ep in enumerate(EPITOPE_ORDER)})
    out["method_order"] = out["Method"].map({method: i for i, method in enumerate(METHOD_ORDER)})
    out = out.sort_values(["epitope_order", "method_order"]).drop(columns=["epitope_order", "method_order"]).reset_index(drop=True)
    out = out.rename(
        columns={
            "epitope_short": "Epitope",
            "precision": "Precision",
            "recall": "Recall",
            "f1": "F1",
        }
    )
    out = out[["Epitope", "Method", "Precision", "Recall", "F1"]]
    return out


def markdown_table(df: pd.DataFrame) -> str:
    rounded = df.copy()
    for col in ["Precision", "Recall", "F1"]:
        rounded[col] = rounded[col].map(lambda value: f"{float(value):.3f}")

    lines = [
        "| Epitope | Method | Precision | Recall | F1 |",
        "|---|---|---:|---:|---:|",
    ]
    for _, row in rounded.iterrows():
        lines.append(f"| {row['Epitope']} | {row['Method']} | {row['Precision']} | {row['Recall']} | {row['F1']} |")
    return "\n".join(lines) + "\n"


def parse_args() -> argparse.Namespace:
    root = find_repo_root()
    result_root = root / "results" / "vdjdb_method_benchmark"
    parser = argparse.ArgumentParser(description="Generate compact paper table for the TRB VDJdb comparison figure.")
    parser.add_argument(
        "--benchmark-detailed",
        type=Path,
        default=result_root / "metrics_detailed.csv",
        help="Detailed baseline metrics CSV.",
    )
    parser.add_argument(
        "--redcea-tuned-detailed",
        type=Path,
        default=result_root / "redcea_trb_param_grid_metrics_detailed.csv",
        help="Detailed tuned RedCEA TRB parameter-grid metrics CSV.",
    )
    parser.add_argument(
        "--output-csv",
        type=Path,
        default=result_root / "paper_table_vdjdb_trb_comparison.csv",
        help="Output CSV path.",
    )
    parser.add_argument(
        "--output-md",
        type=Path,
        default=result_root / "paper_table_vdjdb_trb_comparison.md",
        help="Output Markdown path.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    table = build_table(args.benchmark_detailed, args.redcea_tuned_detailed)

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    args.output_md.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.output_csv, index=False)
    args.output_md.write_text(markdown_table(table), encoding="utf-8")

    print(f"Saved {args.output_csv}")
    print(f"Saved {args.output_md}")
    print(table.to_string(index=False))


if __name__ == "__main__":
    main()
