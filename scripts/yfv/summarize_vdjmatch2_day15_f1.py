from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Summarize YFV day15/F1 vdjmatch2 outputs into one overlap matrix."
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=Path("/projects/immunestatus/pogorelyy/airr_format"),
    )
    parser.add_argument(
        "--matches-dir",
        type=Path,
        default=Path("/projects/immunestatus/pogorelyy/vdjmatch2_yfv/day15_f1_matches"),
    )
    parser.add_argument(
        "--output-prefix",
        type=Path,
        default=Path("results") / "yfv_vdjmatch2" / "pairwise_overlap" / "yfv_day15_f1_vdjmatch2",
    )
    parser.add_argument(
        "--sample-pattern",
        default="*_15_F1.txt",
    )
    parser.add_argument(
        "--chain",
        default="beta",
    )
    return parser.parse_args()


def sample_id_from_path(path: Path) -> str:
    return path.name.removesuffix(".txt").removesuffix(".txt.gz")


def donor_id_from_sample_id(sample_id: str) -> str:
    return sample_id.split("_", 1)[0]


def load_input_metadata(input_dir: Path, sample_pattern: str, chain: str) -> pd.DataFrame:
    records: list[dict[str, object]] = []
    for sample_path in sorted(input_dir.glob(sample_pattern)):
        sample_id = sample_id_from_path(sample_path)
        df = pd.read_csv(sample_path, sep="\t", usecols=["locus"])
        n_rows = int(df["locus"].fillna("").astype(str).str.strip().str.lower().eq(chain.lower()).sum())
        records.append(
            {
                "sample_id": sample_id,
                "donor": donor_id_from_sample_id(sample_id),
                "sample_path": str(sample_path),
                "n_input_rows": n_rows,
            }
        )
    if not records:
        raise FileNotFoundError(f"No sample files matched '*_15_F1.txt' in {input_dir}")
    return pd.DataFrame.from_records(records).sort_values("sample_id").reset_index(drop=True)


def pick_row_column(columns: list[str], side: str) -> str:
    expected = f"{side}_row"
    for column in columns:
        if column.lower() == expected:
            return column
    raise ValueError(f"Could not find {expected} column in vdjmatch2 output: {columns}")


def count_unique_query_rows(match_path: Path) -> tuple[int, int]:
    try:
        df = pd.read_csv(match_path, sep="\t")
    except pd.errors.EmptyDataError:
        return 0, 0
    if df.empty:
        return 0, 0

    query_col = pick_row_column(list(df.columns), "query")
    target_col = pick_row_column(list(df.columns), "target")

    pairs = (
        df[[query_col, target_col]]
        .fillna("")
        .astype(str)
        .rename(columns={query_col: "query_row", target_col: "target_row"})
    )
    pairs = pairs.loc[
        pairs["query_row"].str.strip().ne("") & pairs["target_row"].str.strip().ne("")
    ].drop_duplicates()
    return int(pairs["query_row"].nunique()), int(len(pairs))


def build_pairwise_table(metadata: pd.DataFrame, matches_dir: Path) -> pd.DataFrame:
    sample_ids = metadata["sample_id"].tolist()
    n_rows_by_sample = dict(zip(metadata["sample_id"], metadata["n_input_rows"]))
    records: list[dict[str, object]] = []

    for sample_x in sample_ids:
        n_x = int(n_rows_by_sample[sample_x])
        for sample_y in sample_ids:
            n_y = int(n_rows_by_sample[sample_y])
            if sample_x == sample_y:
                n_matched = n_x
                n_pairs = n_x
            else:
                match_path = matches_dir / f"{sample_x}__vs__{sample_y}.tsv"
                n_matched, n_pairs = count_unique_query_rows(match_path)

            records.append(
                {
                    "sample_x": sample_x,
                    "sample_y": sample_y,
                    "donor_x": donor_id_from_sample_id(sample_x),
                    "donor_y": donor_id_from_sample_id(sample_y),
                    "n_input_rows_x": n_x,
                    "n_input_rows_y": n_y,
                    "n_matched_query_rows_x_in_y": n_matched,
                    "n_unique_query_target_pairs": n_pairs,
                    "share_x_in_y": (n_matched / n_x) if n_x else 0.0,
                }
            )

    return pd.DataFrame.from_records(records)


def main() -> None:
    args = parse_args()
    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)

    metadata = load_input_metadata(args.input_dir, args.sample_pattern, args.chain)
    pairwise = build_pairwise_table(metadata, args.matches_dir)
    matrix = pairwise.pivot(index="sample_x", columns="sample_y", values="share_x_in_y")
    matrix = matrix.sort_index(axis=0).sort_index(axis=1)

    metadata_path = args.output_prefix.with_name(args.output_prefix.name + "_sample_metadata.tsv")
    pairwise_path = args.output_prefix.with_name(args.output_prefix.name + "_pairwise.tsv")
    matrix_path = args.output_prefix.with_name(args.output_prefix.name + "_matrix.tsv")

    metadata.to_csv(metadata_path, sep="\t", index=False)
    pairwise.to_csv(pairwise_path, sep="\t", index=False)
    matrix.to_csv(matrix_path, sep="\t")

    print(f"Saved sample metadata to {metadata_path}")
    print(f"Saved pairwise overlap table to {pairwise_path}")
    print(f"Saved overlap matrix to {matrix_path}")


if __name__ == "__main__":
    main()
