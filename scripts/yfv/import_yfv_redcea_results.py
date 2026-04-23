from __future__ import annotations

import argparse
import filecmp
import shutil
from pathlib import Path


FILE_PATTERNS = (
    "*_summary_tcrempnet.tsv",
    "*_enriched_clonotypes_tcremp.tsv",
    "*_tcremp_clusters.tsv",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Copy YFV RedCEA result files from an external directory into a dedicated "
            "folder inside this project."
        )
    )
    parser.add_argument(
        "--source",
        type=Path,
        required=True,
        help="Directory that contains RedCEA outputs, optionally with nested run subdirectories.",
    )
    parser.add_argument(
        "--dest",
        type=Path,
        default=Path("results") / "yfv_redcea" / "raw",
        help="Destination directory inside the project.",
    )
    parser.add_argument(
        "--prefix",
        default="yfv_",
        help="Only copy files whose basename starts with this prefix.",
    )
    parser.add_argument(
        "--method",
        choices=("leiden", "vdbscan"),
        required=True,
        help=(
            "Clustering method label for this result set. It is appended to copied "
            "filenames so different methods can coexist inside the project."
        ),
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite destination files if they already exist and differ.",
    )
    return parser.parse_args()


def iter_candidate_files(source: Path, prefix: str) -> list[Path]:
    matches: list[Path] = []
    for pattern in FILE_PATTERNS:
        for path in source.rglob(pattern):
            if path.is_file() and path.name.startswith(prefix):
                matches.append(path)
    return sorted(set(matches))


def add_method_to_name(name: str, method: str) -> str:
    suffixes = (
        "_summary_tcrempnet.tsv",
        "_enriched_clonotypes_tcremp.tsv",
        "_tcremp_clusters.tsv",
    )
    for suffix in suffixes:
        if name.endswith(f"_{method}{suffix}"):
            return name
        if name.endswith(suffix):
            return f"{name[:-len(suffix)]}_{method}{suffix}"
    return f"{Path(name).stem}_{method}{Path(name).suffix}"


def copy_file(src: Path, dest_dir: Path, overwrite: bool, method: str) -> str:
    dest = dest_dir / add_method_to_name(src.name, method=method)
    if dest.exists():
        same = filecmp.cmp(src, dest, shallow=False)
        if same:
            return f"skip  {src} -> {dest} (already identical)"
        if not overwrite:
            raise FileExistsError(
                f"Destination file already exists and differs: {dest}. "
                "Re-run with --overwrite to replace it."
            )
    shutil.copy2(src, dest)
    return f"copy  {src} -> {dest}"


def main() -> int:
    args = parse_args()
    source = args.source.resolve()
    dest = args.dest.resolve()

    if not source.exists():
        raise FileNotFoundError(f"Source directory does not exist: {source}")
    if not source.is_dir():
        raise NotADirectoryError(f"Source is not a directory: {source}")

    files = iter_candidate_files(source, args.prefix)
    if not files:
        print(
            "No matching YFV RedCEA result files were found. "
            f"Searched under {source} for prefix {args.prefix!r}."
        )
        return 1

    dest.mkdir(parents=True, exist_ok=True)
    copied = 0
    skipped = 0

    for src in files:
        message = copy_file(src, dest, overwrite=args.overwrite, method=args.method)
        print(message)
        if message.startswith("copy"):
            copied += 1
        else:
            skipped += 1

    print()
    print(f"Imported {copied} {args.method} files into {dest}")
    print(f"Skipped {skipped} already identical files")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
