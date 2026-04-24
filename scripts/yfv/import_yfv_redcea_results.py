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
        "--default-method",
        choices=("dbscan", "vdbscan"),
        default="vdbscan",
        help=(
            "Method label to assign to files that do not already contain an explicit "
            "method token in their basename. In the current server layout, bare "
            "yfv_* files may correspond to legacy dbscan outputs or newer vdbscan outputs."
        ),
    )
    parser.add_argument(
        "--method",
        choices=("dbscan", "leiden", "vdbscan"),
        help="Only import files detected as belonging to the selected clustering method.",
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


def detect_method_from_name(name: str, default_method: str) -> str:
    if "_leiden_" in name:
        return "leiden"
    if "_dbscan_" in name:
        return "dbscan"
    if "_vdbscan_" in name:
        return "vdbscan"
    return default_method


def add_method_to_name(name: str, method: str) -> str:
    suffixes = (
        "_summary_tcrempnet.tsv",
        "_enriched_clonotypes_tcremp.tsv",
        "_tcremp_clusters.tsv",
    )
    for suffix in suffixes:
        if name.endswith(f"_{method}{suffix}") or f"_{method}_" in name[: -len(suffix)]:
            return name
        if name.endswith(suffix):
            return f"{name[:-len(suffix)]}_{method}{suffix}"
    return f"{Path(name).stem}_{method}{Path(name).suffix}"


def copy_file(src: Path, dest_dir: Path, overwrite: bool, default_method: str) -> str:
    method = detect_method_from_name(src.name, default_method=default_method)
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
    methods_seen: dict[str, int] = {}

    for src in files:
        method = detect_method_from_name(src.name, default_method=args.default_method)
        if args.method and method != args.method:
            continue
        message = copy_file(src, dest, overwrite=args.overwrite, default_method=args.default_method)
        print(message)
        methods_seen[method] = methods_seen.get(method, 0) + 1
        if message.startswith("copy"):
            copied += 1
        else:
            skipped += 1

    if not methods_seen:
        print(
            "No matching YFV RedCEA result files were found after method filtering. "
            f"Requested method: {args.method!r}."
        )
        return 1

    print()
    method_summary = ", ".join(f"{method}={count}" for method, count in sorted(methods_seen.items()))
    print(f"Imported {copied} files into {dest}")
    print(f"Detected methods: {method_summary}")
    print(f"Skipped {skipped} already identical files")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
