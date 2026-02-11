#!/usr/bin/env python3
"""
Compare two directories of Python scripts and show changes.

- Matches .py files by relative path.
- Reports added / removed / modified files.
- Prints unified diffs for modified files.
"""

from __future__ import annotations
import argparse
import os
from pathlib import Path
import difflib
import hashlib
from typing import Dict, Set, Tuple


def iter_py_files(root: Path) -> Dict[str, Path]:
    """Return mapping: relative_posix_path -> absolute_path for all .py under root."""
    mapping: Dict[str, Path] = {}
    for p in root.rglob("*.py"):
        if p.is_file():
            rel = p.relative_to(root).as_posix()
            mapping[rel] = p
    return mapping


def file_hash(path: Path) -> str:
    """Fast content hash for change detection."""
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def read_text_lines(path: Path) -> list[str]:
    """Read file as text lines (best-effort)."""
    try:
        text = path.read_text(encoding="utf-8")
    except UnicodeDecodeError:
        # fallback if encoding differs
        text = path.read_text(errors="replace")
    return text.splitlines(keepends=True)


def unified_diff(old_path: Path, new_path: Path, rel: str, context: int) -> str:
    old_lines = read_text_lines(old_path)
    new_lines = read_text_lines(new_path)

    return "".join(
        difflib.unified_diff(
            old_lines,
            new_lines,
            fromfile=f"a/{rel}",
            tofile=f"b/{rel}",
            n=context,
            lineterm="",
        )
    )


def main():
    ap = argparse.ArgumentParser(description="Compare two folders of Python scripts and show diffs.")
    ap.add_argument("old_dir", type=Path, help="Old/base directory")
    ap.add_argument("new_dir", type=Path, help="New/changed directory")
    ap.add_argument("--context", "-c", type=int, default=3, help="Diff context lines (default: 3)")
    ap.add_argument("--only-summary", action="store_true", help="Only print summary, no diffs")
    ap.add_argument("--write", type=Path, default=None, help="Write output to a file")
    ap.add_argument("--ignore-whitespace", action="store_true",
                    help="Ignore whitespace-only changes (approximate)")
    args = ap.parse_args()

    old_dir: Path = args.old_dir.resolve()
    new_dir: Path = args.new_dir.resolve()

    if not old_dir.is_dir() or not new_dir.is_dir():
        raise SystemExit("Both arguments must be directories.")

    old_files = iter_py_files(old_dir)
    new_files = iter_py_files(new_dir)

    old_set: Set[str] = set(old_files.keys())
    new_set: Set[str] = set(new_files.keys())

    added = sorted(new_set - old_set)
    removed = sorted(old_set - new_set)
    common = sorted(old_set & new_set)

    modified: list[str] = []
    unchanged: list[str] = []

    # Detect modified via hash (fast)
    for rel in common:
        if file_hash(old_files[rel]) != file_hash(new_files[rel]):
            modified.append(rel)
        else:
            unchanged.append(rel)

    out_lines: list[str] = []

    out_lines.append(f"OLD: {old_dir}")
    out_lines.append(f"NEW: {new_dir}")
    out_lines.append("")
    out_lines.append("=== Summary ===")
    out_lines.append(f"Added   : {len(added)}")
    out_lines.append(f"Removed : {len(removed)}")
    out_lines.append(f"Modified: {len(modified)}")
    out_lines.append(f"Same    : {len(unchanged)}")
    out_lines.append("")

    if added:
        out_lines.append("=== Added files ===")
        out_lines.extend(f"+ {rel}" for rel in added)
        out_lines.append("")
    if removed:
        out_lines.append("=== Removed files ===")
        out_lines.extend(f"- {rel}" for rel in removed)
        out_lines.append("")
    if modified:
        out_lines.append("=== Modified files ===")
        out_lines.extend(f"* {rel}" for rel in modified)
        out_lines.append("")

    if not args.only_summary and modified:
        out_lines.append("=== Diffs ===")
        for rel in modified:
            old_p = old_files[rel]
            new_p = new_files[rel]

            diff_txt = unified_diff(old_p, new_p, rel, args.context)

            # Optional: ignore whitespace-only changes (approximate)
            if args.ignore_whitespace:
                # If after stripping whitespace the content matches, skip
                old_stripped = [ln.strip() for ln in read_text_lines(old_p)]
                new_stripped = [ln.strip() for ln in read_text_lines(new_p)]
                if old_stripped == new_stripped:
                    out_lines.append(f"(whitespace-only change skipped) {rel}")
                    out_lines.append("")
                    continue

            if diff_txt.strip():
                out_lines.append(diff_txt)
                out_lines.append("")  # blank line between files
            else:
                out_lines.append(f"(no textual diff produced) {rel}")
                out_lines.append("")

    output = "\n".join(out_lines)

    if args.write:
        args.write.parent.mkdir(parents=True, exist_ok=True)
        args.write.write_text(output, encoding="utf-8")
        print(f"Wrote report to: {args.write}")
    else:
        print(output)


if __name__ == "__main__":
    main()
