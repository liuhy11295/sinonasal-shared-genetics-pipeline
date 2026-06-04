#!/usr/bin/env python3
"""Fail fast if a GitHub publication mirror includes non-code project material."""

from __future__ import annotations

import argparse
from pathlib import Path


FORBIDDEN_DIRS = {
    "data",
    "results",
}
SKIPPED_DIRS = {".git", "__pycache__", ".pytest_cache", ".codex_deps"}
FORBIDDEN_NAMES = {
    "authorized_keys",
    "wget-log",
}
FORBIDDEN_NAME_FRAGMENTS = {
    "codex_ssh",
    "tailscale",
    "private_key",
}
FORBIDDEN_SUFFIXES = {
    ".bam",
    ".bai",
    ".bed",
    ".bgen",
    ".cram",
    ".csi",
    ".err",
    ".gz",
    ".h5ad",
    ".log",
    ".out",
    ".parquet",
    ".part",
    ".rds",
    ".tar",
    ".tgz",
    ".tbi",
    ".vcf",
    ".zip",
}
SECRET_MARKERS = {
    "BEGIN " + "OPENSSH PRIVATE KEY",
    "BEGIN " + "RSA PRIVATE KEY",
    "github_" + "pat_",
    "gh" + "p_",
    "Authorization: " + "Bearer ",
}
DEFAULT_MAX_BYTES = 5 * 1024 * 1024


def audit_export(root: Path, max_bytes: int = DEFAULT_MAX_BYTES) -> list[str]:
    failures: list[str] = []
    root = root.resolve()
    if not root.is_dir():
        return [f"missing export directory: {root}"]

    for path in sorted(root.rglob("*")):
        if path.is_symlink():
            failures.append(f"symlink is not publishable: {path.relative_to(root)}")
            continue
        if not path.is_file():
            continue

        relative = path.relative_to(root)
        lowercase_parts = {part.lower() for part in relative.parts}
        lowercase_name = path.name.lower()
        suffixes = "".join(path.suffixes).lower()

        if lowercase_parts & SKIPPED_DIRS:
            continue
        if lowercase_parts & FORBIDDEN_DIRS:
            failures.append(f"forbidden directory: {relative}")
        if lowercase_name in FORBIDDEN_NAMES:
            failures.append(f"forbidden file: {relative}")
        if any(fragment in lowercase_name for fragment in FORBIDDEN_NAME_FRAGMENTS):
            failures.append(f"sensitive filename: {relative}")
        if path.suffix.lower() in FORBIDDEN_SUFFIXES or suffixes.endswith(".tar.gz"):
            failures.append(f"forbidden data/runtime suffix: {relative}")
        if path.stat().st_size > max_bytes:
            failures.append(f"unexpectedly large file: {relative} ({path.stat().st_size} bytes)")

        if path.stat().st_size <= max_bytes:
            try:
                content = path.read_text(encoding="utf-8", errors="ignore")
            except OSError as error:
                failures.append(f"unreadable file: {relative} ({error})")
                continue
            for marker in SECRET_MARKERS:
                if marker in content:
                    failures.append(f"possible credential material in: {relative}")
                    break

    return failures


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Audit the code/documentation GitHub mirror before pushing."
    )
    parser.add_argument("root", type=Path, help="Path to the code-only publication mirror.")
    parser.add_argument("--max-bytes", type=int, default=DEFAULT_MAX_BYTES)
    args = parser.parse_args()

    failures = audit_export(args.root, max_bytes=args.max_bytes)
    if failures:
        for failure in failures:
            print(f"BLOCKED\t{failure}")
        print(f"audit_github_code_export status=FAIL failures={len(failures)}")
        return 1

    files = sum(1 for path in args.root.rglob("*") if path.is_file())
    print(f"audit_github_code_export status=PASS files={files}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
