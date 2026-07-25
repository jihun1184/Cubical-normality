"""Command-line interface for producing a JSON normality certificate."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

from .checker import build_certificate_cli


def _load_json(path: str | None) -> Any:
    if path is None or path == "-":
        return json.load(sys.stdin)
    with Path(path).open("r", encoding="utf-8") as handle:
        return json.load(handle)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Check a nonempty finite set of 3D integer voxels and emit a "
            "schema-versioned cubical-normality certificate."
        )
    )
    parser.add_argument(
        "input",
        nargs="?",
        default="-",
        help="JSON file containing [[i,j,k], ...]; use '-' or omit for stdin.",
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Write JSON to this file instead of stdout.",
    )
    parser.add_argument(
        "--indent",
        type=int,
        default=2,
        help="JSON indentation (default: 2).",
    )
    args = parser.parse_args(argv)

    try:
        raw = _load_json(args.input)
    except (OSError, json.JSONDecodeError) as exc:
        parser.error(f"could not read JSON input: {exc}")

    # JSON arrays are converted to tuples to satisfy the library's explicit
    # length-3 tuple input contract; malformed entries remain visible to the
    # validator rather than being silently repaired.
    voxels = [tuple(v) if isinstance(v, list) else v for v in raw] if isinstance(raw, list) else raw
    certificate = build_certificate_cli(voxels)
    rendered = json.dumps(certificate, indent=args.indent, sort_keys=True)

    if args.output:
        Path(args.output).write_text(rendered + "\n", encoding="utf-8")
    else:
        print(rendered)

    return 0 if certificate.get("overall") == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
