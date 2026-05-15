#!/usr/bin/env python3
"""
Download KnotInfo's database dump and convert each knot's grid notation
into yui-hfk's OX-form JSON files.

Usage:
    python3 scripts/fetch-knotinfo-data.py [--out DIR] [--in FILE] [--force] [--keep-xls]

    --out DIR      output directory (default: $YUI_DATA_DIR, else the platform
                   user-data dir for yui).
    --in FILE      use this .xls file instead of downloading.
    --force        overwrite existing JSON files.
    --keep-xls     keep the downloaded .xls file (printed at the end).

Output:
    <out>/grid_diagrams/<name>.json with shape [[o-coords], [x-coords]]
    (0-indexed; the JSON shape yui-hfk's `GridDiagram::load` expects).

Requires: Python 3.8+ and xlrd<2.0  (`pip install 'xlrd<2.0'`).
"""

import argparse
import ast
import json
import os
import sys
import tempfile
import urllib.request
from pathlib import Path

try:
    import xlrd
except ImportError:
    sys.exit(
        "error: xlrd is not installed. Install with:\n"
        "    pip install 'xlrd<2.0'\n"
        "\n"
        "If pip refuses with `externally-managed-environment` (PEP 668),\n"
        "use a venv. From the workspace root:\n"
        "    python3 -m venv .venv\n"
        "    .venv/bin/pip install 'xlrd<2.0'\n"
        "    .venv/bin/python scripts/fetch-knotinfo-data.py"
    )


URL = "https://knotinfo.org/knotinfo_data_complete.xls"

COL_NAME = "name"
COL_GRID = "grid_notation"

SUBDIR = "grid_diagrams"


# ---- KnotInfo grid_notation -> OX conversion ----------------------------------
# KnotInfo's grid_notation is a flat list of 2n markers [[c,r], [c,r], ...],
# 1-indexed, consecutive pairs share a column. We need OX-form: two row-
# permutations of length n (0-indexed), where list[i] is the row of the O/X
# mark in column i. The pair labelling (O vs X) is not in the source; we
# recover it by constraint propagation (every row has exactly one O and one X).
# Adapted from g0-batch/scripts/knotinfo_to_ox.py.

def _parse_columns(points):
    assert len(points) % 2 == 0, "expected an even number of markers"
    n = len(points) // 2
    cols = [None] * (n + 1)
    for i in range(n):
        (c1, r1), (c2, r2) = points[2 * i], points[2 * i + 1]
        assert c1 == c2 == i + 1, f"pair {i} must sit in column {i + 1}"
        cols[i + 1] = (n + 1 - r1, n + 1 - r2)  # flip row index (top→bottom)
    return n, cols


def _index_rows(n, cols):
    rows = [[] for _ in range(n + 1)]
    for c in range(1, n + 1):
        r1, r2 = cols[c]
        rows[r1].append(c)
        rows[r2].append(c)
    for r in range(1, n + 1):
        if len(rows[r]) != 2:
            raise ValueError(f"row {r} has {len(rows[r])} markers, expected 2")
    return rows


def _assign_os(n, cols, rows):
    o_of_col = [None] * (n + 1)
    o_of_row = [None] * (n + 1)

    def propagate(c, r):
        stack = [(c, r)]
        while stack:
            c, r = stack.pop()
            if o_of_col[c] is not None or o_of_row[r] is not None:
                if o_of_col[c] != r or o_of_row[r] != c:
                    raise ValueError("inconsistent X/O assignment")
                continue
            o_of_col[c] = r
            o_of_row[r] = c
            other_r = sum(cols[c]) - r
            other_c = sum(rows[other_r]) - c
            stack.append((other_c, other_r))

    for r in range(1, n + 1):
        if o_of_row[r] is None:
            propagate(rows[r][0], r)
    return o_of_col


def grid_notation_to_ox(s):
    """Parse a KnotInfo grid_notation string and return (O-list, X-list), 0-indexed."""
    points = ast.literal_eval(s)
    n, cols = _parse_columns(points)
    rows = _index_rows(n, cols)
    o_of_col = _assign_os(n, cols, rows)
    o_list = [o_of_col[c] - 1 for c in range(1, n + 1)]
    x_list = [(sum(cols[c]) - o_of_col[c]) - 1 for c in range(1, n + 1)]
    return o_list, x_list


# ---- output dir resolution (mirrors yui_core::util::data_dir) -----------------

def platform_data_dir() -> Path:
    if sys.platform == "darwin":
        return Path.home() / "Library" / "Application Support" / "yui"
    if sys.platform.startswith("win"):
        appdata = os.environ.get("APPDATA") or str(Path.home() / "AppData" / "Roaming")
        return Path(appdata) / "yui" / "data"
    xdg = os.environ.get("XDG_DATA_HOME") or str(Path.home() / ".local" / "share")
    return Path(xdg) / "yui"


def resolve_out_dir(cli_value):
    if cli_value:
        return Path(cli_value)
    if env := os.environ.get("YUI_DATA_DIR"):
        return Path(env)
    return platform_data_dir()


def find_col(header, key):
    for i, h in enumerate(header):
        if h == key:
            return i
    raise SystemExit(f"error: column '{key}' not found in header row")


# ---- main flow ----------------------------------------------------------------

def download_xls(dest: Path) -> None:
    print(f"==> downloading {URL}")
    with urllib.request.urlopen(URL) as r, open(dest, "wb") as f:
        while chunk := r.read(1 << 16):
            f.write(chunk)


def convert(xls_path: Path, out_dir: Path, force: bool) -> None:
    print(f"==> reading {xls_path}")
    wb = xlrd.open_workbook(xls_path)
    sh = wb.sheet_by_index(0)
    header = sh.row_values(0)

    i_name = find_col(header, COL_NAME)
    i_grid = find_col(header, COL_GRID)

    grid_dir = out_dir / SUBDIR
    grid_dir.mkdir(parents=True, exist_ok=True)

    n_written = skipped = 0

    # row 0 = machine header, row 1 = human header, row 2+ = data
    for r in range(2, sh.nrows):
        row = sh.row_values(r)
        name = (row[i_name] or "").strip()
        raw = (row[i_grid] or "").strip()
        if not name or not raw:
            skipped += 1
            continue
        try:
            o, x = grid_notation_to_ox(raw)
        except (ValueError, AssertionError, SyntaxError) as e:
            print(f"  skip {name}: {e}", file=sys.stderr)
            skipped += 1
            continue

        path = grid_dir / f"{name}.json"
        if force or not path.exists():
            path.write_text(json.dumps([o, x]))
            n_written += 1

    print(f"==> wrote {n_written} grid diagrams, skipped {skipped}")
    print(f"    output dir: {out_dir}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--out", dest="out", default=None,
                    help="destination directory (overrides $YUI_DATA_DIR and the platform default)")
    ap.add_argument("--in", dest="src", default=None, type=Path,
                    help="use this .xls file instead of downloading")
    ap.add_argument("--force", action="store_true",
                    help="overwrite existing JSON files")
    ap.add_argument("--keep-xls", action="store_true",
                    help="keep the downloaded .xls (only relevant without --in)")
    args = ap.parse_args()

    out_dir = resolve_out_dir(args.out)

    if args.src:
        if not args.src.exists():
            sys.exit(f"error: input file not found: {args.src}")
        convert(args.src, out_dir, args.force)
        return

    tmp_dir = Path(tempfile.mkdtemp(prefix="knotinfo-"))
    xls = tmp_dir / "knotinfo_data_complete.xls"
    try:
        download_xls(xls)
        convert(xls, out_dir, args.force)
    finally:
        if args.keep_xls:
            print(f"    kept .xls at: {xls}")
        else:
            try:
                xls.unlink()
                tmp_dir.rmdir()
            except OSError:
                pass


if __name__ == "__main__":
    main()
