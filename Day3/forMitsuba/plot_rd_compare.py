#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def read_csv(path: Path):
    with path.open("r", encoding="utf-8") as f:
        r = csv.DictReader(f)
        rows = list(r)
    if not rows:
        raise ValueError(f"No data in {path}")
    return rows


def extract(rows, x_key, y_key):
    xs = [float(row[x_key]) for row in rows]
    ys = [float(row[y_key]) for row in rows]
    return xs, ys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rd_no_beard", default="Rd_no_beard.csv")
    ap.add_argument("--rd_beard", default="rd_beard.csv")
    ap.add_argument("--out", default="rd_norm_compare.png")
    ap.add_argument("--x_key", default="r_center_m")
    ap.add_argument("--y_key", default="RdY_norm_int",
                    help="Column name to plot. Ex: RdY_norm_int or two_pi_r_RdY_norm_int or Y")
    ap.add_argument("--title", default="Normalized Rd(r) comparison")
    ap.add_argument("--xlabel", default="r [m]")
    ap.add_argument("--ylabel", default="Rd(r) (normalized)")
    ap.add_argument("--dpi", type=int, default=200)
    args = ap.parse_args()

    p_no = Path(args.rd_no_beard)
    p_be = Path(args.rd_beard)
    if not p_no.exists():
        raise FileNotFoundError(p_no)
    if not p_be.exists():
        raise FileNotFoundError(p_be)

    rows_no = read_csv(p_no)
    rows_be = read_csv(p_be)

    x_no, y_no = extract(rows_no, args.x_key, args.y_key)
    x_be, y_be = extract(rows_be, args.x_key, args.y_key)

    plt.figure(figsize=(6.4, 4.0))
    plt.plot(x_no, y_no, label="no beard")
    plt.plot(x_be, y_be, label="beard")
    plt.title(args.title)
    plt.xlabel(args.xlabel)
    plt.ylabel(args.ylabel)
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()

    out_path = Path(args.out)
    plt.savefig(out_path, dpi=args.dpi)
    print(f"wrote: {out_path.resolve()}")


if __name__ == "__main__":
    main()
