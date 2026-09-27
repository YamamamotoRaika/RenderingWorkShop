import csv
import math
import argparse
from typing import List, Dict, Tuple

def trapz(x: List[float], y: List[float]) -> float:
    """Simple trapezoidal integration (pure python)."""
    s = 0.0
    for i in range(len(x) - 1):
        dx = x[i + 1] - x[i]
        s += 0.5 * (y[i] + y[i + 1]) * dx
    return s

def detect_r_column(fieldnames: List[str]) -> str:
    # Try common names first
    candidates = ["r", "r_m", "r_mm", "r_center_m", "r_center", "radius", "radius_m", "radius_mm"]
    for c in candidates:
        if c in fieldnames:
            return c
    # Fallback: pick first column that looks like r*
    for name in fieldnames:
        low = name.lower()
        if low.startswith("r"):
            return name
    raise RuntimeError(f"Could not detect r column. Columns: {fieldnames}")

def detect_rd_columns(fieldnames: List[str]) -> List[str]:
    # Prefer columns like Rd_skin, Rd_beard, Rd_mix
    rd_cols = [c for c in fieldnames if c.lower().startswith("rd")]
    # Filter out already-normalized columns if present
    rd_cols = [c for c in rd_cols if ("norm" not in c.lower() and "two_pi" not in c.lower())]
    # If none found, fail loudly
    if not rd_cols:
        raise RuntimeError(f"No Rd columns found. Columns: {fieldnames}")
    return rd_cols

def read_csv(path: str) -> Tuple[List[str], List[Dict[str, str]]]:
    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        if reader.fieldnames is None:
            raise RuntimeError("CSV has no header row.")
        rows = list(reader)
        return reader.fieldnames, rows

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in_csv", required=True)
    ap.add_argument("--out_csv", required=True)
    ap.add_argument("--r_col", default=None, help="Override r column name if auto-detection fails")
    ap.add_argument("--rd_cols", default=None,
                    help="Comma-separated Rd columns (e.g., Rd_skin,Rd_beard,Rd_mix). Auto-detect if omitted.")
    ap.add_argument("--assume_r_unit", default="as-is",
                    choices=["as-is", "mm_to_m"],
                    help="If your r is in mm, choose mm_to_m to convert r to meters before integration.")
    args = ap.parse_args()

    fieldnames, rows = read_csv(args.in_csv)
    r_col = args.r_col or detect_r_column(fieldnames)
    rd_cols = args.rd_cols.split(",") if args.rd_cols else detect_rd_columns(fieldnames)

    # Parse r (optionally convert mm->m)
    r = []
    for row in rows:
        rv = float(row[r_col])
        if args.assume_r_unit == "mm_to_m":
            rv *= 1e-3
        r.append(rv)

    # Ensure r is sorted (if not, sort rows together)
    # This also protects trapz.
    order = sorted(range(len(r)), key=lambda i: r[i])
    r = [r[i] for i in order]
    rows = [rows[i] for i in order]

    # Compute integrals and normalized columns
    integrals = {}
    rd_norm = {c: [] for c in rd_cols}
    two_pi_r_rd_norm = {c: [] for c in rd_cols}

    for c in rd_cols:
        yd = [float(row[c]) for row in rows]
        weighted = [yd[i] * 2.0 * math.pi * r[i] for i in range(len(r))]
        I = trapz(r, weighted)
        integrals[c] = I

        if not math.isfinite(I) or I <= 0.0:
            print(f"[WARN] Integral for {c} is non-positive/non-finite: {I}. Writing zeros for normalized columns.")
            for i in range(len(r)):
                rd_norm[c].append(0.0)
                two_pi_r_rd_norm[c].append(0.0)
            continue

        for i in range(len(r)):
            v = float(rows[i][c]) / I
            rd_norm[c].append(v)
            two_pi_r_rd_norm[c].append(2.0 * math.pi * r[i] * v)

    # Write output CSV
    out_fields = [r_col] + rd_cols
    for c in rd_cols:
        out_fields.append(f"{c}_norm_int")
    for c in rd_cols:
        out_fields.append(f"two_pi_r_{c}_norm_int")

    with open(args.out_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=out_fields)
        w.writeheader()
        for i in range(len(r)):
            out = {}
            out[r_col] = r[i]
            for c in rd_cols:
                out[c] = float(rows[i][c])
                out[f"{c}_norm_int"] = rd_norm[c][i]
                out[f"two_pi_r_{c}_norm_int"] = two_pi_r_rd_norm[c][i]
            w.writerow(out)

    print("=== integrals (should become ~1 after normalization check) ===")
    for c in rd_cols:
        print(f"{c}: I = {integrals[c]}")

    # Check normalization: ∫ 2π r Rd_norm dr should be ~1
    print("=== check: integral of normalized (target ~1) ===")
    for c in rd_cols:
        y = rd_norm[c]
        weighted = [y[i] * 2.0 * math.pi * r[i] for i in range(len(r))]
        chk = trapz(r, weighted)
        print(f"{c}: check = {chk}")

    print("wrote:", args.out_csv)

if __name__ == "__main__":
    main()
