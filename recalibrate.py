#!/usr/bin/env python3
"""recalibrate.py - re-derive verdicts from a finished results.csv under the
current (calibrated) thresholds, without re-running the expensive search.

Reuses everything already logged (depth, SNR, vetting metrics, the blend
neighbour note) and re-applies classify() + the calibrated blend rule uniformly,
so the whole file ends up on one consistent calibration. Also repairs the stale
27-column header from early rows (missing `method`) by writing a clean 28-column
`results_calibrated.csv`; the raw `results.csv` is left untouched.

Prints the scoreboard under both blend policies:
  demote   - a blend-flagged candidate becomes BLEND_SUSPECT (removed from hits)
  annotate - it stays a candidate; the blend note is kept in the blend_risk column

Usage:
    python3 recalibrate.py                    # default: demote, writes results_calibrated.csv
    python3 recalibrate.py --blend-policy annotate
"""
import argparse
import csv
import math
import re

# calibrated TESS thresholds (mirror exohunter.py)
SNR_MIN, SNR_STRONG = 8.0, 15.0
DEPTH_MIN, DEPTH_MAX = 50.0, 50000.0
ODD_EVEN_MAX, SEC_FRAC, SEC_MIN, DUR_MAX = 0.5, 0.3, 200.0, 0.15
CONTAM_ARCSEC, MARGIN = 42.0, 0.0

FIELDS = ["timestamp", "target", "toi", "tfopwg_disp", "verdict", "agreement",
          "period", "nasa_period_d", "depth_ppm", "snr", "est_prad_earth",
          "nasa_prad_earth", "duration_hr", "duration_frac", "depth_odd_ppm",
          "depth_even_ppm", "secondary_ppm", "blend_risk", "reasons",
          "star_rad_sun", "star_teff_k", "ra", "dec", "n_points", "power",
          "method", "proc_sec", "exofop_url"]
DET = ("STRONG_CANDIDATE", "CANDIDATE")


def num(x):
    try: return float(x)
    except (TypeError, ValueError): return float("nan")


def classify(depth, snr, durf, do, de, sec):
    if not math.isfinite(depth) or snr < SNR_MIN or depth < DEPTH_MIN:
        return "NO_SIGNAL"
    bad = (durf > DUR_MAX) or (depth > DEPTH_MAX)
    if math.isfinite(do) and math.isfinite(de) and max(do, de) > 0 and abs(do-de)/max(do, de) > ODD_EVEN_MAX:
        bad = True
    if math.isfinite(sec) and sec > SEC_MIN and sec > SEC_FRAC*depth:
        bad = True
    if bad: return "FALSE_POSITIVE"
    return "STRONG_CANDIDATE" if snr >= SNR_STRONG else "CANDIDATE"


def blend_flag(blend_risk):
    m = re.search(r'Tmag([\d.]+)@(\d+)"_vs_target([\d.]+)', blend_risk or "")
    if not m: return False
    nb, sep, tg = float(m.group(1)), float(m.group(2)), float(m.group(3))
    return (nb - tg) < MARGIN and sep <= CONTAM_ARCSEC


def agreement(verdict, disp):
    disp = (disp or "").upper()
    det = verdict in DET
    rej = verdict in ("FALSE_POSITIVE", "BLEND_SUSPECT")
    null = verdict == "NO_SIGNAL"
    if disp in ("CP", "KP"):
        return "RE-DETECTED" if det else "MISSED" if null else "DISAGREE" if rej else ""
    if disp in ("FP", "FA"):
        return "DISAGREE" if det else "FP-CAUGHT" if (rej or null) else ""
    if disp in ("PC", "APC"):
        return "OPEN"
    return ""


def read_rows(path):
    """Yield dict rows, inserting method='bls' for stale 27-col rows."""
    for raw in list(csv.reader(open(path)))[1:]:
        if len(raw) == 28:
            yield dict(zip(FIELDS, raw))
        elif len(raw) == 27:                       # old rows: 'method' missing at idx 25
            yield dict(zip(FIELDS, raw[:25] + ["bls"] + raw[25:]))
        # else: truncated/malformed -> skip


def verdict_for(d, policy):
    base = classify(num(d["depth_ppm"]), num(d["snr"]), num(d["duration_frac"]),
                    num(d["depth_odd_ppm"]), num(d["depth_even_ppm"]), num(d["secondary_ppm"]))
    if policy == "demote" and base in DET and blend_flag(d["blend_risk"]):
        return "BLEND_SUSPECT"
    return base


def scoreboard(rows, policy):
    from collections import Counter
    mix = Counter(); conf_n = conf_det = fp_n = fp_rej = blend = 0
    for d in rows:
        v = verdict_for(d, policy)
        mix[v] += 1
        if v == "BLEND_SUSPECT": blend += 1
        disp = (d["tfopwg_disp"] or "").upper()
        if disp in ("CP", "KP"):
            conf_n += 1; conf_det += (v in DET)
        elif disp in ("FP", "FA"):
            fp_n += 1; fp_rej += (v not in DET)
    p = lambda a, b: f"{100*a/b:.0f}%" if b else "n/a"
    n = sum(mix.values())
    print(f"  [{policy}] verdicts: {dict(mix.most_common())}")
    print(f"  [{policy}] confirmed/known detected: {p(conf_det,conf_n)} (n={conf_n})")
    print(f"  [{policy}] false positives rejected: {p(fp_rej,fp_n)} (n={fp_n})")
    if policy == "demote":
        print(f"  [{policy}] blend-flagged: {p(blend,n)}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", default="results.csv")
    ap.add_argument("--out", default="results_calibrated.csv")
    ap.add_argument("--blend-policy", choices=["demote", "annotate"], default="demote")
    a = ap.parse_args()

    rows = list(read_rows(a.inp))
    print(f"recalibrating {len(rows)} rows from {a.inp}\n")
    print("SCOREBOARD (both policies, for comparison):")
    scoreboard(rows, "demote")
    scoreboard(rows, "annotate")

    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        w.writeheader()
        for d in rows:
            v = verdict_for(d, a.blend_policy)
            d["verdict"] = v
            d["agreement"] = agreement(v, d["tfopwg_disp"])
            w.writerow({k: d.get(k, "") for k in FIELDS})
    print(f"\nwrote {a.out}  (blend policy: {a.blend_policy}, clean 28-col header)")


if __name__ == "__main__":
    main()
