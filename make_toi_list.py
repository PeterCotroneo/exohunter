#!/usr/bin/env python3
"""make_toi_list.py - fetch the live TESS Objects of Interest (TOI) list from the
NASA Exoplanet Archive and write a target file for scan.py.

Each TOI is a candidate NASA flagged from TESS data. This pulls their host-star
TIC ids so you can re-run the transit search yourself and practise vetting them.
By default it keeps only the UNCONFIRMED candidates (disposition PC / APC), which
are the ones actually worth vetting.

Pure standard library + NASA's public TAP service. No API key, no tokens.

Usage:
    python3 make_toi_list.py                 # -> tois.txt (PC + APC)
    python3 make_toi_list.py --disp PC       # planet candidates only
    python3 make_toi_list.py --all           # every TOI regardless of disposition
    python3 make_toi_list.py -o my_tois.txt

Disposition codes (tfopwg_disp):
    PC  = planetary candidate (unconfirmed)      <- default keep
    APC = ambiguous planetary candidate          <- default keep
    CP  = confirmed planet
    KP  = known planet
    FP  = false positive
    FA  = false alarm
"""
import argparse
import csv
import io
import os
import ssl
import sys
import urllib.parse
import urllib.request
from collections import Counter

try:  # macOS system Python often lacks root certs; use certifi's bundle
    import certifi
    SSL_CTX = ssl.create_default_context(cafile=certifi.where())
except Exception:
    SSL_CTX = ssl.create_default_context()

TAP = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync"
QUERY = ("select tid,toi,tfopwg_disp,pl_orbper,pl_rade,st_rad,st_teff from toi")
HERE = os.path.dirname(os.path.abspath(__file__))
META_CSV = os.path.join(HERE, "toi_meta.csv")


def fetch_toi_csv():
    url = TAP + "?" + urllib.parse.urlencode({"query": QUERY, "format": "csv"})
    print(f"Fetching TOI table from NASA Exoplanet Archive ...", flush=True)
    with urllib.request.urlopen(url, timeout=120, context=SSL_CTX) as resp:
        return resp.read().decode("utf-8")


def main():
    p = argparse.ArgumentParser(description="Build a TIC target list from the TOI catalog.")
    p.add_argument("-o", "--output", default=os.path.join(HERE, "tois.txt"))
    p.add_argument("--disp", default="PC,APC",
                   help="comma-separated dispositions to keep (default PC,APC)")
    p.add_argument("--all", action="store_true", help="keep every TOI regardless of disposition")
    a = p.parse_args()

    try:
        text = fetch_toi_csv()
    except Exception as e:
        sys.exit(f"Failed to fetch TOI table: {type(e).__name__}: {e}")

    keep = None if a.all else {d.strip().upper() for d in a.disp.split(",") if d.strip()}
    counts = Counter()
    tics = []
    seen = set()
    meta = {}
    for row in csv.DictReader(io.StringIO(text)):
        disp = (row.get("tfopwg_disp") or "").strip().upper()
        counts[disp or "(blank)"] += 1
        if keep is not None and disp not in keep:
            continue
        tid = (row.get("tid") or "").strip()
        if not tid or tid in seen:
            continue
        seen.add(tid)
        t = int(float(tid))
        tics.append(t)
        meta[t] = {"toi": row.get("toi", ""), "tfopwg_disp": disp,
                   "nasa_period_d": row.get("pl_orbper", ""),
                   "nasa_prad_earth": row.get("pl_rade", ""),
                   "star_rad_sun": row.get("st_rad", ""),
                   "star_teff_k": row.get("st_teff", "")}

    tics.sort()
    with open(a.output, "w") as fh:
        fh.write("# TESS Objects of Interest host stars (TIC ids)\n")
        fh.write(f"# dispositions kept: {'ALL' if a.all else a.disp}\n")
        for t in tics:
            fh.write(f"TIC {t}\n")

    mcols = ["tic", "toi", "tfopwg_disp", "nasa_period_d", "nasa_prad_earth",
             "star_rad_sun", "star_teff_k"]
    with open(META_CSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=mcols)
        w.writeheader()
        for t in tics:
            w.writerow({"tic": t, **meta[t]})

    print("\nTOI dispositions in the catalog:")
    for d, n in counts.most_common():
        print(f"  {d:<8s} {n}")
    print(f"\nWrote {len(tics)} host stars to {a.output}")
    print(f"Wrote metadata (TOI/disposition/NASA period/star radius) to {META_CSV}")
    print(f"Now run:  python3 scan.py --targets {os.path.basename(a.output)} --pmax 30")


if __name__ == "__main__":
    main()
