#!/usr/bin/env python3
"""make_koi_list.py - fetch the Kepler Objects of Interest (KOI) cumulative table
from the NASA Exoplanet Archive and write a LABELLED target file for benchmark.py.

Kepler is the best ground-truth set there is for calibrating a transit search:
thousands of stars each already labelled CONFIRMED / FALSE POSITIVE / CANDIDATE by
the Kepler team. We pull those labels so benchmark.py can measure how often
exohunter agrees, per mission - instead of us guessing whether it "will be accurate".

Pure standard library + NASA's public TAP service. No API key, no tokens.

Usage:
    python3 make_koi_list.py                 # -> koi_targets.txt + koi_meta.csv
    python3 make_koi_list.py --disp CONFIRMED,FALSE POSITIVE   # skip open candidates

Disposition codes (koi_disposition):
    CONFIRMED        - a real planet (exohunter SHOULD detect it)
    FALSE POSITIVE   - not a planet (exohunter SHOULD reject it)
    CANDIDATE        - unresolved (no ground truth -> "OPEN")

One star (kepid) can host several KOIs. We collapse to one row per star, keeping
the strongest label (CONFIRMED > CANDIDATE > FALSE POSITIVE) so a star with any
confirmed planet counts as one we should detect.
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
QUERY = ("select kepid,kepoi_name,kepler_name,koi_disposition,koi_period,"
         "koi_prad,koi_srad,koi_steff from cumulative")
HERE = os.path.dirname(os.path.abspath(__file__))
META_CSV = os.path.join(HERE, "koi_meta.csv")

# strongest label wins when a star hosts several KOIs
RANK = {"CONFIRMED": 3, "CANDIDATE": 2, "FALSE POSITIVE": 1}


def fetch_koi_csv():
    url = TAP + "?" + urllib.parse.urlencode({"query": QUERY, "format": "csv"})
    print("Fetching KOI cumulative table from NASA Exoplanet Archive ...", flush=True)
    with urllib.request.urlopen(url, timeout=180, context=SSL_CTX) as resp:
        return resp.read().decode("utf-8")


def main():
    p = argparse.ArgumentParser(description="Build a labelled KIC target list from the KOI table.")
    p.add_argument("-o", "--output", default=os.path.join(HERE, "koi_targets.txt"))
    p.add_argument("--disp", default="CONFIRMED,CANDIDATE,FALSE POSITIVE",
                   help="comma-separated dispositions to keep")
    a = p.parse_args()

    try:
        text = fetch_koi_csv()
    except Exception as e:
        sys.exit(f"Failed to fetch KOI table: {type(e).__name__}: {e}")

    keep = {d.strip().upper() for d in a.disp.split(",") if d.strip()}
    counts = Counter()
    best = {}   # kepid -> row dict (strongest label kept)
    for row in csv.DictReader(io.StringIO(text)):
        disp = (row.get("koi_disposition") or "").strip().upper()
        counts[disp or "(blank)"] += 1
        if disp not in keep:
            continue
        kepid = (row.get("kepid") or "").strip()
        if not kepid:
            continue
        kid = int(float(kepid))
        cur = best.get(kid)
        if cur is None or RANK.get(disp, 0) > RANK.get(cur["disposition"], 0):
            best[kid] = {
                "target": f"KIC {kid}",
                "disposition": disp,
                "kepler_name": row.get("kepler_name", ""),
                "nasa_period_d": row.get("koi_period", ""),
                "nasa_prad_earth": row.get("koi_prad", ""),
                "star_rad_sun": row.get("koi_srad", ""),
                "star_teff_k": row.get("koi_steff", ""),
            }

    rows = sorted(best.values(), key=lambda r: r["target"])
    with open(a.output, "w") as fh:
        fh.write("# Kepler Objects of Interest host stars (KIC ids), labelled\n")
        for r in rows:
            fh.write(f"{r['target']}\n")

    mcols = ["target", "disposition", "kepler_name", "nasa_period_d",
             "nasa_prad_earth", "star_rad_sun", "star_teff_k"]
    with open(META_CSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=mcols)
        w.writeheader()
        w.writerows(rows)

    label_counts = Counter(r["disposition"] for r in rows)
    print("\nKOI dispositions in the catalog:")
    for d, n in counts.most_common():
        print(f"  {d:<16s} {n}")
    print(f"\nWrote {len(rows)} unique host stars to {a.output}")
    print("  kept by label:", dict(label_counts))
    print(f"Wrote labels + NASA period/radius to {META_CSV}")
    print(f"Now benchmark:  python3 benchmark.py --targets {os.path.basename(a.output)} "
          f"--meta {os.path.basename(META_CSV)} --mission Kepler --sample 150")


if __name__ == "__main__":
    main()
