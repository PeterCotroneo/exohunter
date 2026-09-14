#!/usr/bin/env python3
"""benchmark.py - measure exohunter's accuracy against a LABELLED set of stars.

This is how we answer "will it be accurate on this dataset?" with a number rather
than a hope. Point it at a target list plus a meta CSV that carries a ground-truth
`disposition` column (CONFIRMED / FALSE POSITIVE / CANDIDATE), and it runs the same
transit search + verdict on each star, then scores the verdicts against the labels.

Because thresholds are mission-specific (see MISSION_PARAMS in exohunter.py), run
this separately per mission and read the number for THAT mission - never assume a
TESS result carries over to Kepler/K2.

Usage:
    # build the labelled Kepler set first:
    python3 make_koi_list.py
    # fast number (BLS, one sector/quarter - quick, the everyday config):
    python3 benchmark.py --targets koi_targets.txt --meta koi_meta.csv \
            --mission Kepler --sample 150 --pmax 40 --mode fast
    # thorough number (TLS, all stacked sectors - slower, more sensitive):
    python3 benchmark.py --targets koi_targets.txt --meta koi_meta.csv \
            --mission Kepler --sample 150 --pmax 40 --mode thorough

Run both and report the pair: "fast" is the realistic day-to-day recall; "thorough"
is the ceiling when you spend the compute. They write to separate CSVs so neither
run clobbers the other.

Outputs a scoreboard to the terminal and a per-star CSV (default
benchmark_<mission>_<mode>.csv). Resumable: stars already in the output CSV are
skipped. Token-free (NASA archive only).
"""
import argparse
import csv
import os
import random
import time
import warnings
from collections import Counter, defaultdict

warnings.filterwarnings("ignore")

import exohunter as eh

HERE = os.path.dirname(os.path.abspath(__file__))

# how each ground-truth label expects a verdict to fall
DETECT = ("STRONG_CANDIDATE", "CANDIDATE")
REJECT = ("FALSE_POSITIVE", "BLEND_SUSPECT")


def bucket(verdict, disposition):
    """Score one verdict against its ground-truth label (same buckets as scan.py)."""
    disp = (disposition or "").strip().upper()
    det = verdict in DETECT
    rej = verdict in REJECT
    null = verdict == "NO_SIGNAL"
    if disp == "CONFIRMED":
        return "RE-DETECTED" if det else "MISSED" if null else "DISAGREE"
    if disp in ("FALSE POSITIVE", "FALSE_POSITIVE"):
        return "DISAGREE" if det else "FP-CAUGHT"
    return "OPEN"   # CANDIDATE / unlabelled: no ground truth to score against


def load_meta(path):
    meta = {}
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            meta[row["target"]] = row
    return meta


def already_done(path):
    done = set()
    if os.path.exists(path):
        with open(path, newline="") as fh:
            for row in csv.DictReader(fh):
                done.add(row["target"])
    return done


def sample_targets(targets, meta, per_class, seed):
    """Cap the run to `per_class` stars of each disposition (balanced, reproducible)."""
    if not per_class:
        return targets
    by = defaultdict(list)
    for t in targets:
        by[(meta.get(t, {}).get("disposition") or "").upper()].append(t)
    rng = random.Random(seed)
    out = []
    for disp, lst in by.items():
        rng.shuffle(lst)
        out.extend(lst[:per_class])
    rng.shuffle(out)
    return out


def main():
    p = argparse.ArgumentParser(description="Benchmark exohunter against labelled stars.")
    p.add_argument("--targets", required=True, help="one target per line")
    p.add_argument("--meta", required=True, help="CSV with target + disposition columns")
    p.add_argument("--mission", default="Kepler", choices=["TESS", "Kepler", "K2", ""])
    p.add_argument("--mode", default="fast", choices=["fast", "stacked", "thorough"],
                   help="fast = BLS on one sector/quarter (quick, the everyday number); "
                        "stacked = BLS on all sectors/quarters stitched (recovers faint "
                        "planets a single segment misses; the practical thorough number); "
                        "thorough = TLS on all stacked sectors (most sensitive but "
                        "minutes-to-tens-of-minutes per star, so only for small samples; "
                        "needs `pip install transitleastsquares`).")
    p.add_argument("--pmin", type=float, default=0.5)
    p.add_argument("--pmax", type=float, default=40.0)
    p.add_argument("--sample", type=int, default=0,
                   help="cap per disposition (0 = all). Balanced + reproducible.")
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--out", default="",
                   help="per-star CSV (default benchmark_<mission>_<mode>.csv)")
    p.add_argument("--rescan", action="store_true")
    p.add_argument("--tesscut", action="store_true",
                   help="use full-frame-image cutouts (needed for faint FFI-only stars)")
    a = p.parse_args()

    # the two configurations we report side by side
    method = "tls" if a.mode == "thorough" else "bls"
    all_sectors = a.mode in ("stacked", "thorough")
    out = a.out or os.path.join(HERE, f"benchmark_{a.mission or 'auto'}_{a.mode}.csv")
    meta = load_meta(os.path.join(HERE, a.meta) if not os.path.isabs(a.meta) else a.meta)
    tfile = os.path.join(HERE, a.targets) if not os.path.isabs(a.targets) else a.targets
    with open(tfile) as fh:
        targets = [ln.strip() for ln in fh if ln.strip() and not ln.startswith("#")]

    targets = sample_targets(targets, meta, a.sample, a.seed)
    done = set() if a.rescan else already_done(out)
    todo = [t for t in targets if t not in done]

    fields = ["target", "disposition", "verdict", "bucket", "period",
              "nasa_period_d", "depth_ppm", "snr", "blend_risk", "reasons", "proc_sec"]
    new_file = a.rescan or not os.path.exists(out)
    fh_out = open(out, "w" if a.rescan else "a", newline="")
    writer = csv.DictWriter(fh_out, fieldnames=fields)
    if new_file:
        writer.writeheader()

    print(f"benchmark | mission={a.mission} | mode={a.mode} ({method}"
          f"{', all-sectors' if all_sectors else ', single-sector'}) | {len(todo)} to run "
          f"({len(targets) - len(todo)} already scored) | out -> {out}", flush=True)
    print("-" * 78, flush=True)

    tally = Counter()
    t_start = time.time()
    for i, target in enumerate(todo, 1):
        disp = meta.get(target, {}).get("disposition", "")
        t0 = time.time()
        try:
            m = eh.analyze(target, mission=a.mission, pmin=a.pmin, pmax=a.pmax,
                           all_sectors=all_sectors, tesscut=a.tesscut,
                           make_plot=False, method=method)
            b = bucket(m["verdict"], disp)
            # coerce to plain float first: Kepler flux is masked, and round() has no
            # __round__ for MaskedNDArray (TESS SPOC flux isn't masked, so it only bit here)
            pf, sf, df = float(m["period"]), float(m["snr"]), float(m["depth_ppm"])
            row = {"target": target, "disposition": disp, "verdict": m["verdict"],
                   "bucket": b, "period": round(pf, 4),
                   "nasa_period_d": meta.get(target, {}).get("nasa_period_d", ""),
                   "depth_ppm": round(df, 0) if df == df else "",
                   "snr": round(sf, 1), "blend_risk": m.get("blend_risk", ""),
                   "reasons": m.get("reasons", ""), "proc_sec": round(time.time() - t0, 1)}
            writer.writerow(row); fh_out.flush()
            tally[b] += 1
            print(f"#{i} of {len(todo)} | {target:<16s} | {disp:<14s} | "
                  f"SNR={m['snr']:.1f} | {m['verdict']:<16s} -> {b}", flush=True)
        except Exception as e:
            msg = f"{type(e).__name__}: {str(e)[:70]}"
            transient = any(s in str(e) or s in type(e).__name__ for s in (
                "Connection", "RemoteDisconnected", "Timeout", "timed out",
                "reset by peer", "ServiceUnavailable", "HTTPError"))
            tag = "RETRY-LATER" if transient else "SKIPPED"
            print(f"#{i} of {len(todo)} | {target:<16s} | {tag} ({msg})", flush=True)
            if not transient:
                writer.writerow({"target": target, "disposition": disp,
                                 "verdict": "ERROR", "bucket": "ERROR", "reasons": msg})
                fh_out.flush()
    fh_out.close()

    scoreboard(out, a.mission, a.mode, time.time() - t_start)


def scoreboard(out, mission, mode, elapsed):
    """Print the accuracy numbers from the full output CSV (all runs, resumed included)."""
    rows = list(csv.DictReader(open(out)))
    tally = Counter(r["bucket"] for r in rows)
    conf = [r for r in rows if r["disposition"].upper() == "CONFIRMED"]
    fp = [r for r in rows if r["disposition"].upper() in ("FALSE POSITIVE", "FALSE_POSITIVE")]

    def pct(n, d):
        return f"{100*n/d:.0f}%" if d else "n/a"

    print("-" * 78, flush=True)
    print(f"BENCHMARK SCOREBOARD  |  mission={mission}  |  mode={mode}  |  {len(rows)} stars  "
          f"|  {elapsed/60:.1f} min this run", flush=True)
    print("  buckets:", dict(tally), flush=True)
    if conf:
        det = sum(1 for r in conf if r["verdict"] in DETECT)
        miss = sum(1 for r in conf if r["verdict"] == "NO_SIGNAL")
        blend = sum(1 for r in conf if r["verdict"] == "BLEND_SUSPECT")
        print(f"  CONFIRMED planets (n={len(conf)}): "
              f"detected {pct(det,len(conf))} | missed(NO_SIGNAL) {pct(miss,len(conf))} "
              f"| blend-demoted {pct(blend,len(conf))}", flush=True)
    if fp:
        caught = sum(1 for r in fp if r["verdict"] in REJECT or r["verdict"] == "NO_SIGNAL")
        print(f"  FALSE POSITIVES (n={len(fp)}): correctly rejected {pct(caught,len(fp))} "
              f"| leaked as candidate {pct(len(fp)-caught,len(fp))}", flush=True)
    print("  (detection rate on CONFIRMED is the headline accuracy number for this "
          "mission)", flush=True)


if __name__ == "__main__":
    main()
