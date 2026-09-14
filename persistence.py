#!/usr/bin/env python3
"""persistence.py - does the transit actually recur across multiple sectors?

A real transiting planet shows the SAME dip, at the SAME depth, in every sector
(or Kepler quarter) that covers the transit phase. A one-off event, a systematic,
or a single-transit fluke does not. This generalises the odd/even check from two
groups to every sector, and is a key confidence filter for a candidate flagged on
new data.

Given a target, it downloads the per-sector light curves, finds the ephemeris with
a Box Least Squares search (via exohunter), then folds EACH sector on that
ephemeris and measures the per-sector depth. It reports how many sectors show a
consistent dip and whether the signal persists.

Honest limit: with only one sector of data there is nothing to compare, so
persistence cannot be assessed (the within-sector odd/even check is the proxy
there). This returns persistence_ok = None in that case rather than a false pass.

Usage (standalone test):
    python3 persistence.py "WASP-18" --pmax 3 --max-sectors 3
"""
import argparse
import math
import socket
import warnings

warnings.filterwarnings("ignore")
import numpy as np

import exohunter as eh

MIN_IN_TRANSIT = 3      # a sector needs this many in-transit points to "cover" a transit
PER_SECTOR_SNR = 3.0    # per-sector dip significance to count as "showing the transit"
DEPTH_CONSISTENCY = 3.0  # max/min depth across detecting sectors below this = consistent


def _sector_label(lc, i):
    for key in ("SECTOR", "QUARTER", "CAMPAIGN"):
        v = lc.meta.get(key)
        if v is not None:
            return f"{key.lower()} {v}"
    return f"segment {i}"


def persistence_from_lcs(lcs, period, t0, dur):
    """Fold each per-sector light curve on the ephemeris; measure recurrence.

    Returns dict(n_segments, n_covering, n_with_dip, depths_ppm, consistency,
    persistence_ok, note).
    """
    covering = 0
    dips = []          # depths (ppm) of sectors that show a significant dip
    per_sector = []    # (label, depth_ppm or None, snr, shows_dip)
    for i, lc in enumerate(lcs):
        flat = lc.flatten(window_length=1001, niters=3)
        t = np.asarray(flat.time.value, float)
        f = np.asarray(flat.flux.value, float)
        phase = (t - t0 + 0.5 * period) % period - 0.5 * period
        in_tr = np.abs(phase) < (dur / 2)
        n_in = int(in_tr.sum())
        label = _sector_label(lc, i)
        if n_in < MIN_IN_TRANSIT:
            per_sector.append((label, None, 0.0, False))
            continue
        covering += 1
        oot = ~in_tr
        depth = (1 - np.nanmedian(f[in_tr])) * 1e6
        scatter = np.nanstd(f[oot]) * 1e6
        snr = depth / (scatter / math.sqrt(n_in)) if scatter > 0 else 0.0
        shows = depth > 0 and snr >= PER_SECTOR_SNR
        if shows:
            dips.append(depth)
        per_sector.append((label, round(depth, 0), round(snr, 1), shows))

    n_with_dip = len(dips)
    consistency = (max(dips) / min(dips)) if (dips and min(dips) > 0) else None

    if len(lcs) <= 1 or covering <= 1:
        ok = None
        note = f"only {covering} sector(s) cover the transit - cannot assess persistence"
    else:
        ok = (n_with_dip >= 2 and n_with_dip >= 0.5 * covering
              and consistency is not None and consistency < DEPTH_CONSISTENCY)
        note = (f"transit recurs in {n_with_dip}/{covering} covering sectors, "
                f"depth consistency {consistency:.1f}x" if consistency
                else f"dip in {n_with_dip}/{covering} covering sectors")

    return {"n_segments": len(lcs), "n_covering": covering, "n_with_dip": n_with_dip,
            "depths_ppm": dips, "consistency": consistency, "persistence_ok": ok,
            "per_sector": per_sector, "note": note}


def check_target(target, mission="TESS", pmin=0.5, pmax=15.0, max_sectors=None,
                 period=None, t0=None, dur=None):
    """Download per-sector light curves and assess persistence.

    If period/t0/dur are supplied (the scan already found them), skip the BLS
    re-search and fold on the given ephemeris; otherwise find it with BLS.
    """
    import lightkurve as lk
    socket.setdefaulttimeout(120)
    kwargs = {"mission": mission} if mission else {}
    if mission == "TESS":
        kwargs.update(author="SPOC", cadence="short")
    sr = lk.search_lightcurve(target, **kwargs)
    if len(sr) == 0 and mission == "TESS":
        sr = lk.search_lightcurve(target, mission="TESS")
    if len(sr) == 0:
        raise ValueError("no light curves found")
    if max_sectors:
        sr = sr[:max_sectors]
    col = sr.download_all()
    lcs = [lc.remove_nans().normalize() for lc in col]
    if period is None or t0 is None or dur is None:
        stitched = col.stitch().remove_nans().normalize()
        flat = stitched.flatten(window_length=1001, niters=3)
        period, t0, dur, *_ = eh.detect(flat, "bls", pmin, pmax)
    res = persistence_from_lcs(lcs, period, t0, dur)
    res.update(target=target, period=period, t0=t0, duration_hr=dur * 24)
    return res


def main():
    p = argparse.ArgumentParser(description="Multi-sector transit persistence check.")
    p.add_argument("target")
    p.add_argument("--mission", default="TESS", choices=["TESS", "Kepler", "K2", ""])
    p.add_argument("--pmin", type=float, default=0.5)
    p.add_argument("--pmax", type=float, default=15.0)
    p.add_argument("--max-sectors", type=int, default=None)
    a = p.parse_args()
    r = check_target(a.target, a.mission, a.pmin, a.pmax, a.max_sectors)
    print(f"{r['target']}  P={r['period']:.4f}d  dur={r['duration_hr']:.2f}h")
    print(f"  segments downloaded: {r['n_segments']} | covering the transit: {r['n_covering']}")
    for label, depth, snr, shows in r["per_sector"]:
        d = f"{depth:.0f} ppm" if depth is not None else "no coverage"
        print(f"    {label:<12s}: {d:<12s} snr={snr}  {'DIP' if shows else '-'}")
    print(f"  persistence_ok: {r['persistence_ok']}  ({r['note']})")


if __name__ == "__main__":
    main()
