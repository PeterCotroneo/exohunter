#!/usr/bin/env python3
"""single_transit.py - find lone transit-like dips (long-period / single-transit).

BLS and TLS assume a repeating signal: they fold the light curve on a period and
need two or more transits. A planet whose orbital period is longer than the
observing window transits only once or twice, so the periodicity searches miss it.
This is the biggest blind spot of the mission pipelines and where citizen-science
projects still find planets.

This module does NOT assume periodicity. It slides a transit-width window across
the flattened light curve and looks for individual, isolated, significant dips
that have a real (multi-point, physically-sized) shape - not single-cadence
spikes and not edge/gap artefacts.

Vetting is necessarily different from the recurring case: with one event there is
no odd/even and no recurrence to check, so we require enough in-transit points, a
physical duration, a symmetric local baseline on both sides (isolation), and we
reject events at data edges or gaps.

Usage:
    python3 single_transit.py "TIC 307210830" --mission TESS
    python3 single_transit.py --selftest        # synthetic injection-recovery
"""
import argparse
import warnings

warnings.filterwarnings("ignore")
import numpy as np

# candidate transit durations to scan, in days (grazing-short to long-ingress)
DUR_GRID_D = (0.04, 0.08, 0.15, 0.25, 0.4, 0.6)
MIN_SNR = 7.0          # dip significance floor
MIN_IN_TRANSIT = 4     # need at least this many points inside the dip (no 1-point spikes)
BASELINE_SPAN = 3.0    # local baseline uses +/- this many durations either side
EDGE_PAD_D = 0.3       # ignore events within this many days of a data edge/gap


def _gap_edges(t, max_gap=0.5):
    """Times bordering data gaps larger than max_gap days (plus the global ends)."""
    edges = [t[0], t[-1]]
    dt = np.diff(t)
    for i in np.where(dt > max_gap)[0]:
        edges.append(t[i]); edges.append(t[i + 1])
    return np.array(edges)


def single_transit_search(t, f, dur_grid=DUR_GRID_D, min_snr=MIN_SNR,
                          min_in=MIN_IN_TRANSIT):
    """Scan a flattened light curve for isolated transit-like dips.

    t, f: time (days) and normalised flux (median ~1), finite, sorted by time.
    Returns a list of events sorted by SNR:
        dict(t0, duration_d, depth_ppm, snr, n_in)
    """
    t = np.asarray(t, float); f = np.asarray(f, float)
    m = np.isfinite(t) & np.isfinite(f)
    t, f = t[m], f[m]
    order = np.argsort(t); t, f = t[order], f[order]
    edges = _gap_edges(t)

    events = []
    for dur in dur_grid:
        half = dur / 2.0
        lo = np.searchsorted(t, t - half, side="left")
        hi = np.searchsorted(t, t + half, side="right")
        blo = np.searchsorted(t, t - BASELINE_SPAN * dur, side="left")
        bhi = np.searchsorted(t, t + BASELINE_SPAN * dur, side="right")
        for i in range(len(t)):
            n_in = hi[i] - lo[i]
            if n_in < min_in:
                continue
            # skip near a data edge/gap
            if np.min(np.abs(edges - t[i])) < EDGE_PAD_D:
                continue
            in_flux = f[lo[i]:hi[i]]
            # local baseline = points in the wider window but OUTSIDE the dip
            base = np.concatenate([f[blo[i]:lo[i]], f[hi[i]:bhi[i]]])
            if base.size < min_in:
                continue
            depth = np.median(base) - np.median(in_flux)     # positive = a dip
            scatter = np.std(base)
            if scatter <= 0 or depth <= 0:
                continue
            snr = depth / (scatter / np.sqrt(n_in))
            if snr >= min_snr:
                events.append({"t0": float(t[i]), "duration_d": float(dur),
                               "depth_ppm": float(depth * 1e6), "snr": float(snr),
                               "n_in": int(n_in)})

    # deduplicate: collapse events whose centres fall within one duration, keep best SNR
    events.sort(key=lambda e: -e["snr"])
    kept = []
    for e in events:
        if all(abs(e["t0"] - k["t0"]) > max(e["duration_d"], k["duration_d"]) for k in kept):
            kept.append(e)
    return kept


def check_target(target, mission="TESS", all_sectors=True, tesscut=False):
    import exohunter as eh
    lc = eh.load_lightcurve(target, mission, all_sectors=all_sectors, tesscut=tesscut)
    flat = lc.flatten(window_length=1001, niters=3)
    t = np.asarray(flat.time.value, float)
    f = np.asarray(flat.flux.value, float)
    return single_transit_search(t, f)


def _selftest():
    """Inject one box transit into synthetic noise and confirm recovery."""
    rng = np.random.default_rng(0)
    t = np.sort(rng.uniform(0, 27, 20000))          # ~27-day sector, irregular sampling
    f = 1 + rng.normal(0, 3e-4, t.size)             # 300 ppm scatter
    t0_true, dur_true, depth_true = 13.5, 0.2, 0.004  # single 4000 ppm, 0.2 d dip
    f[np.abs(t - t0_true) < dur_true / 2] -= depth_true
    ev = single_transit_search(t, f)
    print(f"injected: t0={t0_true} dur={dur_true}d depth={depth_true*1e6:.0f}ppm")
    print(f"events found: {len(ev)}")
    for e in ev[:5]:
        print(f"  t0={e['t0']:.2f} dur={e['duration_d']}d depth={e['depth_ppm']:.0f}ppm "
              f"snr={e['snr']:.1f} n_in={e['n_in']}")
    ok = ev and min(abs(e["t0"] - t0_true) for e in ev) < 0.2
    print("RECOVERED at correct time" if ok else "FAILED to recover")
    return ok


def main():
    p = argparse.ArgumentParser(description="Single-transit / long-period dip search.")
    p.add_argument("target", nargs="?")
    p.add_argument("--mission", default="TESS", choices=["TESS", "Kepler", "K2", ""])
    p.add_argument("--tesscut", action="store_true")
    p.add_argument("--single-sector", action="store_true", help="do not stitch all sectors")
    p.add_argument("--selftest", action="store_true")
    a = p.parse_args()
    if a.selftest or not a.target:
        _selftest(); return
    ev = check_target(a.target, a.mission, all_sectors=not a.single_sector, tesscut=a.tesscut)
    print(f"{a.target}: {len(ev)} single-transit candidate event(s)")
    for e in ev[:10]:
        print(f"  t0={e['t0']:.3f}  dur={e['duration_d']}d  depth={e['depth_ppm']:.0f}ppm  "
              f"snr={e['snr']:.1f}  pts={e['n_in']}")


if __name__ == "__main__":
    main()
