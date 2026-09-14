#!/usr/bin/env python3
"""scan.py - batch transit scanner. Runs exohunter over a list of stars.

Runs entirely on your machine against the NASA archive - uses NO Claude / AI and
NO tokens. Safe to run in a plain terminal, all day, every day.

What it does per star, in real time:
  - prints progress: [i/N] TARGET | RA/Dec | points | time | VERDICT
  - appends one row to results.csv (the permanent log of EVERYTHING searched)
  - for survivors (CANDIDATE / STRONG_CANDIDATE) it ALSO:
      * saves the plot + light-curve data (kept for follow-up)
      * appends to candidates.csv and hits.log
      * prints a loud banner and (on macOS) pops a desktop notification
  - for NO_SIGNAL / FALSE_POSITIVE it keeps only the one-line log entry
    (no heavy files saved - "run it, bin the crap, keep going")

Resumable: stars already in results.csv are skipped, so a daily run continues
where it left off.

Usage:
    python3 scan.py --demo                    # built-in sample list
    python3 scan.py --targets my_stars.txt    # one target per line
    python3 scan.py --targets tois.txt --pmax 20
"""
import argparse
import csv
import os
import socket
import subprocess
import sys
import time
import warnings

warnings.filterwarnings("ignore")

# A stalled MAST/astroquery download has no timeout of its own, so a dead socket
# blocks the whole scan forever (the classic "stuck on star N for an hour"). A
# process-wide socket timeout turns any read that goes silent for this many
# seconds into an exception, which the loop below treats as transient
# (RETRY-LATER): the star is skipped for now, retried on the next run, and the
# scan keeps moving. This is inactivity, not total download time - a download
# that keeps sending data never trips it.
socket.setdefaulttimeout(120)

import exohunter as eh
import crosscheck
import persistence

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS_CSV = os.path.join(HERE, "results.csv")
CANDIDATES_CSV = os.path.join(HERE, "candidates.csv")
HITS_LOG = os.path.join(HERE, "hits.log")
META_CSV = os.path.join(HERE, "toi_meta.csv")

FIELDS = ["timestamp", "target", "toi", "tfopwg_disp", "verdict", "agreement",
          "period", "nasa_period_d", "depth_ppm", "snr", "est_prad_earth",
          "nasa_prad_earth", "duration_hr", "duration_frac", "depth_odd_ppm",
          "depth_even_ppm", "secondary_ppm", "blend_risk", "reasons",
          "star_rad_sun", "star_teff_k", "ra", "dec", "n_points", "power",
          "method", "proc_sec", "exofop_url", "known", "persist"]


def agreement(verdict, disp):
    """Compare exohunter's verdict to NASA's TFOPWG disposition.

    Buckets that tell you where to spend time:
      RE-DETECTED  - I flagged it, NASA has it confirmed/known (validation)
      MISSED       - I saw nothing, NASA has it confirmed/known (I'm too strict)
      FP-CAUGHT    - I rejected it, NASA calls it a false positive (agree)
      DISAGREE     - my call conflicts with a NASA confirmed/FP label (investigate)
      OPEN         - NASA hasn't resolved it (PC/APC) - the interesting ones
    """
    disp = (disp or "").upper()
    det = verdict in ("STRONG_CANDIDATE", "CANDIDATE")
    rej = verdict in ("FALSE_POSITIVE", "BLEND_SUSPECT")
    null = verdict == "NO_SIGNAL"
    if disp in ("CP", "KP"):
        return "RE-DETECTED" if det else "MISSED" if null else "DISAGREE" if rej else ""
    if disp in ("FP", "FA"):
        return "DISAGREE" if det else "FP-CAUGHT" if (rej or null) else ""
    if disp in ("PC", "APC"):
        return "OPEN"
    return ""


def load_meta():
    meta = {}
    if os.path.exists(META_CSV):
        with open(META_CSV, newline="") as fh:
            for row in csv.DictReader(fh):
                try:
                    meta[int(row["tic"])] = row
                except (ValueError, KeyError):
                    pass
    return meta


def enrich(m, target, meta):
    """Add ExoFOP link, TOI metadata, and an estimated planet radius to a row."""
    digits = "".join(ch for ch in target if ch.isdigit())
    if digits:
        m["exofop_url"] = f"https://exofop.ipac.caltech.edu/tess/target.php?id={digits}"
        md = meta.get(int(digits))
        if md:
            for k in ("toi", "tfopwg_disp", "nasa_period_d", "nasa_prad_earth",
                      "star_rad_sun", "star_teff_k"):
                m[k] = md.get(k, "")
            # planet size estimate from our depth + the star's radius
            try:
                srad = float(md["star_rad_sun"])
                depth = float(m.get("depth_ppm"))
                if srad > 0 and depth > 0:
                    m["est_prad_earth"] = round((depth / 1e6) ** 0.5 * srad * 109.2, 2)
            except (TypeError, ValueError):
                pass
    m["agreement"] = agreement(m.get("verdict", ""), m.get("tfopwg_disp", ""))
    return m

DEMO_TARGETS = [
    "HD 209458", "HD 189733", "WASP-18", "Pi Mensae", "55 Cancri",
    "Zeta2 Reticuli", "Tau Ceti", "Vega", "TOI-700", "GJ 1214",
]


def already_done():
    done = set()
    if os.path.exists(RESULTS_CSV):
        with open(RESULTS_CSV, newline="") as fh:
            for row in csv.DictReader(fh):
                done.add(row["target"])
    return done


def append_row(path, row, header):
    new = not os.path.exists(path)
    with open(path, "a", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=header)
        if new:
            w.writeheader()
        w.writerow({k: row.get(k, "") for k in header})


def notify(title, message):
    """Best-effort macOS desktop notification; silently ignored elsewhere."""
    if sys.platform == "darwin":
        try:
            subprocess.run(
                ["osascript", "-e",
                 f'display notification "{message}" with title "{title}"'],
                check=False, capture_output=True, timeout=5)
        except Exception:
            pass


def vet_survivor(m, target, mission):
    """Vet one survivor in place: catalogue cross-check + multi-sector persistence.

    Returns (known, persist) as short strings and appends the catalogue match note
    to m['reasons']. Reuses the ephemeris the scan already found (no BLS re-run);
    persistence still needs a per-sector download, so this runs on survivors only.
    """
    known = "?"
    try:
        kc = crosscheck.known_object_check(m.get("ra"), m.get("dec"), m.get("period"))
        if kc.get("known") is True:
            known = "yes"
            if kc.get("note"):
                m["reasons"] = (m["reasons"] + ";" if m["reasons"] else "") + kc["note"]
        elif kc.get("known") is False:
            known = "no"
    except Exception:
        pass
    persist = "?"
    try:
        dur_d = float(m.get("duration_hr", 0)) / 24.0
        pc = persistence.check_target(target, mission=mission, period=m.get("period"),
                                      t0=m.get("t0"), dur=dur_d)
        persist = f"{pc['persistence_ok']}({pc['n_with_dip']}/{pc['n_covering']})"
    except Exception:
        pass
    return known, persist


def fmt_elapsed(seconds):
    """Human-readable elapsed time, e.g. '1h 23m', '4m 12s', '38s'."""
    h, rem = divmod(int(seconds), 3600)
    m, s = divmod(rem, 60)
    if h:
        return f"{h}h {m}m"
    if m:
        return f"{m}m {s}s"
    return f"{s}s"


def beep(times=3):
    """Audible alert: macOS system beep, falling back to the terminal bell."""
    if sys.platform == "darwin":
        try:
            subprocess.run(["osascript", "-e", f"beep {times}"],
                           check=False, capture_output=True, timeout=5)
            return
        except Exception:
            pass
    try:
        sys.stdout.write("\a" * times)
        sys.stdout.flush()
    except Exception:
        pass


def get_region_targets(region, maxmag, limit, minmag=None):
    """Return TIC ids for stars in a sky patch.

    region: "RA DEC RADIUS" in degrees, e.g. "285.7 44.5 0.5".
    Pulls the TESS Input Catalog for that circle, keeps stars with Tmag in
    [minmag, maxmag), and returns up to `limit` of them (brightest first).
    Set minmag (e.g. 13.5) to target the faint stars QLP does not search - the
    genuinely under-analysed population - instead of the already-searched bright ones.
    """
    from astroquery.mast import Catalogs
    parts = region.split()
    if len(parts) != 3:
        sys.exit('--region must be "RA DEC RADIUS" in degrees, e.g. "285.7 44.5 0.5"')
    ra, dec, radius = map(float, parts)
    band = f"Tmag < {maxmag}" if minmag is None else f"{minmag} <= Tmag < {maxmag}"
    print(f"Querying TESS Input Catalog around RA={ra} Dec={dec}, radius={radius} deg "
          f"({band}) ...", flush=True)
    cat = Catalogs.query_region(f"{ra} {dec}", radius=radius, catalog="TIC")
    cat = cat[cat["Tmag"] < maxmag]
    if minmag is not None:
        cat = cat[cat["Tmag"] >= minmag]
    cat.sort("Tmag")
    ids = [f"TIC {int(r)}" for r in cat["ID"][:limit]]
    print(f"  {len(cat)} stars in region with {band}; scanning {len(ids)}.",
          flush=True)
    return ids


def load_targets(args):
    if args.demo:
        return DEMO_TARGETS
    if args.region:
        return get_region_targets(args.region, args.maxmag, args.limit, args.minmag)
    if args.targets:
        with open(args.targets) as fh:
            return [ln.strip() for ln in fh if ln.strip() and not ln.startswith("#")]
    sys.exit('Provide one of: --demo, --region "RA DEC RADIUS", or --targets FILE')


def main():
    p = argparse.ArgumentParser(description="Batch transit scanner (token-free).")
    p.add_argument("--targets", help="file with one target per line")
    p.add_argument("--demo", action="store_true", help="use the built-in sample list")
    p.add_argument("--region", help='sweep a sky patch: "RA DEC RADIUS" in degrees, '
                   'e.g. "285.7 44.5 0.3"')
    p.add_argument("--maxmag", type=float, default=12.0,
                   help="region mode: only stars brighter than this TESS mag (default 12)")
    p.add_argument("--minmag", type=float, default=None,
                   help="region mode: only stars fainter than this TESS mag. Set e.g. 13.5 to "
                        "target the faint stars QLP does not search (the under-analysed population).")
    p.add_argument("--limit", type=int, default=300,
                   help="region mode: max stars to scan (brightest first)")
    p.add_argument("--mission", default="TESS", choices=["TESS", "Kepler", "K2", ""])
    p.add_argument("--pmin", type=float, default=0.5)
    p.add_argument("--pmax", type=float, default=15.0)
    p.add_argument("--rescan", action="store_true", help="ignore results.csv and redo all")
    p.add_argument("--tls", action="store_true",
                   help="use Transit Least Squares instead of BLS (more sensitive, slower)")
    p.add_argument("--tesscut", action="store_true",
                   help="search TESS full-frame-image cutouts (stars with no SPOC light curve) "
                        "- this is how you reach the under-searched FFI population")
    p.add_argument("--all-sectors", action="store_true",
                   help="stitch all available sectors/quarters before searching (thorough, slower)")
    p.add_argument("--vet", action="store_true",
                   help="vet each survivor DURING the run: cross-check catalogues (already known?) "
                        "and multi-sector persistence, so candidates.csv holds only genuinely-new "
                        "survivors. Adds a catalogue lookup + a re-download per survivor (use on hunts).")
    p.add_argument("--beep", action=argparse.BooleanOptionalAction, default=True,
                   help="audible beep on a STRONG_CANDIDATE hit (on by default; --no-beep to silence)")
    args = p.parse_args()
    method = "tls" if args.tls else "bls"

    targets = load_targets(args)
    meta = load_meta()
    done = set() if args.rescan else already_done()
    todo = [t for t in targets if t not in done]

    n = len(todo)
    print(f"exohunter batch scan  |  {n} to process "
          f"({len(targets) - n} already logged)  |  results -> results.csv", flush=True)
    print("-" * 78, flush=True)

    hits = 0
    t_start = time.time()
    for i, target in enumerate(todo, 1):
        t0 = time.time()
        try:
            m = eh.analyze(target, mission=args.mission, pmin=args.pmin,
                           pmax=args.pmax, all_sectors=args.all_sectors,
                           tesscut=args.tesscut, make_plot="survivors", method=method)
            m["proc_sec"] = round(time.time() - t0, 1)
            m["timestamp"] = time.strftime("%Y-%m-%d %H:%M:%S")
            m["reasons"] = m.get("reasons", "")
            enrich(m, target, meta)
            is_survivor = m["verdict"] in ("CANDIDATE", "STRONG_CANDIDATE")
            if args.vet and is_survivor:
                m["known"], m["persist"] = vet_survivor(m, target, args.mission)
            append_row(RESULTS_CSV, m, FIELDS)

            ra = f"{m['ra']:.3f}" if isinstance(m["ra"], float) else m["ra"]
            dec = f"{m['dec']:.3f}" if isinstance(m["dec"], float) else m["dec"]
            line = (f"#{i} of {n} | {target:<18s} | RA={ra} Dec={dec} | "
                    f"{m['n_points']} pts | {m['proc_sec']}s | "
                    f"P={m['period']:.3f}d depth={m['depth_ppm']:.0f}ppm "
                    f"SNR={m['snr']:.1f} | {m['verdict']} | "
                    f"elapsed {fmt_elapsed(time.time() - t_start)}")
            print(line, flush=True)

            # A "new" candidate is a survivor that isn't already catalogued. When
            # --vet is off, m['known'] is unset so every survivor still alerts (old
            # behaviour); when --vet is on, already-known survivors are logged but
            # not alerted, so the hits list is only genuinely-new candidates.
            if is_survivor and m.get("known") != "yes":
                hits += 1
                banner = (f">>> POSSIBLE PLANET: {target}  ({m['verdict']})  "
                          f"P={m['period']:.4f}d depth={m['depth_ppm']:.0f}ppm "
                          f"SNR={m['snr']:.1f} <<<")
                print("    " + "!" * 66, flush=True)
                print("    " + banner, flush=True)
                if args.vet:
                    print(f"    vet: known={m.get('known')}  persist={m.get('persist')}", flush=True)
                print("    plot: " + m.get("plot_path", "(none)"), flush=True)
                print("    " + "!" * 66, flush=True)
                append_row(CANDIDATES_CSV, m, FIELDS)
                with open(HITS_LOG, "a") as fh:
                    fh.write(f"{m['timestamp']}  {banner}  plot={m.get('plot_path','')}\n")
                notify(f"exohunter: {m['verdict']}", f"{target}  P={m['period']:.3f}d SNR={m['snr']:.0f}")
                if args.beep and m["verdict"] == "STRONG_CANDIDATE":
                    beep()

        except Exception as e:
            msg = f"{type(e).__name__}: {str(e)[:80]}"
            transient = any(s in str(e) or s in type(e).__name__ for s in (
                "Connection", "RemoteDisconnected", "Timeout", "timed out",
                "reset by peer", "ServiceUnavailable", "HTTPError", "TimeoutError"))
            el = fmt_elapsed(time.time() - t_start)
            if transient:
                # Don't log it -> a later run will retry this star.
                print(f"#{i} of {n} | {target:<18s} | RETRY-LATER ({msg}) | elapsed {el}", flush=True)
            else:
                print(f"#{i} of {n} | {target:<18s} | SKIPPED ({msg}) | elapsed {el}", flush=True)
                errrow = {"timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
                          "target": target, "verdict": "ERROR", "reasons": msg,
                          "proc_sec": round(time.time() - t0, 1)}
                enrich(errrow, target, meta)
                append_row(RESULTS_CSV, errrow, FIELDS)

    elapsed = time.time() - t_start
    print("-" * 78, flush=True)
    print(f"Done. {n} processed in {elapsed/60:.1f} min | {hits} candidate(s) found "
          f"| full log: results.csv | hits: hits.log / candidates.csv", flush=True)


if __name__ == "__main__":
    main()
