#!/usr/bin/env python3
"""exohunter - transit-search engine for a single star.

Downloads a light curve, flattens it, runs a Box Least Squares (BLS) transit
search, vets the result, and assigns a deterministic verdict. Used directly for
one star, or imported by scan.py for batch runs.

Single-star usage:
    python3 exohunter.py "WASP-18"
    python3 exohunter.py "Kepler-10" --mission Kepler
    python3 exohunter.py "TIC 261136679" --pmin 1 --pmax 20 --plot

Nothing here uses any AI/network beyond the NASA MAST archive; results are
deterministic (same star in, same verdict out).
"""
import argparse
import os
import warnings

warnings.filterwarnings("ignore")
import numpy as np
import lightkurve as lk

HERE = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(HERE, "data")
PLOTS_DIR = os.path.join(HERE, "plots")

# ---- verdict thresholds (tune these against known planets/false positives) ----
SNR_STRONG = 15.0        # >= this and clean -> STRONG_CANDIDATE
SNR_MIN = 8.0            # below this -> NO_SIGNAL
POWER_MIN = 200.0        # BLS peak power floor for a real detection
DEPTH_MIN_PPM = 50.0     # shallower than this at low SNR is noise
DEPTH_MAX_PPM = 50000.0  # deeper than ~5% is usually a star, not a planet
ODD_EVEN_MAX_FRAC = 0.5  # odd/even depth mismatch above this -> binary flag
SECONDARY_FRAC = 0.3     # secondary depth > this fraction of transit -> binary
SECONDARY_MIN_PPM = 200.0
DURATION_MAX_FRAC = 0.15   # transit longer than 15% of the orbit is unphysical
BRIGHTER_MARGIN_MAG = 0.0  # only a neighbor BRIGHTER than the target can plausibly fake it

# ---- mission-specific knobs ----------------------------------------------------
# Some thresholds are NOT portable between missions and must be set per mission.
# The clearest is the blend radius: it is an angular distance, but "how many
# pixels" that is depends entirely on the detector. TESS pixels are ~21"/px, so
# 42" is ~2 pixels; Kepler/K2 pixels are ~4"/px, so the same contamination scale
# (~3 px) is only ~12". SNR floors are kept equal to TESS for now and must be
# re-benchmarked against Kepler/K2 known objects before they are trusted there.
MISSION_PARAMS = {
    "TESS":   {"contam_arcsec": 42.0, "snr_min": SNR_MIN, "snr_strong": SNR_STRONG},
    "Kepler": {"contam_arcsec": 12.0, "snr_min": SNR_MIN, "snr_strong": SNR_STRONG},
    "K2":     {"contam_arcsec": 12.0, "snr_min": SNR_MIN, "snr_strong": SNR_STRONG},
}


def params_for(mission):
    """Mission-specific thresholds, defaulting to TESS for unknown/blank missions."""
    return MISSION_PARAMS.get(mission or "TESS", MISSION_PARAMS["TESS"])


def safe_name(target):
    return "".join(c if c.isalnum() else "_" for c in target).strip("_")


def load_lightcurve(target, mission="TESS", all_sectors=False, tesscut=False):
    """Search MAST and return a cleaned, normalized LightCurve (or raise)."""
    if tesscut:
        sr = lk.search_tesscut(target)
        if len(sr) == 0:
            raise ValueError("no TESScut data")
        tpf = sr[0].download(cutout_size=11)
        return tpf.to_lightcurve(aperture_mask="threshold").remove_nans().normalize()

    kwargs = {}
    if mission:
        kwargs["mission"] = mission
    if mission == "TESS":
        kwargs.update(author="SPOC", cadence="short")
    sr = lk.search_lightcurve(target, **kwargs)
    if len(sr) == 0 and mission == "TESS":
        sr = lk.search_lightcurve(target, mission="TESS")
    if len(sr) == 0:
        sr = lk.search_lightcurve(target)
    if len(sr) == 0:
        raise ValueError("no light curves found")
    if all_sectors:
        return sr.download_all().stitch().remove_nans().normalize()
    return sr[0].download().remove_nans().normalize()


def detect(flat, method, pmin, pmax):
    """Find the best transit period. method='bls' (fast, built into lightkurve)
    or 'tls' (Transit Least Squares - fits a realistic transit shape, more
    sensitive to small planets, but slower). Returns period, t0, duration (days),
    a peak score, and the (periods, power) arrays for plotting.
    """
    if method == "tls":
        from transitleastsquares import transitleastsquares
        t = np.asarray(flat.time.value, float)
        f = np.asarray(flat.flux.value, float)
        m = np.isfinite(t) & np.isfinite(f)
        model = transitleastsquares(t[m], f[m])
        r = model.power(period_min=pmin, period_max=pmax, show_progress_bar=False)
        return (float(r.period), float(r.T0), float(r.duration), float(r.SDE),
                np.asarray(r.periods, float), np.asarray(r.power, float))
    grid = np.linspace(pmin, pmax, 20000)
    bls = flat.to_periodogram(method="bls", period=grid, frequency_factor=500)
    return (float(bls.period_at_max_power.value),
            float(bls.transit_time_at_max_power.value),
            float(bls.duration_at_max_power.value),
            float(bls.max_power.value),
            np.asarray(bls.period.value, float),
            np.asarray(bls.power.value, float))


def analyze(target, mission="TESS", pmin=0.5, pmax=15.0, all_sectors=False,
            tesscut=False, make_plot=False, do_centroid=True, method="bls"):
    """Run the full transit search + vetting. Returns a metrics dict.

    Raises on data/download failure; the caller decides how to log that.
    """
    lc = load_lightcurve(target, mission, all_sectors, tesscut)
    flat = lc.flatten(window_length=1001, niters=3)

    period, t0, dur, power, periods_arr, power_arr = detect(flat, method, pmin, pmax)

    # Work in the flattened time series so we can identify individual transits.
    t = flat.time.value
    f = flat.flux.value
    phase = (t - t0 + 0.5 * period) % period - 0.5 * period
    in_tr = np.abs(phase) < (dur / 2)
    oot = ~in_tr

    depth = (1 - np.nanmedian(f[in_tr])) * 1e6 if in_tr.sum() else np.nan
    scatter = np.nanstd(f[oot]) * 1e6
    n_in = max(int(in_tr.sum()), 1)
    snr = depth / (scatter / np.sqrt(n_in)) if scatter > 0 else 0.0

    # Proper odd/even: assign each in-transit point to a transit number.
    tnum = np.round((t - t0) / period).astype(int)
    odd = in_tr & (tnum % 2 == 1)
    even = in_tr & (tnum % 2 == 0)
    depth_odd = (1 - np.nanmedian(f[odd])) * 1e6 if odd.sum() >= 3 else np.nan
    depth_even = (1 - np.nanmedian(f[even])) * 1e6 if even.sum() >= 3 else np.nan

    # Secondary eclipse: look near phase 0.5.
    sec_phase = ((t - t0) % period) - period / 2  # ~0 at phase 0.5
    near_sec = np.abs(sec_phase) < (dur / 2)
    secondary = (1 - np.nanmedian(f[near_sec])) * 1e6 if near_sec.sum() >= 3 else np.nan

    dur_frac = dur / period

    # coordinates for logging
    ra = lc.meta.get("RA_OBJ", lc.meta.get("RA", np.nan))
    dec = lc.meta.get("DEC_OBJ", lc.meta.get("DEC", np.nan))

    m = dict(
        target=target, ra=ra, dec=dec, n_points=len(f),
        period=period, t0=t0, duration_hr=dur * 24, duration_frac=dur_frac,
        depth_ppm=depth, snr=snr, power=power,
        depth_odd_ppm=depth_odd, depth_even_ppm=depth_even, secondary_ppm=secondary,
    )
    m["verdict"], m["reasons"] = classify(m, mission)

    # Deep vetting: only survivors get the extra blend check (a catalog query).
    # A comparably-bright-or-brighter neighbor can fake the transit by bleeding
    # light into TESS's big pixels, so we demote those to BLEND_SUSPECT - not a
    # clean hit, but kept for review rather than silently dropped.
    m["blend_risk"] = ""
    if do_centroid and m["verdict"] in ("CANDIDATE", "STRONG_CANDIDATE"):
        bc = blend_neighbor_check(target, mission)
        if bc is not None:
            m["blend_risk"] = bc["note"] or "clean"
            if bc["flagged"]:
                extra = "possible_blend:" + bc["note"]
                m["reasons"] = (m["reasons"] + ";" + extra) if m["reasons"] else extra
                m["verdict"] = "BLEND_SUSPECT"

    m["method"] = method
    is_survivor = m["verdict"] in ("CANDIDATE", "STRONG_CANDIDATE", "BLEND_SUSPECT")
    if make_plot is True or (make_plot == "survivors" and is_survivor):
        m["plot_path"] = _make_plot(target, lc, flat, periods_arr, power_arr,
                                    period, t0, dur, method)
        flat.to_pandas()[["flux", "flux_err"]].to_csv(
            os.path.join(DATA_DIR, f"{safe_name(target)}_lightcurve.csv"))
    return m


def blend_neighbor_check(target, mission="TESS"):
    """Is there a bright enough star nearby to fake this transit by contamination?

    Uses only the TESS Input Catalog (always available - no pixel file needed).
    A real planet host is usually the brightest star in its patch; a blend victim
    has a comparably-bright-or-brighter neighbor within a few detector pixels whose
    eclipse/variability can bleed into the aperture. The search radius is
    mission-specific (see MISSION_PARAMS): the same "few pixels" is ~42" on TESS
    but only ~12" on Kepler/K2.

    Returns dict(flagged: bool, note: str) or None if it couldn't be checked.
    """
    contam_arcsec = params_for(mission)["contam_arcsec"]
    try:
        from astroquery.mast import Catalogs
        obj = Catalogs.query_object(target, catalog="TIC", radius=0.003)
        obj.sort("Tmag")
        ra, dec, tmag = float(obj[0]["ra"]), float(obj[0]["dec"]), float(obj[0]["Tmag"])
        cat = Catalogs.query_region(f"{ra} {dec}",
                                    radius=contam_arcsec / 3600.0, catalog="TIC")
        dra = (np.asarray(cat["ra"], float) - ra) * np.cos(np.radians(dec)) * 3600.0
        ddec = (np.asarray(cat["dec"], float) - dec) * 3600.0
        sep = np.hypot(dra, ddec)
        tmags = np.asarray(cat["Tmag"], float)
        others = (sep > 1.0) & np.isfinite(tmags) & (tmags < tmag + BRIGHTER_MARGIN_MAG)
        if others.sum():
            j = np.where(others)[0][np.argmin(tmags[others])]
            return {"flagged": True,
                    "note": f"bright_neighbor(Tmag{tmags[j]:.1f}@{sep[j]:.0f}\"_vs_target{tmag:.1f})"}
        return {"flagged": False, "note": ""}
    except Exception:
        return None


def classify(m, mission="TESS"):
    """Return (verdict, reasons) from the metrics. Pure deterministic rules.

    SNR floors are mission-specific (see MISSION_PARAMS); depth / odd-even /
    secondary / duration checks are physical and apply to every mission.
    """
    p = params_for(mission)
    snr_min, snr_strong = p["snr_min"], p["snr_strong"]
    # BLS "power" is on an arbitrary per-star scale, so it is NOT used as a
    # threshold - SNR and depth are the trustworthy significance measures.
    snr, depth = m["snr"], m["depth_ppm"]
    if not np.isfinite(depth) or snr < snr_min or depth < DEPTH_MIN_PPM:
        return "NO_SIGNAL", ""

    reasons = []
    if m["duration_frac"] > DURATION_MAX_FRAC:
        reasons.append("duration_unphysical")
    if depth > DEPTH_MAX_PPM:
        reasons.append("too_deep(likely_stellar)")
    do, de = m["depth_odd_ppm"], m["depth_even_ppm"]
    if np.isfinite(do) and np.isfinite(de) and max(do, de) > 0:
        if abs(do - de) / max(do, de) > ODD_EVEN_MAX_FRAC:
            reasons.append("odd_even_mismatch")
    sec = m["secondary_ppm"]
    if np.isfinite(sec) and sec > SECONDARY_MIN_PPM and sec > SECONDARY_FRAC * depth:
        reasons.append("secondary_eclipse")

    if reasons:
        return "FALSE_POSITIVE", ";".join(reasons)
    if snr >= snr_strong:
        return "STRONG_CANDIDATE", ""
    return "CANDIDATE", ""


def _make_plot(target, lc, flat, periods_arr, power_arr, period, t0, dur, method="bls"):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    folded = flat.fold(period=period, epoch_time=t0)
    fig, ax = plt.subplots(4, 1, figsize=(9, 13))
    lc.scatter(ax=ax[0], s=1); ax[0].set_title(f"1) {target}: raw normalized light curve")
    ax[1].plot(periods_arr, power_arr, lw=0.8)
    ax[1].axvline(period, color="red", ls="--", alpha=.6, label=f"peak = {period:.4f} d")
    ax[1].set_xlabel("Period [days]"); ax[1].set_ylabel(f"{method.upper()} power")
    ax[1].legend(); ax[1].set_title(f"2) {method.upper()} periodogram (power vs trial period)")
    folded.scatter(ax=ax[2], s=1); ax[2].set_xlim(-3 * dur, 3 * dur); ax[2].axvline(0, color="red", ls=":", alpha=.5)
    ax[2].set_title(f"3) Folded on {period:.4f} d - transit close-up")
    folded.scatter(ax=ax[3], s=1); ax[3].axvline(0, color="red", ls=":", alpha=.4)
    ax[3].set_title("4) Full phase (dip near phase 0.5 = binary flag)")
    plt.tight_layout()
    path = os.path.join(PLOTS_DIR, f"{safe_name(target)}_transit.png")
    fig.savefig(path, dpi=110, bbox_inches="tight"); plt.close(fig)
    return path


def _report(m):
    print("=" * 60)
    print(f"RESULT for {m['target']}   (RA={m['ra']}, Dec={m['dec']})")
    print("=" * 60)
    print(f"  Period      : {m['period']:.5f} d      Duration: {m['duration_hr']:.2f} h "
          f"({m['duration_frac']*100:.1f}% of orbit)")
    print(f"  Depth       : {m['depth_ppm']:.0f} ppm   SNR: {m['snr']:.1f}   "
          f"{m.get('method', 'bls').upper()} score: {m['power']:.0f}")
    print(f"  Odd/even    : {m['depth_odd_ppm']:.0f} / {m['depth_even_ppm']:.0f} ppm   "
          f"Secondary@0.5: {m['secondary_ppm']:.0f} ppm")
    if m.get("blend_risk") not in ("", None):
        print(f"  Blend check : {m['blend_risk']}")
    print(f"  VERDICT     : {m['verdict']}" + (f"  ({m['reasons']})" if m["reasons"] else ""))


def main():
    p = argparse.ArgumentParser(description="Search one star's light curve for transiting planets.")
    p.add_argument("target")
    p.add_argument("--mission", default="TESS", choices=["TESS", "Kepler", "K2", ""])
    p.add_argument("--pmin", type=float, default=0.5)
    p.add_argument("--pmax", type=float, default=15.0)
    p.add_argument("--all-sectors", action="store_true")
    p.add_argument("--tesscut", action="store_true")
    p.add_argument("--plot", action="store_true", help="save the 4-panel plot + data csv")
    p.add_argument("--tls", action="store_true",
                   help="use Transit Least Squares (more sensitive, slower) instead of BLS")
    a = p.parse_args()
    method = "tls" if a.tls else "bls"
    m = analyze(a.target, a.mission, a.pmin, a.pmax, a.all_sectors, a.tesscut,
                make_plot=a.plot, method=method)
    _report(m)
    if a.plot:
        print(f"  plot : {m.get('plot_path')}")


if __name__ == "__main__":
    main()
