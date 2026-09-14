#!/usr/bin/env python3
"""crosscheck.py - is this star (and period) already in a NASA catalogue?

The single most important filter for turning exohunter from "re-detects known
things" into "flags genuinely new candidates". Given a survivor's sky position
(and optionally its detected period), it queries the NASA Exoplanet Archive for:

  - confirmed planets      (pscomppars)
  - TESS Objects of Interest (toi)
  - Kepler Objects of Interest (cumulative / KOI)

within a few arcsec of the position. A positional match means the star already
has a catalogued planet or candidate, so the "discovery" is a re-detection and
should be dropped from a new-candidate list.

Pure standard library + NASA's public TAP service. No API key, no tokens.

Usage (standalone test):
    python3 crosscheck.py <ra_deg> <dec_deg> [period_days]
    python3 crosscheck.py 330.7948 18.8844 3.5247     # HD 209458 b -> KNOWN
"""
import csv
import io
import math
import socket
import ssl
import sys
import urllib.parse
import urllib.request

try:  # macOS system Python often lacks root certs; use certifi's bundle
    import certifi
    SSL_CTX = ssl.create_default_context(cafile=certifi.where())
except Exception:
    SSL_CTX = ssl.create_default_context()

TAP = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync"

# (table, ra col, dec col, period col, label expression) for each catalogue
CATALOGS = [
    ("pscomppars", "ra", "dec", "pl_orbper", "pl_name", "confirmed planet"),
    ("toi",        "ra", "dec", "pl_orbper", "toi",     "TOI"),
    ("cumulative", "ra", "dec", "koi_period", "kepoi_name", "KOI"),
]


def _sep_arcsec(ra1, dec1, ra2, dec2):
    """Angular separation in arcsec (small-angle, fine for a few-arcsec match)."""
    dra = (ra1 - ra2) * math.cos(math.radians((dec1 + dec2) / 2)) * 3600.0
    ddec = (dec1 - dec2) * 3600.0
    return math.hypot(dra, ddec)


def _tap_query(adql, timeout=60):
    url = TAP + "?" + urllib.parse.urlencode({"query": adql, "format": "csv"})
    old = socket.getdefaulttimeout()
    socket.setdefaulttimeout(timeout)
    try:
        with urllib.request.urlopen(url, context=SSL_CTX) as resp:
            return list(csv.DictReader(io.StringIO(resp.read().decode("utf-8"))))
    finally:
        socket.setdefaulttimeout(old)


def _period_matches(p, p_known, tol=0.02):
    """True if p equals p_known, its double, or its half (BLS harmonics), within tol."""
    if p is None or p_known in (None, "", 0):
        return None
    try:
        pk = float(p_known)
    except (TypeError, ValueError):
        return None
    if pk <= 0:
        return None
    for factor in (1.0, 2.0, 0.5):
        if abs(p - pk * factor) / (pk * factor) < tol:
            return True
    return False


def known_object_check(ra, dec, period=None, match_arcsec=6.0):
    """Return dict(known: bool, matches: list[str], note: str) or note the error.

    A match within `match_arcsec` of (ra, dec) in any catalogue counts as known.
    If `period` is given, each match is annotated with whether the period agrees.
    """
    ddeg = match_arcsec / 3600.0
    # RA padding widens with declination; guard near the poles.
    cosd = math.cos(math.radians(min(abs(dec), 89.5)))
    ra_pad = ddeg / cosd if cosd > 1e-6 else 180.0
    ra_lo, ra_hi = ra - ra_pad, ra + ra_pad
    dec_lo, dec_hi = dec - ddeg, dec + ddeg
    wraps = ra_lo < 0 or ra_hi > 360  # RA=0/360 boundary (rare) -> skip RA filter

    matches = []
    for table, racol, deccol, pcol, labelcol, kind in CATALOGS:
        cols = f"{labelcol},{racol},{deccol},{pcol}"
        where = f"{deccol} between {dec_lo} and {dec_hi}"
        if not wraps:
            where += f" and {racol} between {ra_lo} and {ra_hi}"
        adql = f"select {cols} from {table} where {where}"
        try:
            rows = _tap_query(adql)
        except Exception as e:
            return {"known": None, "matches": [],
                    "note": f"lookup failed on {table}: {type(e).__name__}: {str(e)[:80]}"}
        for row in rows:
            try:
                rra, rdec = float(row[racol]), float(row[deccol])
            except (TypeError, ValueError):
                continue
            sep = _sep_arcsec(ra, dec, rra, rdec)
            if sep <= match_arcsec:
                label = str(row.get(labelcol, "?")).strip()
                pm = _period_matches(period, row.get(pcol))
                ptag = "" if pm is None else (" P~match" if pm else " P-differs")
                matches.append(f"{kind} {label} (sep {sep:.1f}\"{ptag})")

    if matches:
        return {"known": True, "matches": matches,
                "note": "already catalogued: " + "; ".join(matches)}
    return {"known": False, "matches": [], "note": "no catalogue match (new)"}


def main():
    if len(sys.argv) < 3:
        sys.exit("usage: python3 crosscheck.py <ra_deg> <dec_deg> [period_days]")
    ra, dec = float(sys.argv[1]), float(sys.argv[2])
    period = float(sys.argv[3]) if len(sys.argv) > 3 else None
    r = known_object_check(ra, dec, period)
    print(f"RA={ra} Dec={dec} P={period}")
    print(f"  known: {r['known']}")
    print(f"  {r['note']}")


if __name__ == "__main__":
    main()
