# exohunter

exohunter is a Python tool that searches NASA's Kepler, K2, and TESS light curves
for the dip a transiting planet makes, vets the result against common
false-positive signatures, and assigns a rule-based verdict. It runs locally
against the public NASA archive, with no API keys, paid services, or AI. Same
star in, same verdict out.

It is a re-detection and triage tool, not a discovery tool. It re-finds known
planets around bright, well-observed stars; it does not find new ones. The bright
stars it works on have already been searched by the mission pipelines (NASA's
SPOC and MIT's QLP), and on the faint stars those pipelines skip, exohunter does
not work (numbers below). So its uses are learning how transit detection
works, independently checking a specific star or candidate, and triaging bright
targets, not discovery.

Measured performance, benchmarked against NASA's own labelled catalogues, per
mission. These are the actual numbers, not estimates:

- TESS, single sector: re-detects 72% of confirmed planets; rejects 41% of known
  false positives.
- Kepler, all quarters stacked: re-detects 65% of confirmed planets; rejects 58%
  of known false positives.
- Faint TESS stars (magnitude 13.5 to 15): re-detects 0% of confirmed planets.
  Faint stars carry too little signal for its simple detrending and
  box-least-squares search. That is also why it cannot usefully search the faint,
  un-analysed population where any new planets would most likely be hiding.

The false-positive rejection rates are low because the test set (NASA's Objects of
Interest) is deliberately hard: those are borderline, transit-shaped signals that
required telescope follow-up to rule out, follow-up exohunter cannot do. On
ordinary stars the false-alarm rate is lower.

How it works: it downloads a light curve, removes slow trends, runs a Box Least
Squares search (or the more sensitive Transit Least Squares as an option), then
checks signal-to-noise, transit depth and duration, odd/even transit consistency,
secondary-eclipse depth, and whether a bright neighbouring star could be causing a
blend. It can also cross-check a candidate against known catalogues and test
whether a transit recurs across multiple sectors.

What it is not: it is not novel science, it does not confirm planets (that needs
pixel-level analysis and telescope follow-up), and it does not work on faint or
previously un-searched stars. It is a working, benchmarked implementation
of standard transit-detection methods, with its limits measured and stated.

## Install

```bash
python3 -m pip install lightkurve astroquery numpy matplotlib certifi
# optional, only for the more sensitive TLS search and thorough benchmark mode:
python3 -m pip install transitleastsquares
```

## Running it

### One star

```bash
python3 exohunter.py "WASP-18" --plot
python3 exohunter.py "Kepler-10" --mission Kepler --all-sectors
```

Options:

- `--mission {TESS,Kepler,K2}` (default TESS)
- `--pmin` / `--pmax` period search range in days (defaults 0.5 / 15)
- `--all-sectors` stitch all sectors/quarters before searching (more sensitive, slower)
- `--tesscut` use a full-frame-image cutout (for stars with no 2-minute light curve)
- `--tls` use Transit Least Squares instead of BLS (more sensitive, much slower)
- `--plot` save the 4-panel figure and the light-curve CSV

It prints period, depth, signal-to-noise, the vetting results, and a verdict:
`STRONG_CANDIDATE`, `CANDIDATE`, `BLEND_SUSPECT`, `FALSE_POSITIVE`, or `NO_SIGNAL`.

### A list of stars, or a sky region

```bash
python3 scan.py --targets stars.txt --pmax 30
python3 scan.py --region "90 -66 2" --limit 200 --pmax 30
python3 scan.py --demo
```

Input: `--targets FILE` (one name per line), `--region "RA DEC RADIUS"` (degrees),
or `--demo`.

Region options: `--maxmag` (bright cutoff, default 12), `--minmag` (faint cutoff;
set e.g. 13.5 to target the faint stars the pipelines skip), `--limit` (max stars,
brightest first, default 300).

Search options: `--mission`, `--pmin`, `--pmax`, `--tls`, `--tesscut`,
`--all-sectors` (same meaning as above).

`--vet` cross-checks each survivor against the catalogues and checks multi-sector
persistence during the run, so `candidates.csv` holds only genuinely-new survivors
(adds a lookup and a re-download per survivor).

`--beep` / `--no-beep` audible alert on a `STRONG_CANDIDATE` hit (on by default).

`--rescan` ignore `results.csv` and redo everything.

Outputs: `results.csv` (every star), `candidates.csv` and `hits.log` (survivors),
`plots/` (a figure per survivor). Scans are resumable: rerun the same command and
it continues.

### Vetting a candidate

```bash
python3 crosscheck.py <ra_deg> <dec_deg> [period_days]   # already catalogued?
python3 persistence.py "TARGET" --mission TESS           # recurs across sectors?
python3 single_transit.py "TARGET" --mission TESS        # long-period/single dip
```

`python3 single_transit.py --selftest` runs a synthetic injection-recovery check.

### Reproducing the benchmarks

```bash
python3 make_toi_list.py                                 # TESS labelled hosts
python3 make_koi_list.py                                 # Kepler labelled hosts
python3 benchmark.py --targets koi_targets.txt --meta koi_meta.csv \
        --mission Kepler --mode stacked --sample 40
python3 recalibrate.py                                   # re-derive verdicts on a finished run
```

`benchmark.py` modes: `fast` (single sector/quarter, BLS), `stacked` (all
sectors/quarters, BLS), `thorough` (all sectors, TLS; minutes per star, small
samples only). `--tesscut` for faint FFI-only stars. `--sample N` caps stars per
disposition.

## How it was benchmarked

The accuracy figures are measured against NASA's own labelled catalogues, not
estimated.

**Ground truth.** For TESS, the TESS Objects of Interest (TOI) table; for Kepler,
the KOI cumulative table. Each host star carries a NASA disposition: confirmed or
known planet, false positive, or unresolved candidate. These lists are built by
`make_toi_list.py` and `make_koi_list.py` from the NASA Exoplanet Archive.

**The test.** `benchmark.py` runs exohunter on each labelled star with no
knowledge of its label, then compares the verdict to NASA's:

- Detection rate = the fraction of NASA's confirmed/known planets that exohunter
  flags as `CANDIDATE` or `STRONG_CANDIDATE`.
- False-positive rejection = the fraction of NASA's known false positives that
  exohunter does not flag as a candidate.

**Modes.** `fast` uses a single sector (TESS) or quarter (Kepler); `stacked`
stitches all sectors/quarters. A TLS-on-all-sectors mode exists but runs minutes
to tens of minutes per star, so it was only spot-checked, not run at scale.

**Measured results:**

| Data set | Mode | Confirmed detected | False positives rejected | Sample (conf / FP) |
|---|---|---|---|---|
| TESS TOIs (bright) | single sector, BLS | 72% | 41% | ~1,270 / ~1,375 |
| Kepler KOIs | single quarter, BLS | 44% | 65% | 150 / 150 |
| Kepler KOIs | all quarters, BLS | 65% | 58% | 40 / 40 |
| TESS faint (T 13.5 to 15) | FFI cutout, single sector | 0% | ~100% (note) | 40 / 39 |

Note: the faint false-positive rejection is an artefact, not a strength. exohunter
detects almost nothing on faint stars, so it trivially "rejects" false positives
by returning `NO_SIGNAL`.

**Caveats:**

- The TOI and KOI sets are adversarial. They are pre-selected transit-shaped
  signals, many of which required telescope follow-up to classify. The
  false-positive rejection rates are therefore a worst case; on ordinary field
  stars the false-alarm rate is lower.
- Sample sizes vary. The TESS bright number comes from the full ~7,800-star TOI
  scan (~1,270 confirmed planets), so it is robust. The Kepler-stacked and
  faint-TESS numbers used about 40 confirmed planets each, so they carry roughly a
  plus-or-minus 15% margin.
- Thresholds are calibrated per mission (the blend-contamination radius is 42
  arcsec for TESS's large pixels, 12 arcsec for Kepler; the detection floor is
  signal-to-noise 8). `recalibrate.py` re-derives every verdict from a finished
  run under the current thresholds.

The headline finding: exohunter re-detects most confirmed planets on bright,
already-searched stars and detects none on faint stars, which is why it is a
re-detection and triage tool rather than a discovery tool.

## Files

- `exohunter.py` — single-star search engine (importable: `analyze`, `classify`)
- `scan.py` — batch scanner over a list or a sky region, with logging and alerts
- `crosscheck.py` — is a star/period already in a NASA catalogue?
- `persistence.py` — does a transit recur consistently across sectors?
- `single_transit.py` — long-period / single-transit dip search
- `make_toi_list.py`, `make_koi_list.py` — build labelled target lists
- `benchmark.py` — measure detection and false-positive rates against the labels
- `recalibrate.py` — re-derive verdicts from a finished run under current thresholds
