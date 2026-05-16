---
name: frame-set-verify
description: Use whenever loading a *range* or *list* of Henrietta FITS frames before analysis. Asserts the set is homogeneous in filter, grism, exposure time, mode (CDS/SUTR), ramp depth, and OBSTYPE so reductions don't silently mix dark/flat/object frames or mismatched configurations.
---

# frame-set-verify

It is a recurring bug to grab a range of frames (e.g. 1755..1760) for a
read-noise or flat-field analysis and discover *afterwards* that one frame had
a different filter, the grism wheel moved mid-sequence, or someone aborted a
ramp and the sample count differs.

**Run this check first. Refuse to proceed if it fails.**

## What to verify

Across the set of frames the user named:

| Header key | Why |
|---|---|
| `FILTER` (or `FILTER1`/`FILTER2`) | dark vs flat vs imaging filter mixed = garbage |
| `GRISM` | imaging vs spectral mixed = garbage |
| `SLIT` | spectral slit width / open slot |
| `EXPTIME` | different exposure times = different signal levels |
| `OBSTYPE` (or `PURPOSE`/`IMAGETYP`) | dark vs flat vs object |
| `MODE` (CDS vs SUTR) | sample structure differs |
| ramp depth (count of `hen{NNNN}_{KKK}r.fits` per NNNN) | aborted ramps |

For SUTR data also confirm every ramp has the **same number of samples**
(usually 23) — a short ramp poisons any per-ramp statistic.

## Reference implementation

```python
from pathlib import Path
from collections import Counter
from astropy.io import fits

SCIDATA = Path("/carnegie/scidata/groups/kallisto/Henrietta/data/images/raw/Commissioning")
HEADER_KEYS = ("FILTER", "GRISM", "SLIT", "EXPTIME", "OBSTYPE", "MODE")

def verify_set(frame_ids, expect_samples=None):
    """Raise AssertionError if frames are inhomogeneous. Returns ramp_depth."""
    seen = {k: set() for k in HEADER_KEYS}
    sample_counts = {}
    for fid in frame_ids:
        # use the fitted/summary frame for header inspection if it exists,
        # otherwise the first SUTR sample
        path = SCIDATA / f"hen{fid}.fits"
        if not path.exists():
            path = SCIDATA / f"hen{fid}_001r.fits"
        hdr = fits.getheader(path)
        for k in HEADER_KEYS:
            if k in hdr:
                seen[k].add(hdr[k])
        # count SUTR samples for this ramp
        samples = sorted(SCIDATA.glob(f"hen{fid}_*r.fits"))
        sample_counts[fid] = len(samples)

    for k, vals in seen.items():
        assert len(vals) <= 1, f"frame set is inhomogeneous in {k}: {vals}"

    depth = Counter(sample_counts.values())
    assert len(depth) == 1, f"ramp depth varies across set: {sample_counts}"
    n_samples = next(iter(depth))
    if expect_samples is not None:
        assert n_samples == expect_samples, f"expected {expect_samples} samples/ramp, got {n_samples}"
    return n_samples

# usage
n = verify_set(range(1755, 1761), expect_samples=23)
```

## When the check fails

Print *which key* differs and across which frame IDs. Then **stop and ask the
user** — don't silently drop the odd frame, and don't proceed with a partial
set. The user may have meant to include it, or there may be a real instrument
problem worth investigating.

## Don't skip this for "obviously homogeneous" runs

The commissioning log has many cases where the operator changed something
mid-sequence and the FITS headers disagree with the log. Trust the headers.
