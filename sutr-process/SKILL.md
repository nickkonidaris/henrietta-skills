---
name: sutr-process
description: Use when loading Henrietta sample-up-the-ramp (SUTR) data and computing signal images — CDS pair differences, simple last-minus-first, or OLS slope fits across an N-sample ramp. Covers the hen{NNNN}_{KKK}r.fits naming, reference-pixel layout, and gain conversion to electrons.
---

# sutr-process

Henrietta's H2RG records non-destructive samples up the ramp. Each ramp is
stored as `hen{NNNN}_{KKK}r.fits` where `NNNN` is the ramp/observation ID and
`KKK` is the 1-indexed sample number (`001r` = first sample, typically ≤ 23
samples per ramp). The pre-existing pipeline also writes a fitted slope
image `hen{NNNN}.fits`.

## Conventions (re-state these in every analysis script)

```python
DATA_DIR = Path("/carnegie/scidata/groups/kallisto/Henrietta/data/images/raw/Commissioning")
GAIN_E_PER_DN = 4.0          # header gain is wrong; use this
N_SAMPLES    = 23            # standard ramp depth; verify with frame-set-verify
EXPTIME_S    = 60.0          # standard, but read from header
# H2RG geometry
SHAPE  = (2048, 2048)
REFPIX = 4                   # reference pixel border on each edge (Y fast, X slow)
N_AMPS = 4                   # 4 amplifier columns, split along X (slow read)
```

Don't read these from the headers unless you have verified the value in the
header. Header gain is known wrong.

## Loading a ramp

```python
import numpy as np
from astropy.io import fits

def load_ramp(ramp_id, n_samples=N_SAMPLES):
    """(n_samples, 2048, 2048) float32 array, raw DN with BZERO applied by astropy."""
    out = np.empty((n_samples, 2048, 2048), dtype=np.float32)
    for k in range(n_samples):
        out[k] = fits.getdata(DATA_DIR / f"hen{ramp_id}_{k+1:03d}r.fits").astype(np.float32)
    return out
```

## Signal estimators (pick by context)

### CDS — last minus first
Cheap, used when you only need a single signal estimate per pixel and don't
care about read-noise rejection.

```python
def cds_signal(ramp):  # (n, 2048, 2048) DN
    return ramp[-1] - ramp[0]
```

### UTR slope — OLS fit across the ramp
Used for read-noise and any photometry where you want the lowest-variance
signal estimate. Slope `α` in DN/s, then signal = `α * EXPTIME`.

```python
def utr_slope(ramp, exptime=EXPTIME_S):
    n = ramp.shape[0]
    t = np.arange(n, dtype=np.float32) * (exptime / (n - 1))
    # demean trick: slope = sum(t-tbar) * (y-ybar) / sum(t-tbar)^2
    tc = t - t.mean()
    denom = (tc**2).sum()
    yc = ramp - ramp.mean(axis=0, keepdims=True)
    slope = (tc[:, None, None] * yc).sum(axis=0) / denom
    return slope  # DN/s, shape (2048, 2048)
```

For a one-shot best-estimate of integrated signal in DN, multiply by
`EXPTIME`.

### Pair differences (for read-noise PTC / CDS noise)

```python
def pair_diffs(ramp):  # adjacent pairs
    return np.diff(ramp, axis=0)  # (n-1, H, W); each is √2 × σ_R
```

## Convert to electrons

Always do this at the end, never mid-pipeline, so DN-domain diagnostics stay
in DN:

```python
signal_e = signal_dn * GAIN_E_PER_DN
```

## H2RG geometry pitfalls

- The 4-pixel border (refpix) is **not science** — exclude it from any image
  statistic. See the `h2rg-edge-mask` snippet below.
- The fast-read direction is **Y** (along columns). Slow-read is **X** (along
  rows). When you see periodic patterns along Y, look for clock pickup. When
  they're along X, look for amp-common signal.
- 4 amplifiers split along X (slow-read): columns `[4:516)`, `[516:1028)`,
  `[1028:1540)`, `[1540:2044)` after refpix masking. Per-amp statistics are
  often required.

```python
def h2rg_science_slice():
    """Slice that excludes the 4-pix reference border."""
    return (slice(4, 2044), slice(4, 2044))
```

## Spectral mode reminder

If `GRISM` is open and a `SLIT` is in: spectral direction is along **columns
(Y)**, cross-dispersion along **rows (X)**. Isospatial lines are tilted and
the tilt varies across the slit. Isowavelength lines are curved. Do not
collapse along an axis without re-checking this.

## When in doubt

Sanity-check by computing the CDS image of one ramp, viewing it (ds9 / matplotlib
imshow with reasonable percentile clipping), and confirming it looks like the
expected image type (dark = noise-only, flat = roughly uniform, object = stars).
A wrong-axis or wrong-units bug is usually visible immediately.
