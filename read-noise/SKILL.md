---
name: read-noise
description: Use when measuring or sanity-checking the read noise of the Henrietta H2RG detector — single-sample (σ_R), CDS, or UTR slope-fit (σ_D). Encodes the canonical formulas, the dark-sequence baseline, per-amp split, and the empirical 1755..1760 reference values from the read-noise investigation.
---

# read-noise

The Henrietta H2RG is read-noise-limited in short SUTR exposures. Three
related quantities show up in conversations and they are routinely confused:

| symbol | meaning | typical value (1755..1760, masked) |
|---|---|---|
| σ_R | std of a single non-destructive sample read, in electrons | **17.5 e⁻** |
| σ_CDS | std of (sample_k − sample_j), one pair | σ_R · √2 ≈ 24.7 e⁻ |
| σ_D | std of the OLS slope-fit integrated signal (DN/s × EXPTIME) | **~24 e⁻** baseline; **~16 e⁻** after the L2 (top-2 amp-PCA modes) phase-invariant correction |

The theoretical UTR scaling is

```
σ_D_theory = σ_R · √(12·(N−1) / (N·(N+1)))
```

For N=23: σ_D_theory ≈ 0.69 · σ_R ≈ 12.1 e⁻. **Measured σ_D is well above
theory** because slope drift across samples (per-pixel + clock-locked
alternating signal) adds variance that white-noise theory ignores.
See `/home/npk/Henrietta/read-noise/SYNTHESIS.md` for the full noise budget.

## Reference dataset

Sequence **1755..1760** on 28 Apr 2026, 6 ramps × 23 samples × 60 s,
`GRISM=closed`. This is the clean-dark baseline. Any new read-noise quote
should be compared against the table above on this same sequence.

Persistence-loaded sequence (use only when you specifically want a stressed
test): **2093..2098**.

## Standard recipe

```python
import numpy as np
from astropy.io import fits

# 0. verify homogeneity (see frame-set-verify skill) — do not skip
# 1. load ramps
ramps = np.stack([load_ramp(rid) for rid in range(1755, 1761)])  # (R, N, H, W)
# 2. mask refpix and known bad pixels
SCI = (slice(4, 2044), slice(4, 2044))
ramps = ramps[:, :, SCI[0], SCI[1]]
bpm = np.load("/path/to/bpm.npy")[SCI[0], SCI[1]]  # True = bad

# 3. σ_R per pixel: std of adjacent-sample differences / √2
diffs = np.diff(ramps, axis=1)                   # (R, N-1, H, W)
sigma_R_pix = diffs.std(axis=(0, 1)) / np.sqrt(2)  # DN per pixel
sigma_R_pix_e = sigma_R_pix * 4.0                 # electrons
# masked median is the headline number
print("σ_R median:", np.median(sigma_R_pix_e[~bpm]))

# 4. σ_D from OLS slope fits
slopes = np.stack([utr_slope(r) for r in ramps]) * 60.0  # DN integrated over EXPTIME
sigma_D_pix = slopes.std(axis=0) * 4.0                    # electrons
print("σ_D median:", np.median(sigma_D_pix[~bpm]))
```

## Per-amp split

H2RG has 4 amplifiers along the slow-read (X) direction. Read noise often
differs amp-to-amp; quote a 4-value vector, not a single global number, when
characterizing the detector:

```python
amp_edges = [4, 516, 1028, 1540, 2044]  # 4-pix refpix already excluded by slicing
for i in range(4):
    sl = (slice(0, None), slice(amp_edges[i] - 4, amp_edges[i+1] - 4))
    print(f"amp {i}: σ_R = {np.median(sigma_R_pix_e[sl][~bpm[sl]]):.2f} e-")
```

## Phase-invariant L2 correction (the σ_D reducer)

When σ_D is the metric of interest, apply the two-template amp-PCA correction
from `exp_036`. Templates live at:

```
/home/npk/Henrietta/read-noise/experiments/exp_036_phase_invariant/results.npz
```

Cross-sample test: this reduced σ_D on `2093..2098` from 339 → 36.5 e⁻ (81%
reduction) and on the in-sample baseline `1755..1760` from 24 → 16 e⁻. It is
phase-invariant — does not require tracking the alternating signal's phase
between sequences.

## Sanity floors and ceilings

- Single sample σ_R < 12 e⁻: suspicious, double-check refpix and bpm masking.
- σ_D / σ_D_theory >> 3: something is contaminating the ramps — persistence,
  light leak, a target in the field. Inspect a single ramp's median image.
- σ_R amp imbalance > 30%: investigate amplifier or bias supply.
