---
name: load-frames
description: Use whenever the user asks to open, read, or analyze a set of Henrietta FITS frames identified by frame-ID range or list. Identifies the on-disk paths, confirms every file exists and the set is homogeneous in filter/grism/exptime/ramp-depth, then loads the data into numpy arrays and returns the matching headers. Handles both fitted slope frames (hen{NNNN}.fits) and SUTR ramp samples (hen{NNNN}_{KKK}r.fits).
---

# load-frames

The single canonical way to go from "frames 1755..1760" (or any list of frame
IDs) to in-memory numpy arrays + headers. Always use this skill before
running statistics or visualization on a set of Henrietta frames — it
enforces the homogeneity check from `CLAUDE.md` ("when you open a set of
frames in a range, please double check that they are homogenous") so you
don't quietly mix a dark into a flat sequence.

## What it does, in three steps

1. **Identify** — Expand a frame-ID range or list into concrete file paths
   under the resolved data directory (see *Setup* below — defaults to the
   `~/Henrietta/data` symlink or `$HENRIETTA_DATA_DIR`). Fitted-frame paths
   look like `hen{NNNN}.fits`; ramp-sample paths look like
   `hen{NNNN}_{KKK}r.fits` (1-indexed `KKK`, three digits).
2. **Confirm** — Verify every expected file is on disk, then verify the set
   is homogeneous across the headers that matter for analysis: `FILTER`,
   `GRISM`, `SLIT`, `EXPTIME`, `OBSTYPE`, `MODE`, and ramp depth (for SUTR).
   Raise immediately on any mismatch and tell the caller *which* key differs
   and across which frame IDs.
3. **Read** — Load the data as a numpy array and return the matching headers
   as a numpy object array of `astropy.io.fits.Header` instances. Caller
   gets back a single `(data, headers)` tuple.

## Setup — where does the data live?

The scidata mount path varies per machine (`/carnegie/scidata/...` on Carnegie
hosts, somewhere else on your laptop or a collaborator's box). The skill
resolves the data directory in this order, falling through until one works:

1. Explicit `data_dir=Path(...)` argument to `load_frames(...)`.
2. `$HENRIETTA_DATA_DIR` environment variable.
3. The team-convention symlink `~/Henrietta/data`.

**Pick one and do it once per machine.** Recommended: the symlink, because
analysis scripts then need no per-machine configuration.

```bash
# On a Carnegie host where scidata is mounted at the usual path:
ln -s /carnegie/scidata/groups/kallisto/Henrietta/data/images/raw/Commissioning \
      ~/Henrietta/data

# On a different host where it's mounted elsewhere:
ln -s /mnt/whatever/Commissioning ~/Henrietta/data

# Or skip the symlink and use the env var (e.g. in ~/.bashrc):
export HENRIETTA_DATA_DIR=/mnt/whatever/Commissioning
```

If none of those resolve, `load_frames` raises with a message telling you the
three options. The skill never writes into the resolved directory — it only
reads (see project-wide `/scidata` read-only convention).

## Reference implementation

Drop this into the experiment's helper module, or run it ad-hoc. It is small
on purpose; don't extend it without a reason.

```python
"""load_frames — identify, confirm, read."""
import os
from pathlib import Path
from collections.abc import Iterable
import numpy as np
from astropy.io import fits

# Header keys that must agree across a homogeneous set. Missing keys in some
# headers are tolerated; conflicting values are not.
HOMOGENEITY_KEYS = ("FILTER", "GRISM", "SLIT", "EXPTIME", "OBSTYPE", "MODE")

# Where to find the raw Commissioning frames. Resolution order:
#   1. data_dir= argument to load_frames()
#   2. $HENRIETTA_DATA_DIR environment variable
#   3. ~/Henrietta/data  (team convention: symlink to wherever scidata is
#      mounted on this machine)
# Mount points vary by machine — never hardcode /carnegie/scidata/... in
# analysis scripts; rely on this resolver instead.
DEFAULT_SYMLINK = Path.home() / "Henrietta" / "data"


def _resolve_data_dir(data_dir=None):
    if data_dir is not None:
        return Path(data_dir)
    env = os.environ.get("HENRIETTA_DATA_DIR")
    if env:
        return Path(env)
    if DEFAULT_SYMLINK.exists():
        return DEFAULT_SYMLINK
    raise FileNotFoundError(
        "Cannot locate Henrietta data. Do one of:\n"
        f"  - symlink {DEFAULT_SYMLINK} -> wherever scidata is mounted (e.g. "
        "/carnegie/scidata/groups/kallisto/Henrietta/data/images/raw/Commissioning)\n"
        "  - export HENRIETTA_DATA_DIR=/path/to/Commissioning\n"
        "  - pass data_dir=Path(...) explicitly to load_frames()"
    )


def _expand_ids(ids):
    """Accept a range, list, or single int. Returns a sorted list of ints."""
    if isinstance(ids, int):
        return [ids]
    if isinstance(ids, range):
        return list(ids)
    if isinstance(ids, Iterable):
        return sorted(int(i) for i in ids)
    raise TypeError(f"frame_ids must be int, range, or iterable of ints; got {type(ids).__name__}")


def _paths_for(frame_id, kind, data_dir):
    if kind == "fitted":
        return [data_dir / f"hen{frame_id}.fits"]
    if kind == "ramp":
        paths = sorted(data_dir.glob(f"hen{frame_id}_*r.fits"))
        if not paths:
            raise FileNotFoundError(f"no ramp samples found for hen{frame_id}_*r.fits in {data_dir}")
        return paths
    raise ValueError(f"kind must be 'fitted' or 'ramp'; got {kind!r}")


def _check_exist(paths):
    missing = [p for p in paths if not p.exists()]
    if missing:
        sample = ", ".join(str(p.name) for p in missing[:5])
        more = f" (+{len(missing) - 5} more)" if len(missing) > 5 else ""
        raise FileNotFoundError(f"{len(missing)} expected file(s) missing: {sample}{more}")


def _check_homogeneous(headers_flat, frame_ids_flat):
    """headers_flat and frame_ids_flat are parallel 1-D sequences."""
    by_key = {k: {} for k in HOMOGENEITY_KEYS}
    for hdr, fid in zip(headers_flat, frame_ids_flat):
        for k in HOMOGENEITY_KEYS:
            if k in hdr:
                by_key[k].setdefault(hdr[k], []).append(fid)
    for k, vals in by_key.items():
        if len(vals) > 1:
            summary = "; ".join(f"{v!r}: {sorted(set(ids))}" for v, ids in vals.items())
            raise ValueError(f"frame set is inhomogeneous in {k}: {summary}")


def load_frames(frame_ids, kind="fitted", data_dir=None):
    """Identify, confirm, and read a set of Henrietta frames.

    Parameters
    ----------
    frame_ids : int | range | iterable of int
        Frame IDs to load. e.g. ``range(1755, 1761)`` or ``[1755, 1760]``.
    kind : {"fitted", "ramp"}, default ``"fitted"``
        ``"fitted"`` -> one ``hen{NNNN}.fits`` per ID, data shape ``(N, H, W)``.
        ``"ramp"``   -> every ``hen{NNNN}_*r.fits`` sample on disk, data shape
        ``(N, n_samples, H, W)``. ``n_samples`` is discovered from disk and
        must agree across every frame ID in the set.
    data_dir : Path, optional
        Where the raw frames live. If omitted, resolves from
        ``$HENRIETTA_DATA_DIR`` or the ``~/Henrietta/data`` symlink. See
        ``_resolve_data_dir`` and the README for setup.

    Returns
    -------
    data : np.ndarray
        Float32 array of DN. Shape ``(N, H, W)`` (fitted) or
        ``(N, n_samples, H, W)`` (ramp).
    headers : np.ndarray
        Object array of ``astropy.io.fits.Header``, same leading shape as
        ``data`` minus the spatial axes. For ramps: shape ``(N, n_samples)``.

    Raises
    ------
    FileNotFoundError : any expected path is missing.
    ValueError       : the set is inhomogeneous in any HOMOGENEITY_KEYS field,
                       or for ramps the on-disk sample count disagrees across IDs.
    """
    ids = _expand_ids(frame_ids)
    data_dir = _resolve_data_dir(data_dir)

    # 1. identify (for kind="ramp", this also discovers n_samples per ID)
    paths_per_id = [_paths_for(fid, kind, data_dir) for fid in ids]
    flat_paths = [p for ps in paths_per_id for p in ps]

    # 2. confirm — files exist
    _check_exist(flat_paths)

    # 2. confirm — ramp depth agrees across IDs
    if kind == "ramp":
        counts = {fid: len(ps) for fid, ps in zip(ids, paths_per_id)}
        unique = set(counts.values())
        if len(unique) > 1:
            by_count = {}
            for fid, c in counts.items():
                by_count.setdefault(c, []).append(fid)
            summary = "; ".join(f"{c} samples: {fids}" for c, fids in sorted(by_count.items()))
            raise ValueError(f"ramp depth disagrees across set: {summary}")
        n_samples = unique.pop()

    # 2. confirm — read headers, check homogeneity
    headers_flat = [fits.getheader(p) for p in flat_paths]
    frame_ids_flat = [fid for fid, ps in zip(ids, paths_per_id) for _ in ps]
    _check_homogeneous(headers_flat, frame_ids_flat)

    # 3. read
    arrays = [fits.getdata(p).astype(np.float32) for p in flat_paths]
    if kind == "fitted":
        data = np.stack(arrays, axis=0)                  # (N, H, W)
        headers = np.array(headers_flat, dtype=object)   # (N,)
    else:
        n = len(ids)
        h, w = arrays[0].shape
        data = np.stack(arrays, axis=0).reshape(n, n_samples, h, w)
        headers = np.array(headers_flat, dtype=object).reshape(n, n_samples)

    return data, headers
```

## Usage

```python
# Six fitted slope frames, returns (6, 2048, 2048) + (6,) header array
data, hdr = load_frames(range(1755, 1761))

# Same six ramps as SUTR, returns (6, 23, 2048, 2048) + (6, 23) header array
data, hdr = load_frames(range(1755, 1761), kind="ramp", n_samples=23)

# A non-contiguous list works too
data, hdr = load_frames([1755, 1757, 1759])

# Inspect headers like a numpy array
print(hdr.shape, hdr[0]["FILTER"], hdr[0]["EXPTIME"])
```

## When the confirm step fails

`load_frames` raises with the specific key and the frame IDs grouped by
value, e.g.:

```
ValueError: frame set is inhomogeneous in GRISM: 'closed': [1755, 1756, 1757]; 'R-J': [1758, 1759, 1760]
```

That message is what you want — surface it verbatim to the user. **Do not**
silently drop the odd frame and re-try; the user needs to decide whether to
split the set, fix the input, or investigate the instrument log.

## What this skill intentionally does NOT do

- Convert to electrons (gain is ≈ 4.0 e⁻/DN — apply at the analysis layer).
- Mask reference pixels or bad pixels (do that at the analysis layer too).
- Compute slopes, CDS, or any signal estimate.
- Write anything to disk (and never anywhere under `/scidata`).

Keep this skill as the I/O front door. Build analysis on top of it.
