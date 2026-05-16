---
name: load-frames
description: Use whenever the user asks to open, read, or analyze a set of Henrietta FITS frames — either by explicit frame-ID range/list, or by target/purpose looked up in the commissioning observing-log CSV. Provides two functions: find_frames_for(target=..., purpose=...) returns a sorted list of frame IDs from the log; load_frames(ids, kind=...) confirms the on-disk set is homogeneous in filter/grism/exptime/ramp-depth and returns numpy arrays + matching headers. Handles both fitted slope frames (hen{NNNN}.fits) and SUTR ramp samples (hen{NNNN}_{KKK}r.fits).
---

# load-frames

The canonical way to go from "frames 1755..1760" — or "every HD120270 frame
in the log" — to in-memory numpy arrays + headers. Use this skill before
running statistics or visualization on a set of Henrietta frames; it enforces
the homogeneity check from `CLAUDE.md` ("when you open a set of frames in a
range, please double check that they are homogenous") so you don't quietly
mix a dark into a flat sequence.

The skill exposes two functions, kept separate on purpose:

- **`find_frames_for(target=..., purpose=...)`** — looks up frame IDs in the
  observing-log CSV. Pure metadata lookup; touches no FITS files.
- **`load_frames(frame_ids, kind=...)`** — confirms the on-disk set is
  homogeneous and returns `(data, headers)` as numpy arrays.

Pipe one into the other:

```python
ids = find_frames_for(target="HD120270")
data, hdr = load_frames(ids)
```

## What `load_frames` does, in three steps

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

## Setup — where do the data and the log live?

Two paths vary per machine. Each resolves in the same order: explicit
argument → env var → team-convention symlink.

### Raw frames

1. `data_dir=Path(...)` argument to `load_frames(...)`.
2. `$HENRIETTA_DATA_DIR` environment variable.
3. The symlink `~/Henrietta/data`.

```bash
ln -s /carnegie/scidata/groups/kallisto/Henrietta/data/images/raw/Commissioning \
      ~/Henrietta/data
# or: export HENRIETTA_DATA_DIR=/mnt/.../Commissioning
```

### Observing log CSV

1. `log_csv=Path(...)` argument to `find_frames_for(...)`.
2. `$HENRIETTA_LOG_CSV` environment variable.
3. The symlink `~/Henrietta/observing-log.csv`.

```bash
ln -s "/path/to/Commissioning Run (2026-04-14 to 2026-05-04) Observing Log - Sheet1.csv" \
      ~/Henrietta/observing-log.csv
# or: export HENRIETTA_LOG_CSV=/path/to/that.csv
```

If a resolver finds nothing, the call raises with a message listing the three
options. The skill never writes into either resolved location — read-only,
per the project-wide `/scidata` rule.

## Reference implementation

Drop this into the experiment's helper module, or run it ad-hoc. It is small
on purpose; don't extend it without a reason.

```python
"""load_frames — find_frames_for + load_frames."""
import os
import re
from pathlib import Path
from collections.abc import Iterable
import numpy as np
from astropy.io import fits

# Header keys that must agree across a homogeneous set. Missing keys in some
# headers are tolerated; conflicting values are not.
HOMOGENEITY_KEYS = ("FILTER", "GRISM", "SLIT", "EXPTIME", "OBSTYPE", "MODE")

# Resolution order for both paths: explicit arg → env var → team symlink.
# Mount points vary by machine — never hardcode /carnegie/scidata/... in
# analysis scripts; rely on these resolvers instead.
DATA_SYMLINK = Path.home() / "Henrietta" / "data"
LOG_SYMLINK  = Path.home() / "Henrietta" / "observing-log.csv"


def _resolve_path(explicit, env_var, default_symlink, kind_label):
    if explicit is not None:
        return Path(explicit)
    env = os.environ.get(env_var)
    if env:
        return Path(env)
    if default_symlink.exists():
        return default_symlink
    raise FileNotFoundError(
        f"Cannot locate Henrietta {kind_label}. Do one of:\n"
        f"  - symlink {default_symlink} -> wherever it lives on this host\n"
        f"  - export {env_var}=/path/to/{default_symlink.name}\n"
        f"  - pass the path explicitly as an argument"
    )


def _resolve_data_dir(data_dir=None):
    return _resolve_path(data_dir, "HENRIETTA_DATA_DIR", DATA_SYMLINK, "data directory")


def _resolve_log_csv(log_csv=None):
    return _resolve_path(log_csv, "HENRIETTA_LOG_CSV", LOG_SYMLINK, "observing-log CSV")


def _expand_ids(ids):
    """Accept a range, list, or single int. Returns a sorted list of ints."""
    if isinstance(ids, int):
        return [ids]
    if isinstance(ids, range):
        return list(ids)
    if isinstance(ids, Iterable):
        return sorted(int(i) for i in ids)
    raise TypeError(f"frame_ids must be int, range, or iterable of ints; got {type(ids).__name__}")


# Pattern for a "File #" cell in the observing log. Each cell is either a
# single int ("123") or an inclusive range separated by one of: `-`, `..`,
# `...` (operators have used all three over the run). Whitespace tolerated.
_RANGE_RE = re.compile(r"^\s*(\d+)\s*(?:-|\.{2,3})\s*(\d+)\s*$")
_INT_RE   = re.compile(r"^\s*(\d+)\s*$")


def _expand_csv_cell(cell):
    """Return the list of frame IDs encoded in one File # cell, or []."""
    if cell is None:
        return []
    s = str(cell)
    m = _RANGE_RE.match(s)
    if m:
        a, b = int(m.group(1)), int(m.group(2))
        lo, hi = (a, b) if a <= b else (b, a)
        return list(range(lo, hi + 1))
    m = _INT_RE.match(s)
    if m:
        return [int(m.group(1))]
    return []


def find_frames_for(target=None, purpose=None, log_csv=None):
    """Look up frame IDs in the commissioning observing-log CSV.

    Pure metadata lookup — does not open any FITS files. Matching is
    case-insensitive substring against the named columns; pass both
    arguments to require both to match (AND, not OR).

    Parameters
    ----------
    target : str, optional
        Substring to match against the ``Target (RA/Dec)`` column.
    purpose : str, optional
        Substring to match against the ``Purpose`` column.
    log_csv : Path, optional
        Path to the observing-log CSV. If omitted, resolves via
        ``$HENRIETTA_LOG_CSV`` or the ``~/Henrietta/observing-log.csv``
        symlink.

    Returns
    -------
    list[int]
        Sorted unique frame IDs from the matching rows' ``File #`` cells.

    Notes
    -----
    The CSV log captures *intent* (what the operator typed). Headers are
    authoritative for what was *actually* recorded. After loading the
    matching frames you may want to cross-check ``hdr['OBJECT']`` (or the
    equivalent) against ``target`` and surface disagreements to the user.
    """
    if target is None and purpose is None:
        raise ValueError("provide target=, purpose=, or both")

    # Lazy import — pandas isn't needed unless you actually call this.
    import pandas as pd

    path = _resolve_log_csv(log_csv)
    df = pd.read_csv(path)
    df.columns = [c.strip() for c in df.columns]
    df["Target (RA/Dec)"] = df["Target (RA/Dec)"].astype(str)
    df["Purpose"] = df["Purpose"].astype(str)

    mask = np.ones(len(df), dtype=bool)
    if target is not None:
        mask &= df["Target (RA/Dec)"].str.contains(target, case=False, na=False, regex=False)
    if purpose is not None:
        mask &= df["Purpose"].str.contains(purpose, case=False, na=False, regex=False)

    ids = set()
    for cell in df.loc[mask, "File #"]:
        ids.update(_expand_csv_cell(cell))
    return sorted(ids)


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
        # NB: np.array(list_of_Header, dtype=object) would descend into each
        # Header (it's iterable) and produce a 2-D array. Allocate-then-assign
        # keeps the Header objects atomic.
        headers = np.empty(len(headers_flat), dtype=object)
        headers[:] = headers_flat                        # (N,)
    else:
        n = len(ids)
        h, w = arrays[0].shape
        data = np.stack(arrays, axis=0).reshape(n, n_samples, h, w)
        headers = np.empty(n * n_samples, dtype=object)
        headers[:] = headers_flat
        headers = headers.reshape(n, n_samples)

    return data, headers
```

## Usage

```python
# Explicit frame range
data, hdr = load_frames(range(1755, 1759))           # (4, 2048, 2048) + (4,)

# Same range as SUTR samples — n_samples is discovered from disk
data, hdr = load_frames(range(1755, 1759), kind="ramp")
# data.shape -> (4, 23, 2048, 2048); hdr.shape -> (4, 23)

# Look up by target via the CSV log
ids = find_frames_for(target="HD120270")
print(f"{len(ids)} frames in log for HD120270")
data, hdr = load_frames(ids)                          # may raise on inhomogeneity

# Narrow to just the science (Purpose contains 'spectrum' or 'science')
ids = find_frames_for(target="HD120270", purpose="science")
data, hdr = load_frames(ids)

# Inspect headers like a numpy array
print(hdr.shape, hdr[0]["FILTER"], hdr[0]["EXPTIME"])
```

## When `find_frames_for` returns a long mixed list

Operator entries are messy — a single target often spans acquisition images,
spectra, multiple nights, nods, and "extra" steps. The matching list will
almost always be inhomogeneous in some header, which is what you want: it
forces you to narrow with `purpose=` (or post-filter the returned IDs) until
the set is scientifically meaningful before passing it to `load_frames`.

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
