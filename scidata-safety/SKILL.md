---
name: scidata-safety
description: Use whenever about to read, write, edit, or move files anywhere under /scidata or /carnegie/scidata (Henrietta raw data store). Enforces the copy-before-modify rule so raw images and pipeline outputs in /scidata are never mutated in place.
---

# scidata-safety

`/scidata` (mounted at `/carnegie/scidata/...`) is the **authoritative,
shared, read-only** store for Henrietta raw frames, fitted slope frames, and
pipeline outputs. Multiple people and pipelines depend on it.

## The rule

**Never modify, overwrite, rename, delete, or write into anything under
`/scidata` or `/carnegie/scidata`.** That includes:

- Editing FITS headers in place
- Running pipelines that write outputs under the raw-data tree
- Creating sidecar files (`.bak`, `.csv`, `.npz`) alongside raw frames
- `chmod`, `chown`, `touch`, anything that changes mtime

If a tool *might* write under `/scidata`, route its output to the local
project workdir instead.

## Allowed operations

- Read: `fits.open(...)`, `astropy.io.fits.getdata(...)`, `np.fromfile(...)`
- List/glob: `Path("/carnegie/scidata/...").glob("hen0*.fits")`
- Anything else: **copy first**, then modify the copy.

## Copy-to-workdir pattern

When a task requires modifying a frame (e.g. cosmetic flagging, re-headering,
re-fitting), copy the inputs to a per-experiment subdirectory **first**:

```python
from pathlib import Path
import shutil

SCIDATA = Path("/carnegie/scidata/groups/kallisto/Henrietta/data/images/raw/Commissioning")
WORK = Path("data/copies")           # under the experiment dir, relative to cwd
WORK.mkdir(parents=True, exist_ok=True)

for frame_id in (1755, 1756, 1757):
    src = SCIDATA / f"hen{frame_id}.fits"
    dst = WORK / src.name
    if not dst.exists():
        shutil.copy2(src, dst)       # copy2 preserves mtime
```

For ramp samples (`hen{NNNN}_{KKK}r.fits`), copy the whole set in one loop and
verify count before proceeding.

## If a write to /scidata is really needed

Stop. Surface it to the user — they will decide whether to run it themselves
with elevated permissions, or to add the result to a side directory. Do not
"just try it."

## Verifying

Before any pipeline run that takes a path argument, check the resolved output
path:

```python
out_dir = Path(args.out).resolve()
assert "/scidata" not in str(out_dir), f"refusing to write under /scidata: {out_dir}"
```

A one-line guard at pipeline entry is cheap insurance.
