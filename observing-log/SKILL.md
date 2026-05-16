---
name: observing-log
description: Use when the user asks to look up frames, targets, conditions, or purposes from the Henrietta commissioning observing log. Knows the CSV's path, column schema, and common query patterns (by file range, filter, grism, purpose, seeing).
---

# observing-log

The commissioning run (2026-04-14 to 2026-05-04) observing log is a Google
Sheets export that lives at:

```
/home/npk/Henrietta/2026A/Commissioning Run (2026-04-14 to 2026-05-04) Observing Log - Sheet1.csv
```

(Path contains spaces and parentheses — always quote it.)

This CSV is the human record. The FITS headers are authoritative for what was
*actually* recorded; the log captures *intent* and *commentary*, which often
disagrees with the headers (operator notes, mid-sequence changes). When the
two disagree, headers win for analysis; the log explains why.

## Schema

| column | meaning |
|---|---|
| `File #` | frame ID (matches `hen{NNNN}.fits`). May be a single int or a range like `123-145`. |
| `Target (RA/Dec)` | target name or coordinates; `None` for cals |
| `Oops` | error flag / operator note about what went wrong |
| `Exp. Time` | seconds |
| `mode` | `CDS` or `SUTR` |
| `Filter` | filter wheel position |
| `Grism` | grism wheel position (`closed` = dark or imaging path) |
| `Slit` | slit width (e.g. `2''`) or empty |
| `Rotator` | rotator angle |
| `Pupil Diffuser` | `None` for science |
| `Slide` | spectrograph slide position |
| `AM` | airmass |
| `Purpose` | high-level intent (`Warm Frame`, `Dark`, `Flat`, target name, etc.) |
| `Note:` | free-form operator notes |
| `Guider X`, `Guider Y`, `Box X`, `Box Y` | guider state |
| `Seeing` | seeing estimate in arcsec |

## Reading the log

```python
import pandas as pd

LOG_PATH = "/home/npk/Henrietta/2026A/Commissioning Run (2026-04-14 to 2026-05-04) Observing Log - Sheet1.csv"
df = pd.read_csv(LOG_PATH)
df.columns = [c.strip() for c in df.columns]   # the trailing colons in "Note:" matter — keep them
```

## File # is sometimes a range

Some rows in the `File #` column are written like `1755-1760` to denote a
sequence taken as one block. Always expand these before filtering by integer
ID:

```python
def expand_ids(s):
    s = str(s).strip()
    if "-" in s:
        a, b = s.split("-")
        return list(range(int(a), int(b) + 1))
    try:
        return [int(s)]
    except ValueError:
        return []

df["ids"] = df["File #"].apply(expand_ids)
df_expanded = df.explode("ids").rename(columns={"ids": "id"})
df_expanded["id"] = pd.to_numeric(df_expanded["id"], errors="coerce")
df_expanded = df_expanded.dropna(subset=["id"]).astype({"id": int})
```

## Common queries

```python
# What was frame 1755?
df_expanded.query("id == 1755")[["Purpose", "Filter", "Grism", "Exp. Time", "Note:"]]

# All R-J flats taken before May 1
flats = df_expanded[df_expanded["Purpose"].str.contains("flat", case=False, na=False) &
                    (df_expanded["Filter"] == "RJ")]

# Sequences flagged Oops by the operator
df[df["Oops"].notna() & (df["Oops"] != "None")][["File #", "Purpose", "Oops"]]
```

## Cross-checking against FITS headers

When the user asks "what was the filter on frame 1755", read it from the
FITS header, not the log:

```python
from astropy.io import fits
hdr = fits.getheader("/carnegie/scidata/groups/kallisto/Henrietta/data/images/raw/Commissioning/hen1755.fits")
print(hdr["FILTER"], hdr["GRISM"], hdr["EXPTIME"])
```

If the header disagrees with the log, **report both**. The log's `Note:` field
often has the explanation.
