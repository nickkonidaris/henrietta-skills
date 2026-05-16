# henrietta-skills

Shared Claude Code skills for the Henrietta instrument (Swope 40", H2RG SUTR
imager/spectrograph) — data reduction, instrument analysis, and observing-night
ops.

These are *skills* in the Claude Code sense: each lives in its own directory
with a `SKILL.md`, and Claude loads them on demand when its `description:`
matches what you're doing. They are not Python packages and don't import code
from this repo at runtime.

## Skills

| Skill | What it's for |
|---|---|
| [`scidata-safety`](scidata-safety/SKILL.md) | Never modify files under `/scidata`/`/carnegie/scidata`. Copy-to-workdir pattern. |
| [`frame-set-verify`](frame-set-verify/SKILL.md) | When opening a range of Henrietta frames, assert filter/exptime/ramp-depth homogeneity before computing anything. |
| [`sutr-process`](sutr-process/SKILL.md) | Load `hen{NNNN}_{KKK}r.fits` ramp samples, compute CDS/UTR signal, slope-fit. |
| [`read-noise`](read-noise/SKILL.md) | Standard read-noise measurement on H2RG SUTR data, with UTR/CDS scaling. |
| [`observing-log`](observing-log/SKILL.md) | Parse and query the commissioning observing log CSV. |

## Install

Clone, then run `install.sh` to symlink each skill into `~/.claude/skills/`:

```bash
git clone <your-fork-url> ~/src/henrietta-skills
cd ~/src/henrietta-skills
./install.sh
```

After install, `git pull` is enough to pick up updates — symlinks resolve to
the new content automatically. Run `./install.sh --uninstall` to remove the
symlinks.

To verify, start `claude` in any directory and ask "what henrietta skills do
you have?" — it should list the five above.

## Adding a skill

1. `mkdir my-skill/`
2. Drop a `SKILL.md` in it with this frontmatter:
   ```markdown
   ---
   name: my-skill
   description: One sentence describing when Claude should use this. Be specific — Claude reads this to decide whether to invoke the skill.
   ---
   ```
3. The body is plain markdown — recipes, code snippets, gotchas.
4. `./install.sh` to re-link.

Keep skills small and single-purpose. If a skill grows past ~150 lines, split it.

## Conventions

- All skills assume Python 3 and `astropy`/`numpy` available in the active venv.
- Data live under `/carnegie/scidata/groups/kallisto/Henrietta/data/images/raw/Commissioning/`.
- Ramp samples are named `hen{NNNN}_{KKK}r.fits` (1-indexed KKK, 3 digits). Fitted slope frames are `hen{NNNN}.fits`.
- Gain ≈ 4.0 e⁻/DN (header value is wrong). Plate scale ≈ 0.776"/pix.
- H2RG: 2048×2048, 4-pixel reference border on each side, 4 amplifiers split along the slow-read (X) direction.
