# henrietta-skills

Shared Claude Code skills for the Henrietta instrument (Swope 40", H2RG SUTR
imager/spectrograph).

Each skill lives in its own directory with a `SKILL.md`. Claude loads it on
demand when the `description:` matches what you're doing.

## Skills

| Skill | What it's for |
|---|---|
| [`load-frames`](load-frames/SKILL.md) | Identify, confirm (existence + homogeneity), and read a set of Henrietta FITS frames into numpy arrays plus a header array. Handles fitted slope frames and SUTR ramp samples. |

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
you have?".

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
