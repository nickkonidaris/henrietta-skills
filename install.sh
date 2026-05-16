#!/usr/bin/env bash
# Symlink each skill directory in this repo into ~/.claude/skills/.
# Run again after `git pull` is not required — symlinks pick up changes.
# Use `./install.sh --uninstall` to remove the symlinks.

set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEST_DIR="${HOME}/.claude/skills"
UNINSTALL=0

for arg in "$@"; do
    case "$arg" in
        --uninstall|-u) UNINSTALL=1 ;;
        -h|--help)
            echo "Usage: $0 [--uninstall]"
            exit 0
            ;;
        *) echo "Unknown arg: $arg" >&2; exit 2 ;;
    esac
done

mkdir -p "$DEST_DIR"

# Every immediate subdir that contains a SKILL.md is a skill.
shopt -s nullglob
for skill_dir in "$REPO_DIR"/*/; do
    [[ -f "${skill_dir}SKILL.md" ]] || continue
    name="$(basename "$skill_dir")"
    link="${DEST_DIR}/${name}"

    if (( UNINSTALL )); then
        if [[ -L "$link" ]]; then
            rm "$link"
            echo "removed  $link"
        fi
        continue
    fi

    if [[ -e "$link" && ! -L "$link" ]]; then
        echo "SKIP     $link (exists and is not a symlink — move it aside first)" >&2
        continue
    fi

    ln -sfn "${skill_dir%/}" "$link"
    echo "linked   $link -> ${skill_dir%/}"
done

if (( UNINSTALL )); then
    echo "Uninstall complete."
else
    echo
    echo "Done. Skills are now available in ~/.claude/skills/."
fi
