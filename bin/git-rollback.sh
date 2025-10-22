#!/usr/bin/env bash
# git-rollabck.sh
# Version: v1.1.0
# Safe rollback helper with guardrails.
#
# Features:
#   --reset <N>                 Hard reset back N commits (HEAD~N)
#   --restore <file> --to <N>   Restore one file from N commits ago (commits it)
#   --revert <N>                Revert the commit N commits ago (N=1 → HEAD)
#   --revert-commit <SHA>       Revert a specific commit SHA
#   --discard                   Discard tracked changes only (safe)  == --discard-tracked
#   --discard-tracked           Discard tracked changes only (safe)
#   --discard-untracked         DANGEROUS: delete untracked files/dirs (needs FORCE=YES)
#   --preview-untracked         Preview what --discard-untracked would delete
#   --help                      Show usage
#
# Notes:
# - Creates a safety tag before changes: git reset --hard <tag> to recover
# - After ops, syncs submodules & LFS if present (best-effort)
# - Reverting a MERGE may need: git revert -m 1 <merge-sha>

set -euo pipefail

usage() {
	cat <<'USAGE'
Usage:
  git-rollabck.sh --reset <N>
  git-rollabck.sh --restore <file> --to <N>
  git-rollabck.sh --revert <N>
  git-rollabck.sh --revert-commit <SHA>
  git-rollabck.sh --discard
  git-rollabck.sh --discard-tracked
  git-rollabck.sh --discard-untracked        # requires FORCE=YES
  git-rollabck.sh --preview-untracked
  git-rollabck.sh --help
USAGE
	exit 1
}

MODE=""
N=0
FILE=""
TO=0
REVERT_SHA=""

while [[ $# -gt 0 ]]; do
	case "$1" in
	--reset)
		MODE="reset"
		N="${2:-0}"
		shift 2
		;;
	--restore)
		MODE="restore"
		FILE="${2:-}"
		shift 2
		;;
	--to)
		TO="${2:-0}"
		shift 2
		;;
	--revert)
		MODE="revert"
		N="${2:-0}"
		shift 2
		;;
	--revert-commit)
		MODE="revert-commit"
		REVERT_SHA="${2:-}"
		shift 2
		;;
	--discard)
		MODE="discard-tracked"
		shift
		;; # alias
	--discard-tracked)
		MODE="discard-tracked"
		shift
		;;
	--discard-untracked)
		MODE="discard-untracked"
		shift
		;;
	--preview-untracked)
		MODE="preview-untracked"
		shift
		;;
	--help | -h) usage ;;
	*)
		echo "Unknown arg: $1"
		usage
		;;
	esac
done

git rev-parse --is-inside-work-tree >/dev/null 2>&1 || {
	echo "Not in a git repo."
	exit 1
}
[[ -n "$MODE" ]] || usage

# Safety tag
stamp="$(date +%Y%m%d-%H%M%S)"
safety_tag="rollback-safety-${stamp}"
git tag "$safety_tag" >/dev/null 2>&1 || true
echo "Created safety tag: $safety_tag"

sync_after() {
	echo "Sync submodules (if any)…"
	git submodule update --init --recursive --force || true
	if git lfs version >/dev/null 2>&1; then
		echo "Sync LFS…"
		git lfs pull || true
	fi
	echo "✅ Done. Current HEAD:"
	git --no-pager log -1 --oneline || true
	echo "Safety tag: $safety_tag"
}

preview_untracked() {
	echo "Preview (untracked files/dirs that would be removed):"
	git clean -nd || true
}

discard_tracked() {
	echo "Discarding changes to TRACKED files only (keeping untracked/ignored)…"
	# Equivalent to: git reset --hard HEAD (does not move HEAD)
	git reset --hard
}

discard_untracked() {
	if [[ "${FORCE:-}" != "YES" ]]; then
		echo "Refusing to delete UNTRACKED files/dirs."
		echo "Preview first:     git-rollabck.sh --preview-untracked"
		echo "To proceed:        FORCE=YES git-rollabck.sh --discard-untracked"
		exit 1
	fi
	echo "DANGEROUS: Deleting UNTRACKED files/dirs (ignored are preserved)…"
	git clean -fd
}

case "$MODE" in
reset)
	[[ "$N" -gt 0 ]] || {
		echo "Need N > 0 for --reset"
		exit 1
	}
	echo "Hard resetting back $N commit(s)…"
	git reset --hard "HEAD~$N"
	sync_after
	;;
restore)
	[[ -n "$FILE" && "$TO" -gt 0 ]] || {
		echo "Need --restore <file> --to <N>"
		exit 1
	}
	echo "Restoring $FILE from HEAD~$TO…"
	git restore --source "HEAD~$TO" -- "$FILE"
	git add "$FILE"
	git commit -m "Restore $FILE from HEAD~$TO"
	sync_after
	;;
revert)
	[[ "$N" -gt 0 ]] || {
		echo "Need N > 0 for --revert"
		exit 1
	}
	target="HEAD~$((N - 1))" # N=1 -> HEAD
	echo "Reverting commit at $target…"
	git revert --no-edit "$target"
	sync_after
	;;
revert-commit)
	[[ -n "$REVERT_SHA" ]] || {
		echo "Need a SHA for --revert-commit"
		exit 1
	}
	echo "Reverting commit $REVERT_SHA…"
	git revert --no-edit "$REVERT_SHA"
	sync_after
	;;
discard-tracked)
	discard_tracked
	sync_after
	;;
discard-untracked)
	preview_untracked
	discard_untracked
	sync_after
	;;
preview-untracked)
	preview_untracked
	echo "No changes made."
	;;
*)
	usage
	;;
esac
