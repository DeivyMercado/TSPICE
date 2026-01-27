#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

##################################################################
# tSPICE release script
#
# - Updates the version in:
#   - setup.py
# - Builds (python -m build), validates (twine check), and uploads (twine upload)
# - If anything fails, automatically rolls back to the previous version
#
# NOTE: Unlike MultiNEAs, this script DOES NOT update source files
# (src/tspice/__init__.py) to modify version numbers, as per project rules.
##################################################################

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

SETUP_PY="$ROOT_DIR/setup.py"

BACKUP_DIR=""
ROLLED_BACK=0

log() { printf '%s\n' "$*"; }
err() { printf '%s\n' "$*" >&2; }
die() { err "ERROR: $*"; exit 1; }

have() { command -v "$1" >/dev/null 2>&1; }

python_bin() {
  if have python; then echo "python"; return 0; fi
  if have python3; then echo "python3"; return 0; fi
  return 1
}

PY="$(python_bin)" || die "Could not find 'python' nor 'python3' in PATH."
have git || die "Could not find 'git' in PATH."

rollback() {
  local exit_code=$?

  # Avoid double rollback if multiple traps are triggered.
  if [[ $ROLLED_BACK -eq 1 ]]; then
    exit "$exit_code"
  fi
  ROLLED_BACK=1

  err
  err "Release failed (exit code $exit_code). Rolling back version changes..."

  if [[ -n "${BACKUP_DIR}" && -d "${BACKUP_DIR}" ]]; then
    # Restore modified files
    if [[ -f "${BACKUP_DIR}/setup.py" ]]; then cp -f "${BACKUP_DIR}/setup.py" "$SETUP_PY" || true; fi
  fi

  # Clean build artifacts (do not touch anything else).
  rm -rf "$ROOT_DIR/dist" "$ROOT_DIR/build" "$ROOT_DIR"/*.egg-info "$ROOT_DIR/src"/*.egg-info 2>/dev/null || true

  err "Rollback complete. Repository restored to the previous version (locally)."
  exit "$exit_code"
}

trap rollback ERR INT TERM

usage() {
  cat <<'EOF'
Usage:
  bash bin/release.sh <test|release> <version> [--dry-run|--no-upload]

Examples:
  bash bin/release.sh test 0.1.1
  bash bin/release.sh release 0.2.0
  bash bin/release.sh test 0.2.0 --dry-run

Notes:
  - Requires a clean working tree (no uncommitted changes).
  - If anything fails in build/check/upload, it restores the previous version.
  - With --dry-run/--no-upload, it runs build + twine check, but does NOT upload to PyPI.
EOF
}

TYPE="${1:-}"
VERSION_NEW="${2:-}"
NO_UPLOAD=0

for arg in "${@:3}"; do
  case "$arg" in
    --dry-run|--no-upload) NO_UPLOAD=1 ;;
    *) die "Unknown argument: $arg (use --help)" ;;
  esac
done

if [[ -z "$TYPE" || "$TYPE" == "-h" || "$TYPE" == "--help" ]]; then
  usage
  exit 0
fi

case "$TYPE" in
  test|release) ;;
  *) die "Unknown release type '$TYPE'. Must be 'test' or 'release'." ;;
esac

if [[ ! -f "$SETUP_PY" ]]; then die "Missing file: $SETUP_PY"; fi

# Require a clean repo so rollback is safe.
if ! git diff --quiet || ! git diff --cached --quiet; then
  die "Working tree is not clean. Commit or stash your changes before releasing."
fi

CURRENT_VERSION="$("$PY" - <<'PY'
import pathlib, re, sys

root = pathlib.Path(".").resolve()
setup_path = root / "setup.py"

def extract(path: pathlib.Path) -> str:
    text = path.read_text(encoding="utf-8")
    m = re.search(r"^\s*version\s*=\s*['\"]([^'\"]+)['\"]\s*,?\s*$", text, flags=re.M)
    if not m:
        raise SystemExit(f"Could not find 'version=' in {path}")
    return m.group(1)

print(extract(setup_path))
PY
)"

if [[ -z "$VERSION_NEW" ]]; then
  err "Current version: $CURRENT_VERSION"
  die "You must specify a new version."
fi

# Minimal validation (accepts 0.2.0, 0.2.0.2, 1.0.0rc1, etc.)
if ! printf '%s' "$VERSION_NEW" | grep -Eq '^[0-9]+(\.[0-9]+)+([a-zA-Z0-9\.\-]+)?$'; then
  die "Invalid version '$VERSION_NEW'. Use something like 0.2.0 or 0.2.0.2."
fi

if [[ "$VERSION_NEW" == "$CURRENT_VERSION" ]]; then
  die "New version ($VERSION_NEW) must be different from current ($CURRENT_VERSION)."
fi

log "Releasing tspice $VERSION_NEW (current: $CURRENT_VERSION) in '$TYPE' mode..."

BACKUP_DIR="$(mktemp -d -t tspice-release.XXXXXX)"
cp -f "$SETUP_PY" "${BACKUP_DIR}/setup.py"

update_file() {
  local path="$1"
  local pattern="$2"
  local replacement="$3"

  UPDATE_PATH="$path" PATTERN="$pattern" REPLACEMENT="$replacement" "$PY" - <<'PY'
import os
import pathlib
import re

p = pathlib.Path(os.environ["UPDATE_PATH"])
pattern = os.environ["PATTERN"]
replacement = os.environ["REPLACEMENT"]

text = p.read_text(encoding="utf-8")
new_text, n = re.subn(pattern, replacement, text, count=1, flags=re.M)
if n != 1:
    raise SystemExit(f"Could not update {p} (replacements: {n})")
p.write_text(new_text, encoding="utf-8")
PY
}

log "Updating version in setup.py..."
# Only updating setup.py, not source files
update_file "$SETUP_PY" '^(\s*version\s*=\s*)["'\'']([^"'\'']+)["'\''](\s*,?\s*)$' "\\1'${VERSION_NEW}'\\3"

log "Cleaning previous build artifacts..."
rm -rf "$ROOT_DIR/dist" "$ROOT_DIR/build" "$ROOT_DIR"/*.egg-info "$ROOT_DIR/src"/*.egg-info 2>/dev/null || true

log "Building distributions (sdist/wheel)..."
"$PY" -m build

log "Validating distributions (twine check)..."
"$PY" -m twine check dist/*

if [[ $NO_UPLOAD -eq 1 ]]; then
  log "Dry-run: skipping upload. (build + twine check OK)"
  trap - ERR INT TERM
  rm -rf "$BACKUP_DIR" 2>/dev/null || true
  log "Done. Release (no upload) completed: $VERSION_NEW"
  exit 0
fi

if [[ "$TYPE" == "test" ]]; then
  log "Uploading to TestPyPI..."
  "$PY" -m twine upload --repository testpypi dist/* --verbose
else
  log "Uploading to PyPI..."
  "$PY" -m twine upload dist/* --verbose
fi

# If we got here, everything succeeded.
trap - ERR INT TERM
#rm -rf "$BACKUP_DIR" 2>/dev/null || true
log "Done. Release completed: $VERSION_NEW"