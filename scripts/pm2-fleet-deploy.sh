#!/usr/bin/env bash
# (Re)start the spring fleet from ONE repo. Use this instead of `pm2 reload`
# when processes were started from old paths (subnet-niome-429, subnet-niome-win, …).
set -euo pipefail

# Default: script lives in <repo>/scripts/ — do not trust NIOME_REPO from ~/.bashrc
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$SCRIPT_DIR/.." && pwd)"
PYTHON="${NIOME_PYTHON:-$REPO/venv/bin/python}"
# pang8512: ecosystem.multi.config.js   spring wallet: ecosystem.spring.config.js
CONFIG="${NIOME_PM2_CONFIG:-$REPO/ecosystem.multi.config.js}"

cd "$REPO"
# Clear placeholder docs path if set in shell (breaks ecosystem.multi.config.js)
case "${NIOME_REPO:-}" in
  *path/to*|*PATH/TO*) unset NIOME_REPO ;;
esac
export NIOME_REPO="$REPO"

if [[ ! -x "$PYTHON" ]]; then
  echo "ERROR: Python not found: $PYTHON" >&2
  echo "Create venv: cd $REPO && python3 -m venv venv && source venv/bin/activate && pip install -e ." >&2
  exit 1
fi

if [[ ! -f "$CONFIG" ]]; then
  echo "ERROR: PM2 config not found: $CONFIG" >&2
  exit 1
fi

echo "==> Repo:   $REPO"
echo "==> Python: $PYTHON"
echo "==> Config: $CONFIG"

"$PYTHON" -c "from niome_subnet.genomics import task_strategy; print('import ok')"

# Delete apps defined in the config (miner-* or spring*)
if [[ "$CONFIG" == *spring* ]]; then
  NAMES=(springhot spring01 spring02 spring03 spring04 spring05 spring06 spring07 spring08 spring09)
else
  NAMES=(miner-1 miner-2 miner-3 miner-4 miner-5)
fi
echo "==> Removing stale PM2 apps (if any)…"
for n in "${NAMES[@]}"; do
  pm2 delete "$n" 2>/dev/null || true
done

echo "==> Starting fleet…"
export NIOME_PYTHON="$PYTHON"
pm2 start "$CONFIG" --update-env
pm2 save

pm2 status

echo ""
echo "==> spring01 paths (must be under $REPO):"
pm2 describe spring01 2>/dev/null | grep -E "exec cwd|script path|NIOME_STRATEGY" || true

echo ""
echo "Tail logs: pm2 logs spring01 --lines 50 --nostream | tail -20"
