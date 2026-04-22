#!/usr/bin/env bash
# Pandrosion lean-incremental wrapper with pre-run container cleanup.
#
# Why: `docker compose run --rm lean-incremental` normally cleans up its own
# container via --rm, but when a run is interrupted (Ctrl+C, broken pipe,
# stale claude subprocess, etc.) the container can be left alive. After a few
# interrupted runs, several poussiere-lean-incremental-run-* containers
# accumulate. Each one holds the shared .lake volume and CPU — a build that
# normally takes ~15 s can slow to 4+ minutes because the Docker host is
# juggling 7 concurrent lean processes.
#
# This wrapper kills any stale poussiere-lean-incremental* containers before
# invoking docker compose, guaranteeing a clean slate.
#
# Usage (from repo root):
#   bash scripts/lean-incremental.sh

set -eo pipefail

# 1. Kill any lingering poussiere-lean-incremental-* containers.
STALE=$(docker ps -q --filter "name=poussiere-lean-incremental-run-" 2>/dev/null || true)
if [ -n "$STALE" ]; then
  echo "⚠️  Killing $(echo "$STALE" | wc -l | tr -d ' ') stale lean-incremental container(s)..."
  echo "$STALE" | xargs -r docker kill >/dev/null 2>&1 || true
fi

# 2. Remove any stopped poussiere-lean-incremental-* containers.
STOPPED=$(docker ps -aq --filter "name=poussiere-lean-incremental-run-" --filter "status=exited" 2>/dev/null || true)
if [ -n "$STOPPED" ]; then
  echo "$STOPPED" | xargs -r docker rm -f >/dev/null 2>&1 || true
fi

# 3. Invoke the compose service.
exec docker compose run --rm lean-incremental "$@"
