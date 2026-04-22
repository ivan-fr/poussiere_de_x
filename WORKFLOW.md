# Pandrosion — Claude Build Workflow

## ⛔ NEVER run `lake build` directly from macOS

macOS iCloud/Time Machine sync corrupts `.lake/build/lib/Pandrosion/Core/`
and `.lake/build/ir/Pandrosion/Core/` by creating phantom duplicate
directories (e.g. `Core 2`, `Legacy 3`). When this happens, the Lean
compiler fails with:

```
failed to write './.lake/build/lib/Pandrosion/Core/<Module>.olean': 2 No such file or directory
```

**Cause:** macOS filesystem sync races with Lake's incremental build.
**Rule:** all Lean compilation goes through the Docker services in
`docker-compose.yml`. The Linux container has no macOS sync layer.

## ⛔ NEVER run `lean-check` (reserved for the user)

**Claude's authorised build service: `lean-incremental` ONLY.**

`lean-check` (clean-build + axiom audit) is **too slow** (10–15 min)
for iterative Claude work and is reserved for the user's manual
final-gate runs. Claude must never invoke it.

If the user explicitly asks Claude to run the strict check, still
only run `lean-incremental` and tell the user to run `lean-check`
themselves from their terminal.

## 🔁 Standard iteration cycle (Claude-only)

### Step 1 — Edit `.lean` files on macOS

Use the normal file-editing tools. Do **not** run `lake build` locally.

### Step 2 — Pre-build macOS sync cleanup

Before every Docker build, purge macOS-sync phantom directories:

```bash
find /Users/ivanbesevic/Documents/poussiere/lean/.lake/build \
  -maxdepth 5 \( -name "Core [0-9]*" -o -name "Legacy [0-9]*" \
    -o -name "* [0-9]" \) -exec rm -rf {} + 2>/dev/null
mkdir -p /Users/ivanbesevic/Documents/poussiere/lean/.lake/build/lib/Pandrosion/Core
mkdir -p /Users/ivanbesevic/Documents/poussiere/lean/.lake/build/lib/Pandrosion/Legacy
mkdir -p /Users/ivanbesevic/Documents/poussiere/lean/.lake/build/ir/Pandrosion/Core
mkdir -p /Users/ivanbesevic/Documents/poussiere/lean/.lake/build/ir/Pandrosion/Legacy
```

### Step 3 — Trust mode (fast, Claude's only build path)

```bash
cd /Users/ivanbesevic/Documents/poussiere
docker compose run --rm lean-incremental
```

**What it does:**
- Reuses `.lake/` cache (no clean rebuild).
- Compiles only modules with dirty dependencies.
- Reports `✅ INCREMENTAL OK — N modules compiled` or `❌ FAILED modules: <list>`.
- **Does NOT** run the axiom audit — that is the user's `lean-check`
  responsibility.

**When `lean-incremental` reports `❌ FAILED modules`:**
1. Check the error log (Docker prints it in the `error: stderr:` section).
2. Fix the offending module.
3. **Re-run step 2 + step 3** (cleanup + incremental).

### Step 4 — Hand-off to the user

When `lean-incremental` reports ✅ green, Claude's job is done on the
build side. Tell the user:

> Ready for strict check. Run from your terminal:
>
> ```bash
> cd /Users/ivanbesevic/Documents/poussiere
> docker compose run --rm lean-check
> ```

The user will then decide whether to run the authoritative gate
(which enforces axiom whitelist and zero warnings/sorries).

## 🚫 Never-do list

- Never `lake build …` directly (macOS sync corruption).
- Never run `docker compose run --rm lean-check` from Claude.
- Never run `docker compose run --rm lean-build` from Claude
  (also a clean rebuild, same reason).
- Never skip Step 2 (cleanup) before Docker builds.
- Never add a new module to `Pandrosion.lean` until it compiles solo
  via `lean-incremental` (test solo via the `lean` service — see
  module-add checklist).

## 🩹 Emergency: if the build state is corrupted

**⚠️ IMPORTANT path:** the actual `.lake/build` lives under
`lean/.lake/build`, **not** at the repo root. A naive `rm -rf .lake/build`
from the repo root does nothing — it silently no-ops.

**Claude's emergency action: purge + incremental.**

```bash
# Full .lake/build purge (correct path):
rm -rf /Users/ivanbesevic/Documents/poussiere/lean/.lake/build

# Then from the repo root:
cd /Users/ivanbesevic/Documents/poussiere
docker compose run --rm lean-incremental   # Claude's full rebuild
```

Note: `lean-incremental` with no cache rebuilds everything from
scratch, same effective cost as `lean-check` but without the strict
audit. This is Claude's correct emergency reset path.

**If incremental rebuild fails with `Directory not empty` during
Lake's own writes:** macOS sync left phantoms. Purge them first:

```bash
find /Users/ivanbesevic/Documents/poussiere/lean/.lake/build \
  -maxdepth 5 \( -name "Core [0-9]*" -o -name "Legacy [0-9]*" \
    -o -name "* [0-9]" \) -exec rm -rf {} + 2>/dev/null
```

Then retry `docker compose run --rm lean-incremental`.

## 📜 Module-add checklist

When adding a new `Pandrosion/Core/<NewModule>.lean`:

1. Write the module. Don't touch `Pandrosion.lean` yet.
2. Cleanup (step 2) + `docker compose run --rm lean-incremental`.
   - Note: `lean-incremental` won't build modules unreachable from
     `Pandrosion.lean`. To test the new module solo:

     ```bash
     cd /Users/ivanbesevic/Documents/poussiere
     docker compose run --rm lean \
       bash -c "lake build Pandrosion.Core.<NewModule>"
     ```

3. Once clean, add `import Pandrosion.Core.<NewModule>` to
   `Pandrosion.lean`.
4. Cleanup (step 2) + `docker compose run --rm lean-incremental` again.
5. **Stop here.** Hand off to the user for `lean-check` (Step 4 above).
