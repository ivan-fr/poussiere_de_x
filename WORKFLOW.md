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

## 🔁 Standard iteration cycle

### Step 1 — Edit `.lean` files on macOS

Use the normal file-editing tools. Do **not** run `lake build` locally.

### Step 2 — Pre-build macOS sync cleanup

Before every Docker build, purge macOS-sync phantom directories:

```bash
find /Users/ivanbesevic/Documents/poussiere/lean/.lake/build \
  -maxdepth 5 -name "*[ ][0-9]*" 2>/dev/null | xargs -I {} rm -rf {} 2>/dev/null
mkdir -p /Users/ivanbesevic/Documents/poussiere/lean/.lake/build/lib/Pandrosion/Core
mkdir -p /Users/ivanbesevic/Documents/poussiere/lean/.lake/build/lib/Pandrosion/Legacy
mkdir -p /Users/ivanbesevic/Documents/poussiere/lean/.lake/build/ir/Pandrosion/Core
mkdir -p /Users/ivanbesevic/Documents/poussiere/lean/.lake/build/ir/Pandrosion/Legacy
```

### Step 3 — Trust mode (fast, for iteration)

```bash
cd /Users/ivanbesevic/Documents/poussiere
docker compose run --rm lean-incremental
```

**What it does:**
- Reuses `.lake/` cache (no clean rebuild).
- Compiles only modules with dirty dependencies.
- Reports `✅ INCREMENTAL OK — N modules compiled` or `❌ FAILED modules: <list>`.
- **Does NOT** run the axiom audit — use `lean-check` for that.

**When `lean-incremental` reports `❌ FAILED modules`:**
1. Check the error log (Docker prints it in the `error: stderr:` section).
2. Fix the offending module.
3. **Re-run step 2 + step 3** (cleanup + incremental).

### Step 4 — Authoritative gate (strict, slow, once everything is green)

After `lean-incremental` passes clean:

```bash
cd /Users/ivanbesevic/Documents/poussiere
docker compose run --rm lean-check
```

**What it does:**
- Clean build of all modules.
- Enforces `-Kwarning.level=error -KlinterSorries=true`.
- Runs the axiom audit: every Pandrosion.* theorem's axiom dependencies
  must be in the whitelist `{propext, Classical.choice, Quot.sound}`.
- Fails on: missing oleans, errors, warnings, `sorry`, `admit`, raw
  `axiom`, `sorryAx` hits, off-whitelist axioms.

**Only ship when `lean-check` prints:**

```
✅ MAX-STRICT CHECK PASSED — N modules, M theorems/lemmas, …
   0 sorry, 0 admit, 0 raw axiom, 0 warning, 0 error, 0 sorryAx, 0 off-whitelist axiom
```

## 🚫 Never-do list

- Never `lake build …` directly (macOS sync corruption).
- Never skip Step 2 (cleanup) before Docker builds.
- Never add a new module to `Pandrosion.lean` until it compiles solo
  via `lean-incremental`.
- Never ship without `lean-check` ✅.

## 🩹 Emergency: if the build state is corrupted

**⚠️ IMPORTANT path:** the actual `.lake/build` lives under
`lean/.lake/build`, **not** at the repo root. A naive `rm -rf .lake/build`
from the repo root does nothing — it silently no-ops.

```bash
# Full nuclear reset (correct path):
rm -rf /Users/ivanbesevic/Documents/poussiere/lean/.lake/build

# Then from the repo root:
cd /Users/ivanbesevic/Documents/poussiere
docker compose run --rm lean-check   # full clean rebuild
```

This forces Lake to rebuild everything from scratch inside the Linux
container. Takes ~10–15 min but guarantees a known-good state.

**If the strict_check.sh clean step fails with `Directory not empty`:**
macOS sync left phantoms. Purge them first:

```bash
find /Users/ivanbesevic/Documents/poussiere/lean/.lake/build \
  -maxdepth 5 \( -name "Core [0-9]*" -o -name "Legacy [0-9]*" \
    -o -name "* [0-9]" \) -exec rm -rf {} + 2>/dev/null
```

Then retry the `docker compose run --rm lean-check`.

## 📜 Module-add checklist

When adding a new `Pandrosion/Core/<NewModule>.lean`:

1. Write the module. Don't touch `Pandrosion.lean` yet.
2. Cleanup (step 2) + `docker compose run --rm lean-incremental`.
   - Note: incremental won't build modules unreachable from
     `Pandrosion.lean`. To test the new module solo:

     ```bash
     cd /Users/ivanbesevic/Documents/poussiere
     docker compose run --rm lean \
       bash -c "lake build Pandrosion.Core.<NewModule>"
     ```

3. Once clean, add `import Pandrosion.Core.<NewModule>` to
   `Pandrosion.lean`.
4. Cleanup (step 2) + `docker compose run --rm lean-incremental` again.
5. Once clean, `docker compose run --rm lean-check` for the gate.
