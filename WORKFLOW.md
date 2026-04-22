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

### Step 2 — Build (single command)

```bash
cd /Users/ivanbesevic/Documents/poussiere
docker compose run --rm lean-incremental
```

The `lean-incremental` service already handles lock-purge and phantom-dir
cleanup itself. No pre-step is required.

**Expected duration:** a few seconds per touched module (typically
**under 50 seconds** for incremental edits).

**If duration > 50 seconds, investigate immediately:**
- Docker daemon stuck? → `docker ps -a` and restart Docker Desktop.
- Stale Lake locks? → The service clears them on start, but a prior
  crash may have left `.lake/lakefile.olean.lock` held by a dead pid.
- macOS-sync phantom dirs (`Core 2`, `Legacy 3`, `* [0-9]`)?
  → Purge them:

  ```bash
  find /Users/ivanbesevic/Documents/poussiere/lean/.lake/build \
    -maxdepth 5 \( -name "Core [0-9]*" -o -name "Legacy [0-9]*" \
      -o -name "* [0-9]" \) -exec rm -rf {} + 2>/dev/null
  ```

- Full cache corruption? → Emergency reset (see below).

### Step 3 — Hand-off to the user

When `lean-incremental` reports `✅ INCREMENTAL OK — N modules compiled`,
Claude's job is done on the build side. Tell the user:

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
- Never add a new module to `Pandrosion.lean` until it compiles solo
  via `lean-incremental` (test solo via the `lean` service — see
  module-add checklist).

## 🩹 Emergency: if the build state is corrupted

**⚠️ IMPORTANT path:** the actual `.lake/build` lives under
`lean/.lake/build`, **not** at the repo root. A naive `rm -rf .lake/build`
from the repo root does nothing — it silently no-ops.

**Claude's emergency action: full purge + incremental.**

```bash
rm -rf /Users/ivanbesevic/Documents/poussiere/lean/.lake/build
cd /Users/ivanbesevic/Documents/poussiere
docker compose run --rm lean-incremental
```

Note: with no cache, this rebuilds everything from scratch — same
effective cost as `lean-check` but without the strict audit.

## 🔬 Empirical exploration with Python (before formalization)

**Strong recommendation:** before attempting to formalize an ambitious
property in Lean, **test it numerically first** with Python. This has
repeatedly caught dead-ends (false conjectures) and revealed the
*minimal* structural claim that needs formalizing (e.g., reducing a
"measure zero" claim to a "finite set" claim).

### Python environment

A venv with `numpy`, `scipy`, `matplotlib` is pre-provisioned at
`/tmp/pandros_venv`. Recreate it if missing:

```bash
python3 -m venv /tmp/pandros_venv
/tmp/pandros_venv/bin/pip install --quiet numpy scipy matplotlib
```

Run scripts via:

```bash
/tmp/pandros_venv/bin/python3 /tmp/your_script.py
```

Prefer `/tmp/` for exploratory scripts — they stay out of the repo.
Promote a script to `articles/` or similar only if it becomes a
reproducible artefact referenced from a module header.

### When to use each library

- **`numpy`** — vectorised scans of initial conditions (2D grids of
  `z₀ ∈ ℂ`), iteration of `h`, `σ`, basin classification by distance
  to the known roots. Foundation for "does the conjecture even look
  true?".
- **`scipy`** — root-finding (`scipy.optimize.brentq`, `newton`),
  period-k cycle search, solving systems to locate repelling
  periodic points. Useful when you suspect a Julia-set component
  other than the super-attracting fixed points.
- **`matplotlib`** — visualise basins of attraction, escape-time
  fractals, level sets. Save figures to `/tmp/` and check them by
  reading with the `Read` tool (it renders PNGs inline).

### Typical exploration pattern

1. **Formulate a candidate Lean theorem** (e.g., "the Julia set at
   x=2, p=3 has Lebesgue measure zero").
2. **Translate into a falsifiable numeric test**:
   - Grid-scan a large window (say 600×600 on `[−3, 3]²`).
   - Classify each `z₀` by final behaviour after `max_iter = 500`.
   - Count divergent/cycling starts.
3. **Interpret**:
   - If the count is `~0` except at isolated points → the structural
     obstacle is *finite* (or countable) — aim for a `.Finite` or
     `.Countable` lemma in Lean.
   - If the count scales with grid density → the obstacle has
     positive measure; the conjecture is likely false as stated.
4. **Refine the Lean target** based on what the numerics show is
   achievable. Formalize only the *minimal* unconditional piece,
   isolate the rest as a named `Prop` hypothesis.

### Example (§61 `JuliaNullX2P3`)

The path that worked:
1. Hypothesis: Julia set of σ at `x=2, p=3` has measure zero.
2. Python scan `[−3, 3]² @ 600×600, max_iter=500`: **zero
   divergent starts** outside the two Sp-zeros `(−1 ± i√3)/2`.
3. Zoom `radius ∈ {10⁻¹, …, 10⁻⁹}` around each Sp-zero: still zero
   divergence except *exactly* on the singularity.
4. scipy period-2 cycle search: **zero candidates**.
5. Conclusion: the "bad set" is (conjecturally) the **countable**
   orbit of Sp-zeros. Formalize the unconditional piece — the
   Sp-zero set is `{ω, ω²}` (finite, measure zero).

That reduced a full-McMullen-strength claim to a 200-line algebraic
module.

## 📜 Module-add checklist

When adding a new `Pandrosion/Core/<NewModule>.lean`:

1. Write the module. Don't touch `Pandrosion.lean` yet.
2. `docker compose run --rm lean-incremental`.
   - Note: `lean-incremental` won't build modules unreachable from
     `Pandrosion.lean`. To test the new module solo:

     ```bash
     cd /Users/ivanbesevic/Documents/poussiere
     docker compose run --rm lean \
       bash -c "lake build Pandrosion.Core.<NewModule>"
     ```

3. Once clean, add `import Pandrosion.Core.<NewModule>` to
   `Pandrosion.lean`.
4. `docker compose run --rm lean-incremental` again.
5. **Stop here.** Hand off to the user for `lean-check` (Step 3 above).
