#!/usr/bin/env bash
# Strict, verbose Lean build + integrity check for the Pandrosion corpus.
set -eo pipefail

cd "$(dirname "$0")/.."

echo "=== Pandrosion Lean Check (strict, verbose) ==="
echo "Toolchain: $(cat lean-toolchain)"
echo "Lean:      $(lean --version)"
echo "Lake:      $(lake --version | head -1)"
echo

echo "--- Build (verbose, warnings = errors) ---"
LOG=$(mktemp)
lake -v build Pandrosion \
  -Kwarning.level=error \
  -KlinterSorries=true 2>&1 | tee "$LOG"

echo
echo "--- Raw lean diagnostics scan ---"
ERRORS=$(grep -cE "^([^:]+\.lean:[0-9]+:[0-9]+: )?error:" "$LOG" || true)
WARNS=$(grep -cE "^([^:]+\.lean:[0-9]+:[0-9]+: )?warning:" "$LOG" || true)
SORRIES=$(grep -cE "(sorry|declaration uses .sorry.)" "$LOG" || true)
AXIOM_USES=$(grep -cE "uses axiom" "$LOG" || true)
echo "  errors in log:   $ERRORS"
echo "  warnings in log: $WARNS"
echo "  sorry mentions:  $SORRIES"
echo "  'uses axiom':    $AXIOM_USES"

echo
echo "--- Source-tree scan (comments stripped) ---"
SRC=$(find Pandrosion/ -name "*.lean" | wc -l | tr -d ' ')
OLE=$(find .lake/build/lib/Pandrosion -name "*.olean" 2>/dev/null | wc -l | tr -d ' ')
THM=$(grep -rE "^(theorem|lemma) " Pandrosion/ --include="*.lean" | wc -l | tr -d ' ')
DEF=$(grep -rE "^(def|noncomputable def) " Pandrosion/ --include="*.lean" | wc -l | tr -d ' ')

SRC_SORRY=$(mktemp)
SRC_ADMIT=$(mktemp)
SRC_AXIOM=$(mktemp)

# Strip /- ... -/ block comments (preserving line count) and -- line comments,
# then scan for real sorry / admit / raw axiom declarations.
while IFS= read -r f; do
  CLEAN=$(perl -0777 -pe 's{(/-.*?-/)}{"\n" x (($1 =~ tr/\n//))}ges; s{--[^\n]*}{}g' "$f")
  printf '%s\n' "$CLEAN" | grep -nP '(?<![-\w])sorry(?![-\w])' | sed "s|^|$f:|" >> "$SRC_SORRY" || true
  printf '%s\n' "$CLEAN" | grep -nP '(?<![-\w])admit(?![-\w])' | sed "s|^|$f:|" >> "$SRC_ADMIT" || true
  printf '%s\n' "$CLEAN" | grep -nE '^axiom '                   | sed "s|^|$f:|" >> "$SRC_AXIOM" || true
done < <(find Pandrosion/ -name "*.lean")

SORRY_CT=$(wc -l < "$SRC_SORRY" | tr -d ' ')
ADMIT_CT=$(wc -l < "$SRC_ADMIT" | tr -d ' ')
AXIOM_CT=$(wc -l < "$SRC_AXIOM" | tr -d ' ')

echo "  .lean files:     $SRC"
echo "  .olean built:    $OLE"
echo "  theorems+lemmas: $THM"
echo "  defs:            $DEF"
echo "  sorry in source: $SORRY_CT"
[ "$SORRY_CT" -gt 0 ] && cat "$SRC_SORRY"
echo "  admit in source: $ADMIT_CT"
[ "$ADMIT_CT" -gt 0 ] && cat "$SRC_ADMIT"
echo "  axiom in source: $AXIOM_CT"
[ "$AXIOM_CT" -gt 0 ] && cat "$SRC_AXIOM"

echo
echo "--- Missing oleans ---"
FAILED=$(comm -23 \
  <(find Pandrosion/ -name "*.lean" | sed 's|/|.|g;s|\.lean||' | sort) \
  <(find .lake/build/lib/Pandrosion -name "*.olean" 2>/dev/null | sed 's|.lake/build/lib/||;s|\.olean||;s|/|.|g' | sort))
if [ -n "$FAILED" ]; then echo "$FAILED"; else echo "  (none)"; fi

FAIL=0
[ -n "$FAILED" ]      && { echo "❌ Missing oleans";                FAIL=1; }
[ "$ERRORS"  -gt 0 ]  && { echo "❌ Errors in build log";           FAIL=1; }
[ "$WARNS"   -gt 0 ]  && { echo "❌ Warnings in build log (strict)"; FAIL=1; }
[ "$SORRY_CT" -gt 0 ] && { echo "❌ sorry present in source";       FAIL=1; }
[ "$ADMIT_CT" -gt 0 ] && { echo "❌ admit present in source";       FAIL=1; }
[ "$AXIOM_CT" -gt 0 ] && { echo "❌ raw axiom present in source";   FAIL=1; }
[ "$OLE" -ne "$SRC" ] && { echo "❌ olean count != src count ($OLE/$SRC)"; FAIL=1; }

echo
if [ "$FAIL" -eq 0 ]; then
  echo "✅ STRICT CHECK PASSED — $SRC modules, $THM theorems/lemmas, $DEF defs, 0 sorry, 0 admit, 0 raw axiom, 0 warning, 0 error"
else
  echo "✗ STRICT CHECK FAILED"
  exit 1
fi
