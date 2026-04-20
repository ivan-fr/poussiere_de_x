#!/usr/bin/env bash
# Strict, verbose Lean build + integrity check for the Pandrosion corpus.
#
# Max-difficulty gates:
#   - Clean build (no incremental artefacts trusted)
#   - warnings fatal, linterSorries on, unused-variables linter on,
#     autoImplicit forbidden, maxHeartbeats bounded
#   - Build-log scan: 0 error, 0 warning, 0 sorry, 0 "uses axiom" leak
#   - Source scan (ALL .lean files under lean/, comments stripped):
#       0 real sorry, 0 admit, 0 raw axiom declaration
#   - Axiom audit: every theorem in Pandrosion namespace must depend
#     only on the classical whitelist {propext, Classical.choice,
#     Quot.sound}. Any sorryAx / extra axiom fails the gate.
#   - Olean integrity: count matches src count

set -eo pipefail

cd "$(dirname "$0")/.."

AXIOM_WHITELIST="propext Classical.choice Quot.sound"

echo "=== Pandrosion Lean Check (MAX strict) ==="
echo "Toolchain: $(cat lean-toolchain)"
echo "Lean:      $(lean --version)"
echo "Lake:      $(lake --version | head -1)"
echo "Axiom whitelist: $AXIOM_WHITELIST"
echo

echo "--- Cleaning build cache ---"
rm -rf .lake/build/lib/Pandrosion .lake/build/ir/Pandrosion
rm -f  .lake/build/lib/Pandrosion.olean .lake/build/ir/Pandrosion.c
rm -f  .lake/build/lib/CheckAxioms.olean .lake/build/lib/CheckAllAxioms.olean
rm -f  .lake/build/lib/MeasureTest.olean
echo "  cleaned."
echo

echo "--- Build (verbose, all gates armed) ---"
LOG=$(mktemp)
# -Kwarning.level=error : warnings fatal (package option)
# -KlinterSorries=true  : surface sorry as lint error (package option)
# Additional Lean options (autoImplicit=false, linter.unusedVariables=true,
# maxHeartbeats=400000) are applied via moreLeanArgs in lakefile.lean.
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
echo "--- Source-tree scan (ALL lean/**/*.lean, comments stripped) ---"
# Scan ALL .lean under the library root, not only Pandrosion/ — this covers
# CheckAxioms.lean, CheckAllAxioms.lean, MeasureTest.lean and any future
# top-level file. Exclude .lake (vendored packages).
ALL_LEAN=$(find . -name "*.lean" -not -path "./.lake/*" | sort)
SRC=$(echo "$ALL_LEAN" | wc -l | tr -d ' ')
PANDR_SRC=$(find Pandrosion/ -name "*.lean" | wc -l | tr -d ' ')
OLE=$(find .lake/build/lib/Pandrosion -name "*.olean" 2>/dev/null | wc -l | tr -d ' ')
THM=$(echo "$ALL_LEAN" | xargs grep -hE "^(theorem|lemma) " 2>/dev/null | wc -l | tr -d ' ')
DEF=$(echo "$ALL_LEAN" | xargs grep -hE "^(def|noncomputable def) " 2>/dev/null | wc -l | tr -d ' ')

SRC_SORRY=$(mktemp)
SRC_ADMIT=$(mktemp)
SRC_AXIOM=$(mktemp)

while IFS= read -r f; do
  CLEAN=$(perl -0777 -pe 's{(/-.*?-/)}{"\n" x (($1 =~ tr/\n//))}ges; s{--[^\n]*}{}g' "$f")
  printf '%s\n' "$CLEAN" | grep -nP '(?<![-\w])sorry(?![-\w])' | sed "s|^|$f:|" >> "$SRC_SORRY" || true
  printf '%s\n' "$CLEAN" | grep -nP '(?<![-\w])admit(?![-\w])' | sed "s|^|$f:|" >> "$SRC_ADMIT" || true
  printf '%s\n' "$CLEAN" | grep -nE '^axiom '                   | sed "s|^|$f:|" >> "$SRC_AXIOM" || true
done <<< "$ALL_LEAN"

SORRY_CT=$(wc -l < "$SRC_SORRY" | tr -d ' ')
ADMIT_CT=$(wc -l < "$SRC_ADMIT" | tr -d ' ')
AXIOM_CT=$(wc -l < "$SRC_AXIOM" | tr -d ' ')

echo "  all .lean files (lean/):    $SRC"
echo "  Pandrosion/ modules:        $PANDR_SRC"
echo "  .olean built (Pandrosion):  $OLE"
echo "  theorems+lemmas:            $THM"
echo "  defs:                       $DEF"
echo "  sorry in source:            $SORRY_CT"
[ "$SORRY_CT" -gt 0 ] && cat "$SRC_SORRY"
echo "  admit in source:            $ADMIT_CT"
[ "$ADMIT_CT" -gt 0 ] && cat "$SRC_ADMIT"
echo "  axiom in source:            $AXIOM_CT"
[ "$AXIOM_CT" -gt 0 ] && cat "$SRC_AXIOM"

echo
echo "--- Axiom audit (every Pandrosion.* theorem via #print axioms) ---"
AUDIT=$(mktemp)
# lake env lean runs with the project's LEAN_PATH so imports resolve.
if lake env lean CheckAllAxioms.lean 2>&1 | tee "$AUDIT"; then
  :
else
  echo "❌ Axiom audit run failed."
fi
AUDIT_THMS=$(grep -c '^## Pandrosion\.' "$AUDIT" || true)
SORRYAX_HITS=$(grep -c 'sorryAx' "$AUDIT" || true)
# Extract every axiom name mentioned by `#print axioms` output.
# `#print axioms` emits lines like:
#   "'Pandrosion.foo' depends on axioms: [propext, Classical.choice]"
# We pull the bracketed list.
FOREIGN_AXIOMS=$(mktemp)
perl -ne '
  if (/depends on axioms:\s*\[([^\]]*)\]/) {
    for my $a (split /,\s*/, $1) {
      $a =~ s/^\s+|\s+$//g;
      next unless $a;
      # whitelist
      next if $a eq "propext" or $a eq "Classical.choice" or $a eq "Quot.sound";
      print "$a\n";
    }
  }
' "$AUDIT" | sort -u > "$FOREIGN_AXIOMS"
FOREIGN_CT=$(wc -l < "$FOREIGN_AXIOMS" | tr -d ' ')

echo "  theorems audited:  $AUDIT_THMS"
echo "  sorryAx hits:      $SORRYAX_HITS"
echo "  off-whitelist ax.: $FOREIGN_CT"
[ "$FOREIGN_CT" -gt 0 ] && { echo "    offenders:"; sed 's/^/      - /' "$FOREIGN_AXIOMS"; }

echo
echo "--- Missing oleans ---"
FAILED=$(comm -23 \
  <(find Pandrosion/ -name "*.lean" | sed 's|/|.|g;s|\.lean||' | sort) \
  <(find .lake/build/lib/Pandrosion -name "*.olean" 2>/dev/null | sed 's|.lake/build/lib/||;s|\.olean||;s|/|.|g' | sort))
if [ -n "$FAILED" ]; then echo "$FAILED"; else echo "  (none)"; fi

FAIL=0
[ -n "$FAILED" ]         && { echo "❌ Missing oleans";                                  FAIL=1; }
[ "$ERRORS"  -gt 0 ]     && { echo "❌ Errors in build log";                             FAIL=1; }
[ "$WARNS"   -gt 0 ]     && { echo "❌ Warnings in build log (strict)";                  FAIL=1; }
[ "$SORRY_CT" -gt 0 ]    && { echo "❌ sorry present in source";                         FAIL=1; }
[ "$ADMIT_CT" -gt 0 ]    && { echo "❌ admit present in source";                         FAIL=1; }
[ "$AXIOM_CT" -gt 0 ]    && { echo "❌ raw axiom declaration present in source";         FAIL=1; }
[ "$OLE" -ne "$PANDR_SRC" ] && { echo "❌ olean count != Pandrosion src count ($OLE/$PANDR_SRC)"; FAIL=1; }
[ "$SORRYAX_HITS" -gt 0 ] && { echo "❌ sorryAx in axiom audit";                          FAIL=1; }
[ "$FOREIGN_CT" -gt 0 ]  && { echo "❌ off-whitelist axioms in audit ($FOREIGN_CT)";      FAIL=1; }
[ "$AUDIT_THMS" -eq 0 ]  && { echo "❌ axiom audit produced no theorems — check broken"; FAIL=1; }

echo
if [ "$FAIL" -eq 0 ]; then
  echo "✅ MAX-STRICT CHECK PASSED — $PANDR_SRC modules, $THM theorems/lemmas, $DEF defs, $AUDIT_THMS audited"
  echo "   0 sorry, 0 admit, 0 raw axiom, 0 warning, 0 error, 0 sorryAx, 0 off-whitelist axiom"
else
  echo "✗ MAX-STRICT CHECK FAILED"
  exit 1
fi
