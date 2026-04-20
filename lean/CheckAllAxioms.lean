/-
  Strict axiom gate (generated section below).

  This file is NOT hand-maintained beyond this header. `scripts/strict_check.sh`
  generates the `#print axioms Pandrosion.<name>` block at build time by
  scanning every `theorem`/`lemma` declaration in `Pandrosion/*.lean`.

  The resulting `#print axioms` output is parsed to enforce:
    - zero `sorryAx`
    - every axiom outside {propext, Classical.choice, Quot.sound} is flagged.
-/
import Pandrosion

-- BEGIN GENERATED AXIOM AUDIT
-- END GENERATED AXIOM AUDIT
