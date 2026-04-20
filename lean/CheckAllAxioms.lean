/-
  Strict axiom gate: enumerates every theorem/lemma in the Pandrosion
  namespace and prints its transitive axiom dependencies via
  `#print axioms`. The build script parses the resulting output and
  fails if any declaration depends on `sorryAx` or on an axiom outside
  the classical whitelist {propext, Classical.choice, Quot.sound}.
-/
import Pandrosion
import Lean

open Lean Elab Command

/-- Emit `#print axioms` for every theorem in the Pandrosion namespace. -/
elab "#check_all_pandrosion_axioms" : command => do
  let env ← getEnv
  let mut names : Array Name := #[]
  for (n, ci) in env.constants.toList do
    -- Only real declarations (not internal / private / auto-generated)
    if n.isInternal then continue
    -- Scope to the Pandrosion namespace
    unless (`Pandrosion).isPrefixOf n do continue
    -- Only propositions: theorems/lemmas (axioms are also allowed; we want
    -- the gate to blow up if one exists).
    let isProp := match ci with
      | .thmInfo _   => true
      | .axiomInfo _ => true
      | _ => false
    if isProp then
      names := names.push n
  -- Sort for deterministic output
  let names := names.qsort (fun a b => a.toString < b.toString)
  IO.println s!"--- AXIOM AUDIT: {names.size} theorems ---"
  for n in names do
    IO.println s!"## {n}"
    let stx ← `(command| #print axioms $(mkIdent n))
    elabCommand stx

#check_all_pandrosion_axioms
