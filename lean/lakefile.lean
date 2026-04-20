import Lake
open Lake DSL

package «pandrosion» {
  -- Max-strict Lean options applied to every module in the library.
  moreLeanArgs := #[
    "-Dlinter.unusedVariables=true",
    "-DautoImplicit=false",
    "-DmaxHeartbeats=400000"
  ]
}

lean_lib «Pandrosion» {
  -- add library configuration options here
}

require mathlib from git
  "https://github.com/leanprover-community/mathlib4.git" @ "v4.7.0"
