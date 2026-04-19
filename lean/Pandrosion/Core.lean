-- Pandrosion Core: Minimal self-contained definitions
-- No heavy Mathlib imports to avoid cache issues

namespace Pandrosion

/-- A Pandrosion configuration: degree d, radius R -/
structure Config where
  d : ℕ
  hd : d ≥ 2
  R : Float
  hR : R > 0

end Pandrosion
