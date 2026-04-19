import Mathlib.Analysis.Complex.Polynomial

open Complex

namespace Pandrosion

/-- The Pandrosion Evaluation Ratio.
For a given polynomial `P`, evaluate the ratio at the target `z` and anchor `a`.
-/
noncomputable def pandrosionRatio (P : Polynomial ℂ) (z a : ℂ) : ℂ :=
  (P.eval z) / (P.eval a)

/-- The Cauchy Circular Evaluation Generator.
Generates `d` equidistant elements on a Cauchy circle of radius `R`.
-/
noncomputable def cauchyRoots (R : ℝ) (d : ℕ) (s : ℕ) : ℂ :=
  R * exp (I * (2 * Real.pi * (s : ℝ) / (d : ℝ)))

/-- Discrete Fourier Transform of the Pandrosion evaluated vector. -/
noncomputable def pandrosionDft (P : Polynomial ℂ) (R : ℝ) (d : ℕ) (k : ℕ) : ℂ :=
  (1 / (d : ℂ)) * ∑' (s : ℕ), 
    if s < d then 
      (pandrosionRatio P 
        (cauchyRoots R d s * exp (I * Real.pi / d)) 
        (cauchyRoots R d s)
      ) * exp (-I * 2 * Real.pi * s * k / d)
    else 0

end Pandrosion
