/-
  Universitas Pandrosion — §59. **Steffensen loglog convergence
  for `p = 3` on ℂ (conditional).**

  Chaîne le Steffensen global loglog §34.2 avec la conditionnelle
  McMullen p=3 §41.4 pour obtenir :

      *Modulo `ComplexJuliaP3MeasureZero x α`, pour Lebesgue-presque
      tout `z₀ ∈ ℂ` et tout `ε > 0`, il existe un nombre d'itérations
      `N(z₀, ε)` au-delà duquel l'itéré de Pandrosion-Steffensen est
      à distance ≤ ε d'une racine cyclotomique cubique.*

  C'est la *version loglog* (quadratique, super-attractive) de
  §41.5's `steffensen_p3_solves_complex_mod_julia_null` qui ne donnait
  que la convergence qualitative.

  Contents.

    §59.1  `steffensen_p3_global_loglog_mod_julia_null` — loglog
           global conditionnel à p=3.
-/

import Pandrosion.Core.SteffensenGlobalLoglog
import Pandrosion.Core.ComplexMcMullenP3Conditional

namespace Pandrosion

open Filter Topology MeasureTheory Complex

/-- **★★★★★ Steffensen global loglog for p=3 (conditional).**

    Sous `ComplexJuliaP3MeasureZero x α` (ensemble de Julia de mesure
    nulle à p=3, isolé en §41.4), pour a.e. `z₀ ∈ ℂ` et tout `ε > 0`,
    il existe `N : ℕ` tel que pour tout `k ≥ N`, l'itéré vérifie

        `∃ s : Fin 3, ‖(S^k)(z₀) − γ_s‖ ≤ ε`.

    Le tail `N − k₀(z₀)` est uniforme loglog (O(log log(1/ε))),
    le temps d'entrée `k₀(z₀)` est fini a.e.

    **Obtenu par chaînage direct** :
      §34.2 global loglog a.e. (mod McMullen)
    + §41.4 conditional `McMullenAEEntry 3`. -/
theorem steffensen_p3_global_loglog_mod_julia_null
    (x : ℂ) (hx_ne : x ≠ 0)
    (α : ℂ) (hα_ne_zero : α ≠ 0) (hα_pow : α ^ 3 = 1 / x)
    (hSp : ∀ s : Fin 3, Sp_C 3 (cycAnchor α 3 s) ≠ 0)
    (hJulia : ComplexJuliaP3MeasureZero x α) :
    ∀ᵐ z₀ : ℂ ∂volume, ∀ ε > 0, ∃ N : ℕ,
      ∀ k ≥ N, ∃ s : Fin 3,
        ‖(steffensen_step_C x 3)^[k] z₀ - cycAnchor α 3 s‖ ≤ ε :=
  steffensen_global_loglog_ae_mod_mcmullen 3 (by norm_num)
    x hx_ne α hα_ne_zero hα_pow hSp
    (mcmullen_p3_complex_mod_julia_null x α hJulia)

end Pandrosion
