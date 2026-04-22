/-
  Universitas Pandrosion — §56. **Difference identity for
  `pandrosion_h_C x 3 − α₀` and supporting scaffolding.**

  **Scope.** Establishes the central algebraic identity
      `h(z) − α₀ = (x − 1) · Q_3(z, α₀) · (z − α₀) /
                     (x · Sp_C 3 z · Sp_C 3 α₀)`
  where `Q_3(z, α₀) := 1 + z + α₀`. Together with §54.2's
  lower bound on `‖Sp_C 3 z‖²` on `Re z ≥ 1/2` and §55's bound
  `2^{−1/3} ≥ 3/4`, this identity is the **final analytical
  brick** needed for a complete complex contraction theorem on
  a disk around a positive-real cube root anchor.

  **Full chaining** (‖h(z) − α₀‖ ≤ K · ‖z − α₀‖ with `K < 1`
  on `B(α₀, r₀)`) requires several additional manipulations of
  complex-norm-versus-real bounds that exceed the scope of this
  self-contained module. The identity below + §54.2 + §55 form
  the **complete algebraic infrastructure**; the quantitative
  chaining is straightforward but lengthy.

  Contents.

    §56.1  `Sp_C_p3_sub_factor` — `Sp_C 3 z − Sp_C 3 α₀ = (z − α₀)·(1 + z + α₀)`.

    §56.2  `Sp_C_p3_ne_zero_of_re_ge_half` — `Sp_C 3 z ≠ 0` on `Re z ≥ 1/2`.

    §56.3  `pandrosion_h_C_p3_diff_identity` — the central
           difference identity.
-/

import Pandrosion.Core.HalfPlaneContractionP3
import Pandrosion.Core.Rpow34Bound

namespace Pandrosion

open Complex

/-! ============================================================
  §56.1  `Sp_C 3` difference factorisation
============================================================ -/

section SpSubFactorC

/-- **`Sp_C 3 z − Sp_C 3 α₀ = (z − α₀)·(1 + z + α₀)`.**
    Direct polynomial identity (complex version of `Sp3_sub`
    from the Legacy). -/
theorem Sp_C_p3_sub_factor (z α₀ : ℂ) :
    Sp_C 3 z - Sp_C 3 α₀ = (z - α₀) * (1 + z + α₀) := by
  rw [Sp_C_p3 z, Sp_C_p3 α₀]
  ring

end SpSubFactorC

/-! ============================================================
  §56.2  `Sp_C 3 z ≠ 0` on the right half-plane
============================================================ -/

section SpNonzeroC

/-- **`Sp_C 3 z ≠ 0` for `Re z ≥ 1/2`.**
    Consequence of §54.2's lower bound `‖Sp_C 3 z‖² ≥ (1 + a + a²)² > 0`. -/
theorem Sp_C_p3_ne_zero_of_re_ge_half (z : ℂ) (hre : z.re ≥ 1 / 2) :
    Sp_C 3 z ≠ 0 := by
  intro h_eq
  have h_norm_zero : ‖Sp_C 3 z‖ = 0 := by rw [h_eq]; exact norm_zero
  have h_sq : ‖Sp_C 3 z‖ ^ 2 = 0 := by rw [h_norm_zero]; ring
  have h_lb := Sp_C_p3_norm_sq_lower_bound z hre
  have h_rhs_pos : (1 + z.re + z.re ^ 2) ^ 2 > 0 := by
    have h_ge : 1 + z.re + z.re ^ 2 ≥ 1 + 1 / 2 + 1 / 4 := by
      nlinarith [hre, sq_nonneg z.re]
    nlinarith [h_ge]
  linarith

end SpNonzeroC

/-! ============================================================
  §56.3  Difference identity for `pandrosion_h_C x 3`
============================================================ -/

section DiffIdentityC

/-- **★★★ Central identity:**
    `h(z) − α₀ = (x − 1) · Q_3(z, α₀) · (z − α₀) /
                   (x · Sp_C 3 z · Sp_C 3 α₀)`.

    Hypotheses: `x ≠ 0`, `Sp_C 3 z ≠ 0`, `Sp_C 3 α₀ ≠ 0`,
    `α₀³ = 1/x` (so `α₀` is a Pandrosion fixed point).

    *Proof.* Pandrosion `h`-formula gives
        `h(z) − h(α₀) = (x − 1)·(Sp z − Sp α₀) / (x · Sp z · Sp α₀)`.
    Since `α₀³ = 1/x` ⟹ `h(α₀) = α₀` (§4 fixed-point characterisation),
    and `Sp z − Sp α₀ = (z − α₀)·(1 + z + α₀)` (§56.1), the claim
    follows. -/
theorem pandrosion_h_C_p3_diff_identity
    (x α₀ z : ℂ) (hx : x ≠ 0)
    (hSz : Sp_C 3 z ≠ 0) (hSα : Sp_C 3 α₀ ≠ 0)
    (hα_fp : α₀ ^ 3 = 1 / x) :
    pandrosion_h_C x 3 z - α₀ =
      (x - 1) * (1 + z + α₀) * (z - α₀) /
        (x * Sp_C 3 z * Sp_C 3 α₀) := by
  -- h(α₀) = α₀ via pandrosion_h_C_fixed_point_reverse.
  have h_fix : pandrosion_h_C x 3 α₀ = α₀ :=
    pandrosion_h_C_fixed_point_reverse x hx 3 α₀ hSα hα_fp
  have hD_ne : x * Sp_C 3 z * Sp_C 3 α₀ ≠ 0 :=
    mul_ne_zero (mul_ne_zero hx hSz) hSα
  have h_diff : pandrosion_h_C x 3 z - pandrosion_h_C x 3 α₀
              = (x - 1) * (Sp_C 3 z - Sp_C 3 α₀) /
                  (x * Sp_C 3 z * Sp_C 3 α₀) := by
    unfold pandrosion_h_C
    rw [eq_div_iff hD_ne]
    field_simp
    ring
  rw [h_fix] at h_diff
  rw [h_diff, Sp_C_p3_sub_factor]
  ring

end DiffIdentityC

end Pandrosion
