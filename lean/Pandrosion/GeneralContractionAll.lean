/-
  Universitas Pandrosion — Lean 4 Formalization
  DEEP IX: THE GENERAL CONTRACTION THEOREM — ALL p ≥ 2

  The crown jewel of the Pandrosion theory:
    |h(s) - s*| ≤ (x-1)/x · |s - s*| < |s - s*|
  for ALL degrees p ≥ 2, ALL x > 1, ALL s ≥ 0.

  Proof architecture:
  1. Factor: Sp(s) - Sp(t) = (s-t) · Qp(s,t)
  2. Bound:  Qp(s,t) ≤ Sp(s) · Sp(t)  [subset of monomials]
  3. Conclude: contraction ratio ≤ (x-1)/x < 1

  Reference: pandrosion_master.tex, Theorem 1729 (Universal)
-/
import Mathlib.Data.Real.Basic
import Mathlib.Algebra.GeomSum
import Mathlib.Tactic

open Finset BigOperators

namespace Pandrosion

/-! ## §70. The Divided Difference Qp -/

/-- The divided difference: Qp(s,t) = Σ_{k<p} Σ_{j<k} s^j · t^{k-1-j}. -/
noncomputable def Qp (p : ℕ) (s t : ℝ) : ℝ :=
  ∑ k in range p, ∑ j in range k, s ^ j * t ^ (k - 1 - j)

/-! ## §71. The Factoring Identity (General p) -/

/-- **Factoring: Sp(s) - Sp(t) = (s-t) · Qp(s,t).**
    Uses Mathlib's `geom_sum₂_mul` for each s^k - t^k factor. -/
theorem Sp_sub_factor (p : ℕ) (s t : ℝ) :
    Sp p s - Sp p t = (s - t) * Qp p s t := by
  unfold Sp Qp
  rw [← Finset.sum_sub_distrib]
  rw [Finset.mul_sum]
  congr 1; ext k
  exact pow_sub_factor s t k

/-! ## §72. The Key Bound: Qp ≤ Sp · Sp

Each monomial s^a · t^b in Qp has a+b ≤ p-2.
Each such monomial also appears in Sp(s)·Sp(t).
Since all monomials are non-negative for s,t ≥ 0,
the partial sum ≤ the full sum. -/

/-- **Qp is non-negative for s, t ≥ 0.** -/
theorem Qp_nonneg (p : ℕ) (s t : ℝ) (hs : s ≥ 0) (ht : t ≥ 0) :
    Qp p s t ≥ 0 := by
  unfold Qp
  apply Finset.sum_nonneg; intro k _
  apply Finset.sum_nonneg; intro j _
  exact mul_nonneg (pow_nonneg hs _) (pow_nonneg ht _)

/-- **Sp(s) · Sp(t) ≥ Qp(s,t) for s, t ≥ 0.**

    The p(p-1)/2 monomials of Qp are a subset of the p² monomials
    of Sp·Sp. All monomials are non-negative. Hence Qp ≤ Sp·Sp.

    Technical proof: we convert the nested sum to a sum over
    Finset.sigma, inject it into the product Finset,
    and apply monotonicity of summation. -/
theorem Qp_le_Sp_mul (p : ℕ) (s t : ℝ) (hs : s ≥ 0) (ht : t ≥ 0) :
    Qp p s t ≤ Sp p s * Sp p t := by
  unfold Sp Qp
  rw [Finset.sum_mul_sum]
  -- Goal: Σ_{k<p} Σ_{j<k} s^j*t^{k-1-j} ≤ Σ_{i<p} Σ_{j<p} s^i*t^j
  -- Strategy: for each k, bound Σ_{j<k} s^j*t^{k-1-j} ≤ Σ_{j<p} s^k_stuff
  -- Actually, use the product-filter approach:
  -- Define T = {(a,b) ∈ range p × range p | a+b ≤ p-2}
  -- Show Qp = Σ_(a,b)∈T s^a*t^b (by re-indexing)
  -- Show T ⊆ range p × range p
  -- Apply sum_le_sum_of_subset_of_nonneg
  --
  -- Instead, we use a direct comparison of the two double sums.
  -- For each k < p and j < k, the term s^j*t^{k-1-j} also appears
  -- in the RHS as the (i=j, l=k-1-j) term (since j<k≤p and k-1-j<k≤p).
  --
  -- We prove this using Finset.sum_le_sum on the outer k,
  -- then for each k, we show Σ_{j<k} s^j*t^{k-1-j} ≤ Σ_{l<p} s^k_appropriate.
  -- But matching the inner sums across different outer indices is hard.
  --
  -- Cleanest approach: filter the product sum
  calc ∑ k in range p, ∑ j in range k, s ^ j * t ^ (k - 1 - j)
      = ∑ ab in (range p ×ˢ range p).filter (fun ab : ℕ × ℕ => ab.1 + ab.2 + 2 ≤ p),
          s ^ ab.1 * t ^ ab.2 := by
        -- Prove by induction on p
        induction p with
        | zero =>
          simp only [Finset.range_zero, Finset.sum_empty]
          rw [Finset.empty_product]
          simp
        | succ n ihn =>
          rw [Finset.sum_range_succ]
          -- Split the filter on range(n+1) ×ˢ range(n+1)
          -- = filter on range(n) ×ˢ range(n) ∪ new pairs with k=n or l=n
          -- But the condition a+b+2 ≤ n+1 means a+b ≤ n-1
          -- For the new layer (k=n): Σ_{j<n} s^j * t^{n-1-j}
          --   corresponds to pairs (j, n-1-j) with j+n-1-j = n-1, so a+b+2 = n+1 ≤ n+1 ✓
          -- We use a direct computation
          conv_rhs =>
            rw [show (range (n + 1) ×ˢ range (n + 1)).filter (fun ab : ℕ × ℕ => ab.1 + ab.2 + 2 ≤ n + 1)
              = (range n ×ˢ range n).filter (fun ab : ℕ × ℕ => ab.1 + ab.2 + 2 ≤ n)
              ∪ (range n).map ⟨fun j => (j, n - 1 - j), fun a b h => by simp [Prod.ext_iff] at h; omega⟩
              from by
                ext ⟨a, b⟩
                simp only [Finset.mem_filter, Finset.mem_product, Finset.mem_range,
                           Finset.mem_union, Finset.mem_map, Function.Embedding.coeFn_mk]
                constructor
                · intro ⟨⟨ha, hb⟩, hab⟩
                  by_cases h : a + b + 2 ≤ n
                  · left; exact ⟨⟨by omega, by omega⟩, h⟩
                  · right
                    refine ⟨a, by omega, ?_⟩
                    exact Prod.ext rfl (by omega)
                · intro h
                  rcases h with ⟨⟨ha, hb⟩, hab⟩ | ⟨j, hj, hprod⟩
                  · exact ⟨⟨by omega, by omega⟩, by omega⟩
                  · simp only [Prod.mk.injEq] at hprod
                    exact ⟨⟨by omega, by omega⟩, by omega⟩]
          rw [Finset.sum_union]
          · -- Goal: LHS_old + layer = Σ filter_old + Σ map
            conv_rhs => rw [Finset.sum_map]
            rw [← ihn]
            congr 1
          · -- Disjointness
            rw [Finset.disjoint_left]
            intro ⟨a, b⟩ hmem1 hmem2
            simp only [Finset.mem_filter, Finset.mem_product, Finset.mem_range] at hmem1
            simp only [Finset.mem_map, Function.Embedding.coeFn_mk] at hmem2
            obtain ⟨j, _, hprod⟩ := hmem2
            simp only [Prod.mk.injEq] at hprod
            omega
    _ ≤ ∑ ab in (range p ×ˢ range p),
          s ^ ab.1 * t ^ ab.2 := by
        apply Finset.sum_le_sum_of_subset_of_nonneg (Finset.filter_subset _ _)
        intro i _ _
        exact mul_nonneg (pow_nonneg hs _) (pow_nonneg ht _)
    _ = _ := by
        rw [Finset.sum_product]

/-! ## §73. The General Contraction Theorem -/

/-- **THE GENERAL CONTRACTION THEOREM.**
    For ALL p ≥ 2, x > 1, s ≥ 0, s ≠ s*:
    |h(s) - s*| ≤ (x-1)/x · |s - s*| < |s - s*|.

    Contraction ratio: (x-1)/x < 1. Global convergence. -/
theorem contraction_general (p : ℕ) (hp : p ≥ 2) (x sstar : ℝ) (hx : x > 1)
    (hss_pos : sstar > 0) (hss_lt : sstar < 1)
    (hss_eq : sstar ^ p = 1 / x)
    (s : ℝ) (hs : s ≥ 0) (_hs_ne : s ≠ sstar) :
    |pandrosion_h x p s - sstar| ≤ (x - 1) / x * |s - sstar| := by
  have hx_pos : x > 0 := by linarith
  have hSs : Sp p s > 0 := Sp_pos p (by omega) s hs
  have hSt : Sp p sstar > 0 := Sp_pos p (by omega) sstar (le_of_lt hss_pos)
  have hfix : pandrosion_h x p sstar = sstar :=
    (fixed_point_iff x hx_pos p (by omega) sstar (le_of_lt hss_pos) (ne_of_lt hss_lt)).mpr hss_eq
  have hD_ne : x * Sp p s * Sp p sstar ≠ 0 := by positivity
  -- h(s) - h(s*) = (x-1)(Sp s - Sp s*)/(x·Sp s·Sp s*)
  have hdiff : pandrosion_h x p s - pandrosion_h x p sstar =
      (x - 1) * (Sp p s - Sp p sstar) / (x * Sp p s * Sp p sstar) := by
    unfold pandrosion_h; rw [eq_div_iff hD_ne]; field_simp; ring
  rw [← hfix, hdiff, Sp_sub_factor]
  -- |(x-1) * ((s-s*) * Qp) / denom| ≤ (x-1)/x * |s-s*|
  have hQnn : Qp p s sstar ≥ 0 := Qp_nonneg p s sstar hs (le_of_lt hss_pos)
  have hQbound : Qp p s sstar ≤ Sp p s * Sp p sstar :=
    Qp_le_Sp_mul p s sstar hs (le_of_lt hss_pos)
  -- |(x-1) * (s-s*) * Qp / (x * Sp * Sp*)| ≤ (x-1)/x * |s-s*|
  rw [show (x - 1) * ((s - sstar) * Qp p s sstar) =
      (x - 1) * Qp p s sstar * (s - sstar) from by ring]
  rw [abs_div, abs_mul, abs_mul,
      abs_of_pos (by linarith : x - 1 > 0),
      abs_of_nonneg hQnn,
      abs_of_pos (by positivity : x * Sp p s * Sp p sstar > 0)]
  -- (x-1) * Qp * |s-s*| / (x * Sp * Sp*) ≤ (x-1)/x * |s-s*|
  rw [hfix]
  rw [div_le_iff (by positivity : x * Sp p s * Sp p sstar > 0)]
  -- (x-1)/x * |s-s*| * (x * Sp * Sp*) = (x-1) * |s-s*| * Sp * Sp*
  have hsimp : (x - 1) / x * |s - sstar| * (x * Sp p s * Sp p sstar)
      = (x - 1) * |s - sstar| * (Sp p s * Sp p sstar) := by
    field_simp; ring
  rw [hsimp]
  -- (x-1) * Qp * |s-s*| ≤ (x-1) * |s-s*| * (Sp * Sp*)
  have hnn : (x - 1) * |s - sstar| ≥ 0 := mul_nonneg (by linarith) (abs_nonneg _)
  calc (x - 1) * Qp p s sstar * |s - sstar|
      = (x - 1) * |s - sstar| * Qp p s sstar := by ring
    _ ≤ (x - 1) * |s - sstar| * (Sp p s * Sp p sstar) :=
        mul_le_mul_of_nonneg_left hQbound hnn

/-- **Corollary: strict distance decrease for ALL p ≥ 2.** -/
theorem distance_decreases_general (p : ℕ) (hp : p ≥ 2) (x sstar : ℝ) (hx : x > 1)
    (hss_pos : sstar > 0) (hss_lt : sstar < 1)
    (hss_eq : sstar ^ p = 1 / x)
    (s : ℝ) (hs : s ≥ 0) (hs_ne : s ≠ sstar) :
    |pandrosion_h x p s - sstar| < |s - sstar| := by
  have hbound := contraction_general p hp x sstar hx hss_pos hss_lt hss_eq s hs hs_ne
  have hfrac : (x - 1) / x < 1 := by rw [div_lt_one (by linarith)]; linarith
  linarith [mul_lt_mul_of_pos_right hfrac (abs_pos.mpr (sub_ne_zero.mpr hs_ne))]

end Pandrosion
