/-
  Universitas Pandrosion — Axiom audit.

  This file exists so the CI axiom-audit pipeline can enumerate the
  trusted base of every theorem in the corpus by running

      `#print axioms <theorem_name>`

  for each public theorem in `Pandrosion.Core.*`. The expected output
  for almost every theorem is the standard Lean 4 / Mathlib base:

      `Classical.choice  Quot.sound  propext`

  The only non-standard axiom admitted on purpose is
  `Pandrosion.kung_traub_conjecture` (an explicit, documented axiom
  in `Pandrosion/Core/KungTraub.lean`), which propagates into
  theorems that depend on Kung-Traub optimality.

  Organized in the dependency order of `Pandrosion.lean`:
    §1  Foundations
    §2  MultiStartBasins
    §3  QuadraticRate
    §4  QuadraticComplexity
    §5  VoronoiMeasure
    §6  CyclotomicMcMullen
    §7  KungTraub
    §8  UniformComplexity
    §9  SuperGrandMaster
    §10 UniformContractionRate
    §11 ConcreteIteration
    §12 MasterAbsolu
    §13 ComplexMultiplier
    §14 LocalAttraction
    §15 DynamicalConvergence
-/

import Pandrosion

open Pandrosion

/-! ================================================================
  §1  Foundations
================================================================ -/

#print axioms Sp_mul_one_sub
#print axioms Sp_at_one
#print axioms Sp_pos
#print axioms pandrosion_fixed_point_iff
#print axioms pandrosion_quasi_newton
#print axioms scaling_power
#print axioms scaling_factorization
#print axioms pandrosion_s0_opt_eq_h_one
#print axioms offset_monotone
#print axioms scaled_s0_opt_closer_to_one
#print axioms double_preconditioning
#print axioms steffensen_step_fixed_point_of_h
#print axioms steffensen_step_fixed_point
#print axioms steffensen_rate_ne_one
#print axioms Sp_C_mul_one_sub
#print axioms Sp_C_at_one
#print axioms pandrosion_h_C_quasi_newton
#print axioms pandrosion_h_C_fixed_point_forward
#print axioms pandrosion_h_C_fixed_point_reverse
#print axioms pandrosion_h_C_fixed_point_iff
#print axioms steffensen_step_C_fixed_point_of_h
#print axioms steffensen_step_C_fixed_point
#print axioms scaling_power_C
#print axioms scaling_factorization_C
#print axioms scaling_steffensen_C
#print axioms pandrosion_s0_opt_C_eq_h_one
#print axioms optimal_scaling_steffensen_C
#print axioms cauchy_radius_ge_two
#print axioms cauchy_radius_pos
#print axioms voronoi_nearest_exists
#print axioms multi_start_grand_master
#print axioms voronoi_unique_off_boundary
#print axioms Sp_gt_p_of_gt_one
#print axioms lambda_closed_pos
#print axioms lambda_closed_lt_one
#print axioms lambda_closed_at_one
#print axioms lambda_closed_lt_uniform_bound
#print axioms one_sub_lambda_closed_ne_zero
#print axioms steffensen_constant_denom_ne_zero
#print axioms K_steffensen_nonneg

/-! ================================================================
  §2  MultiStartBasins
================================================================ -/

#print axioms fin_injective_iff_bijective
#print axioms multi_start_sigma_bijective
#print axioms generic_convergence_conditional
#print axioms basin_eq_iInter
#print axioms halfPlane_convex
#print axioms basin_convex
#print axioms basin_self_mem
#print axioms basin_nonempty
#print axioms basin_isConnected
#print axioms basin_isPreconnected
#print axioms basin_isClosed
#print axioms basin_frontier_on_bisector
#print axioms basin_boundary_finite_union
#print axioms basin_anchor_mem_interior
#print axioms basin_interior_nonempty
#print axioms basins_cover
#print axioms VoronoiBoundary_subset_bisectors
#print axioms VoronoiBoundary_eq_bisectors_filtered
#print axioms constructive_mcmullen

/-! ================================================================
  §3  QuadraticRate
================================================================ -/

#print axioms pandrosion_h_sub_mul
#print axioms Qpoly_mul_sub
#print axioms pow_sub_pow_factor
#print axioms Qpoly_diag
#print axioms pandrosion_h_sub_fp_factor
#print axioms pandrosion_h_sub_fp_quotient
#print axioms lambda_fp_eq_quotient
#print axioms pandrosion_linear_rate_ratio
#print axioms pandrosion_ratio_at_fp
#print axioms quadratic_deviation_at_fp
#print axioms Sp_sub_Sp_factor
#print axioms Qpoly_sub_Qpoly_factor
#print axioms lambda_fp_mul_Sp
#print axioms pandrosion_h_sub_fp_factor_no_x
#print axioms pandrosion_h_quadratic_polynomial
#print axioms pandrosion_h_quadratic_quotient
#print axioms pandrosion_h_quadratic_bound_of
#print axioms Sp_reflect
#print axioms lambda_fp_eq_lambda_closed

/-! ================================================================
  §4  QuadraticComplexity
================================================================ -/

#print axioms quadratic_tower_bound
#print axioms quadratic_loglog_complexity
#print axioms pandrosion_steffensen_loglog_complexity

/-! ================================================================
  §5  VoronoiMeasure
================================================================ -/

#print axioms perpBisectorForm_apply
#print axioms perpBisectorForm_ne_zero_of_ne
#print axioms perpBisector_volume_zero
#print axioms voronoiBoundary_volume_zero
#print axioms constructive_mcmullen_ae
#print axioms voronoi_selector_unique_ae
#print axioms constructive_mcmullen_ae_proper

/-! ================================================================
  §6  CyclotomicMcMullen
================================================================ -/

#print axioms cycAnchor_pow
#print axioms cycAnchor_injective
#print axioms unconditional_mcmullen

/-! ================================================================
  §7  KungTraub

  Note: theorems in this module legitimately depend on the explicit
  axiom `Pandrosion.kung_traub_conjecture`, which is the load-bearing
  number-theoretic assumption quoted in the paper.
================================================================ -/

#print axioms efficiencyIndex_two_two
#print axioms pandrosion_const_lt_newton
#print axioms pandrosion_const_pos
#print axioms pandrosion_attains_kung_traub
#print axioms pandrosion_kung_traub_optimal
#print axioms pandrosion_kung_traub_optimal_bound

/-! ================================================================
  §8  UniformComplexity
================================================================ -/

#print axioms uniform_quadratic_loglog
#print axioms uniform_quadratic_loglog_scaled
#print axioms pandrosion_uniform_loglog

/-! ================================================================
  §9  SuperGrandMaster
================================================================ -/

#print axioms linear_decay_bound
#print axioms scaled_basin_entry_time
#print axioms two_phase_complexity
#print axioms super_grand_master_uniform

/-! ================================================================
  §10  UniformContractionRate
================================================================ -/

#print axioms lamInfty_lt_one_of_one_lt
#print axioms lamInfty_pos_of_one_lt
#print axioms lamInfty_mem_Ioo
#print axioms uniform_contraction_from_rate_convergence
#print axioms midpoint_lt_one
#print axioms midpoint_pos_of_pos
#print axioms super_grand_master_uniform_bounded_regime
#print axioms alphaP_pos
#print axioms lamModel_sub_lamInfty
#print axioms lamModel_residue_bound
#print axioms alphaP_times_x_minus_one_lb
#print axioms lamModel_taylor_bound
#print axioms super_grand_master_uniform_closed_form
#print axioms lamModel_eq_lambda_fp
#print axioms lamModel_eq_lambda_closed

/-! ================================================================
  §11  ConcreteIteration
================================================================ -/

#print axioms real_iterate_scaled_basin_invariant
#print axioms real_iterate_linear_from_quadratic
#print axioms pandrosion_produces_admissible_sequences
#print axioms pandrosion_concrete_super_grand_master

/-! ================================================================
  §12  MasterAbsolu
================================================================ -/

#print axioms master_absolu

/-! ================================================================
  §13  ComplexMultiplier
================================================================ -/

#print axioms Sp_C_hasDerivAt
#print axioms SpDerivC_mul_one_sub
#print axioms pandrosion_h_C_hasDerivAt
#print axioms hMultC_eq_lambdaClosedC_at_fp
#print axioms pandrosion_h_C_hasDerivAt_at_fp

/-! ================================================================
  §14  LocalAttraction
================================================================ -/

#print axioms local_contraction_iterate_converges
#print axioms local_quadratic_iterate_converges
#print axioms steffensen_local_attraction

/-! ================================================================
  §15  DynamicalConvergence
================================================================ -/

#print axioms tendsto_iterate_of_tail
#print axioms steffensen_dynamical_convergence_pointwise
#print axioms steffensen_dynamical_convergence_ae
