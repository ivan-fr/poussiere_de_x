/-
  Universitas Pandrosion — Axiom audit.

  This file exists so the CI axiom-audit pipeline (`scripts/strict_check.sh`)
  can enumerate the trusted base of every theorem in the corpus by
  running `#print axioms <name>` for each theorem / lemma declared in
  `Pandrosion/**/*.lean`.

  The CI script autogenerates the body of this file by:

    1. Scanning every `theorem` / `lemma` in `Pandrosion/**/*.lean`.
    2. Replacing everything between the markers

           `-- BEGIN GENERATED AXIOM AUDIT`
           `-- END GENERATED AXIOM AUDIT`

       with one `#print axioms Pandrosion.<name>` line per declaration.
    3. Compiling the result with `lake env lean CheckAllAxioms.lean` and
       checking that no theorem depends on any axiom outside the
       classical whitelist `{propext, Classical.choice, Quot.sound}`.

  The hand-curated placeholder block below is also a working audit (it
  lists every theorem in dependency order), so this file is useful on
  its own when the CI script is not invoked. Running it produces one
  `'foo' depends on axioms: [...]` line per theorem.
-/

import Pandrosion

-- BEGIN GENERATED AXIOM AUDIT

/-! ================================================================
  §1  Foundations
================================================================ -/

#print axioms Pandrosion.Sp_mul_one_sub
#print axioms Pandrosion.Sp_at_one
#print axioms Pandrosion.Sp_pos
#print axioms Pandrosion.pandrosion_fixed_point_iff
#print axioms Pandrosion.pandrosion_quasi_newton
#print axioms Pandrosion.scaling_power
#print axioms Pandrosion.scaling_factorization
#print axioms Pandrosion.pandrosion_s0_opt_eq_h_one
#print axioms Pandrosion.offset_monotone
#print axioms Pandrosion.scaled_s0_opt_closer_to_one
#print axioms Pandrosion.double_preconditioning
#print axioms Pandrosion.steffensen_step_fixed_point_of_h
#print axioms Pandrosion.steffensen_step_fixed_point
#print axioms Pandrosion.steffensen_rate_ne_one
#print axioms Pandrosion.Sp_C_mul_one_sub
#print axioms Pandrosion.Sp_C_at_one
#print axioms Pandrosion.pandrosion_h_C_quasi_newton
#print axioms Pandrosion.pandrosion_h_C_fixed_point_forward
#print axioms Pandrosion.pandrosion_h_C_fixed_point_reverse
#print axioms Pandrosion.pandrosion_h_C_fixed_point_iff
#print axioms Pandrosion.steffensen_step_C_fixed_point_of_h
#print axioms Pandrosion.steffensen_step_C_fixed_point
#print axioms Pandrosion.scaling_power_C
#print axioms Pandrosion.scaling_factorization_C
#print axioms Pandrosion.scaling_steffensen_C
#print axioms Pandrosion.pandrosion_s0_opt_C_eq_h_one
#print axioms Pandrosion.optimal_scaling_steffensen_C
#print axioms Pandrosion.cauchy_radius_ge_two
#print axioms Pandrosion.cauchy_radius_pos
#print axioms Pandrosion.voronoi_nearest_exists
#print axioms Pandrosion.multi_start_grand_master
#print axioms Pandrosion.voronoi_unique_off_boundary
#print axioms Pandrosion.Sp_gt_p_of_gt_one
#print axioms Pandrosion.lambda_closed_pos
#print axioms Pandrosion.lambda_closed_lt_one
#print axioms Pandrosion.lambda_closed_at_one
#print axioms Pandrosion.lambda_closed_lt_uniform_bound
#print axioms Pandrosion.one_sub_lambda_closed_ne_zero
#print axioms Pandrosion.steffensen_constant_denom_ne_zero
#print axioms Pandrosion.K_steffensen_nonneg

/-! ================================================================
  §2  MultiStartBasins
================================================================ -/

#print axioms Pandrosion.fin_injective_iff_bijective
#print axioms Pandrosion.multi_start_sigma_bijective
#print axioms Pandrosion.generic_convergence_conditional
#print axioms Pandrosion.basin_eq_iInter
#print axioms Pandrosion.halfPlane_convex
#print axioms Pandrosion.basin_convex
#print axioms Pandrosion.basin_self_mem
#print axioms Pandrosion.basin_nonempty
#print axioms Pandrosion.basin_isConnected
#print axioms Pandrosion.basin_isPreconnected
#print axioms Pandrosion.basin_isClosed
#print axioms Pandrosion.basin_frontier_on_bisector
#print axioms Pandrosion.basin_boundary_finite_union
#print axioms Pandrosion.basin_anchor_mem_interior
#print axioms Pandrosion.basin_interior_nonempty
#print axioms Pandrosion.basins_cover
#print axioms Pandrosion.VoronoiBoundary_subset_bisectors
#print axioms Pandrosion.VoronoiBoundary_eq_bisectors_filtered
#print axioms Pandrosion.constructive_mcmullen

/-! ================================================================
  §3  QuadraticRate
================================================================ -/

#print axioms Pandrosion.pandrosion_h_sub_mul
#print axioms Pandrosion.Qpoly_mul_sub
#print axioms Pandrosion.pow_sub_pow_factor
#print axioms Pandrosion.Qpoly_diag
#print axioms Pandrosion.pandrosion_h_sub_fp_factor
#print axioms Pandrosion.pandrosion_h_sub_fp_quotient
#print axioms Pandrosion.lambda_fp_eq_quotient
#print axioms Pandrosion.pandrosion_linear_rate_ratio
#print axioms Pandrosion.pandrosion_ratio_at_fp
#print axioms Pandrosion.quadratic_deviation_at_fp
#print axioms Pandrosion.Sp_sub_Sp_factor
#print axioms Pandrosion.Qpoly_sub_Qpoly_factor
#print axioms Pandrosion.lambda_fp_mul_Sp
#print axioms Pandrosion.pandrosion_h_sub_fp_factor_no_x
#print axioms Pandrosion.pandrosion_h_quadratic_polynomial
#print axioms Pandrosion.pandrosion_h_quadratic_quotient
#print axioms Pandrosion.pandrosion_h_quadratic_bound_of
#print axioms Pandrosion.Sp_reflect
#print axioms Pandrosion.lambda_fp_eq_lambda_closed

/-! ================================================================
  §4  QuadraticComplexity
================================================================ -/

#print axioms Pandrosion.quadratic_tower_bound
#print axioms Pandrosion.quadratic_loglog_complexity
#print axioms Pandrosion.pandrosion_steffensen_loglog_complexity

/-! ================================================================
  §5  VoronoiMeasure
================================================================ -/

#print axioms Pandrosion.perpBisectorForm_apply
#print axioms Pandrosion.perpBisectorForm_ne_zero_of_ne
#print axioms Pandrosion.perpBisector_volume_zero
#print axioms Pandrosion.voronoiBoundary_volume_zero
#print axioms Pandrosion.constructive_mcmullen_ae
#print axioms Pandrosion.voronoi_selector_unique_ae
#print axioms Pandrosion.constructive_mcmullen_ae_proper

/-! ================================================================
  §6  CyclotomicMcMullen
================================================================ -/

#print axioms Pandrosion.cycAnchor_pow
#print axioms Pandrosion.cycAnchor_injective
#print axioms Pandrosion.unconditional_mcmullen

/-! ================================================================
  §7  KungTraub
================================================================ -/

#print axioms Pandrosion.efficiencyIndex_two_two
#print axioms Pandrosion.pandrosion_const_lt_newton
#print axioms Pandrosion.pandrosion_const_pos
#print axioms Pandrosion.pandrosion_attains_kung_traub
#print axioms Pandrosion.pandrosion_kung_traub_optimal
#print axioms Pandrosion.pandrosion_kung_traub_optimal_bound

/-! ================================================================
  §8  UniformComplexity
================================================================ -/

#print axioms Pandrosion.uniform_quadratic_loglog
#print axioms Pandrosion.uniform_quadratic_loglog_scaled
#print axioms Pandrosion.pandrosion_uniform_loglog

/-! ================================================================
  §9  SuperGrandMaster
================================================================ -/

#print axioms Pandrosion.linear_decay_bound
#print axioms Pandrosion.scaled_basin_entry_time
#print axioms Pandrosion.two_phase_complexity
#print axioms Pandrosion.super_grand_master_uniform

/-! ================================================================
  §10  UniformContractionRate
================================================================ -/

#print axioms Pandrosion.lamInfty_lt_one_of_one_lt
#print axioms Pandrosion.lamInfty_pos_of_one_lt
#print axioms Pandrosion.lamInfty_mem_Ioo
#print axioms Pandrosion.uniform_contraction_from_rate_convergence
#print axioms Pandrosion.midpoint_lt_one
#print axioms Pandrosion.midpoint_pos_of_pos
#print axioms Pandrosion.super_grand_master_uniform_bounded_regime
#print axioms Pandrosion.alphaP_pos
#print axioms Pandrosion.lamModel_sub_lamInfty
#print axioms Pandrosion.lamModel_residue_bound
#print axioms Pandrosion.alphaP_times_x_minus_one_lb
#print axioms Pandrosion.lamModel_taylor_bound
#print axioms Pandrosion.super_grand_master_uniform_closed_form
#print axioms Pandrosion.lamModel_eq_lambda_fp
#print axioms Pandrosion.lamModel_eq_lambda_closed

/-! ================================================================
  §11  ConcreteIteration
================================================================ -/

#print axioms Pandrosion.real_iterate_scaled_basin_invariant
#print axioms Pandrosion.real_iterate_linear_from_quadratic
#print axioms Pandrosion.pandrosion_produces_admissible_sequences
#print axioms Pandrosion.pandrosion_concrete_super_grand_master

/-! ================================================================
  §12  MasterAbsolu
================================================================ -/

#print axioms Pandrosion.master_absolu

/-! ================================================================
  §13  ComplexMultiplier
================================================================ -/

#print axioms Pandrosion.Sp_C_hasDerivAt
#print axioms Pandrosion.SpDerivC_mul_one_sub
#print axioms Pandrosion.pandrosion_h_C_hasDerivAt
#print axioms Pandrosion.hMultC_eq_lambdaClosedC_at_fp
#print axioms Pandrosion.pandrosion_h_C_hasDerivAt_at_fp

/-! ================================================================
  §14  LocalAttraction
================================================================ -/

#print axioms Pandrosion.local_contraction_iterate_converges
#print axioms Pandrosion.local_quadratic_iterate_converges
#print axioms Pandrosion.steffensen_local_attraction

/-! ================================================================
  §15  DynamicalConvergence
================================================================ -/

#print axioms Pandrosion.tendsto_iterate_of_tail
#print axioms Pandrosion.steffensen_dynamical_convergence_pointwise
#print axioms Pandrosion.steffensen_dynamical_convergence_ae

/-! ================================================================
  §16  LinearAsymptotics
================================================================ -/

#print axioms Pandrosion.lambdaClosedC_ne_one_at_fp
#print axioms Pandrosion.cycAnchor_ne_zero
#print axioms Pandrosion.cycAnchor_lambdaClosedC_ne_one
#print axioms Pandrosion.pandrosion_h_C_locally_lipschitz_at_fp
#print axioms Pandrosion.pandrosion_h_C_strict_contraction_when_abs_lt_one

/-! ================================================================
  §17  QuadraticBound
================================================================ -/

#print axioms Pandrosion.pandrosion_h_C_factorization
#print axioms Pandrosion.hFactorAux_at_alpha
#print axioms Pandrosion.hFactor_at_alpha
#print axioms Pandrosion.Sp_C_analyticAt
#print axioms Pandrosion.hFactorAux_analyticAt
#print axioms Pandrosion.hFactor_analyticAt
#print axioms Pandrosion.pandrosion_h_C_analyticAt
#print axioms Pandrosion.hFactor_isBigO_sub_at_fp
#print axioms Pandrosion.pandrosion_h_C_quadratic_bound
#print axioms Pandrosion.pandrosion_h_C_quadratic_bound_with_small

-- END GENERATED AXIOM AUDIT
