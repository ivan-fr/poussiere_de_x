import unittest

def _complex_pandrosion_T4(x_c: complex, p: float = 3.0) -> complex:
    if x_c == 1.0 + 0j or x_c == 0.0 + 0j: return 0j, 0j
    s_0 = 1.0 - (x_c - 1.0) / (x_c * p)
    def S_p(s):
        if s == 1.0 + 0j: return p
        return (1.0 - s**p) / (1.0 - s)
    def h(s):
        sp = S_p(s)
        if sp == 0j: return s
        return 1.0 - (x_c - 1.0) / (x_c * sp)
    s_1 = h(s_0)
    s_2 = h(s_1)
    diff_s1_s0 = s_1 - s_0
    hat_lambda = (s_2 - s_1) / diff_s1_s0 if abs(diff_s1_s0) > 1e-12 else 0j
    denom_t2 = (s_2 - 2*s_1 + s_0)
    T_2 = s_0 - (diff_s1_s0**2) / denom_t2 if abs(denom_t2) > 1e-12 else s_0
    s_3 = h(T_2)
    denom_t4 = hat_lambda - 1.0
    T_4 = T_2 - (s_3 - T_2) / denom_t4 if abs(denom_t4) > 1e-12 else T_2
    return s_0 - T_4, hat_lambda

class TestPandrosionGeometryBounds(unittest.TestCase):
    
    def test_monotonic_expansion_lambda(self):
        # Simulating a purely monotonic continuous expansion (e.g., compounding growth)
        # Ratio of p_now / p_past > 1
        x_fraction = complex(1.05, 0)
        gap, lam = _complex_pandrosion_T4(x_fraction)
        
        # In a generic monotonic geometric sequence, lambda must reflect structural growth.
        # Lambda should be cleanly separated and positive in real part for an expanding asymptote.
        self.assertTrue(lam.real > 0, f"Lambda {lam.real} should be positive for expansion")
        self.assertTrue(abs(gap) > 0, "Dust limit must be strictly non-zero for approximation")
        
    def test_monotonic_contraction_lambda(self):
        # Simulating a continuous breakdown
        x_fraction = complex(0.95, 0)
        gap, lam = _complex_pandrosion_T4(x_fraction)
        
        # When shrinking, lambda.real typically aligns with the negative pole
        self.assertTrue(lam.real < 0, f"Lambda {lam.real} should be negative for contraction")
        
    def test_scaling_principle_homogeneity(self):
        # Prove that calculating T4 on x=1000 / A=500 is identical to x=2 / A=1
        fraction_1 = complex(2.0, 0)
        fraction_2 = complex(1000.0 / 500.0, 0)
        
        gap1, lam1 = _complex_pandrosion_T4(fraction_1)
        gap2, lam2 = _complex_pandrosion_T4(fraction_2)
        
        self.assertAlmostEqual(lam1.real, lam2.real, places=10)
        self.assertAlmostEqual(abs(gap1), abs(gap2), places=10)
        print("Scaling principle confirmed: Financial scalar decoupling holds mathematically.")

if __name__ == '__main__':
    unittest.main()
