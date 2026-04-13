import math
import scipy.integrate as integrate
import unittest

class PoussiereContinue:
    """
    Implémentation de la Formule Continue de Pandrosion.
    Plutôt que d'utiliser des suites récursives, on utilise le "Flot de Newton" continu.
    """
    
    @staticmethod
    def formule_analytique(x: float, N_resolution: float) -> float:
        """
        Donne la poussière directement via une formule explicite.
        Pour une dimension géométrique d'étude de cube (p=3, comme Pandrosion).
        poussière = x * (1 - (1 - e^{-N})^(1/3))
        """
        # (1 - e^-N)^(1/3)
        terme = (1.0 - math.exp(-N_resolution)) ** (1.0 / 3.0)
        return x * (1.0 - terme)
        
    @staticmethod
    def formule_integrale(x: float, N_resolution: float) -> float:
        """
        Donne la poussière de Pandrosion calculée sous forme intégrale stricte :
        poussière = (x / 3) * ∫_{0}^{e^{-N}} (1 - u)^{-2/3} du
        Cette intégrale exprime géométriquement l'accumulation de matière le long de la résolution N.
        """
        def integrand(u):
            # Pour éviter la division par zéro stricte à 1 (bien que notre borne soit e^-N < 1)
            return 1.0 / ((1.0 - u) ** (2.0 / 3.0))
            
        borne_sup = math.exp(-N_resolution)
        
        resultat, erreur = integrate.quad(integrand, 0, borne_sup)
        return (x / 3.0) * resultat


class TestPoussiereContinue(unittest.TestCase):
    
    def test_equivalence_formules(self):
        """
        Teste que l'intégrale et l'analytique donnent exactement la même dimension de poussière.
        """
        x_cible = 42.0
        # À différentes résolutions continues N
        for N in [0.1, 1.0, 5.0, 10.0]:
            p_ana = PoussiereContinue.formule_analytique(x_cible, N)
            p_int = PoussiereContinue.formule_integrale(x_cible, N)
            self.assertAlmostEqual(p_ana, p_int, places=7, 
                                   msg=f"Les formules diffèrent pour N={N}")

    def test_convergence_vers_zero(self):
        """
        Vérifie qu'à une résolution infinie (N très grand), la poussière tend vers 0.
        """
        x_cible = 100.0
        p_petite_res = PoussiereContinue.formule_analytique(x_cible, 1.0)
        p_grande_res = PoussiereContinue.formule_analytique(x_cible, 20.0)
        
        # La poussière doit diminuer
        self.assertTrue(p_grande_res < p_petite_res)
        # La grande résolution doit donner une poussière très faible
        self.assertAlmostEqual(p_grande_res, 0.0, places=5)
        
        # Si on dessine ce segment, l'erreur d'observation respecte la poussière
        nombre_physique = x_cible - p_grande_res
        self.assertAlmostEqual(nombre_physique, x_cible, places=3)

if __name__ == '__main__':
    unittest.main()
