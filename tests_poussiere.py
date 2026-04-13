import unittest

class QuantitePoussiere:
    @staticmethod
    def pandrosion(n_iterations: int) -> float:
        """
        Calcule la poussière selon la méthode de Pandrosion.
        Suite récursive : a_{n+1} = (2 * a_n + 1) / (a_n + 2), avec a_0 = 0.
        Limite = 1. La poussière est l'écart à 1 (1 - a_n).
        """
        a = 0.0
        for _ in range(n_iterations):
            a = (2 * a + 1) / (a + 2)
        return 1.0 - a

    @staticmethod
    def thalescopique(n_iterations: int) -> float:
        """
        Calcule la poussière selon la méthode Thaléscopique.
        Somme de k=1 à n de 1/(k(k+1)) = 1 - 1/(n+1).
        La poussière est 1/(n+1).
        """
        return 1.0 / (n_iterations + 1)

    @staticmethod
    def pandrosion_racine_cubique(n_iterations: int) -> float:
        """
        Suite à convergence quadratique expliquée dans l'étude de Pandrosion 
        pour la racine cubique de 2 (Méthode de Newton / Héron).
        u_{n+1} = (2 * u_n^3 + 2) / (3 * u_n^2).
        u_0 = 1. La limite est la vraie racine cubique de 2.
        La "poussière" est alors l'écart.
        """
        u = 1.0
        for _ in range(n_iterations):
            u = (2 * (u**3) + 2) / (3 * (u**2))
        return abs((2.0**(1.0/3.0)) - u)

    @staticmethod
    def theorie_generale(x: float, granularite: int) -> float:
        """
        Pour chaque nombre x, sa poussière dépend d'une observation à la granularité N.
        La formulation générale stipule que l'imprécision inéluctable physique suit 
        l'approximation numérique (par exemple série de Taylor locale).
        Ici, poussière_generale = x * (1 / (10 ** granularite)) 
        (La poussière au niveau quantique N correspond à son résidu décimal de rang N)
        """
        # La valeur x a une incertitude proportionnelle à 10^-N.
        return abs(x) * (1.0 / (10 ** granularite))

class NombreAvecPoussiere:
    def __init__(self, valeur: float, granularity: int = 10, methode: str = "theorie_generale"):
        self.valeur = valeur
        self.granularity = granularity
        self.methode = methode
        
        if methode == "thalescopique":
            self.poussiere = QuantitePoussiere.thalescopique(granularity)
        elif methode == "pandrosion":
            self.poussiere = QuantitePoussiere.pandrosion(granularity)
        elif methode == "pandrosion_racine_cubique":
            self.poussiere = QuantitePoussiere.pandrosion_racine_cubique(granularity)
        elif methode == "theorie_generale":
            self.poussiere = QuantitePoussiere.theorie_generale(valeur, granularity)
        else:
            raise ValueError("Méthode inconnue")
            
    def resolution_physique(self) -> float:
        if self.methode == "pandrosion_racine_cubique":
             # Dans la methode cubique de pandrosion, on approxime 2^(1/3).
             # La valeur physique est la vrai valeur + la poussière absolue.
             return self.valeur + self.poussiere
        elif self.methode == "theorie_generale":
             return self.valeur + self.poussiere
        else:
            return self.valeur + (self.valeur * self.poussiere)


class TestTheoriePoussiere(unittest.TestCase):
    def test_pandrosion_converge_vers_1(self):
        self.assertAlmostEqual(QuantitePoussiere.pandrosion(100), 0.0, places=5)

    def test_pandrosion_cubique(self):
        p0 = QuantitePoussiere.pandrosion_racine_cubique(0)
        p1 = QuantitePoussiere.pandrosion_racine_cubique(1)
        self.assertTrue(p1 < p0)
        # Convergence quadratique très rapide
        self.assertAlmostEqual(QuantitePoussiere.pandrosion_racine_cubique(5), 0.0, places=7)

    def test_theorie_generale(self):
        n = NombreAvecPoussiere(42.0, granularity=3, methode="theorie_generale")
        # granularite 3 -> poussiere est 42 * 10^-3 = 0.042
        self.assertAlmostEqual(n.resolution_physique(), 42.042)

if __name__ == '__main__':
    unittest.main()
