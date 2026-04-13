import os
from fpdf import FPDF

class PDFArticle(FPDF):
    def header(self):
        self.set_font('helvetica', 'B', 12)
        # Title
        self.cell(0, 10, 'La Theorie de la Poussiere Mathematique', border=0, align='C')
        self.ln(20)
        
    def footer(self):
        self.set_y(-15)
        self.set_font('helvetica', 'I', 8)
        # Page number
        self.cell(0, 10, f'Page {self.page_no()}', 0, 0, 'C')

def generer_article():
    pdf = PDFArticle()
    pdf.add_page()
    
    # Title
    pdf.set_font("helvetica", "B", 18)
    pdf.cell(0, 10, "Une Notion de Poussiere pour Chaque Nombre", ln=True, align='C')
    pdf.set_font("helvetica", "I", 12)
    pdf.cell(0, 10, "Application de la Methode de Pandrosion et de la Somme Thalescopique", ln=True, align='C')
    pdf.ln(10)
    
    # Abstract
    pdf.set_font("helvetica", "B", 14)
    pdf.cell(0, 10, "Resume", ln=True)
    pdf.set_font("helvetica", "", 12)
    abstract = (
        "Cet article explore une idee fondamentale : les entites mathematiques pures (comme le point, ou le nombre 1) "
        "possedent, dans le monde physique, une 'granularite' ineluctable. Par exemple, dessiner un segment de 1 cm "
        "n'est jamais exact : il mesure '1 + une petite poussiere'. Nous formalisons cette intuition a l'aide de deux "
        "methodes d'approximation."
    )
    pdf.multi_cell(0, 8, abstract)
    pdf.ln(5)
    
    # Section 1
    pdf.set_font("helvetica", "B", 14)
    pdf.cell(0, 10, "1. La Methode Thalescopique", ln=True)
    pdf.set_font("helvetica", "", 12)
    sec1 = (
        "Inspirée par la décomposition projective, nous définissons la poussière à travers une somme "
        "télescopique (rebaptisée Thaléscopique en hommage au théorème de Thalès). "
        "Le nombre 1 peut être vu comme la somme infinie de découpes toujours plus fines :\n\n"
        "   1 = Somme(de k=1 a l'infini) [ 1/k - 1/(k+1) ]\n\n"
        "Si notre résolution spatiale s'arrête à l'itération N, la somme partielle nous donne :\n"
        "   T_N = 1 - 1/(N+1)\n\n"
        "Ainsi, la 'Poussiere' definissant l'ecart à l'idealite est epsilon_t(N) = 1/(N+1). "
        "Chaque nombre physique x porte une granularite proportionnelle : P(x) = x + x * epsilon_t(N)."
    )
    pdf.multi_cell(0, 8, txt=sec1.encode("latin-1", "replace").decode("latin-1"))
    pdf.ln(5)

    # Section 2
    pdf.set_font("helvetica", "B", 14)
    pdf.cell(0, 10, "2. La Methode de Pandrosion pour le nombre 1", ln=True)
    pdf.set_font("helvetica", "", 12)
    sec2 = (
        "La mathématicienne classique Pandrosion développa des méthodes itératives pour approcher des nombres continus. "
        "Au lieu de chercher la racine cubique, reprenons son idée pour approcher le nombre 1 lui-meme par la "
        "moyenne harmonique-arithmétique. Nous definissons la suite :\n\n"
        "   a_{0} = 0\n"
        "   a_{n+1} = (2 * a_n + 1) / (a_n + 2)\n\n"
        "La limite de cette suite infinie est exactement 1. Dans un monde de resolution finie N, "
        "l'approche s'arrete à l'etape N. La Poussiere est alors la difference entre la perfection (1) "
        "et la realite atteinte (a_N) :\n\n"
        "   epsilon_p(N) = 1 - a_N\n\n"
        "Ce epsilon est non-lineaire, decroissant de plus en plus lentement, refletant mieux "
        "la lutte physique contre les dernieres impuretes de la matiere (la poussiere ireductible)."
    )
    pdf.multi_cell(0, 8, txt=sec2.encode("latin-1", "replace").decode("latin-1"))
    pdf.ln(5)
    
    # Section 3
    pdf.set_font("helvetica", "B", 14)
    pdf.cell(0, 10, "3. Definition Formelle de la Poussiere", ln=True)
    pdf.set_font("helvetica", "", 12)
    sec3 = (
        "Soit un nombre reel ideal x. Sa projection dans le monde denombrable a une "
        "resolution N est note P(x, N).\n"
        "   P(x, N) = x * (1 + epsilon(N))\n\n"
        "L'ajout d'une telle theorie est fondamental. Quand un ingenieur prevoit 1 cm, "
        "la theorie mathematique de la poussiere assure que sa mesure contiendra "
        "exactement la quantite d'imprecision liee a l'observation du quantique et du continu.\n"
        "Des tests unitaires Python implementant cette theorie ont ete codes (""tests_poussiere.py"") "
        "pour prouver cette theorie."
    )
    pdf.multi_cell(0, 8, txt=sec3.encode("latin-1", "replace").decode("latin-1"))
    
    pdf_path = "/Users/ivanbesevic/Documents/poussiere/Theorie_de_la_Poussiere.pdf"
    pdf.output(pdf_path)
    print(f"Article g\u00e9n\u00e9r\u00e9 avec succ\u00e8s : {pdf_path}")

if __name__ == "__main__":
    generer_article()
