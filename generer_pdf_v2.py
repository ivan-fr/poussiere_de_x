import os
import matplotlib.pyplot as plt
from fpdf import FPDF
from fpdf.enums import XPos, YPos

def create_math_image(formula, filename, fontsize=24):
    """ Rend la formule mathématique en PNG avec fond transparent grace à Matplotlib """
    fig = plt.figure(figsize=(8, 1))
    fig.text(0.5, 0.5, f"${formula}$", fontsize=fontsize, ha='center', va='center')
    plt.axis('off')
    plt.savefig(filename, dpi=300, bbox_inches='tight', transparent=True)
    plt.close()

class PDFArticle(FPDF):
    def header(self):
        self.set_font('helvetica', 'B', 12)
        self.cell(0, 10, 'La Theorie de la Poussiere Mathematique', border=0, align='C', new_x=XPos.LMARGIN, new_y=YPos.NEXT)
        self.ln(5)
        
    def footer(self):
        self.set_y(-15)
        self.set_font('helvetica', 'I', 8)
        self.cell(0, 10, f'Page {self.page_no()}', 0, 0, 'C')

def generer_article():
    pdf = PDFArticle()
    pdf.add_page()
    
    # Titre
    pdf.set_font("helvetica", "B", 18)
    pdf.cell(0, 10, "Une Notion de Poussiere pour Chaque Nombre", align='C', new_x=XPos.LMARGIN, new_y=YPos.NEXT)
    pdf.set_font("helvetica", "I", 12)
    pdf.cell(0, 10, "Racine Cubique, Methode de Pandrosion et Generalisation", align='C', new_x=XPos.LMARGIN, new_y=YPos.NEXT)
    pdf.ln(10)

    # 1. Methode Thalescopique
    pdf.set_font("helvetica", "B", 14)
    pdf.cell(0, 10, "1. La Methode Thalescopique", new_x=XPos.LMARGIN, new_y=YPos.NEXT)
    pdf.set_font("helvetica", "", 12)
    sec1 = "Le nombre 1 peut etre vu comme la somme infinie de decoupes projetees grace au theoreme de Thales :"
    pdf.multi_cell(0, 8, text=sec1)
    
    # Image Math 1
    y1 = pdf.get_y()
    create_math_image(r"1 = \sum_{k=1}^{\infty} \left( \frac{1}{k} - \frac{1}{k+1} \right)", "eq1.png")
    pdf.image("eq1.png", x=50, y=y1+2, w=100)
    pdf.set_y(y1 + 25)
    
    pdf.multi_cell(0, 8, text="Si notre resolution s'arrete à l'etape N, on obtient la poussiere Thalescopique :")
    y2 = pdf.get_y()
    create_math_image(r"\varepsilon_T(N) = \frac{1}{N+1}", "eq2.png", fontsize=28)
    pdf.image("eq2.png", x=70, y=y2+2, w=50)
    pdf.set_y(y2 + 20)

    # 2. Methode de Pandrosion
    pdf.set_font("helvetica", "B", 14)
    pdf.cell(0, 10, "2. Poussiere de la Racine Cubique (Pandrosion Quadratique)", new_x=XPos.LMARGIN, new_y=YPos.NEXT)
    pdf.set_font("helvetica", "", 12)
    sec2 = "Pour doubler le cube comme Pandrosion, nous utilisons l'approximation de la racine cubique de 2. La suite a convergence quadratique s'ecrit :"
    pdf.multi_cell(0, 8, text=sec2)
    
    # Image Math 3
    y3 = pdf.get_y()
    create_math_image(r"u_{n+1} = \frac{2 u_n^3 + 2}{3 u_n^2}", "eq3.png", fontsize=28)
    pdf.image("eq3.png", x=70, y=y3+2, w=50)
    pdf.set_y(y3 + 25)
    
    sec3 = "La limite etant cbrt(2), la quantite de 'poussiere' accrochee au nombre au stade physique decroit quadratiquement :"
    pdf.multi_cell(0, 8, text=sec3)
    
    y4 = pdf.get_y()
    create_math_image(r"\varepsilon(N) = \left| \sqrt[3]{2} - u_N \right|", "eq4.png")
    pdf.image("eq4.png", x=60, y=y4+2, w=70)
    pdf.set_y(y4 + 20)

    # 3. Formule de Generalisation
    pdf.set_font("helvetica", "B", 14)
    pdf.cell(0, 10, "3. Formule Universelle de la Poussiere", new_x=XPos.LMARGIN, new_y=YPos.NEXT)
    pdf.set_font("helvetica", "", 12)
    sec4 = "On peut maintenant generaliser ! Chaque nombre physique x possede une granularite dependante de la resolution de l'observation macroscopique N. La quantite d'imprecision de n'importe quel point numerique est decrite par son developpement en base 10 vu le theoreme d'incompletude physique :"
    pdf.multi_cell(0, 8, text=sec4)
    
    y5 = pdf.get_y()
    create_math_image(r"\mathcal{D}(x, N) = x \cdot 10^{-N}", "eq5.png", fontsize=28)
    pdf.image("eq5.png", x=70, y=y5+2, w=50)
    pdf.set_y(y5 + 25)
    
    sec5 = "Ainsi, si vous tracez un point de largeur theorethique 'x', il emportera toujours avec lui une poussiere de fond, resolvant du meme coup les mysteres de l'erreur d'arrondi de notre univers quantique.\n(Voir 'tests_poussiere.py' pour implémentation)."
    pdf.multi_cell(0, 8, text=sec5)
    
    pdf_path = "/Users/ivanbesevic/Documents/poussiere/Theorie_de_la_Poussiere.pdf"
    pdf.output(pdf_path)
    print(f"Article g\u00e9n\u00e9r\u00e9 avec succ\u00e8s et formules Latex : {pdf_path}")
    
    # Nettoyage des images
    for f in ["eq1.png", "eq2.png", "eq3.png", "eq4.png", "eq5.png"]:
        if os.path.exists(f):
            os.remove(f)

if __name__ == "__main__":
    generer_article()
