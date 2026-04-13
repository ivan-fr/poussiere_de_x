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
    pdf.cell(0, 10, "Pandrosion, Sommes Thalescopiques et Flot de Newton Continu", align='C', new_x=XPos.LMARGIN, new_y=YPos.NEXT)
    pdf.ln(10)

    # 1. Methode Thalescopique
    pdf.set_font("helvetica", "B", 14)
    pdf.cell(0, 10, "1. La Methode Thalescopique", new_x=XPos.LMARGIN, new_y=YPos.NEXT)
    pdf.set_font("helvetica", "", 12)
    sec1 = "Le nombre 1 peut etre vu comme la somme infinie de decoupes grace au theoreme de Thales. Si notre resolution s'arrete à l'etape N, on observe la poussiere 'Thalescopique' :"
    pdf.multi_cell(0, 8, text=sec1)
    
    y1 = pdf.get_y()
    create_math_image(r"1 = \sum_{k=1}^{\infty} \left( \frac{1}{k} - \frac{1}{k+1} \right) \Rightarrow \varepsilon_T(N) = \frac{1}{N+1}", "eq1.png", fontsize=20)
    pdf.image("eq1.png", x=20, y=y1+2, w=150)
    pdf.set_y(y1 + 25)

    # 2. Pandrosion - Discret
    pdf.set_font("helvetica", "B", 14)
    pdf.cell(0, 10, "2. Poussiere Geometrique Discrete (Pandrosion & Newton)", new_x=XPos.LMARGIN, new_y=YPos.NEXT)
    pdf.set_font("helvetica", "", 12)
    sec2 = "Historiquement, la mathematicienne Pandrosion chercha a approcher la duplication du cube. L'equilibre geometrique de ce processus conduit à la suite de convergence quadratique definie par :"
    pdf.multi_cell(0, 8, text=sec2)
    
    y2 = pdf.get_y()
    create_math_image(r"u_{n+1} = \frac{2 u_n^3 + 2}{3 u_n^2} \Rightarrow \varepsilon_{\sqrt[3]{2}}(N) = \left| \sqrt[3]{2} - u_N \right|", "eq2.png", fontsize=20)
    pdf.image("eq2.png", x=20, y=y2+2, w=150)
    pdf.set_y(y2 + 25)

    # 3. Flot Continu de Pandrosion
    pdf.set_font("helvetica", "B", 14)
    pdf.cell(0, 10, "3. La Generalisation Integrale (Flot de Newton Continu)", new_x=XPos.LMARGIN, new_y=YPos.NEXT)
    pdf.set_font("helvetica", "", 12)
    sec3 = "La poussiere quantique n'est pas uniquement un 'saut'. En transformant la recurrence de Pandrosion en son flot differentiel, nous trouvons que pour un nombre purement objectif x et un degre geometrique p=3, la poussiere analytique restante a la resolution N est exactement :"
    pdf.multi_cell(0, 8, text=sec3)
    
    y3 = pdf.get_y()
    create_math_image(r"\varepsilon_{Pand}(x, N) = x \left( 1 - \sqrt[3]{1 - e^{-N}} \right)", "eq3.png", fontsize=24)
    pdf.image("eq3.png", x=50, y=y3+2, w=100)
    pdf.set_y(y3 + 22)
    
    sec4 = "Ce meme flot peut merveilleusement s'ecrire de maniere integrale pure pour repondre a la vision d'une addition infinie de la 'scorie' : "
    pdf.multi_cell(0, 8, text=sec4)
    
    y4 = pdf.get_y()
    create_math_image(r"\varepsilon_{Pand}(x, N) = \frac{x}{3} \int_{0}^{e^{-N}} \frac{1}{(1-u)^{2/3}} du", "eq4.png", fontsize=24)
    pdf.image("eq4.png", x=50, y=y4+2, w=100)
    pdf.set_y(y4 + 28)

    # 4. Formule d'Imprecision Finale
    pdf.set_font("helvetica", "B", 14)
    pdf.cell(0, 10, "4. Theoreme Final de la Poussiere Universelle", new_x=XPos.LMARGIN, new_y=YPos.NEXT)
    pdf.set_font("helvetica", "", 12)
    sec5 = "Tout point representant un nombre physique x entraine une incertitude liee a la fraction d'ombre dictee par son environnement. Dans son expression la plus empirique quantique (resolution N en base 10), cela est decrit par "
    pdf.multi_cell(0, 8, text=sec5)
    
    y5 = pdf.get_y()
    create_math_image(r"\mathcal{D}(x, N) = x \cdot 10^{-N}", "eq5.png", fontsize=28)
    pdf.image("eq5.png", x=70, y=y5+2, w=50)
    pdf.set_y(y5 + 25)

    sec6 = "Ainsi, l'espace physique n'est jamais vide autour d'une valeur exacte : la poussiere de Pandrosion l'enrobe comme le reflet d'une dimension inachevee."
    pdf.multi_cell(0, 8, text=sec6)
    
    pdf_path = "/Users/ivanbesevic/Documents/poussiere/Article_Final_Theorie_Poussiere.pdf"
    pdf.output(pdf_path)
    print(f"Article g\u00e9n\u00e9r\u00e9 avec succ\u00e8s : {pdf_path}")
    
    # Nettoyage
    for f in ["eq1.png", "eq2.png", "eq3.png", "eq4.png", "eq5.png"]:
        if os.path.exists(f):
            os.remove(f)

if __name__ == "__main__":
    generer_article()
