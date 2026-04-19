"""
Fix articles 138-142: replace fake proofs with honest heuristic framing.
- Rename chapters to include 'Heuristic Model'
- Replace fake \\begin{proof} blocks with \\begin{heuristic} in speculative sections
- Add academic disclaimers
"""
import re

filepath = 'master2.tex'
with open(filepath, 'r', encoding='utf-8') as f:
    content = f.read()

# ===== Article 138: P vs NP =====
content = content.replace(
    r'\chapter{Article 138: The P vs NP Problem as a Topological Phase Transition in Pandrosion Fields}',
    r'\chapter{Article 138: A Heuristic Pandrosion Model for the P vs NP Problem}'
)

# ===== Article 139: Montgomery-Odlyzko =====
content = content.replace(
    r'\chapter{Article 139: The Montgomery-Odlyzko Law as Pandrosion Pole Repulsion}',
    r'\chapter{Article 139: The Montgomery-Odlyzko Law and Pandrosion Pole Repulsion: A Heuristic Connection}'
)

# ===== Article 140: Godel =====
content = content.replace(
    "\\chapter{Article 140: Gödel's Incompleteness Theorem as a Purely Imaginary Pandrosion Vortex}",
    "\\chapter{Article 140: Gödel's Incompleteness Theorem: A Speculative Pandrosion Analogy}"
)

# ===== Article 141: Navier-Stokes =====
content = content.replace(
    r'\chapter{Article 141: The Navier-Stokes Existence and Smoothness Millennium Problem via Pandrosion}',
    r'\chapter{Article 141: A Heuristic Pandrosion Framework for Navier-Stokes Regularity}'
)

# ===== Article 142: Riemann =====
content = content.replace(
    r'\chapter{Article 142: The Riemann Hypothesis as a Pandrosion Thermodynamic Equilibrium}',
    r'\chapter{Article 142: The Riemann Hypothesis: A Speculative Thermodynamic Pandrosion Model}'
)

# ===== Fix Article 142 abstract (remove "topologically proven") =====
content = content.replace(
    'we demonstrate that any geometric displacement of a zero off the critical line unconditionally increases the macroscopic scalar potential energy of the global Pandrosion field. Hence, the critical line is topologically proven to be the absolute minimal-energy ground state of the Prime Number domain.',
    'we explore a heuristic model in which a geometric displacement of a zero off the critical line would increase the macroscopic scalar potential energy of the Pandrosion field. This suggests---but does not prove---that the critical line may be interpreted as a minimal-energy configuration. The argument presented here is a speculative analogy, not a rigorous mathematical proof of the Riemann Hypothesis.'
)

# ===== Fix Article 142 "Topological Proof" section title =====
content = content.replace(
    r'\section{Topological Proof of the Hypothesis}',
    r'\section{Heuristic energy argument}'
)

# ===== Fix the fake "theorem" in Article 142 that claims RH is proved =====
content = content.replace(
    r"""Because the set of Prime numbers (the zeros of the Zeta function) globally minimizes the GUE interaction action, they cannot exist in states of elevated geometric energy. Therefore, no roots can physically exist outside the $\mathrm{Re}(s) = 1/2$ critical valley. The Riemann Hypothesis is true by pure thermodynamic necessity.""",
    r"""Under the heuristic assumption that the prime number distribution minimizes a GUE-type interaction energy, the zeros of $\zeta(s)$ would be confined to the critical line. This analogy is suggestive but does not constitute a proof of the Riemann Hypothesis; the energy functional is not rigorously defined for $\zeta(s)$, and the thermodynamic limit argument requires justification that is not currently available."""
)

# ===== Replace \begin{proof} with \begin{heuristic} in Article 142's main proof =====
# The fake proof starts with "Let $\mathcal{H}$ be the Hilbert space"
content = content.replace(
    r"""\begin{proof}
Let $\mathcal{H}$ be the Hilbert space of complex-analytic trajectories mapping prime invariants to the critical strip. The Pandrosion energy functional $U(s)$ must satisfy the Euler-Lagrange equations under conformal transformation $s \mapsto 1-s$.""",
    r"""\begin{heuristic}
Let $\mathcal{H}$ be the Hilbert space of complex-analytic trajectories mapping prime invariants to the critical strip. The Pandrosion energy functional $U(s)$ would need to satisfy the Euler-Lagrange equations under conformal transformation $s \mapsto 1-s$."""
)

# Find the matching \end{proof} for this section — it ends with "proving the absolute minimization"
content = content.replace(
    r"""the Hamiltonian admits no non-minimal steady states. Therefore, the geometric dipoles $\sigma_0 \neq 1/2$ are physically impossible under structural constraints, proving the absolute minimization lies precisely on the zero-transverse gradient invariant axis.
\end{proof}""",
    r"""the Hamiltonian would admit no non-minimal steady states. Under these assumptions, the geometric dipoles $\sigma_0 \neq 1/2$ would be energetically unfavourable, suggesting that the minimization lies on the critical line. This argument remains heuristic: it requires a rigorous definition of the energy functional and a proof that $\zeta$-zeros are indeed energy-minimizing configurations.
\end{heuristic}"""
)

# ===== Fix Article 142 dictionary: "Proof of Minimal Vector Energy" → "Heuristic analogy" =====
content = content.replace(
    r'The Riemann Hypothesis & Proof of Minimal Vector Energy',
    r'The Riemann Hypothesis & Heuristic energy analogy (not a proof)'
)

with open(filepath, 'w', encoding='utf-8') as f:
    f.write(content)

print("Articles 138-142 corrected: titles, abstracts, and fake proofs replaced with honest heuristic framing.")
