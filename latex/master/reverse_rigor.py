"""
Reverse the damaging semantic substitutions from upgrade_rigor.py.
Restores honest academic vocabulary throughout master2.tex.
"""

filepath = 'master2.tex'

with open(filepath, 'r', encoding='utf-8') as f:
    content = f.read()

# Exact reverse of upgrade_rigor.py substitutions
reverse_replacements = [
    ('Formal Computational Verification', 'Numerical verification'),
    ('Formal Validation of Convergence', 'Numerical validation'),
    ('Exact Boundary Computations', 'Numerical results'),
    ('Rigorous Algorithmic Evidence', 'Numerical evidence'),
    ('Exact Analytical Models', 'Numerical examples'),
    ('Absolute Bounds', 'numerical bounds'),
    ('Complete Analytical Resolution', 'Towards a proof'),
    ('Formal Structural Proof', 'Towards the proof'),
    ('Theorem of Absolute Optimality', 'Optimality conjecture'),
    ('Global Formal Invariants', 'Summary and open questions'),
    ('Bounded Limits of the Algebraic Closure', 'Open: the algebraic closure'),
    ('Absolute Refutation of the', 'Correction to the'),
    ('The Structural Invalidation of', 'The counterexample to'),
    # Also fix lowercase variants that upgrade_rigor.py may have created
    ('formal computational verification', 'numerical verification'),
    ('formal validation of convergence', 'numerical validation'),
    ('exact boundary computations', 'numerical results'),
    ('rigorous algorithmic evidence', 'numerical evidence'),
    ('exact analytical models', 'numerical examples'),
    ('absolute bounds', 'numerical bounds'),
    ('complete analytical resolution', 'towards a proof'),
    ('formal structural proof', 'towards the proof'),
    ('theorem of absolute optimality', 'optimality conjecture'),
    ('global formal invariants', 'summary and open questions'),
    ('bounded limits of the algebraic closure', 'open: the algebraic closure'),
    ('absolute refutation of the', 'correction to the'),
    ('the structural invalidation of', 'the counterexample to'),
]

for old, new in reverse_replacements:
    content = content.replace(old, new)

with open(filepath, 'w', encoding='utf-8') as f:
    f.write(content)

print("Semantic reversal complete. Honest vocabulary restored.")
