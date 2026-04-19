import os
import shutil
import re

repodir = '/Users/ivanbesevic/Documents/poussiere'
latex_dir = os.path.join(repodir, 'latex')
figures_dir = os.path.join(repodir, 'figures')
verification_dir = os.path.join(repodir, 'verification')
scripts_dir = os.path.join(repodir, 'scripts')

# Ensure target directories exist
for d in [figures_dir, verification_dir, scripts_dir]:
    os.makedirs(d, exist_ok=True)

print("1. Consolidating figures...")
# Move all remaining images to the root figures/ folder
for root, dirs, files in os.walk(latex_dir):
    for f in files:
        if f.endswith('.png') or f.endswith('.pdf'):
            if root != figures_dir: # avoid moving from figures to figures
                src = os.path.join(root, f)
                dst = os.path.join(figures_dir, f)
                if not os.path.exists(dst):
                    print(f"  Moving {f} to figures/")
                    shutil.move(src, dst)

print("\n2. Consolidating LaTeX...")
old_master = os.path.join(latex_dir, 'master', 'master2.tex')
new_master = os.path.join(latex_dir, 'pandrosion_master.tex')
if os.path.exists(old_master):
    print("  Moving master2.tex to pandrosion_master.tex")
    shutil.move(old_master, new_master)

print("\n3. Updating LaTeX image paths...")
if os.path.exists(new_master):
    with open(new_master, 'r') as f:
        content = f.read()

    # Change all graphicspath combinations to just {{../figures/}}
    content = re.sub(r'\\graphicspath\{\{.*\}\}', r'\\graphicspath{{../figures/}}', content)

    # Change absolute or subdirectory image imports to just the filename
    # e.g. \includegraphics{figures_conclusion/spec_spectral.png} -> \includegraphics{spec_spectral.png}
    def replace_includegraphics(match):
        options = match.group(1) if match.group(1) else ""
        path = match.group(2)
        basename = os.path.basename(path)
        return f"\\includegraphics{options}{{{basename}}}"

    content = re.sub(r'\\includegraphics(\[.*?\])?\{(.+?)\}', replace_includegraphics, content)

    with open(new_master, 'w') as f:
        f.write(content)
    print("  Updated all \\includegraphics and \\graphicspath")

print("\n4. Distributing Python scripts...")
# Look for any .py files in latex_dir
for root, dirs, files in os.walk(latex_dir):
    for f in files:
        if f.endswith('.py'):
            src = os.path.join(root, f)
            # Strategy: if it plots/generates stuff, goes to scripts. 
            # If it mathematically proves or scans, goes to verification.
            if f.startswith('generate_') or f.startswith('plot_') or f.startswith('figures_'):
                target = os.path.join(scripts_dir, f)
            else:
                target = os.path.join(verification_dir, f)
            
            # Don't overwrite existing
            if not os.path.exists(target):
                print(f"  Moving {f} to {os.path.basename(os.path.dirname(target))}/")
                shutil.move(src, target)

print("\n5. Cleaning up redundant LaTeX directories & files...")
# Delete all remaining .tex files natively inside latex/ EXCEPT pandrosion_master.tex
for f in os.listdir(latex_dir):
    item = os.path.join(latex_dir, f)
    if os.path.isfile(item) and item.endswith('.tex') and f != 'pandrosion_master.tex':
        print(f"  Deleting redundant file: {f}")
        os.remove(item)

# Delete subdirectories of latex!
for f in os.listdir(latex_dir):
    item = os.path.join(latex_dir, f)
    if os.path.isdir(item):
        print(f"  Deleting directory: latex/{f}")
        shutil.rmtree(item)

print("\nRestructuring complete!")
