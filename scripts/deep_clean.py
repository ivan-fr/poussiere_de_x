import os
import re

# Repo root
repodir = '/Users/ivanbesevic/Documents/poussiere'
master2_path = os.path.join(repodir, 'latex', 'master', 'master2.tex')

print("1. Parsing master2.tex for \\includegraphics...")
with open(master2_path, 'r') as f:
    latex_content = f.read()

# Find all \includegraphics instances
# E.g., \includegraphics[width=\textwidth]{figures_conclusion/spec_spectral.png}
# We just want the basename
matches = re.findall(r'\\includegraphics(?:\[.*?\])?\{(.+?)\}', latex_content)

used_basenames = set()
for m in matches:
    basename = os.path.basename(m)
    used_basenames.add(basename)
    
print(f"Found {len(used_basenames)} uniquely referenced figures:")
for name in sorted(used_basenames):
    print(f"  - {name}")

print("\n2. Sweeping for unused figures...")
# Directories that contain figures
figure_dirs = [
    os.path.join(repodir, 'figures'),
    os.path.join(repodir, 'latex', 'smale', 'figures'),
    os.path.join(repodir, 'latex', 'mcmullen') # sometimes images are here
]

deleted_count = 0
for d in figure_dirs:
    if not os.path.exists(d):
        continue
    for root, dirs, files in os.walk(d):
        for file in files:
            if file.endswith('.png') or file.endswith('.pdf'):
                if file not in used_basenames:
                    filepath = os.path.join(root, file)
                    print(f"Deleting unused asset: {filepath}")
                    os.remove(filepath)
                    deleted_count += 1

print(f"\nFinished. Deleted {deleted_count} unused figures.")
