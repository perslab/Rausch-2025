"""
Lepr/Glp1r co-expression across the whole mouse brain
using the Allen Brain Cell Atlas (WMB-10Xv3).

Downloads each brain division h5ad one at a time, extracts only
Lepr and Glp1r expression, then deletes the h5ad to save disk.
Joins with pre-downloaded cell annotations for brain region labels.

Usage:
    python paper/abc_coexpression.py
"""

import anndata
import pandas as pd
import numpy as np
import scipy.sparse
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Paths
cache_dir = os.path.join(os.path.dirname(__file__), '..', 'external_data', 'abc_atlas')
cache_dir = os.path.abspath(cache_dir)
out_dir = os.path.join(os.path.dirname(__file__), '..', 'output_figures', 'paper', 'extended_figures', 'abc_coexpression')
os.makedirs(out_dir, exist_ok=True)

# Ensembl IDs
LEPR_ID = 'ENSMUSG00000057722'
GLP1R_ID = 'ENSMUSG00000024027'

# Brain divisions in WMB-10Xv3
DIVISIONS = ['CB', 'CTXsp', 'HPF', 'HY', 'Isocortex-1', 'Isocortex-2',
             'MB', 'MY', 'OLF', 'P', 'PAL', 'STR', 'TH']

# Full names for plotting
REGION_NAMES = {
    'CB': 'Cerebellum', 'CTXsp': 'Cortical subplate', 'HPF': 'Hippocampus',
    'HY': 'Hypothalamus', 'Isocortex-1': 'Isocortex (1)', 'Isocortex-2': 'Isocortex (2)',
    'MB': 'Midbrain', 'MY': 'Medulla', 'OLF': 'Olfactory', 'P': 'Pons',
    'PAL': 'Pallidum', 'STR': 'Striatum', 'TH': 'Thalamus'
}

# ── Load cell annotations ──────────────────────────────────────────
annot_file = os.path.join(out_dir, 'cell_annotations.csv.gz')
if os.path.exists(annot_file):
    cells = pd.read_csv(annot_file)
    print(f"Loaded {len(cells):,} cell annotations")
else:
    print("ERROR: Run the annotation step first (cell_annotations.csv.gz not found)")
    exit(1)

# ── Extract Lepr/Glp1r per division ────────────────────────────────
from abc_atlas_access.abc_atlas_cache.abc_project_cache import AbcProjectCache
abc = AbcProjectCache.from_s3_cache(cache_dir)

all_results = []

for div in DIVISIONS:
    file_key = f"WMB-10Xv3-{div}/log2"
    print(f"Processing {div}...", end=' ', flush=True)

    try:
        path = abc.get_data_path(directory="WMB-10Xv3", file_name=file_key)
    except Exception as e:
        print(f"SKIP ({e})")
        continue

    adata = anndata.read_h5ad(path, backed='r')

    lepr_idx = list(adata.var_names).index(LEPR_ID)
    glp1r_idx = list(adata.var_names).index(GLP1R_ID)

    # Extract just these 2 columns from the sparse matrix
    X_sub = adata.X[:, [lepr_idx, glp1r_idx]]
    if scipy.sparse.issparse(X_sub):
        arr = np.asarray(X_sub.todense())
    else:
        arr = np.asarray(X_sub)

    df = pd.DataFrame({
        'cell_label': list(adata.obs_names),
        'Lepr': arr[:, 0],
        'Glp1r': arr[:, 1],
        'region': div
    })

    n = len(df)
    n_lepr = int((df['Lepr'] > 0).sum())
    n_glp1r = int((df['Glp1r'] > 0).sum())
    n_both = int(((df['Lepr'] > 0) & (df['Glp1r'] > 0)).sum())
    print(f"{n:,} cells | Lepr+={n_lepr:,} ({n_lepr/n*100:.1f}%) | "
          f"Glp1r+={n_glp1r:,} ({n_glp1r/n*100:.1f}%) | Both={n_both:,} ({n_both/n*100:.2f}%)")

    all_results.append(df)
    adata.file.close()

    # Delete the cached h5ad to save disk
    if os.path.exists(path) and div != 'HY':  # keep HY for potential further analysis
        os.remove(path)

# ── Combine and annotate ───────────────────────────────────────────
combined = pd.concat(all_results, ignore_index=True)
merged = combined.merge(cells[['cell_label', 'brain_class', 'subclass']], on='cell_label', how='left')

total_both = int(((merged['Lepr'] > 0) & (merged['Glp1r'] > 0)).sum())
print(f"\nTotal cells: {len(merged):,}")
print(f"Total co-expressing: {total_both:,}")

# ── Summary tables ─────────────────────────────────────────────────
summary = []
for region, grp in merged.groupby('region'):
    both = int(((grp['Lepr'] > 0) & (grp['Glp1r'] > 0)).sum())
    summary.append({
        'region': region, 'region_name': REGION_NAMES.get(region, region),
        'total_cells': len(grp),
        'lepr_pos': int((grp['Lepr'] > 0).sum()),
        'glp1r_pos': int((grp['Glp1r'] > 0).sum()),
        'both_pos': both,
        'pct_lepr': round((grp['Lepr'] > 0).mean() * 100, 2),
        'pct_glp1r': round((grp['Glp1r'] > 0).mean() * 100, 2),
        'pct_both': round(both / len(grp) * 100, 3)
    })
summary_df = pd.DataFrame(summary).sort_values('both_pos', ascending=False)
summary_df.to_csv(os.path.join(out_dir, 'coexpression_by_region.csv'), index=False)

print("\n=== Co-expression by region ===")
print(summary_df.to_string(index=False))

# Subclass breakdown of co-expressing cells
both_mask = (merged['Lepr'] > 0) & (merged['Glp1r'] > 0)
coexpr_sub = merged[both_mask]['subclass'].value_counts().reset_index()
coexpr_sub.columns = ['subclass', 'n_coexpressing']
coexpr_sub.to_csv(os.path.join(out_dir, 'coexpression_by_subclass.csv'), index=False)

print(f"\n=== Top 20 subclasses with Lepr/Glp1r co-expression ===")
print(coexpr_sub.head(20).to_string(index=False))

# % of co-expressing cells in hypothalamus
hy_both = int(((merged[merged['region'] == 'HY']['Lepr'] > 0) &
               (merged[merged['region'] == 'HY']['Glp1r'] > 0)).sum())
print(f"\n{hy_both:,} / {total_both:,} ({hy_both/total_both*100:.1f}%) of co-expressing cells are in the hypothalamus")

# ── Plots ──────────────────────────────────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(14, 5))

# 1. % cells expressing each gene by region
ax = axes[0]
plot_df = summary_df.sort_values('pct_both', ascending=True)
y = range(len(plot_df))
ax.barh(y, plot_df['pct_lepr'], height=0.3, label='Lepr+', color='#2980B9', alpha=0.7)
ax.barh([i + 0.3 for i in y], plot_df['pct_glp1r'], height=0.3, label='Glp1r+', color='#E74C3C', alpha=0.7)
ax.set_yticks([i + 0.15 for i in y])
ax.set_yticklabels(plot_df['region_name'], fontsize=7)
ax.set_xlabel('% cells expressing', fontsize=8)
ax.legend(fontsize=7)
ax.set_title('Expression by region', fontsize=9)

# 2. % co-expressing by region
ax = axes[1]
ax.barh(y, plot_df['pct_both'], color='#8E44AD')
ax.set_yticks(y)
ax.set_yticklabels(plot_df['region_name'], fontsize=7)
ax.set_xlabel('% cells co-expressing', fontsize=8)
ax.set_title('Lepr/Glp1r co-expression', fontsize=9)

# 3. Absolute number of co-expressing cells
ax = axes[2]
ax.barh(y, plot_df['both_pos'], color='#27AE60')
ax.set_yticks(y)
ax.set_yticklabels(plot_df['region_name'], fontsize=7)
ax.set_xlabel('# co-expressing cells', fontsize=8)
ax.set_title('Co-expressing cell count', fontsize=9)

plt.tight_layout()
plt.savefig(os.path.join(out_dir, 'lepr_glp1r_coexpression_brain.pdf'), dpi=300)
plt.savefig(os.path.join(out_dir, 'lepr_glp1r_coexpression_brain.png'), dpi=300)
print(f"\nSaved plots to {out_dir}")
