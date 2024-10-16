import palantir
import scanpy as sc
import numpy as np
import pandas as pd
import os
from harmony import harmonize
# Plotting
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
sns.set_style("ticks")
matplotlib.rcParams["figure.figsize"] = [4, 4]
matplotlib.rcParams["figure.dpi"] = 100
matplotlib.rcParams["image.cmap"] = "Spectral_r"

# warnings
import warnings
from numba.core.errors import NumbaDeprecationWarning

warnings.filterwarnings(action="ignore", module="matplotlib", message="findfont")
warnings.filterwarnings(action="ignore", category=NumbaDeprecationWarning)
warnings.filterwarnings(
    action="ignore", module="scanpy", message="No data for colormapping"
)

#Read in h5ad file
adata = sc.read_h5ad('data.h5ad')
#Choose clusters relevant to the current analysis
clusters_of_choice = ['ASPC 1', 'ASPC 2', 'ASPC 3', 'ASPC 4', 'VA1', 'VA2', 'VA3', 'VA4', 'VA5','VA6','VA7','VA8', 'Mac 1', 'Mac 2']

adata = adata[adata.obs.new_detailed.isin(clusters_of_choice)]
sc.pp.normalize_per_cell(adata)
palantir.preprocess.log_transform(adata)

#Change n_top_genes to 2500 if running only ASPC-VA, all the rest 1500

sc.pp.highly_variable_genes(adata, n_top_genes=1500, flavor="cell_ranger")
sc.pp.pca(adata)
#Batch correction
Z = harmonize(adata.obsm['X_pca'], adata.obs, batch_key = 'orig.ident')
adata.obsm['X_harmony'] = Z

dm_res = palantir.utils.run_diffusion_maps(adata, pca_key='X_harmony', n_components=5)
ms_data = palantir.utils.determine_multiscale_space(adata)

sc.pp.neighbors(adata, use_rep='X_harmony')
sc.tl.umap(adata)
sc.tl.draw_graph(adata)
graph = pd.DataFrame(adata.obsm['X_draw_graph_fa'], index=adata.obs_names)
#Plot clusters map on FA graph
sc.pl.draw_graph(adata, color='new_detailed', title='',frameon=False,show=False)
plt.savefig('clusters.png', dpi=300, bbox_inches='tight')
plt.clf()

#Plot gene trends based on MAGIC imputation, change 'color' parameter for different genes
imputed_X = palantir.utils.run_magic_imputation(adata)

sc.pl.embedding(
    adata,
    basis="draw_graph_fa",
    layer="MAGIC_imputed_data",
    color=["ADIPOQ", "PDGFRA", "PPARG"],
    frameon=False, show=False
)
plt.savefig('marker_genes.png', dpi=300)
plt.clf()

n_cells = adata.n_obs #Used to later adjust dot size
cluster_name = "ASPC 1" #Cluster to start the analysis from
cells_in_cluster = adata.obs_names[adata.obs['new_detailed'] == cluster_name]

# Get the expression values of F3 and DPP4 for cells in the cluster
expression_f3 = adata[cells_in_cluster, 'F3'].X.toarray().flatten()
expression_dpp4 = adata[cells_in_cluster, 'DPP4'].X.toarray().flatten()

# Find cells with the lowest expression of F3 (or expression value equal to 0)
lowest_f3_expression = np.min(expression_f3)
lowest_f3_expression_cells = cells_in_cluster[expression_f3 == lowest_f3_expression]

# If there are multiple cells with expression 0, choose the one with the highest expression of DPP4
if len(lowest_f3_expression_cells) > 1:
    # Find the indices of the 100 highest DPP4 expressing cells within the Low F3 expression population
    indices_sorted_by_dpp4 = np.argsort(expression_dpp4[expression_f3 == lowest_f3_expression])[::-1][:100]
    # Get the cells with the 100 highest DPP4 expression
    cells_highest_dpp4_expression = lowest_f3_expression_cells[indices_sorted_by_dpp4]
    # Choose a random cell from the 100 highest DPP4 expressing cells
    start_cell_barcode = np.random.choice(cells_highest_dpp4_expression)
else:
    start_cell_barcode = lowest_f3_expression_cells[0]

start_cell_index = adata.obs_names.get_loc(start_cell_barcode)
adata.uns['iroot'] = start_cell_index

#After setting start cell, run palantir
pr_res = palantir.core.run_palantir(adata, start_cell_barcode, num_waypoints=500)
palantir.plot.plot_palantir_results(graph, pr_res, s=25000/n_cells)
plt.savefig('palantir_pseudotime.png', dpi=300)
plt.clf()

#Run dpt for validation
sc.tl.dpt(adata)
sc.pl.draw_graph(adata, color=['dpt_pseudotime'], size=50000/n_cells,frameon=False, show=False)
plt.savefig('dpt_comparison_pseudotime.png', dpi=300)
plt.clf()
