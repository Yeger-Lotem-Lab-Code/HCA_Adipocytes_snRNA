import os
import anndata as ad
import pandas as pd
import ktplotspy as kpy
import matplotlib.pyplot as plt
import scanpy as sc
from scipy import io
from cellphonedb.src.core.methods import cpdb_statistical_analysis_method

## Make from your h5ad data compatible files for CellphoneDB analysis

adata = ad.read_h5ad('vis_data.h5ad')

with open('data/barcodes.tsv', 'w') as f:
    for item in adata.obs_names:
        f.write(item + '\n')

with open('data/features.tsv', 'w') as f:
    for item in ['\t'.join([x, x, 'Gene Expression']) for x in adata.var_names]:
        f.write(item + '\n')

io.mmwrite('data/matrix', adata.X.T)

adata.obs.to_csv('metadata.csv')

metadata = pd.read_csv('metadata.csv', index_col=0)

#CellphoneDB requires a metadata with the barcodes as index and one column of metadata - relevant clusters annotation

metadata = metadata['new_detailed']
metadata.to_csv('metadata_2columns.csv')

##Run CellphoneDB analysis - input path should be a folder (!!!) with mtx, features, barcodes files only. You should also supply path to the database itself and the metadata

from cellphonedb.src.core.methods import cpdb_statistical_analysis_method
deconvoluted, means, pvalues, significant_means = cpdb_statistical_analysis_method.call(
        cpdb_file_path = 'cellphonedb.zip',
        meta_file_path = 'metadata_2columns.csv',
        counts_file_path = 'data',
        counts_data = 'hgnc_symbol',
        output_path = 'results')


## Next step will be creating the heatmap table, we will use the adata and the pvalues object created by the cpdb analysis

a = kpy.plot_cpdb_heatmap(
        adata=adata,
        pvals=pvals,
        celltype_key="new_detailed",
        figsize = (5,5),
        title = "Sum of significant interactions",
        return_tables=True
    )

data = pd.DataFrame.from_dict(a['count_network'])

#Now plotting the triangle heatmaps using the 'df' object

def get_custom_color_palette_hash():
    return LinearSegmentedColormap.from_list("", [
        '#104e8b', '#ffdab9', '#8b0a50'
    ])



for i in range(len(data)):
    for j in range(i+1, len(data)):
        data.iloc[j,i] += data.iloc[i,j]
        data.iloc[i,j] = 0

df_lt = data.where(np.tril(np.ones(data.shape)).astype(np.bool_))
cmap = get_custom_color_palette_hash()
hmap = sns.heatmap(df_lt, cmap=cmap)
hmap.collections[0].set_edgecolor("white")
hmap.collections[0].set_linewidth(0.5)
hmap.set(ylabel=None)
plt.title('VAT')
plt.tight_layout()
plt.savefig(fname='triangle_heatmap.png', dpi=300)

