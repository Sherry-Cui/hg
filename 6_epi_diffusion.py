#!/usr/bin/env python
# coding: utf-8

# fundus and antrum 
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as pl
adata = sc.read_h5ad("epirmovegland.h5ad")

sc.pl.embedding(adata,'X_diffmap',color='cell.type')
adata.obs['cell.type'].value_counts()
adata.uns['iroot'] = np.flatnonzero((adata.obs['cell.type'] == 'Fundus#1') | (adata.obs['cell.type'] == 'Antrum#1'))[0]
sc.tl.dpt(adata)
sc.pl.embedding(adata,'X_diffmap',color=['dpt_pseudotime','cell.type'])
pseudotime=adata.obs['dpt_pseudotime']
pseudotime.to_csv('epiremovegland_pseudotime.csv',index=True)

# antrum 
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as pl
adata = sc.read_h5ad("antral.h5ad")

sc.pp.neighbors(adata, n_neighbors=30, n_pcs=15, random_state= 42)
sc.tl.diffmap(adata,n_comps=20)
sc.pl.diffmap(adata, color='cell.type', legend_loc='on data')
sc.pl.umap(adata, color=['day','cell.type'], legend_loc='on data')
sc.tl.dpt(adata)
adata.obs['cell.type'].unique()
adata.obs['cell.type'] = adata.obs['cell.type'].astype('category')
new_cluster_names = ["Antrum#2","Antrum#1","Antrum#3"]
adata.rename_categories('cell.type', new_cluster_names)
sc.pl.embedding(adata,'X_diffmap',color='cell.type')
adata.uns['iroot'] = np.flatnonzero(adata.obs['cell.type']  == 'Antrum#1')[0]
sc.tl.dpt(adata)
sc.pl.embedding(adata,'X_diffmap',color=['dpt_pseudotime','cell.type'])
pseudotime=adata.obs['dpt_pseudotime']
pseudotime.to_csv('antral_pseudotime.csv',index=True)

# fundus 
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as pl
adata = sc.read_h5ad("fudic.h5ad")

adata.obs['cell.type'].value_counts()
adata.uns['iroot'] = np.flatnonzero(adata.obs['cell.type']  == 'Fundus#1')[0]
sc.tl.dpt(adata)
sc.pl.embedding(adata,'X_diffmap',color=['dpt_pseudotime','cell.type'])
pseudotime=adata.obs['dpt_pseudotime']
pseudotime.to_csv('fundic_pseudotime.csv',index=True)







