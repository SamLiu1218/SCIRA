import numpy as np
import pandas as pd
import scanpy as sc


import warnings
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

def row_norm(df, mode = 'max'):
    if mode == 'max':
        return (df.T/df.T.max()).T
    if mode == 'sum':
        return (df.T/df.T.sum()).T
        
from sklearn.neighbors import NearestNeighbors
def Edf(
    st_data, 
    pairs_df, 
    exp = False, 
    fill = -1,
    sop = True, 
    radius = 15, 
    autocrine = False, 
    col = 'label_fine',
):

    xy = st_data.obsm['spatial']
    neigh = NearestNeighbors(radius=radius, algorithm='auto')
    # Fit the model
    neigh.fit(xy)
    indices = neigh.radius_neighbors(xy, return_distance=False)
    if autocrine:
        El = [[cell_idx, v]  for cell_idx, neighbors in enumerate(indices) for v in neighbors]
    else:
        El = [[cell_idx, v]  for cell_idx, neighbors in enumerate(indices) for v in neighbors if v != cell_idx]

    E_df = pd.DataFrame(El, columns = ['Cell_id_R','Cell_id_S'])
    E_df[pairs_df['RL']] = fill
    for rl_i in range(len(pairs_df)):
        Rgene, Lgene, RL = pairs_df.loc[rl_i, ['R','L','RL']].values
        Rs = st_data.X[E_df.Cell_id_R, st_data.var_names == Rgene]
        Ls = st_data.X[E_df.Cell_id_S, st_data.var_names == Lgene]

        if sop:
            # E_df[RL] = np.sqrt(Ls*Rs)
            E_df.loc[Rs > 0, RL] = np.sqrt(Ls[Rs > 0] * Rs[Rs > 0])
        else:
            E_df.loc[Rs > 0, RL] = Ls[Rs > 0]
            
    E_df[['Cell_id_R','Cell_id_S']] = E_df[['Cell_id_R','Cell_id_S']]

    if exp:
        E_df[pairs_df['RL']] = np.exp(E_df[pairs_df['RL']]) - 1

    E_df['Cell_idx_R'] = st_data.obs.index.values[E_df['Cell_id_R']]
    E_df['Cell_idx_S'] = st_data.obs.index.values[E_df['Cell_id_S']]

    E_df['Sender'] = st_data.obs.loc[E_df['Cell_idx_S'], col].astype(str).values
    E_df['Receiver'] = st_data.obs.loc[E_df['Cell_idx_R'], col].astype(str).values
    E_df['SR'] = E_df['Sender'] + ' to ' + E_df['Receiver']

    Ls = pairs_df['L'].unique()
    E_df[Ls] = st_data.to_df().loc[E_df['Cell_idx_S'], Ls].values

    return E_df

def RA(E_df, RL_list):
    RA_ls = []
    for RL in RL_list:
        df = E_df[['Cell_idx_R', 'Receiver', RL]].copy()
        df.loc[df[RL] < 0, RL] = 0
        ra = df.groupby([ 'Cell_idx_R', 'Receiver'])[RL].sum().fillna(0).reset_index().groupby(['Receiver'])[RL].mean().fillna(0)
        RA_ls.append(ra)
    RA = pd.concat(RA_ls, axis = 1)
    return RA
    
def SR_df(E_df, RL_list):
    E_df = E_df.copy()
    sr_ls = []
    for RL in RL_list:
        E_df.loc[E_df[RL] < 0, RL] = 0
        
        ra = E_df.groupby([ 'Cell_idx_R', 'Receiver'])[RL].sum().fillna(0).reset_index().groupby(['Receiver'])[RL].mean().fillna(0)
        R, L = RL.split('_')
        sl = E_df.groupby(['Cell_idx_S', 'Receiver', 'Sender'])[L].sum().reset_index().groupby(['Receiver', 'Sender'])[L].mean().unstack().fillna(0)
        
        # sns.heatmap((ztl.row_norm(sl, mode = 'sum').T*ra).T)
        
        sr_ls.append((row_norm(sl, mode = 'sum').T*ra).T.stack().rename(RL))

    sr_df = pd.concat(sr_ls, axis = 1)
    return sr_df