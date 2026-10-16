# SCIRA: Spatial Cellular Interaction and Receptor Activation

**Spatial Cellular Interaction and Receptor Activation (SCIRA)** quantifies receptor–ligand signaling between neighboring cells in spatial single-cell data, with results at the cell, receiver-cell-type, and sender–receiver-cell-type levels.

This repository provides the Python implementation (`SCIRA.py`). A small example dataset and Jupyter notebook (`example_scira.ipynb`) are also provided as a tutorial. For a detailed explanation of the method, see the accompanying manuscript.

## Setup

Use a Python environment with the following packages:

```bash
python -m pip install numpy pandas scikit-learn scanpy matplotlib seaborn
```

Note that matplotlib and seaborn are used for plotting in the example notebook, not required by SCIRA itself. These packages are also available via Conda.

To use, simply copy `SCIRA.py` to a place that's accessible by your Python environment.

## Quick start

For the full workflow, including plots, see [`example_scira.ipynb`](example_scira.ipynb).

### Prepare the data

The required inputs include an AnnData object with gene expression, spatial coordinates, and cell-type labels, together with a curated list of potential receptor–ligand pairs.

The receptor and ligand genes in the pair list must be present in `st_data.var_names`. Spatial coordinates should use the same distance units as `radius`.

SCIRA expect spatial coordiantes at `st_data.obsm['spatial']`.

### Calculate cell–cell interactions

```python
E_df = SCIRA.Edf(
    st_data,
    pairs_df=pairs_df,
    radius=15,
    col='label',
    autocrine=False,
    sop=True,
    fill=-1,
    exp=False,
)
```

`E_df` contains neighboring cell pairs, their sender and receiver identities, and interaction values for each receptor–ligand pair. The `fill` value (default: `-1`) indicates no receptor expression and `0` indicates receptor expression without ligand expression in the neighboring sender cell, for those cases where the distinction matters.

**Common settings**

| Parameter | What to change it for |
| --- | --- |
| `radius` | Maximum center-to-center distance for neighboring cells; default = `15` . |
| `col` | Cell-type annotation column in `st_data.obs`. |
| `autocrine` | Whether to include a cell's interaction with itself. |
| `sop` | If `True`, use the geometric mean of receptor and ligand expression for a cell pair; if `False`, use ligand expression when the receptor is expressed. |
| `fill` | Value assigned when the receiver does not express the receptor (default = -1). |
| `exp` | Whether to transform interaction values using `exp(value) - 1`. For oddly scaled data. |

### Summarize receptor activation

To sum interaction values over the neighbors of each receiving cell:

```python
RL_df = E_df.groupby('Cell_idx_R')[RL_list].sum()
```

`RL_df` can be used to visualize spatial receptor activation or make a dot plot, as shown in the notebook. To calculate mean receptor activation for each receiver cell type:

```python
RA = SCIRA.RA(E_df, RL_list=pairs_df['RL'])
```

`RA` is a DataFrame with receiver cell types as rows and receptor–ligand pairs as columns.

### Summarize sender–receiver interactions

```python
sr_df = SCIRA.SR_df(E_df, RL_list=pairs_df['RL'])
```

`sr_df` contains interaction scores for sender–receiver cell-type pairs, with one column per receptor–ligand pair. Its row index is `(Receiver, Sender)`, in that order. The example notebook illustrates one way to visualize selected interactions as a heatmap.

For visualization, SCIRA.row_norm() can normalize scores by row maximum (mode='max') or row sum (mode='sum').

## Citation

If you use SCIRA, please cite the accompanying manuscript.
