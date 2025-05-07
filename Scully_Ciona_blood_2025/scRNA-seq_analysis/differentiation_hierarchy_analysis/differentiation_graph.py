import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import networkx as nx
import sys
import os
from tqdm import tqdm

# Change this path to point to folder containing helper functions scripts
path_to_repo_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(os.path.join(path_to_repo_dir, 'helper_functions'))
import scrna_helper_functions as hf

plt.style.use('tal_paper_spine')

# Set random seed
np.random.seed(seed = 0)

# ============================================================================
# IMPORT DATA

# Import counts matrix
adata = sc.read_h5ad(hf.path.Crob_adata_file)

# Set up output folder
out_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        'differentiation_graph_output')
if not os.path.exists(out_path): os.mkdir(out_path)

# ============================================================================
# DEFINE CYCLING PROGENITOR CELL STATES
# Cell states are defined to be cycling if their median expression of MKI67
# is above 0.

percentile = 50 # 50th percentile = median
percentile_dict = {}
for cell_type in adata.obs['cell_type'].cat.categories:
    adata_subset = adata[adata.obs['cell_type'] == cell_type]
    percentile_dict[cell_type] = np.percentile(
        np.array(adata_subset.raw[:, hf.human2ciona['MKI67']].X.todense()),
        percentile)

with plt.style.context('tal_light'):
    f = plt.figure(figsize=(10, 8))
    ax = plt.subplot(2, 2, 1)
    sc.pl.umap(adata, color=hf.human2ciona['MKI67'], show=False, ax=ax)
    ax = plt.subplot(2, 2, 2)
    sc.pl.umap(adata, color='cell_type', show=False, ax=ax,
               legend_loc='on data')

    ax = plt.subplot(2, 1, 2)
    plt.bar([c for c in percentile_dict],
            [percentile_dict[c] for c in percentile_dict])
    plt.xticks(rotation=90)
    plt.ylabel(f'{percentile}th percentile expression')

    plt.tight_layout()
    plt.savefig(os.path.join(out_path, 'mki67_expr_per_cell_type.png'),
                dpi=200)
    plt.close()

# ============================================================================
# BUILD A COARSE-GRAIN WEIGHTED GRAPH FROM THE SINGLE-CELL GRAPH

# From Wagner et al. 2018 methods section for building coarse-grain graph:
"""
Construction of Coarse-Grained Graphs
A coarse-graining procedure to abstract the major features of the single-cell
graph was performed as follows. First, single-cell nodes belonging to the same
annotated tSNE cluster ID were collapsed into a single state node.  Edges
between each pair of state nodes were then weighted by calculating the norm
index of original shared single-cell edges (i.e. the ratio of shared single-
cell edges to the total number of outgoing edges for that node pair). State
edges were then discarded if they received a norm index weight < 0.01.
Finally, a spanning tree was traced through the weighted edges as follows.
Beginning with the final timepoint, edges for all nodes were ranked according
to weight. Edges then were then removed recursively, starting with the weakest
edges, unless doing so would increase the total number of graph connected
components. This process was then repeated for each timepoint. The resulting
spanning tree connects all nodes to a single 4hpf “root” node defined by all
cells of the first timepoint.
"""

# ------------------------------------
# Get cell-cell edges that we care about

# Convert the single-cell k-nearest neighbors graph sparse matrix to a
# NetworkX graph
G = nx.from_scipy_sparse_matrix(adata.obsp['connectivities'])

# Filter cell edges to retain only those connected to EXACTLY 1 progenitor
# cell_edges = []
# for i in selected_cells:
#     neighbors = list(G.neighbors(i))
#     for j in neighbors:
#         # if j not in selected_cells:
#         cell_edges.append((i, j))

cell_edges = list(G.edges)

# ------------------------------------
# Count cell edges to create a weighted cell type graph

cell_type_list = adata.obs['cell_type'].cat.categories
cLRP_list = [c for c in cell_type_list if 'cLRP' in c]
prog_list = [c for c in cell_type_list if c[0] == 'c']
other_list = [c for c in cell_type_list if c not in prog_list]
prog_to_mature = pd.DataFrame(0, index=other_list, columns=prog_list)
all_to_all = pd.DataFrame(0, index=cell_type_list, columns=cell_type_list)

for (i, j) in tqdm(cell_edges):
    cell1 = adata.obs.iloc[i]['cell_type']
    cell2 = adata.obs.iloc[j]['cell_type']
    if cell1 != cell2:
        # Edges between progenitor and non-progenitor cell types
        if (cell1 in prog_list) and (cell2 not in prog_list):
            prog_to_mature.loc[cell2, cell1] += 1
        elif (cell2 in prog_list) and (cell1 not in prog_list):
            prog_to_mature.loc[cell1, cell2] += 1
        # Edges between cLRP and other cell types
        else:
            all_to_all.loc[cell1, cell2] += 1
            all_to_all.loc[cell2, cell1] += 1

# Calculate normalized edge counts
prog_to_mature_norm = prog_to_mature.copy()
total_edges_row = prog_to_mature.sum(axis=1)
total_edges_col = prog_to_mature.sum(axis=0)
for row in prog_to_mature_norm.index:
    for col in prog_to_mature_norm.columns:
        prog_to_mature_norm.loc[row, col] = (
            prog_to_mature.loc[row, col]
            / np.sqrt(total_edges_row[row] * total_edges_col[col])
        )
all_to_all_norm = all_to_all.copy()
total_edges_row = all_to_all.sum(axis=1)
total_edges_col = all_to_all.sum(axis=0)
for row in all_to_all_norm.index:
    for col in all_to_all_norm.columns:
        all_to_all_norm.loc[row, col] = (
            all_to_all.loc[row, col]
            / np.sqrt(total_edges_row[row] * total_edges_col[col])
        )

# ------------------------------------
# Plot as heatmap

graph_dict = {
    'prog_to_mature': prog_to_mature_norm,
    'all_to_all': all_to_all_norm,
}

for graph_name in graph_dict:
    norm_index = graph_dict[graph_name]

    # Reorder for plotting
    marginals_col = norm_index.max(axis=0)
    prog_order = norm_index.columns[np.argsort(marginals_col)[::-1]]
    marginals_row = norm_index.max(axis=1)
    other_order = norm_index.index[np.argsort(marginals_row)[::-1]]

    f = plt.figure(figsize=(1.126 + 0.1*norm_index.shape[1],
                            0.6297 + 0.1*norm_index.shape[0]))
    sns.heatmap(norm_index.loc[other_order, prog_order],
                cmap='rocket_r', xticklabels=True, yticklabels=True,
                cbar_kws={'label': 'Normalized edge counts'})
    plt.tight_layout()
    plt.savefig(os.path.join(out_path, f'connections_heatmap_{graph_name}.pdf'))
    plt.close()

# ============================================================================
# REMOVING LOW-WEIGHT EDGES AND GET MAXIMUM SPANNING TREE

# Graph of cell types
G = nx.Graph()
G.add_nodes_from(cell_type_list)
G.nodes['cMPP']['group'] = 'cMPP'
for n in cLRP_list:
    G.nodes[n]['group'] = 'cLRP'
for n in other_list:
    G.nodes[n]['group'] = 'mature'

# Set edge weight threshold
f = plt.figure(figsize=(3, 4))
percentage = 30
for i, graph_name in enumerate(graph_dict):
    norm_index = graph_dict[graph_name]
    threshold = (percentage/100) * norm_index.max().max()

    freq, bins = np.histogram(norm_index.values.flatten(), bins=10)
    ax = plt.subplot(2, 1, i+1)
    plt.bar(bins[:-1], freq, width=np.diff(bins), color='#888888')
    plt.axvline(x=threshold,label=f'Threshold = {percentage}%\nof max weight')
    plt.yscale('log')
    plt.ylabel('Frequency')
    plt.legend()
    plt.title(graph_name)
# Finish plot
plt.xlabel('Edge weight')
plt.tight_layout()
plt.savefig(os.path.join(out_path, 'edge_weight_distribution.pdf'))
plt.close()

# Add edges to grpah object
graph_edges = []
for node1 in G.nodes:
    for node2 in G.nodes:
        if node1 != node2:
            # Weights from each graph, set to zero if below threshold
            weight1 = all_to_all_norm.loc[node1, node2]
            if weight1 < ((percentage/100) * all_to_all_norm.max().max()):
                weight1 = 0
            if (node2 in prog_list) and (node1 not in prog_list):
                weight2 = prog_to_mature_norm.loc[node1, node2]
                if weight2 < ((percentage/100) * prog_to_mature_norm.max().max()):
                    weight2 = 0
            elif (node1 in prog_list) and (node2 not in prog_list):
                weight2 = prog_to_mature_norm.loc[node2, node1]
                if weight2 < ((percentage/100) * prog_to_mature_norm.max().max()):
                    weight2 = 0
            else:
                weight2 = 0
            
            # Take max of weights from each graph
            weight = max(weight1, weight2)
            if weight != 0:
                if (node2, node1, {'weight': weight}) not in graph_edges:
                    graph_edges.append((node1, node2, {'weight': weight}))
G.add_edges_from(graph_edges)

# Save nodes without edges, then remove them
with open(os.path.join(out_path, 'nodes_without_edges.txt'), 'w') as f:
    f.write('\n'.join(list(nx.isolates(G))))
G.remove_nodes_from(list(nx.isolates(G)))

# Take the maximum spanning tree
G = nx.maximum_spanning_tree(G, weight='weight')

# ============================================================================
# VISUALIZE COARSE-GRAIN TREE

fig_len = {
    0: 2.25,
    1: 0.5,
    2: 0.7,
}
for i, node_set in enumerate(list(nx.connected_components(G))):
    subG = G.subgraph(node_set).copy()

    # Plot as a tree
    # (Have to rename nodes to not get an error)
    mapping = {node: f"node_{i}" for i, node in enumerate(subG.nodes)}
    mapping_reverse = {mapping[node]: node for node in mapping}
    subG = nx.relabel_nodes(subG, mapping)
    pos = nx.drawing.nx_pydot.graphviz_layout(subG, prog="neato")
    subG = nx.relabel_nodes(subG, mapping_reverse)
    pos_renamed = {mapping_reverse[node]: pos[node] for node in mapping_reverse}

    # Step 3: Draw the graph
    edge_weights = nx.get_edge_attributes(subG, 'weight')
    color_map = {'cMPP': '#d94a34', 'cLRP': '#c5a801', 'mature': '#7db9d1'}
    node_colors = [color_map[subG.nodes[node]['group']] for node in subG.nodes()]
    f = plt.figure(figsize=(fig_len[i], fig_len[i]))
    nx.draw(
        subG,
        pos_renamed,
        with_labels=True,
        node_size=350,
        node_color=node_colors,
        edge_color="gray",
        font_size=7,
        width=[7*weight for weight in edge_weights.values()],
        labels={c: c.replace(' (', '\n(') for c in subG.nodes},
        alpha=1,
    )
    f.savefig(os.path.join(out_path, f'cell_type_graph_{i+1}.pdf'))
    plt.close()

# ------------------------------------
# Plot as heatmap

df = pd.DataFrame(nx.adjacency_matrix(G).todense(), index=G.nodes(),
                columns=G.nodes())

# Reorder for plotting
marginals_col = df.max(axis=0)
prog_order = df.columns[np.argsort(marginals_col)[::-1]]
marginals_row = df.max(axis=1)
other_order = df.index[np.argsort(marginals_row)[::-1]]

f = plt.figure(figsize=(1.126 + 0.1*df.shape[1], 0.6297 + 0.1*df.shape[0]))
sns.heatmap(df.loc[other_order, prog_order],
            cmap='rocket_r', xticklabels=True, yticklabels=True,
            cbar_kws={'label': 'norm index'})
plt.tight_layout()
plt.savefig(os.path.join(out_path, f'cell_type_graph_adjacency_matrix.pdf'))
plt.close()
