The scripts in this folder were used to generate a hypothesized hematopoietic hierarchy for _C. robsuta_, shown in **Fig. 4c** of the paper. See methods section "Differentiation Hierarchy Analysis".

The files are:
- `differentiation_graph.py` is the python script that uses the scRNA-seq data to generate a hypothesized hierarchy. The script outputs a folder `differentiation_graph_output/` containing plots of the inferred hierarchy (graph, adjacency matrix), a text file listing cell states not connected to other nodes in the graph, as well as miscellaneous supporting plots.
- `differentiation_graph_env.txt` contains information on the conda environment used, saved to a text file (by running `conda list --explicit > differentiation_graph_env.txt`).
