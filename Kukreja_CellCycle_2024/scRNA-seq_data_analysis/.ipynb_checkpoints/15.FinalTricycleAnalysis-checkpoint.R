# Script containing all of the code necessary for the work done for Kalki
# Kukreja's paper on the perturbation of cell cycle in zebrafish.

# Command clearing all variables
rm(list=ls())

# STEPS TO CREATE CONDA ENVIRONMENT TO RUN THIS SCRIPT:
# conda create -n cc_perturb r-base=4.3.3
# conda install bioconda::bioconductor-zellkonverter
# conda install bioconda::bioconductor-tricycle

# STEPS TO ANALYZE THE PERTURBED CELL CYCLE DATA AND GENERATE A FIGURE:
# 1. Replace the following path with the corect path to the perturbed data set.
K_INPUT_PATH = "/Users/sean/Dropbox (HMS)/_klein_lab/cell_cycle_heterogeneity_shared/adata24_perturbations.h5ad"
# 2. Within the `cc_perturb` conda environment, navigate to the path where the
#     script is located.
# 3. Type `R` in the terminal.
# 3. Within R, type `source("FinalTricycleAnalysis.R")`


# `tricycle`` enables the SingleCellExperiment package function to work, such as
# `colData`.
library(tricycle)
#  `zellkonverter` enables the `readH5AD` function to work.
library(zellkonverter)


## FUNCTIONS ###################################################################
# 1.____________________________________________________________________________
LoadCellCycleHetData <- function() {
  # This function loads the dataset and renames the column name associated with
  # the cell states "cluster". This is inherited from when applying these
  # functions to multiple datasets where all datasets were processed to have
  # the same column information.
  sce <- (readH5AD(K_INPUT_PATH))
  colnames(colData(sce))[which(colnames(colData(sce)) == "state_v3p2")] <- "cluster"
  return(sce)
}


# 2.____________________________________________________________________________
LogNormalizeData <- function(sce, pseudocount=1) {
  # Function that performs library-size normalization, log transformation, and
  # optionally adds a pseduocount. This is currently performed identically to
  # how it is performed within the Tricycle paper, which uses the
  # NormalizeCounts function from the scRNA-seq analysis package scuttle.
  # 1. Get the counts matrix.
  counts <- assay(sce)
  # Normalize by mean library size. 
  libsizes <- colSums(counts)
  size.factors <- libsizes/mean(libsizes)
  # Make new assay that is log2-transform of data + (the optional)
  # pseudocount.
  logcounts(sce) <- log2(t(t(counts)/size.factors) + pseudocount)
  return(sce)
}


# 3.____________________________________________________________________________
GetCompatibleNeuroRef <- function(sce) {
  # 1.  Load the mouse orthology table, and remove the right-most columne (which 
  #     is generated due to a processing error with the raw table downloaded from
  #     ZFIN.
  ortho_table <- read.table(
    "mouse_orthos_2024.04.30.txt", sep="\t",
    header=TRUE, row.names=NULL, fill=TRUE, skip=1, quote=""
  )
  # 2.  Subset the dataframe to only include rows that are based on AA homology,
  #     to remove the last two columns of data, and only unique rows. The columns
  #     being dropped are an empty column and a column listing the published
  #     reference associated with each row. In cases where multiple references
  #     support the same orthology annotation, this leads to redundant rows
  #     which are collapsed using the "unique" function.
  ortho_table <- unique(
    ortho_table[which(ortho_table$Evidence == "AA"), 1:(ncol(ortho_table) - 2)]
  )
  #SYMBOLS in ZFIN TABLE:
  # AA: Amino-acid sequence composition
  # CL: Conserved genome location (synteny)
  # PT: Phylogenetic tree
  # CE: Coincident expression
  # FC: Functional complementation
  # NT: Nucleotide sequence comparison
  # OT: Other
  # 3.  Load the neuroRef table from the tricycle package.
  data(neuroRef)
  rownames(neuroRef) <- neuroRef$symbol
  # 4.  Subset the orthology table to contain only genes for which the 
  #     gene symbol of the zebrafish ortholog is found in the scRNA-seq
  #     experiment object, and the gene symbol of the mouse ortholog is found in
  #     the neuroRef projection matrix provided by the Tricycle package.
  ortho_row_inds_use <- which(
    ((ortho_table$ZFIN.Symbol %in% rownames(sce)) &       # gene in fish
     (ortho_table$Mouse.Symbol %in% rownames(neuroRef)))  # gene in mouse
  )
  ortho_table <- ortho_table[ortho_row_inds_use, ]
  # 5.  Subset the neuroRef using the orthology table, and rename the rows to be
  ##    compatible with the scRNA-seq dataset.
  neuroRef <- neuroRef[ortho_table$Mouse.Symbol, ]
  rownames(neuroRef) <- ortho_table$ZFIN.Symbol
  return(neuroRef)
}


# 4. ___________________________________________________________________________
AssignTricycleAngle <- function(sce) {
    #   This function returns a `SingleCellExperiment` object in which both the
    # projection vectors and the tricycle angle (in radians) has been added.
    #
    #	sce : The `SingleCellExperiment` object for which the projected angle is
    #       added.
    # 1. PROCESS SCRNA-SEQ DATA. #################################################
    # 1a. Load the scRNA-seq data. _______________________________________________
    if (names(assays(sce))[1] == "X") {
        assayNames(sce)[1] <- "counts"
    }
    # 1b. Library normalize and log+1-transform. _______________________________
    if (!("logcounts" %in% names(assays(sce)))) {
        sce <- LogNormalizeData(sce)
    }
    # 2. PROCESS REFERENCE. ####################################################
    # 2a. Load the projection data reference.
    neuroRef <- GetCompatibleNeuroRef(sce)
    # 3. EXCLUDE GENES BASED ON OVERLAP BETWEEN SCRNA-SEQ AND REFERENCE. #######
    # 4. ORDER GENES AND FORMAT REFERENCE. #####################################
    sce <- sce[rownames(neuroRef), ]
    neuroRef <- as.matrix(neuroRef[, 1:2])
    # 5. PROJECT THE DATA INTO THE FORMATTED EMBEDDING. ########################
    sce_proj <- project_cycle_space(sce, ref.m=neuroRef)
    sce_proj <-  estimate_cycle_position(sce_proj)
    # 6. Explicitly bring the pc components of each cell out of the
    # `tricycleEmbedding` that it is put in by the tricycle function.
    colData(sce_proj)$pc1 <- reducedDim(sce_proj, "tricycleEmbedding")[, 1]
    colData(sce_proj)$pc2 <- reducedDim(sce_proj, "tricycleEmbedding")[, 2]
    return(sce_proj)
}


# 5. ___________________________________________________________________________
ClassifyTricyclePhases <- function(values) {
  # 1.  Define the break points based on the description in the following website:
  # https://bioconductor.org/packages/release/bioc/vignettes/tricycle/inst/doc/tricycle.html#alternative-infer-cell-cycle-stages
  # "The estimated cell cycle position is bound between 0 and 2pi. Note that we
  #  strive to get high resolution of cell cycle state, and we think the
  #  continuous position is more appropriate when describing the cell cycle.
  #  However, to help users understand the position variable, we also note that
  #  users can approximately relate 0.5pi to be the start of S stage, pi to be
  #  the start of G2M stage, 1.5pi to be the middle of M stage, and
  #  1.75pi-0.25pi to be G1/G0 stage."
  breakpoints <- c(0,       0.25 * pi,       0.5 * pi,    pi,       1.5 * pi,    1.75 * pi,      2 * pi)
  labels <-      c(   "M/G1",          "G1/S",         "S",   "G2/M",         "M",          "M/G1")
  # 2. Convert values outside of the 0–2pi range using modulo operator.
  values_mod <- (values %% (2 * pi))
  # 3. Classify each value based on the range it falls into
  cut_points <- cut(values_mod, breaks = breakpoints, labels = labels, include.lowest = TRUE)  
  return(as.character(cut_points))
}


# 9. ___________________________________________________________________________
MakeAngleAndPhaseTable <- function(sce, tissue="") {
    # Can't I pass this in a larger loop?
    # sce_zf <- LoadCellCycleHetData()
    # Down-sample to one tissue.
    sce_use <- sce[, colData(sce)$cluster == tissue]
    # Assign the cell-cycle angle using the tricylce package.
    sce_p <- AssignTricycleAngle(sce_use)
    # Add a categorical label of the cell cycle phase.
    colData(sce_p)$tricycleLabel <- ClassifyTricyclePhases(colData(sce_p)$tricyclePosition)

    output <- colData(sce_p)[, c("cluster", "treatment", "tricyclePosition", "tricycleLabel", "pc1", "pc2")]
    colnames(output) <- c("cluster", "Treatment", "Tricycle angle", "Tricycle phase", "PC1", "PC2")
    return(output)
}


# 10. __________________________________________________________________________
MakeAllPhaseAndAngleTables <- function(outputfile=FALSE) {
  # 1. Pre-allocate an tempty dataframe:
  full_output_df <- data.frame(matrix(ncol = 0, nrow = 0))
  # 2. Load the scRNA-seq data to pass through the loop and to get each of the
  # tissues/cell types/ cell states for the loop,
  sce <- LoadCellCycleHetData()
  tissues <- unique(colData(sce)$cluster)
  # 3.  Iterate over the loop, generating a matrix in which each row is a cell
  #     and the columns provide the tricycle Position, the label, and the the
  #     value of the two cell-cycle projection components.
  for (tissue in tissues) {
    output <- MakeAngleAndPhaseTable(sce, tissue=tissue)
    full_output_df <- rbind(full_output_df, data.frame(output))
  }
  if (outputfile) {
      write.csv(full_output_df, "final_values.csv", row.names=TRUE)
  }
  return(full_output_df)
}

tricycle_df <- MakeAllPhaseAndAngleTables(outputfile=TRUE)
print(head(tricycle_df))
