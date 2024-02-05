#' Generate seeds for reproducibility when sampling from the starting dataset to perform bootstrapping
#'
#' This function will generate seeds for bootstrapping for each modeled gene.
#' It takes in the desired number of bootstrap samples and the list of genes as input
#' and returns a list, named by inputted genes, of lists of seeds.
#'
#' @param num_bootstraps Number of bootstrap samples
#' @param gene_names Modeled gene names
#' @param dims.list A list containing the dimensions for each reduction to use
#'
#' @return List named by inputted gene names, of lists of seeds, each of length num_bootstraps
#' @export
#'
make_bootstrap_sequences <- function(num_bootstraps, gene_names) {
	bootstrap_sequence <- sample((-2^31 + 1):(2^31 - 1), length(gene_names))
	bootstrap_sequence_input_lists <- lapply(bootstrap_sequence, function(x) {
                                        	set.seed(x)
                                        	return(sample((-2^31 + 1):(2^31 - 1), num_bootstraps))
                                  	  })
	names(bootstrap_sequence_input_lists) <- gene_names

	return(bootstrap_sequence_input_lists)
}

#' Select differentially expressed genes for regression
#'
#' This is a wrapper function for the Seurat FindMarkers function,
#' taking as input a Seurat object and other parameters to call differentially
#' expressed genes between two cell groups, then removing genes that are not 
#' in the Seurat's object genome annotation, as CREs to model cannot be selected
#' without gene coordinates. 
#' 
#' @param seurat A Seurat object
#' @param ident.1 Character or character vector of cell identities in group 1
#' @param ident.2 Character or character vector of cell identities in group 2
#' @param group.by Name of metadata column to use to group cells by, or NULL if using stored cell identities
#' @param min.pct min.pct for FindMarkers
#' @param assay Expression assay to use
#' @param slot Expression assay slot to use
#' @param p_val_cutoff Adjusted p value cutoff to filter DEG list
#' @param peak_assay Assay that contains genomic annotation
#' @return List named by inputted gene names, of lists of seeds, each of length num_bootstraps
#'
#' @return Differential gene expression data frame
#' @export
#'
prepare_degs <- function(seurat,
                         ident.1,
                         ident.2,
                         group.by = NULL,
                         min.pct = 0.1,
                         assay = "SCT",
                         slot = "data",
                         p_val_cutoff = 0.05,
                         peak_assay = "peaks", ...) {
	require(Seurat)
	require(Signac)

	if (!is.null(group.by)) Idents(seurat) <- seurat@meta.data[, group.by]
	de_genes <- FindMarkers(seurat, ident.1 = ident.1, ident.2 = ident.2,
									group.by = group.by, min.pct = min.pct,
									assay = assay, slot = slot, ...)
	de_genes <- de_genes[which(de_genes$p_val_adj < p_val_cutoff), ]
	DefaultAssay(seurat) <- peak_assay
	gene_coords <- Signac:::CollapseToLongestTranscript(Annotation(seurat))
	de_genes <- de_genes[which(rownames(de_genes) %in% gene_coords$gene_name), ]

	return(de_genes)
}

#' Generate pseudocell matrix
#'
#' In order to reduce sparsity and noise, you can elect to generate a pseudocell dataset.
#' This function takes as input a Seurat object, the desired assay and slot to use for 
#' pseudocell creation, and the number of cells per partition to target. It pools individual
#' cells into pseudocells with a modified version of VISION's microclustering algorithm to 
#' work on a multiomic WNN graph. Each pseudocell's expression/accessibility is the average of
#' its constituent individual cells' expression/accessibility. The output is a matrix with 
#' a row for each pseudocell, and column for each gene/peak in the specified assay.
#'
#' @param seurat A Seurat object
#' @param assay Expression or accessibility assay to use
#' @param slot Expression or accessibility assay slot to use
#' @param cells_per_partition Targeted number of cells per partition
#' @param find_neighbors Boolean designating whether to run Seurat's FindMultiModalNeighbors. Set to FALSE if this has already been run
#' @param reduction1 Reduction 1 to use for FindMultiModalNeighbors, if find_neighbors is TRUE
#' @param reduction2 Reduction 2 to use for FindMultiModalNeighbors, if find_neighbors is TRUE
#' @param dim_list List of dimensions to use for FindMultiModalNeighbors, if find_neighbors is TRUE
#' @param k.nn Number of nearest neighbors to use for KNN algorithm, if find_neighbors is TRUE
#' @param seed Seed for reproducibility, set to NULL if not needed
#'
#' @return Pseudocell expression/accessibility as a data frame
#' @export
#'
prepare_pseudocell_matrix <- function(seurat,
                                      assay,
                                      slot = "data",
                                      cells_per_partition = 100,
                                      find_neighbors = FALSE,
                                      reduction1 = "harmony",
                                      reduction2 = "harmony_peaks",
                                      non_multiome = FALSE,
                                      umap_name = "WNN.UMAP",
                                      neighbor_name = "weighted.nn",
                                      dims = list(1:50, 1:50),
                                      k.nn = 5,
                                      seed = 489284) {
	require(Seurat)
	require(Signac)
	require(tidyverse)
	require(Matrix)

	if (!is.null(seed)) set.seed(seed)

	if (is.null(reduction2) && non_multiome) {
		require(VISION)
		data <- GetAssayData(seurat, slot = slot, assay = assay)
		mpools <- applyMicroClustering(data, cellsPerPartition = cells_per_partition) # from VISION

		mpool_column <- vector("numeric", dim(seurat)[2])
		names(mpool_column) <- colnames(seurat)
		for (i in 1:length(mpools)) {
		    mpool_column[mpools[[i]]] <- rep(i, length(mpools[[i]]))
		}
		seurat <- AddMetaData(seurat, metadata = as.data.frame(mpool_column))

		pseudocell_matrix <- as_tibble(t(sapply(mpools,
		                                        function(x) {
		                                                    if (length(x) > 1) {
		                                                        rowSums(data[, x]) / length(x);
		                                                    } else {
		                                                        data[, x]
		                                                    }
		                                                }, 
		                                        simplify = TRUE)),
		                               rownames = "pseudocell_IDs")
		pseudocell_mat <- pseudocell_matrix %>% dplyr::select(-pseudocell_IDs) %>% base::as.data.frame()

		return(pseudocell_mat)
	}



	pools <- apply_multi_micro_clustering(seurat,
										  cells_per_partition,
										  find_neighbors,
										  reduction1,
										  reduction2,
										  umap_name,
										  assay,
										  neighbor_name,
										  dims,
										  k.nn)
	pool_column <- vector("numeric", dim(seurat)[2])
	names(pool_column) <- colnames(seurat)
	for (i in 1:length(pools)) {
		pool_column[pools[[i]]] <- rep(i, length(pools[[i]]))
	}
	seurat <- AddMetaData(seurat, metadata = as.data.frame(pool_column))

	pseudocell_mats <- vector("list", length = length(assay))
	names(pseudocell_mats) <- assay
	for (a in assay) {
		mat <- GetAssayData(seurat, assay = a, slot = slot)
		pseudocell_mat <- as_tibble(t(sapply(pools, 
									    	 function(x) {
									    	 		if (length(x) > 1) {
									    	 			Matrix::rowSums(mat[, x]) / length(x)
									    	 		} else {
									    	 			mat[, x]
									    	 		}
									    	 	},
									    	 simplify = TRUE)),
									rownames = "pseudocell_IDs")
		pseudocell_mat <- pseudocell_mat %>% select(-pseudocell_IDs) %>% as.data.frame
		pseudocell_mats[[a]] <- pseudocell_mat
	}

	if (length(assay) == 1) return(pseudocell_mat)

	return(pseudocell_mats)
}

#' Perform microclustering on dataset for pseudocell generation
#'
#' This function takes as input a Seurat object and generates a 
#' weighted nearest neighbor (WNN) graph if needed. It then performs
#' commnunity detection with Louvain clustering and returns a list of pools
#' denoting the individual cells that comprise each pseudocell
#' 
#' @param seurat A Seurat object
#' @param cells_per_partition Targeted number of cells per partition
#' @param find_neighbors Boolean designating whether to run Seurat's FindMultiModalNeighbors. Set to FALSE if this has already been run
#' @param reduction1 Reduction 1 to use for FindMultiModalNeighbors, if find_neighbors is TRUE
#' @param reduction2 Reduction 2 to use for FindMultiModalNeighbors, if find_neighbors is TRUE
#' @param dim_list List of dimensions to use for FindMultiModalNeighbors, if find_neighbors is TRUE
#' @param k.nn Number of nearest neighbors to use for KNN algorithm, if find_neighbors is TRUE
#' @return List of individual cell communities that make up each called pseudocell
#'
#' @return List of individual cell barcodes for each pseudocell
#' @export
#'
apply_multi_micro_clustering <- function(seurat,
                                         cells_per_partition = 100,
                                         find_neighbors = FALSE,
                                         reduction1 = "harmony",
                                         reduction2 = "harmony_peaks",
                                         umap_name = "WNN.UMAP",
                                         assay_name = "SCT",
                                         neighbor_name = "weighted.nn",
                                         dims = list(1:50, 1:50),
                                         k.nn = 5) {
	require(Seurat)
	require(SeuratWrappers)
	require(VISION)
	require(SingleCellExperiment)

	if (find_neighbors || is.null(neighbor_name) || is.null(seurat@neighbors[[neighbor_name]])) {
		if (is.null(reduction2)) {
			seurat <- FindNeighbors(seurat, dims = dims, reduction = reduction1, return.neighbor = TRUE, k.param = k.nn)
		} else {
			seurat <- FindMultiModalNeighbors(seurat, reduction.list = list(reduction1, reduction2), k.nn = k.nn, dims.list = dims)
		}
		neighbor_name <- names(seurat@neighbors)[1]
	}

	kn <- list(seurat@neighbors[[neighbor_name]]@nn.idx, seurat@neighbors[[neighbor_name]]@nn.dist)

	seurat.cds <- as.cell_data_set(seurat, assay = assay_name)
	res <- reducedDim(seurat.cds, type = umap_name, withDimnames = TRUE)

    message("Performing initial coarse-clustering...")
    cl <- VISION:::louvainCluster(kn, res)

    message("Further partitioning coarse clusters...")
    pools <- VISION:::readjust_clusters(cl, res, cellsPerPartition = cells_per_partition)

    # Rename clusters
    cn <- paste0("microcluster_", 1:length(pools))
    names(pools) <- cn

    message(
        sprintf("Micro-pooling completed reducing %i cells into %i pools",
            nrow(res), length(pools))
        )

    return(pools)
}