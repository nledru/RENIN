#' Plot motif enrichment volcano plot
#'
#' This function will plot a one-sided volcano plot of motif
#' enrichment in the provided CRE/peak set. Axes are fold enrichment
#' along the x-axis and -log(p-value) along the y-axis. Optional
#' parameters control labeling of specific motifs. The function returns
#' a ggplot, in which the x-axis can be flipped if desired, which
#' may be useful for comparing motif enrichment between two sets of CREs.
#'
#' @param seurat A Seurat object
#' @param cres CRE/peak set to perform motif enrichment on
#' @param motifs_to_label Names of motifs to label on plot
#' @param num_top_label Number of top motifs to label on plot if motifs_to_label is NULL
#' @param label_bottom Boolean to label the least enriched motifs 
#' @param num_bottom_label If labeling least enriched motifs, number of motifs to label
#' @param color Color of points on plot
#' @param flip Boolean to flip the x axis. May be useful to plot motif enrichment of two sets of CREs side-by-side
#'
#' @return ggplot of motif enrichment, fold enrichment vs -log(p-value)
#' @export
#'
plot_motif_enrichment <- function(seurat,
                                  cres,
                                  motifs_to_label = NULL,
                                  num_top_label = 10,
                                  label_bottom = TRUE,
                                  num_bottom_label = 5,
                                  color = "#5862AD",
                                  flip = FALSE) {
	require(Seurat)
	require(Signac)
	require(ggplot2)
	require(ggrepel)

	motifs <- FindMotifs(seurat, features = cres)
	motifs$logp <- -log(motifs$pvalue, 10)
	max_y <- max(motifs$logp) * 1.05
	max_x <- max(motifs$fold.enrichment) * 1.05
	if (is.null(motifs_to_label)) {
		if (label_bottom) {
			motifs_to_label <- c(head(motifs$motif.name, n = num_top_label), tail(motifs$motif.name, n = num_bottom_label))
		} else {
			motifs_to_label <- head(motifs$motif.name, n = num_top_label)
		}
	}
	motifs$label <- ""
	motifs$label[which(motifs$motif.name %in% motifs_to_label)] <- motifs$motif.name[which(motifs$motif.name %in% motifs_to_label)]
	motifs$border <- 0
	motifs$border[which(motifs$label != "")] <- 1

	if (flip) {
		xlimits <- c(-max_x, 0)
	} else {
		xlimits <- c(0, max_x)
	}
	if (flip) motifs$fold.enrichment <- motifs$fold.enrichment * -1

	g <- ggplot(motifs, aes(x = fold.enrichment, y = logp, label = label)) + 
				geom_point(fill = color, color = "black", stroke = motifs$border, pch = 21, size = 1.5) + 
                theme_classic() + geom_text_repel(max.overlaps = 500) + 
                ylab("-log(p_value)") + xlab("Fold enrichment") + ylim(c(0, max_y)) + xlim(xlimits) +
                theme(text = element_text(size = 14), axis.text = element_text(size = 14))
    
    g
}

#' Plot Tn5 integration enrichment between two sets of CREs/peaks
#'
#' This function is modified from Signac's PlotFootprint function. 
#' It plots Tn5 enrichment around inputted motifs in two sets of 
#' peaks, stored in assay.1 and assay.2. Functionality remains
#' preserved, with the exception of plotting Tn5 enrichment
#' for two distinct peak sets, returning plots of Tn5 integration
#' enrichment for each motif
#'
#' @param object A Seurat object
#' @param assay.1 Assay containing first set of peaks
#' @param assay.2 Assay containing second set of peaks
#' @param features TF motifs to plot adjacent Tn5 enrichment
#' @param group.by Grouping variable to use
#' @param idents Identities to group cells by
#' @param label Boolean whether to label groups
#' @param repel Boolean whether to repel labels
#' @param show.expected Boolean whether to plot expected Tn5 integration
#' @param normalization Method of normalization of Tn5 integration bias
#' @param label.top Number of top groups to label
#' @param label.idents Names of identities to label
#' @param colors Colors for each plotted CRE set
#'
#' @return ggplot of Tn5 enrichment plots for each motif
#' @export
#'
plot_footprinting <- function(object,
                              assay.1,
                              assay.2,
                              features,
                              group.by = NULL,
                              idents = NULL,
                              label = TRUE,
                              repel = TRUE,
                              show.expected = TRUE,
                              normalization = "subtract", 
                              label.top = 3,
                              label.idents = NULL,
                              colors = c("#39B54A", "#5862AD")) {
    require(Seurat)
    require(Signac)
    require(ggplot2)
    require(ggrepel)

    colors <- unlist(colors)
    color.1 <- colors[1]
    color.2 <- colors[2]

    plot.data <- GetFootprintData(object = object,
    							  features = features, 
                                  assay = assay.1,
                                  group.by = group.by,
                                  idents = idents)
    plot.data.2 <- GetFootprintData(object = object,
    								 features = features,
                                     assay = assay.2,
                                     group.by = group.by,
                                     idents = idents)
    motif.sizes <- Signac:::GetMotifSize(object = object,
    									 features = features, 
                                         assay = assay.1)

    obs <- plot.data[plot.data$class == "Observed", ]
    obs.2 <- plot.data.2[plot.data.2$class == "Observed",]

    expect <- plot.data[plot.data$class == "Expected", ]
    expect.2 <- plot.data[plot.data$class == "Expected", ]

    base <- ceiling(motif.sizes/2)

    obs$flanks <- sapply(X = seq_len(length.out = nrow(x = obs)), 
                        	function(x) {
	                            pos <- abs(obs[x, "position"])
	                            size <- base[[obs[x, "feature"]]]
	                            return((pos > size) & (pos < (size + 50)))
	                         })
    obs.2$flanks <- sapply(X = seq_len(length.out = nrow(x = obs.2)), 
                            	function(x) {
	                                pos <- abs(obs.2[x, "position"])
	                                size <- base[[obs.2[x, "feature"]]]
	                                return((pos > size) & (pos < (size + 50)))
	                            })

    if (!is.null(normalization)) {
        correction.vec <- expect$norm.value
        correction.vec.2 <- expect.2$norm.value

        names(correction.vec) <- paste(expect$position, expect$feature)
        names(correction.vec.2) <- paste(expect.2$position, expect.2$feature)

        if (normalization == "subtract") {
            obs$norm.value <- obs$norm.value - correction.vec[paste(obs$position, 
                                                                    obs$feature)]
            obs.2$norm.value <- obs.2$norm.value - correction.vec.2[paste(obs.2$position,
                                                                             obs.2$feature)]
        }
        else if (normalization == "divide") {
            obs$norm.value <- obs$norm.value/correction.vec[paste(obs$position, 
                                                                  obs$feature)]
            obs.2$norm.value <- obs.2$norm.value/correction.vec[paste(obs.2$position,
                                                                        obs.2$feature)]
        }
        else {
            stop("Unknown normalization method requested")
        }
    }

    flanks <- obs[obs$flanks, ]
    flanks.2 <- obs.2[obs.2$flanks, ]

    flanks <- group_by(.data = flanks, feature, group)
    flanks.2 <- group_by(.data = flanks.2, feature, group)

    flankmeans <- summarize(.data = flanks, mn = mean(x = norm.value))
    flankmeans.2 <- summarize(.data = flanks.2, mn = mean(x = norm.value))

    topmean <- top_n(x = flankmeans, n = label.top, wt = mn)
    topmean.2 <- top_n(x = flankmeans.2, n = label.top, wt = mn)

    ymax <- top_n(x = flankmeans, n = 1, wt = mn)
    ymax$mn <- ifelse(top_n(x = flankmeans, n = 1, wt = mn)$mn > top_n(x = flankmeans.2, n = 1, wt = mn)$mn,
                      top_n(x = flankmeans, n = 1, wt = mn)$mn,
                      top_n(x = flankmeans.2, n = 1, wt = mn)$mn)

    ymin <- top_n(x = flankmeans, n = 1, wt = -mn)
    ymin$mn <- ifelse(top_n(x = flankmeans, n = 1, wt = -mn)$mn < top_n(x = flankmeans.2, n = 1, wt = -mn)$mn,
                      top_n(x = flankmeans, n = 1, wt = -mn)$mn,
                      top_n(x = flankmeans.2, n = 1, wt = -mn)$mn)

    label.df <- data.frame()
    label.df.2 <- data.frame()

    sub <- obs[obs$position == 75, ]
    sub.2 <- obs.2[obs.2$positio == 75, ]

    for (i in seq_along(along.with = features)) {
        if (is.null(x = label.idents)) {
            groups.use <- topmean[topmean$feature == features[[i]], 
            ]$group
        }
        else {
            groups.use <- label.idents
        }
        df.sub <- sub[(sub$feature == features[[i]]) & (sub$group %in% 
                                                            groups.use), ]
        df.sub.2 <- sub.2[(sub.2$feature == features[[i]]) & (sub.2$group %in%
                                                                     groups.use), ]
        label.df <- rbind(label.df, df.sub)
        label.df.2 <- rbind(label.df.2, df.sub.2)
    }

    obs$label <- NA
    obs.2$label <- NA
    label.df$label <- label.df$group
    label.df.2$label <- label.df.2$group
    obs <- rbind(obs, label.df)
    obs.2 <- rbind(obs.2, label.df.2)

    plotlist <- list()
    for (i in seq_along(along.with = features)) {
        df <- obs[obs$feature == features[[i]], ]
        df.2 <- obs.2[obs.2$feature == features[[i]], ]

        min.use <- -.5
        axis.min <- min(min.use, ymin[ymin$feature == features[[i]], 
        ]$mn)
        axis.max <- ymax[ymax$feature == features[[i]], ]$mn + 0.5
        # p <- ggplot(data = df, mapping = aes(x = position, y = norm.value, 
        #     color = group, label = "H"))
        p <- ggplot()
        p <- p + geom_line(data = df, mapping = aes(x = position, y = norm.value,
                                                    color = "H", label = "H"), color = color.1, size = 0.2)
        p <- p + geom_line(data = df.2, mapping = aes(x = position, y = norm.value,
                                                       color = "FR", label = "FR"), color = color.2, size = 0.2)
        p <- p + xlab("Distance from motif") + 
            ylab(label = "Tn5 insertion\nenrichment") + theme_classic() + 
            ggtitle(label = features[[i]]) + ylim(c(axis.min, axis.max)) + guides(color = guide_legend(override.aes = list(size = 1)))
        #if (label) {
        #    if (repel) {
        p <- p + geom_label_repel(box.padding = 0.5, 
                                  show.legend = TRUE)
        #    }
        #    else {
        #        p <- p + geom_label(show.legend = FALSE)
        #    }
        # }
        #if (show.expected) {
        df <- expect[expect$feature == features[[i]], ]
        p1 <- ggplot(data = df, mapping = aes(x = position, 
                                              y = norm.value)) + geom_line(size = 0.2) + xlab("Distance from motif") + 
            ylab(label = "Expected\nTn5 enrichment") + theme_classic()
        p <- p + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), 
                       axis.line.x.bottom = element_blank(), axis.ticks.x.bottom = element_blank())
        p <- p + p1 + patchwork::plot_layout(ncol = 1, heights = c(3, 1))
        plotlist[[i]] <- p
        #}
    }
    plots <- patchwork::wrap_plots(plotlist)

    return(plots)
}

#' Function to generate graph analysis plots
#'
#' This function will produce plots of graph analysis of predicted
#' gene regulatory networks. It takes in the output of run_tf_aen
#' and generates graph and ranking plots of top TFs. It calculates
#' two centrality measures, betweenness and PageRank, to determine
#' most central TFs to plot. In order to prioritize TFs that are
#' central, rather than targets, for PageRank, edges are reversed.
#' The function returns a graph and ranking plot for each measure
#' of centrality.
#'
#' @param tf_results List output returned by run_tf_aen with trained model information
#' @param regulator_tf_names List of TFs with motif and expression information in dataset
#' @param top_n_graph Number of most central TFs to plot on gene regulatory network graph plots
#' @param top_n_ranking Number of most central TFs to plot on ranking plots
#' @param layout_algorithm Algorithm used to construct gene regulatory network graphs. Options are star, circle, gem, dh, graphopt, grid, mds, randomly, fr, drl, lgl, and kk
#' @param color Color of points on ranking plots
#'
#' @return list of plots, two graphs with top TFs by betweenness or PageRank labeled and two ranking plots of top TFs by betweenness or PageRank
#' @export
#'
plot_graph_rankings <- function(tf_results,
                                regulator_tf_names,
                                seurat,
                                top_n_graph = 20, # >20 can be a little too cluttered, so split this out
                                top_n_ranking = 25,
                                layout_algorithm = "kk",
                                color = "#5862AD") {
	require(Seurat)
	require(igraph)
	require(ggraph)

	if (is.null(regulator_tf_names)) {
		regulator_tf_names <- unique(unlist(lapply(tf_results, function(x) rownames(x[[4]]))))
		regulator_tf_names <- regulator_tf_names[which(regulator_tf_names != "(Intercept)")]
	}

	mean_expr <- Matrix::rowSums(GetAssayData(seurat, assay = "SCT", slot = "data")) / dim(seurat)[2]

	edge_df <- lapply(tf_results, function(x) {
	                    coef_df <- x[[4]]
	                    start <- rownames(coef_df)[which(coef_df[, "coef_if_kept"] != 0)]
	                    start <- start[which(start != "(Intercept)")]
	                    target <- rep(x[[1]], length(start))
	                    coef <- coef_df[start, "coef_if_kept"]
	                    abscoef <- abs(coef)
	                    coefmean <- coef * mean_expr[start]
	                    abscoefmean <- abs(coefmean)
	                    return(data.frame(start = start, target = target, abscoef = abscoef))  
	                 }) %>% bind_rows() %>% as.data.frame

	graph_tf <- graph_from_data_frame(edge_df[, c("start", "target", "abscoef")], directed = TRUE)
	graph_tf_only <- graph_from_data_frame(edge_df[which(edge_df$target %in% regulator_tf_names), c("start", "target", "abscoef") ])
	graph_centrality <- graph_from_data_frame(edge_df[, c("target", "start", "abscoef")], directed = TRUE) # reverse order to prioritize central TFs rather than target genes

	V(graph_tf)$betweenness <- betweenness(graph_tf)
	betweenness_values <- V(graph_tf)$betweenness
	names(betweenness_values) <- names(V(graph_tf))
	pagerank_values <- page_rank(graph_centrality)[[1]]
	V(graph_tf_only)$betweenness <- betweenness_values[names(V(graph_tf_only))] # plotting tf only network for visual purposes, but using calculated betweenness values using any edges to de.genes as well
	V(graph_tf_only)$pagerank <- page_rank(graph_centrality)[[1]][names(V(graph_tf_only))] # same for pagerank
	label_values <- rep(NA, length(names(V(graph_tf_only))))
	names(label_values) <- names(V(graph_tf_only))
	label_pagerank_values <- label_values
	
	tf_only_betweenness <- betweenness_values[names(V(graph_tf_only))]
	tf_only_pagerank <- page_rank(graph_centrality)[[1]][names(V(graph_tf_only))]
	label_values[names(tail(sort(tf_only_betweenness), top_n_graph))] <- names(tail(sort(tf_only_betweenness), top_n_graph))
	label_pagerank_values[names(tail(sort(tf_only_pagerank), top_n_graph))] <- names(tail(sort(tf_only_pagerank), top_n_graph))
	fill_values <- rep("black", length(names(V(graph_tf_only))))
	names(fill_values) <- names(V(graph_tf_only))
	fill_values[names(tail(sort(tf_only_betweenness), top_n_graph))] <- "white"
	V(graph_tf_only)$label <- label_values
	V(graph_tf_only)$labelpage <- label_pagerank_values

	betweenness_graph <- ggraph(graph_tf_only, layout = layout_algorithm) + 
	            geom_edge_link(aes(color = "red"), alpha = .1, color = color) + 
	            geom_node_point(size = V(graph_tf_only)$betweenness / max(V(graph_tf_only)$betweenness) * 10, fill = fill_values, colour = "black", stroke = 1, shape = 21) +
	            geom_node_label(aes(label = label), repel = TRUE, max.overlaps = 50) +
	            theme_void()

	fill_values <- rep("black", length(names(V(graph_tf_only))))
	names(fill_values) <- names(V(graph_tf_only))
	fill_values[names(tail(sort(tf_only_pagerank), top_n_graph))] <- "white"
	pagerank_graph <- ggraph(graph_tf_only, layout = layout_algorithm) + 
	            geom_edge_link(aes(color = "red"), alpha = .1, color = color) + 
	            geom_node_point(size = V(graph_tf_only)$pagerank / max(V(graph_tf_only)$pagerank) * 10, fill = fill_values, colour = "black", stroke = 1, shape = 21) +
	            geom_node_label(aes(label = labelpage), repel = TRUE, max.overlaps = 50) +
	            theme_void()

	betweenness_rank <- ggplot() + geom_point(aes(y = top_n_ranking:1, x = log(sort(V(graph_tf_only)$betweenness, decreasing = TRUE)[1:top_n_ranking])),
	                                           color = color, size = 2.5) + 
	                                xlab("Log(Betweenness)") + ylab(element_blank()) +
	                                scale_y_continuous(breaks = seq(1,top_n_ranking,1), labels = names(V(graph_tf_only))[order(V(graph_tf_only)$betweenness, decreasing = TRUE)][top_n_ranking:1]) + 
	                                theme_classic() + 
	                                theme(text = element_text(size = 14), axis.text = element_text(size = 14)) + NoLegend()

	pagerank_rank <- ggplot() + geom_point(aes(y = top_n_ranking:1, x = sort(V(graph_tf_only)$pagerank, decreasing = TRUE)[1:top_n_ranking],
	                                           fill = names(V(graph_tf_only))[order(V(graph_tf_only)$pagerank, decreasing = TRUE)][1:top_n_ranking]),
	                                           color = color, size = 2.5) + 
	                                xlab("PageRank") + ylab(element_blank()) +
	                                scale_y_continuous(breaks = seq(1,top_n_ranking,1), labels = names(V(graph_tf_only))[order(V(graph_tf_only)$pagerank, decreasing = TRUE)][top_n_ranking:1]) + 
	                                theme_classic() + 
	                                theme(text = element_text(size = 14), axis.text = element_text(size = 14)) + NoLegend()
	results <- list(betweenness_graph, pagerank_graph, betweenness_rank, pagerank_rank)
	names(results) <- c("Betweenness graph", "PageRank graph", "Betweenness ranking", "PageRank ranking")
	return(results)
}

#' Plot TFs ordered by regulatory score
#'
#' This function inputs the output of rank_tfs to plot TF
#' by regulatory scores, ordered along the x-axis. TFs with
#' estimated scores that do not differ from 0 based on 
#' p_value_cutoff or do not meet the score_cutoff threshold 
#' are trimmed to trim the middle of the plot of less relevant
#' TFs. The function returns a ggplot of ordered rankings.
#'
#' @param results_df Data frame output of rank_tfs, with regulatory score and standard error for each TF
#' @param tfs_to_label Names of TFs to label on plot
#' @param p_value_cutoff P value cutoff below which TFs are removed from plot. Set to 1 to exclude no TFs
#' @param score_cutoff Cutoff for absolute value of regulatory score below which TFs are removed from plot. Set to 0 to exclude no TFs
#' @param two_tailed Boolean whether P value cutoff is for a two-tailed test
#' @param top_n_to_label Number of top positive and negative TFs to label on plot if tfs_to_label is NULL
#' @param label_tfs Boolean whether to label TFs on plot
#' @param colors A vector of two color values: color1 for negative and color2 for positive scores
#' @param ident1 Label for negative scores, default FR
#' @param ident2 Label for positive scores, default H
#'
#' @return ggplot of motif enrichment, fold enrichment vs -log(p-value)
#' @export
#'
plot_tf_rankings <- function(results_df,
							 tfs_to_label = NULL,
							 p_value_cutoff = 0.05,
							 score_cutoff = 0.1,
							 two_tailed = TRUE,
							 top_n_to_label = 5,
							 label_tfs = TRUE,
							 colors = c("#5862AD", "#39B54A"), 
							 ident1 = "FR", ident2 = "H") {
	require(ggplot2)
	require(ggrepel)

	results_df$CI <- results_df$SE * qnorm((1 - p_value_cutoff) / ifelse(two_tailed, 1, 2))
	results_df <- results_df[which(abs(results_df$Score) - results_df$CI > 0), ]# remove nonsig TFs by p_value_cutoff
	results_df <- results_df[which(abs(results_df$Score) > score_cutoff), ] # trim middle for plotting with score_cutoff
	results_df <- results_df[order(results_df$Score, decreasing = TRUE), ]

	if (is.null(tfs_to_label)) {
		tfs_to_label <- c(rownames(head(results_df, n = top_n_to_label)),
						  rownames(tail(results_df, n = top_n_to_label)))
	}
	
	results_df$label <- rep("", dim(results_df)[1])

	if (label_tfs) { 
		results_df$label[which(results_df$TF_name %in% tfs_to_label)] <- results_df$TF_name[which(results_df$TF_name %in% tfs_to_label)]
	}

	# allow change the legend label for ident1 and ident2. 
	results_df$comp <- ifelse(results_df$Score > 0, ident2, ident1)
	# fix the order of color with color 1 corresponding to ident 1, and color 2 corresponding to ident 2.
	results_df$comp = factor(results_df$comp, levels = c(ident1, ident2))
	
	if (length(which(results_df$Score < 0)) == 0) {
		colors <- colors[1]
	}

	results_df$axis <- results_df$axis <- nrow(results_df):1
	g <- ggplot(results_df, aes(x=axis, y=Score, fill=comp)) + 
		    geom_bar(stat = "identity", color = "black", size = 0.0, width = 1, alpha = .8) + 
		    scale_fill_manual(values = colors) +
	    	theme_classic() + xlab("TF") + ylab("Predicted regulatory influence") +
		    geom_text_repel(aes(label = label), max.overlaps = 100, size = 5) + 
		    theme(text = element_text(size = 14), axis.text = element_text(size = 14))
	return(g)
}

#' Function to simulate the results of perturbing a set of TFs in Seurat dataset
#'
#' This function inputs a Seurat object, samples cells from specified 
#' cell type, target_celltype, perturbs specified TFs up or down, then
#' simulates the results of those perturbations using the model data frame
#' that contains the estimated regulatory coefficients of each TF for each
#' modeled gene. The function merges the original Seurat object with the 
#' perturbed and simulated cells for plotting.
#' 
#' @param seurat A Seurat object
#' @param target_celltype Starting cell type(s) to sample from, in order to perturb and simulate
#' @param end_celltype Ending cell type(s) to sample from, useful to plot with simulated cells to assess trajectory
#' @param model Data frame with regulatory coefficients of each TF for each gene
#' @param genes Gene names that are present in model
#' @param tfs TF names that are present in model
#' @param perturbed_tfs_up Names of TFs to upregulate
#' @param perturbed_tfs_down Names of TFs to knock down
#' @param mean_var_df Data frame with mean expression and standard deviation for modeled genes
#' @param how_much Vector of multiples to multiply standard deviation by to determine perturbation amount up or down
#' @param num_starting_cells Number of starting cells to sample from target_celltype group to simulate
#' @param set.seed Boolean whether to set seed for reproducibility during testing
#'
#' @return A Seurat object with the perturbed and simulated cells added
#' @export
#'
simulate_dataset <- function(seurat,
                             target_celltype = c(),
                             end_celltype = c(),
                             model = aen_model,
                             genes = rownames(aen_model),
                             tfs = colnames(aen_model),
                             perturbed_tfs_up = c(),
                             perturbed_tfs_down = c(),
                             mean_var_df = mean_var,
                             how_much = c(1, 5, 10),
                             num_starting_cells = 500,
                             set.seed = TRUE) {
    require(Seurat)
    require(tidyverse)

    if (set.seed) set.seed(12345)
    cells <- sample(Cells(seurat)[which(Idents(seurat) %in% end_celltype)], num_starting_cells) %>% as.character()
    cells_expr_ending <- GetAssayData(subset(seurat, cells = cells), assay = "SCT", slot = "data")
    seu_metadata <- seurat@meta.data %>% as.data.frame
    colnames(seu_metadata) <- colnames(seurat@meta.data)
    rownames(seu_metadata) <- Cells(seurat) %>% as.character()
    metadata <- seu_metadata[which(rownames(seu_metadata) %in% cells), ]
    colnames(metadata) <- colnames(seu_metadata)
    metadata$sim <- "ending sample"

    all_genes <- rownames(cells_expr_ending)
    cells <- sample(Cells(seurat)[which(Idents(seurat) %in% target_celltype)], num_starting_cells)
    metadata_add <- seu_metadata[which(rownames(seu_metadata) %in% cells), ]
    cells <- subset(seurat, cells = cells)
    cells_expr <- GetAssayData(cells, assay = "SCT", slot = "data")

    cells_expr_list <- vector("list", dim(cells_expr)[2])
    names(cells_expr_list) <- c(1:length(cells_expr_list))
    for (i in 1:dim(cells_expr)[2]) {
        cells_expr_list[[i]] <- cells_expr[, i]
    }
    seurat_list <- vector("list", length = length(how_much) + 2)
    
    metadata_add$sim <- "starting sample"
    metadata <- rbind(metadata_add, metadata)
    seurat_list[[1]] <- (exp(cells_expr) - 1) %>% as.matrix
    seurat_list[[2]] <- (exp(cells_expr_ending) - 1) %>% as.matrix
    print("initialize done")

    for (i in how_much) {
        sim_expr <- lapply(cells_expr_list, function(x) perturb_cells(x, perturbed_tfs_up, perturbed_tfs_down, mean_var_df, i))
        sim_expr <- lapply(sim_expr, function(x) simulate_cell(x, model, genes, tfs))
        sim_expr_df <- sim_expr %>% bind_cols %>% as.matrix
        rownames(sim_expr_df) <- names(sim_expr[[1]])
        seurat_list[[i + 2]] <- exp(sim_expr_df) - 1
        metadata_add <- seu_metadata[which(Cells(seurat) %in% Cells(cells)), ]
        metadata_add$sim <- paste("sim", i)
        metadata <- rbind(metadata, metadata_add)
        print(paste(i, "std sim done"))
    }

    sim <- bind_cols(seurat_list) %>% as.matrix
    rownames(sim) <- all_genes
    sim <- CreateSeuratObject(counts = sim)
    rownames(metadata) <- Cells(sim)
    sim <- AddMetaData(sim, metadata, col.name = colnames(metadata))
    seurat$sim <- "original"
    
    return(merge(seurat, sim))
}

#' Utility function to perturb specific TFs in an idnvidual cell
#'
#' This function takes in the expression of a cell and the names
#' of TFs to either upregulate or knock down. The magnitude of 
#' perturbation is measured in multiples of the standard deviation
#' of that TF's expression in the dataset. In order to mimic
#' long term TF upregulation or knock down, higher multiples can be
#' inputted. Negative gene expression after knock down is set to 0.
#' The function returns the expression values of the cell with 
#' altered expression for the specified perturbed TFs.
#' 
#' @param cell Expression of an individual cell
#' @param perturbed_tfs_up Names of TFs to upregulate
#' @param perturbed_tfs_down Names of TFs to knock down
#' @param mean_var_df Data frame with mean expression and standard deviation for modeled genes
#' @param how_much Vector of multiples to multiply standard deviation by to determine perturbation amount up or down
#'
#' @return A cell with updated, perturbed expression
#' @export
#'
perturb_cells <- function(cell,
                          perturbed_tfs_up,
                          perturbed_tfs_down,
                          mean_var_df,
                          how_much = 1) {
    cell[perturbed_tfs_down] <- cell[perturbed_tfs_down] - (how_much * mean_var_df["sd", perturbed_tfs_down])
    cell[which(cell < 0)] <- 0

    cell[perturbed_tfs_up] <- cell[perturbed_tfs_up] + (how_much * mean_var_df["sd", perturbed_tfs_up])
    rownames <- names(cell)
    
    cell <- as.numeric(cell)
    names(cell) <- rownames

    cell
}

#' Utility function to simulate a cell based on current expression and trained gene regulatory network model
#'
#' This function takes in the expression of a cell and a data frame
#' consisting of the estimated regulatory coefficients of each TF
#' for each gene. It simulates new expression values as a result
#' of perturbing a set of TFs.
#' 
#' @param cell Expression of an individual cell
#' @param model Data frame with regulatory coefficients of each TF for each gene
#' @param genes Gene names that are present in model
#' @param tfs TF names that are present in model
#'
#' @return A cell with updated, simulated expression
#' @export
#'
simulate_cell <- function(cell, 
                          model,
                          genes,
                          tfs) {
    cell_tfs <- cell[tfs]
    cell_tfs[1] <- 1
    cell[genes] <- model[genes, ] %*% cell_tfs

    cell
}
