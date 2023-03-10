% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{simulate_dataset}
\alias{simulate_dataset}
\title{Function to simulate the results of perturbing a set of TFs in Seurat dataset}
\usage{
simulate_dataset(
  seurat,
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
  set.seed = TRUE
)
}
\arguments{
\item{seurat}{A Seurat object}

\item{target_celltype}{Starting cell type(s) to sample from, in order to perturb and simulate}

\item{end_celltype}{Ending cell type(s) to sample from, useful to plot with simulated cells to assess trajectory}

\item{model}{Data frame with regulatory coefficients of each TF for each gene}

\item{genes}{Gene names that are present in model}

\item{tfs}{TF names that are present in model}

\item{perturbed_tfs_up}{Names of TFs to upregulate}

\item{perturbed_tfs_down}{Names of TFs to knock down}

\item{mean_var_df}{Data frame with mean expression and standard deviation for modeled genes}

\item{how_much}{Vector of multiples to multiply standard deviation by to determine perturbation amount up or down}

\item{num_starting_cells}{Number of starting cells to sample from target_celltype group to simulate}

\item{set.seed}{Boolean whether to set seed for reproducibility during testing}
}
\value{
A Seurat object with the perturbed and simulated cells added
}
\description{
This function inputs a Seurat object, samples cells from specified
cell type, target_celltype, perturbs specified TFs up or down, then
simulates the results of those perturbations using the model data frame
that contains the estimated regulatory coefficients of each TF for each
modeled gene. The function merges the original Seurat object with the
perturbed and simulated cells for plotting.
}
