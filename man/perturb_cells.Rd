% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{perturb_cells}
\alias{perturb_cells}
\title{Utility function to perturb specific TFs in an idnvidual cell}
\usage{
perturb_cells(
  cell,
  perturbed_tfs_up,
  perturbed_tfs_down,
  mean_var_df,
  how_much = 1
)
}
\arguments{
\item{cell}{Expression of an individual cell}

\item{perturbed_tfs_up}{Names of TFs to upregulate}

\item{perturbed_tfs_down}{Names of TFs to knock down}

\item{mean_var_df}{Data frame with mean expression and standard deviation for modeled genes}

\item{how_much}{Vector of multiples to multiply standard deviation by to determine perturbation amount up or down}
}
\value{
A cell with updated, perturbed expression
}
\description{
This function takes in the expression of a cell and the names
of TFs to either upregulate or knock down. The magnitude of
perturbation is measured in multiples of the standard deviation
of that TF's expression in the dataset. In order to mimic
long term TF upregulation or knock down, higher multiples can be
inputted. Negative gene expression after knock down is set to 0.
The function returns the expression values of the cell with
altered expression for the specified perturbed TFs.
}
