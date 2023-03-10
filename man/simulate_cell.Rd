% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{simulate_cell}
\alias{simulate_cell}
\title{Utility function to simulate a cell based on current expression and trained gene regulatory network model}
\usage{
simulate_cell(cell, model, genes, tfs)
}
\arguments{
\item{cell}{Expression of an individual cell}

\item{model}{Data frame with regulatory coefficients of each TF for each gene}

\item{genes}{Gene names that are present in model}

\item{tfs}{TF names that are present in model}
}
\value{
A cell with updated, simulated expression
}
\description{
This function takes in the expression of a cell and a data frame
consisting of the estimated regulatory coefficients of each TF
for each gene. It simulates new expression values as a result
of perturbing a set of TFs.
}
