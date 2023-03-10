% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/modeling.R
\name{make_model_df}
\alias{make_model_df}
\title{Utility function to convert list of lists run_aen_tf output to data frame}
\usage{
make_model_df(tf_results)
}
\arguments{
\item{tf_results}{Output of run_tf_aen, in list of lists format}
}
\value{
Data frame of regulatory interactions, with TFs by row and modeled genes by column
}
\description{
This function extracts predicted non-zero regulatory interactions
between TFs and genes from the output of run_aen_tf, and stores
the information in a data frame.
}
