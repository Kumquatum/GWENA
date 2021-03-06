% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_preprocessing.R
\name{filter_low_var}
\alias{filter_low_var}
\title{Filtering genes with low variability}
\usage{
filter_low_var(data_expr, pct = 0.8, type = c("mean", "median", "mad"))
}
\arguments{
\item{data_expr}{matrix or data.frame or SummarizedExperiment, table of
expression values (either microarray or RNA-seq), with genes as column and
samples as row}

\item{pct}{float, percentage of gene to keep, value must be in ]0;1[}

\item{type}{string, function name used for filtration. Should be either
"mean", "median", or "mad"}
}
\value{
A data.frame of filtered genes
}
\description{
Remove low variating genes based on the percentage given and the type of
variation specified.
}
\examples{
df <- matrix(abs(rnorm(15*45)), 15)
colnames(df) <- paste0("gene_", seq_len(ncol(df)))
rownames(df) <- paste0("sample_", seq_len(nrow(df)))
df_filtered <- filter_low_var(df)

}
