% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{kuehne_expr}
\alias{kuehne_expr}
\title{Transcriptomic data from the Kuehne et al. publication}
\format{
A data frame with 48 rows (samples) and 15801 columns (genes).
}
\source{
\url{https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-
017-3547-3}
}
\usage{
kuehne_expr
}
\description{
A dataset containing the expression levels collapsed to the gene level.
Obtained from script provided in additional data n°10 runned on GSE85358 and
reduced from probe to gene by WGCNA::collapseRows with median as fucntion.
}
\keyword{datasets}
