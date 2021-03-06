% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/netrep_functions.R
\name{matchsub}
\alias{matchsub}
\alias{\%nin\%}
\alias{\%sub_in\%}
\alias{\%sub_nin\%}
\title{Functions from NetRep packages}
\usage{
x \%nin\% table

x \%sub_in\% table

x \%sub_nin\% table
}
\arguments{
\item{x}{vector or \code{NULL}: the values to be matched. \link{Long vectors}
are supported.}

\item{table}{vector or \code{NULL}: the values to be matched against.
\link{Long vectors} are not supported.}
}
\value{
A vector of the same length as x.
}
\description{
Sources: https://github.com/sritchie73/NetRep
CRAN page: https://cran.r-project.org/web/packages/NetRep
Package licence: GPL-2
}
\details{
Reason of the copy: the contingency function weren't exported from NetRep
package and neither CRAN or Bioconductor allows to use un-exported function
through the `:::` operator.

Note: functions description have been modified from the original work
Value Matching and Subsetting

This set of functions provides shortcuts for value matching and subsetting,
on top of the functionality provided by \code{\link{\%in\%}}.

\code{\%nin\%} returns a logical vector indicating if elements of \code{x}
are not in \code{table}, This is the opposite of \code{\%in\%}.

\code{\%sub_in\%} returns the elements \code{x} that are \code{\%in\%}
\code{table} rather than a logical vector.

\code{\%sub_nin\%} returns the elements \code{x} that are \code{\%nin\%}
\code{table} rather than a logical vector.
}
\keyword{internal}
