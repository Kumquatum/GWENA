#' Functions from NetRep packages
#'
#' Sources: https://github.com/sritchie73/NetRep
#' CRAN page: https://cran.r-project.org/web/packages/NetRep
#' Package licence: GPL-2
#'
#' Reason of the copy: the contingency function weren't exported from NetRep package and neither CRAN or Bioconductor
#' allows to use un-exported function through the `:::` operator.


#' Value Matching and Subsetting
#'
#' @description
#' This set of functions provides shortcuts for value matching and subsetting,
#' on top of the functionality provided by \code{\link[base]{\%in\%}}.
#'
#' \code{\%nin\%} returns a logical vector indicating if elements of \code{x}
#' are not in \code{table}, This is the opposite of \code{\%in\%}.
#'
#' \code{\%sub_in\%} returns the elements \code{x} that are \code{\%in\%}
#' \code{table} rather than a logical vector.
#'
#' \code{\%sub_nin\%} returns the elements \code{x} that are \code{\%nin\%}
#' \code{table} rather than a logical vector.
#'
#' @param x vector or \code{NULL}: the values to be matched. \link[base]{Long vectors}
#'   are supported.
#' @param table vector or \code{NULL}: the values to be matched against.
#'   \link[base]{Long vectors} are not supported.
#' @name matchsub
NULL
#' @rdname matchsub
#' @keywords internal
`%nin%` <- function(x, table) !(x %in% table)
#' @rdname matchsub
#' @keywords internal
`%sub_in%` <- function(x, table) x[x %in% table]
#' @rdname matchsub
#' @keywords internal
`%sub_nin%` <- function(x, table) x[x %nin% table]


#' Order the module vector numerically
#'
#' The module assingments may be numeric, but coded as characters.
#'
#' @param vec module vector to order
#'
#' @return the order of the vector
#'
#' @keywords internal
orderAsNumeric <- function(vec) {
  tryCatch({
    order(as.integer(vec))
  }, warning = function(w) {
    order(vec)
  })
}


#' Calculate a contigency table of module overlap between datasets
#'
#' @param modAssignments a list where the first element is the
#'  'moduleAssignments' vector in the discovery dataset, and the second element
#'  is the 'moduleAssignments vector in the test dataset.
#' @param mods the 'modules' vector for the discovery dataset.
#' @param tiNodelist a vector of node IDs in the test dataset.
#'
#' @return
#'  A list containing a contigency table, a vector of the proportion of
#'  nodes present in the test dataset for each module, a vector containing
#'  the number of nodes present in the test dataset for each module, a vector
#'  of the node names present in both the discovery and test datasets, a vector
#'  of modules that are both requested and have nodes present in the test
#'  dataset, and the \code{modAssignments} vector containing only nodes present
#'  in the test dataset.
#'
#' @export

contingencyTable <- function(modAssignments, mods, tiNodelist) {
  # To simplify later function calls, we need to get a vector of module
  # assignments only for (a) modules of interest and (b) the variables
  # present in both datasets for those modules.
  overlapVars <- intersect(names(modAssignments[[1]]), tiNodelist)
  overlapAssignments <- modAssignments[[1]][overlapVars]

  overlapAssignments <- overlapAssignments %sub_in% mods
  overlapModules <- unique(overlapAssignments)
  overlapModules <- overlapModules[orderAsNumeric(overlapModules)]

  # How many variables are present in the test dataset for the modules
  # of interest?
  varsPres <- table(overlapAssignments)

  modulesWithNoOverlap <- mods %sub_nin% overlapAssignments
  varsPres <- c(varsPres, rep(0, length(modulesWithNoOverlap)))
  names(varsPres)[names(varsPres) == ""] <- modulesWithNoOverlap
  varsPres <- varsPres[orderAsNumeric(names(varsPres))]

  if (any(varsPres == 0)) {
    noNodes <- names(varsPres[varsPres == 0])
    warning(
      "None of the nodes in module(s) ",
      paste(paste0('"', noNodes, '"'), collapse = ", "),
      " are present in the test dataset",
      immediate. = TRUE
    )
  }

  # What proportion?
  moduleSizes <- table(modAssignments[[1]])
  moduleSizes <- moduleSizes[names(moduleSizes) %sub_in% mods]
  propVarsPres <- varsPres / moduleSizes[names(varsPres)]

  # Calculate some basic cross-tabulation statistics so we can assess
  # which modules in both datasets map to each other, if module
  # detection has also been performed for the test network
  contingency <- NULL
  if (!is.null(modAssignments[[2]])) {
    # Get total number of nodes from each discovery subset in each test subset
    contingency <- table(
      modAssignments[[1]][overlapVars],
      modAssignments[[2]][overlapVars]
    )

    ## Add in the sizes of each module, and the number of variables present in
    ## the other dataset

    # Get the module sizes irrespective of dataset overlap
    discSizes <- table(modAssignments[[1]])
    testSizes <- table(modAssignments[[2]])

    # The dummy row / column handles modules with no variables present
    discInfo <- cbind(discSizes, 0)
    testInfo <- rbind(testSizes, 0)
    colnames(discInfo) <- c("size", "present")
    rownames(testInfo) <- c("size", "present")

    # Add the number of variables present for each module in the other dataset
    discPresent <- table(modAssignments[[1]][overlapVars])
    testPresent <- table(modAssignments[[2]][overlapVars])
    discInfo[names(discPresent), "present"] <- discPresent
    testInfo["present", names(testPresent)] <- testPresent

    # Add missing modules to the contigency table
    discMiss <- rownames(discInfo) %sub_nin% rownames(contingency)
    testMiss <- colnames(testInfo) %sub_nin% colnames(contingency)
    discMissMat <- matrix(0, nrow = length(discMiss), ncol = ncol(contingency),
                          dimnames = list(discMiss, colnames(contingency)))
    contingency <- rbind(contingency, discMissMat)
    testMissMat <- matrix(0, ncol = length(testMiss), nrow = nrow(contingency),
                          dimnames = list(rownames(contingency), testMiss))
    contingency <- cbind(contingency, testMissMat)

    # Order and filter contingency table prior to adding information
    contingency <- contingency[mods,, drop = FALSE]
    contingency <- contingency[orderAsNumeric(rownames(contingency)),
                               orderAsNumeric(colnames(contingency)),
                               drop = FALSE]

    # Make sure the info tables are in the same order so we can add them
    # correctly
    discInfo <- discInfo[rownames(contingency),, drop = FALSE]
    testInfo <- testInfo[, colnames(contingency), drop = FALSE]

    # Create the empty upper left square
    na_mat <- matrix(NA, nrow = 2, ncol = 2)
    dimnames(na_mat) <- list(c("size", "present"), c("size", "present"))

    # combine
    contingency <- cbind(discInfo, contingency)
    contingency <- rbind(cbind(na_mat, testInfo), contingency)
  }

  return(list(
    contingency = contingency, propVarsPres = propVarsPres, overlapVars = overlapVars,
    varsPres = varsPres, overlapModules = overlapModules,
    overlapAssignments = overlapAssignments
  ))
}
