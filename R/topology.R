# Removing errors about magrittr's placeholder `.` not declared
utils::globalVariables(c("."))

# Removing errors about dplyr data-variables
utils::globalVariables(c("gene_from", "gene_to"))

#' Return graph from squared matrix network
#'
#' Takes a squared matrix containing the pairwise similarity scores for each
#' gene and return a igraph object.
#'
#' @param sq_mat matrix or data.frame, squared matrix representing
#'
#' @importFrom igraph graph_from_data_frame
#' @importFrom dplyr filter
#' @importFrom magrittr %>%
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr pivot_longer
#'
#' @return An igraph object
#'
#' @examples
#' mat <- matrix(runif(40*40), 40)
#' build_graph_from_sq_mat(mat)
#'
#' @export

build_graph_from_sq_mat <- function(sq_mat) {
  # Checks
  if (!(is.matrix(sq_mat) | is.data.frame((sq_mat))))
    stop("sq_mat should be a matrix or a data.frame")
  if (any(is.na(sq_mat)))
    warning("sq_mat should not contain any missing value")
  if ((any(sq_mat > 1) | any(sq_mat < -1)) & !any(is.na(sq_mat)))
    stop("sq_mat should be filled with value in the [-1,1] range")
  if (nrow(sq_mat) != ncol(sq_mat))
    stop("sq_mat must be a squared matrix")

  # From matrix to graph
  sq_mat %>%
    as.data.frame %>%
    tibble::rownames_to_column("gene_from") %>%
    tidyr::pivot_longer(-gene_from,
                        names_to = "gene_to",
                        values_to = "weight") %>%
    # removing pairwise duplicates
    .[!duplicated(t(apply(.[, c("gene_from", "gene_to")], 1, sort))), ] %>%
    as.data.frame %>%
    # removing self edge
    dplyr::filter(gene_from != gene_to) %>%
    igraph::graph_from_data_frame(directed = FALSE)
}


#' Determine hub genes based on connectivity
#'
#' Compute connectivity of each gene by module if provided or for whole network
#' if not, and return the top_n highest connected ones.
#'
#' @param network matrix or data.frame, square table representing connectivity
#' between each genes as returned by build_net. Can be whole network or a
#' single module.
#' @param modules list, modules defined as list of gene vectors. If null,
#' network is supposed to be the whole network or an already split module
#' @param top_n integer, number of genes to be considered as hub genes
#'
#' @return A list of vectors, or single vector of gene names
#'
#' @importFrom magrittr %>%
#'
#' @examples
#' mat <- matrix(runif(40*40), 40)
#' colnames(mat) <- paste0("gene_", seq_len(ncol(mat)))
#' rownames(mat) <- paste0("gene_", seq_len(nrow(mat)))
#' get_hub_high_co(mat)
#'
#' @export

get_hub_high_co <- function(network, modules = NULL, top_n = 5) {
  # Checks
  .check_network(network)
  if (!is.null(modules)) {
    .check_module(modules, is.list(modules))}
  if (length(top_n) > 1) stop("top_n must be a single numeric value")
  if (!is.numeric(top_n)) stop("top_n must be a numeric value")
  if (top_n < 1 | top_n %% 1 != 0) stop("If not NULL, block_size must be a ",
  " whole number >= 1")

  # Highly connected genes selection
  # Considered network is already a single module, or looking for hubs
  # independently of modules split
  if (is.null(modules)) {
    if (top_n > ncol(network)) stop("top_n > to number of genes into network")
    modules <- list(colnames(network))
  }

  hubs <- lapply(modules, function(x){
    net <- network[x,x]
    if (top_n > ncol(net)) {
      warning("top_n > to number of genes into one of the modules. Taking all",
              " genes in this case.")
      top_n <- ncol(net)
    } else {
      hubs_name <- rowSums(net) %>%
        sort(decreasing = TRUE) %>%
        .[seq_len(top_n)]
    }
  })
  return(hubs)
}

#' Determine hub genes based on degree
#'
#' Remove edges from the graph which value is under weight_th then compute
#' degree of each node (gene). Hub gene are genes whose degree value is above
#' average degree value of the thresholded network.
#'
#' @param network matrix or data.frame, square table representing connectivity
#' between each genes as returned by build_net. Can be whole network or a
#' single module.
#' @param modules list, modules defined as list of gene vectors. If null,
#' network is supposed to be the whole network or an already split module
#' @param weight_th decimal, weight threshold under or equal to which edges
#' will be removed
#'
#' @details GWENA natively build networks using WGCNA. These networks are
#' complete in a graph theory sens, meaning all nodes are connected to each
#' other. Therefore a threshold need to be applied so degree of all nodes
#' isn't the same.
#'
#' @importFrom igraph delete.edges degree E
#'
#' @return A list of vectors, or single vector of gene names
#'
#' @examples
#' mat <- matrix(runif(40*40), 40)
#' colnames(mat) <- paste0("gene_", seq_len(ncol(mat)))
#' rownames(mat) <- paste0("gene_", seq_len(nrow(mat)))
#' get_hub_degree(mat)
#'
#' @export

get_hub_degree <- function(network, modules = NULL, weight_th = 0.2) {
  # Checks
  .check_network(network)
  if (!is.null(modules)) {
    .check_module(modules, is.list(modules))}
  if (!is.null(weight_th)) {
    if (length(weight_th) > 1) stop("weight_th must be a single numeric value")
    if (!is.numeric(weight_th)) stop("weight_th must be a numeric value")
    if (weight_th < 0 | weight_th >= 1) stop("weight_th must be a in [0;1[")
  }
  # TODO : check if can avoid transforming to igraph object.
  # It takes a lot of time...

  # Above average degree genes
  # Considered network is already a single module, or looking for hubs
  # independently of modules split
  if (is.null(modules)) {
    modules <- list(colnames(network))
  }
  hubs <- lapply(modules, function(x){
    net_degree <- network[x,x] %>%
      build_graph_from_sq_mat() %>%
      igraph::delete.edges(which(igraph::E(.)$weight <= weight_th)) %>%
      igraph::degree()
    if (table(net_degree) %>% length == 1) {
      x_hubs <- c()
    } else {
      x_hubs <- net_degree[which(net_degree > (net_degree %>% mean))] }

  })

  return(hubs)
}


#' Determine hub genes based on Kleinberg's score
#'
#' Compute Kleinberg's score (defined as the principal eigenvector of A*t(A),
#' where A is the similarity matrix of the graph) of each gene by module if
#' provided or for whole network if not, and return the top_n highest ones.
#'
#' @param network matrix or data.frame, square table representing connectivity
#' between each genes as returned by build_net. Can be whole network or a
#' single module.
#' @param modules list, modules defined as list of gene vectors. If null,
#' network is supposed to be the whole network or an already split module
#' @param top_n integer, number genes to be considered as hub genes
#' @param k_th decimal, Kleinberg's score threshold above or equal to which
#' genes are considered as hubs
#'
#' @details If you provide a top_n value, you can't provide a k_th value and
#' vice versa. If none of them is provided, top_n = 5.
#' For more information on Kleinberg's score, look at
#' \code{\link[igraph]{hub_score}} from igraph.
#'
#' @importFrom igraph hub_score
#'
#' @return A list of vectors, or single vector of gene names
#'
#' @examples
#' mat <- matrix(runif(40*40), 40)
#' colnames(mat) <- paste0("gene_", seq_len(ncol(mat)))
#' rownames(mat) <- paste0("gene_", seq_len(nrow(mat)))
#' get_hub_degree(mat)
#' get_hub_kleinberg(mat, top_n = NULL, k_th = 0.9)
#'
#' @export

get_hub_kleinberg <- function(network, modules = NULL, top_n = 5,
                              k_th = NULL) {
  # Checks
  .check_network(network)
  if (!is.null(modules)) {
    .check_module(modules, is.list(modules))}
  if (!is.null(top_n) & !is.null(k_th)) {
    k_th <- NULL
    warning("Conflict: top_n and k_th value provided.",
    " Keeping only top_n value") }
  if (!is.null(k_th)) {
    if (length(k_th) > 1) stop("k_th must be a single numeric value")
    if (!is.numeric(k_th)) stop("k_th must be a numeric value")
    if (k_th < 0 | k_th >= 1) stop("k_th must be a in [0;1[")
  }
  if (!is.null(top_n)) {
    if (length(top_n) > 1) stop("top_n must be a single numeric value")
    if (!is.numeric(top_n)) stop("top_n must be a numeric value")
    if (top_n < 1 | top_n %% 1 != 0)
      stop("If not NULL, block_size must be a whole number >= 1")
  }

  # Considered network is already a single module, or looking for hubs
  # independently of modules split
  if (is.null(modules)) {
    modules <- list(module = colnames(network))
  }
  hubs <- lapply(modules, function(x){
    net_hub_score <- network[x,x] %>%
      build_graph_from_sq_mat() %>%
      igraph::hub_score(scale = TRUE) %>%
      .$vector

    # With top_n
    if (!is.null(top_n)) {
      x_hubs <- net_hub_score %>% sort(decreasing = TRUE) %>% .[seq_len(top_n)]
    } else { # With k_th
      x_hubs <- net_hub_score[which(net_hub_score >= k_th)]
    }
    return(x_hubs)
  })

  return(hubs)
}


#' Determine hub genes inside each module
#'
#' Return genes considered as hub genes inside each module of a network
#' following the selected method. Method will be lauched with default
#' parameters. If specific parameters desired, please use directly the
#' function \code{get_hub_...} itself.
#'
#' @param network matrix or data.frame, square table representing connectivity
#' between each genes as returned by build_net. Can be whole network or a
#' single module.
#' @param modules list, modules defined as list of gene vectors. If null,
#' network is supposed to be the whole network or an already split module
#' @param method string, name of the method to be used for hub gene detection.
#' See details.
#'
#' @details
#' \describe{
#'   \item{highest connectivity}{Select the top n (n depending on
#'   parameter given) highest connected genes. Similar to
#'   WGCNA::chooseTopHubInEachModule.}
#'   \item{superior degree}{Select genes which degree is greater than average
#'   connection degree of the network. Definition from network theory.}
#'   \item{Kleinberg's score}{Select genes which Kleinberg's score superior to
#'   provided threshold.}
#' }
#'
#' @return A list of vectors representing hub genes, by module
#'
#' @examples
#' mat <- matrix(runif(40*40), 40)
#' colnames(mat) <- paste0("gene_", seq_len(ncol(mat)))
#' rownames(mat) <- paste0("gene_", seq_len(nrow(mat)))
#' get_hub_genes(mat)
#'
#' @export

get_hub_genes <- function(network, modules = NULL, method =
                            c("highest connectivity", "superior degree",
                              "Kleinberg's score")) {
  # Checks
  method = match.arg(method)

  # Call
  if (method == "highest connectivity") {
    hubs <- get_hub_high_co(network, modules)
  } else if (method == "superior degree") {
    hubs <- get_hub_degree(network, modules)
  } else if (method == "Kleinberg's score") {
    hubs <- get_hub_kleinberg(network, modules) }

  return(hubs)
}


# Removing errors about dplyr data-variables
utils::globalVariables(c("vertex.size", "edge.width"))

#' Plot co-expression network
#'
#' Display a graph representing the co-expression network and different
#' informations like hubs, enrichments
#'
#' @param graph_module igraph object, module to plot.
#' @param hubs character vector or numeric vector with names, optionnal,
#' vector of gene names or vector of numeric values named with gene names.
#' @param groups matrix or data.frame, a two cols table with the gene id in the
#' first one, and the group assignation in the second one.
#' @param lower_weight_th,upper_weight_th decimal, weight threshold above
#' lower_weight_th or below upper_weight_th which edges will be removed.
### @param enrichment list, representing a gost object. # Not yet implemented
#' @param title string, main title that will be displayed on the plot.
#' @param degree_node_scaling boolean, indicates if node size should represent
#' the degree of this node.
#' @param node_scaling_min,node_scaling_max integer, if degree_node_scaling is
#' TRUE, it is the min/max size of the node, else it is the exact size of all
#' node.
#' @param edge_scaling_min,edge_scaling_max integer, min/max width of the edge
#' @param nb_row_legend integer, number of levels in the legend.
#' @param layout numeric matrix or function or string, numeric matrix for nodes
#' coordinates, or function for layout, or name of a layout function available
#' in \code{igraph}. Default "auto" will choose the best layout depending on
#' the graph. For more information, see \code{\link{igraph.plotting}}.
#' @param zoom integer, scaling factor by which it's possible to have
#' compact graph (< 1) or larger graph (> 1) display.
#' @param vertex.label.cex,legend_cex float, font size for vertex labels. It is
#' interpreted as a multiplication factor of some device-dependent base font
#' size. If 0, no labels displayed.
#' @param vertex.label.color,edge.color,vertex.frame.color,vertex.color
#' character and/or integer vector , color of the labels.
#' It may either contain integer values, named colors or RGB specified colors
#' with three or four bytes. All strings starting with ‘#’ are assumed to be
#' RGB color specifications. It is possible to mix named color and RGB colors.
#' @param vertex.label.family character, font family to be used for vertex
#' labels.
#' @param vertex.label.dist integer, distance of the label from the center of
#' the vertex. If it is 0 then the label is centered on the vertex. If it is 1
#' then the label is displayed beside the vertex.
#' @param groups_palette character and/or integer vector, vertices group palette
#' of colors for the groups specified. It may either contain integer values, 
#' named colors or RGB specified colors with three or four bytes. All strings 
#' starting with ‘#’ are assumed to be RGB color specifications. It is possible 
#' to mix named color and RGB colors.
#' @param window_x_min decimal, value for the bottom limit of the window.
#' @param window_x_max decimal, value for the top limit of the window.
#' @param window_y_min decimal, value for the left limit of the window.
#' @param window_y_max decimal, value for the right limit of the window.
#' @param legend boolean, indicates if the legend should be plotted.
# #' @param dymamic_plot boolean, does the final plot should be dynamic (plotly
# #' plot) or static (ggplot2 plot).
#' @param ... any other parameter compatible with the
#' \code{\link[igraph]{plot.igraph}} function.
#'
#' @details Take care of you intend to compare modules' graphs, the same size
#' of node will not correspond to the same values because of the scaling.
#'
#' @import igraph
#' @importFrom magrittr %>%
#' @importFrom graphics legend symbols
#' @importFrom utils lsf.str
#'
#' @return matrix, layout of the graph as a two column matrix (x, y)
#'
#' @examples
#' mat <- matrix(runif(40*40), 40)
#' g <- build_graph_from_sq_mat(mat)
#' plot_module(g, lower_weight_th = -0.5, upper_weight_th = 0.5)
#'
#' @export

plot_module <- function(graph_module, hubs = NULL, groups = NULL,
                        lower_weight_th = NULL, upper_weight_th = NULL,
                        # enrichment = NULL,
                        # dynamic_plot = TRUE,
                        title = "Module", degree_node_scaling = TRUE,
                        node_scaling_min = 1, edge_scaling_min = 0.2,
                        node_scaling_max = 6, edge_scaling_max = 1,
                        nb_row_legend = 6, layout = "auto", zoom = 1,
                        vertex.label.cex = 0.7,
                        vertex.label.color = "gray20",
                        vertex.label.family = "Helvetica",
                        edge.color = "gray70", vertex.frame.color = "white",
                        vertex.color = "gray60", vertex.label.dist = 1, 
                        legend_cex = 0.8, groups_palette = NULL,
                        window_x_min = -1, window_x_max = 1,
                        window_y_min = -1, window_y_max = 1,
                        legend = TRUE, ...) {
  # TODO: add a layer arg that could take a list where each named element could
  # be list of genes of interest (like one displaying gene known in aging and
  # another for inflammation)

  # Checking args
  if (!igraph::is.igraph(graph_module))
    stop("graph_module must be an igraph object")
  named_num_vec <- char_vec <- FALSE
  if (!is.null(hubs)) {
    # Named vector with numeric values
    if (is.vector(hubs, "numeric") & !is.null(names(hubs)))
      named_num_vec <- TRUE
    # Vector of characters
    if (is.vector(hubs, "character")) char_vec <- TRUE
    if (!(named_num_vec & char_vec))
      stop("hubs must be a named vector of numeric values, or a vector of",
      " characters")
  }
  if (!is.null(groups)) {
    if (!is.data.frame(groups) & !is.matrix(groups))
      stop("groups must be a matrix or a data.frame")
    if (length(groups) != 2)
      stop("groups must be a two columns matrix with genes id and group",
           " assignation")
  }
  if (!is.null(groups_palette)){
    if (is.null(groups)) 
      stop("groups need to be provided to use groups_palette")
  }
  lapply(c(lower_weight_th, upper_weight_th), function(weight_th){
    if (!is.null(weight_th)) {
      if (length(weight_th) > 1)
        stop("Weight threshold must be a single numeric value")
      if (!is.numeric(weight_th))
        stop("Weight threshold must be a numeric value")
      # if (weight_th <= -1 | weight_th >= 1)
      #   stop("Weight thresholds must be a in ]-1;1[")
    }
  })
  # TODO: Change when implemented
  # if (!is.null(enrichment)) .check_gost(enrichment)
  if (!is.character(layout) & !is.matrix(layout) & !is.function(layout)) {
    stop("layout must be a layout function, its name as a string, or a matrix",
    " giving position of each node ") }
  if (!is.numeric(zoom))
    stop("zoom must be a numeric value superior to 0.")
  if (!is.logical(degree_node_scaling))
    stop("degree_node_scaling must be a boolean value")
  if (isTRUE(degree_node_scaling) & exists("vertex.size"))
    stop("If degree_node_scaling is TRUE, you cannot specify a vertex.size")
  if (!is.numeric(c(node_scaling_max, edge_scaling_max,
                    node_scaling_min, edge_scaling_min)))
      stop("All scaling values must be a numeric value")
  if (any(node_scaling_max < node_scaling_min,
          edge_scaling_max < edge_scaling_min))
      stop("Scalling values must have a min inferior to the corresponding max")
  if (!is.null(edge_scaling_max) & exists("edge.width"))
    stop("If edge_scaling_max is not NULL, you cannot specify a edge.width")

  # Keeping only gene names if hubs is a named numeric vector
  if (named_num_vec) hubs <- names(hubs)
  if (!(all(hubs %in% igraph::V(graph_module)$names)))
    stop("Not all hubs are in graph_module")

  # Removing edges whose weight is below or above thresholds given
  if (!is.null(lower_weight_th) & !is.null(upper_weight_th)) {
    graph_to_plot <- graph_module %>%
      igraph::delete.edges(which(igraph::E(.)$weight > lower_weight_th &
                                   igraph::E(.)$weight < upper_weight_th))
  } else if (!is.null(lower_weight_th)) {
    graph_to_plot <- graph_module %>%
      igraph::delete.edges(which(igraph::E(.)$weight > lower_weight_th))
  } else if (!is.null(upper_weight_th)) {
    graph_to_plot <- graph_module %>%
      igraph::delete.edges(which(igraph::E(.)$weight < upper_weight_th))
  } else graph_to_plot <- graph_module

  # Détermining layout
  if (is.character(layout)) {
    if (layout != "auto") {
      # Checking if layout function name exists
      igraph_layouts <- grep("^layout_\\w[\\w|_]*",
                             utils::lsf.str("package:igraph"),
                             value = TRUE)
      if (!any(layout %in% igraph_layouts))
        stop("layout name provided not found in igraph layout functions")
      l <- igraph::layout_(graph_to_plot, get(layout))
    } else {
      l <- igraph::layout_nicely(graph_to_plot)}
  } else if (is.function(layout)) {
    l <- layout(graph_to_plot)
  } else if (is.matrix(layout)){
    l <- layout
  } else stop("Should never be triggered.")

  # Normalizing boundaries to ensure right plotting
  if (zoom != 1) {
      l <- igraph::norm_coords(l, xmin = window_x_min, xmax = window_x_max,
                          ymin = window_y_min, ymax = window_y_max) }

  # Node scaling on degree or provided value
  if (degree_node_scaling) {
    deg <- igraph::degree(graph_to_plot)
    # Checking deg is the same for all (meaning graph is fully connected)
    delta_deg <- max(deg) - min(deg)
    if (delta_deg == 0) {
      warning("max and min degree of the nodes are equal. ",
      "Consider changing the weight thresholds value.")
      delta_deg <- 1
    }
    vertex_size <- lapply(deg, function(x) {
      (x - min(deg)) / delta_deg *
      (node_scaling_max - node_scaling_min) +
      node_scaling_min }) %>% unlist
  } else {
    if (exists("vertex.size")) {
      vertex_size <- vertex.size
    } else {
      vertex_size <- node_scaling_max
    }
  }

  # Edge scaling if no edge width provided
  if (exists("edge.width")) {
    edge_width <- edge.width
  } else {
    edge <- igraph::E(graph_to_plot)$weight
    edge_width <- lapply(edge, function(x) {
      (x - min(edge)) / (max(edge) - min(edge)) *
      (edge_scaling_max - edge_scaling_min) +
      edge_scaling_min}) %>% unlist
  }
  
  
  # FIXME: hotfix to manage no label
  if (vertex.label.cex == 0) vertex_label = NA else vertex_label = NULL

  # Plotting static or dynamic graph
  # TODO: implementer
  # if (dynamic_plot) {
  # 
  # 
  # } else {
    # Vertex color if groups provided
    if (!is.null(groups)) {
      # No palette specified
      if (is.null(groups_palette)) {
        graph_to_plot$palette <- RColorBrewer::brewer.pal(
          length(unique(groups[,2])), "Paired")
      # Palette name specified
      } else if (is.character(groups_palette) & length(groups_palette) == 1) {
        graph_to_plot$palette <- RColorBrewer::brewer.pal(
          length(unique(groups[,2])), groups_palette)
      # Palette of custom colors specified
      } else {
        if (length(unique(groups[,2])) > length(groups_palette)) {
          warning("Number of groups superior to number of colors in the",
                  "palette. Using ggplot2 like colors instead.")
          graph_to_plot$palette <- gg_palette(length(unique(groups[,2])))
        } else {
          graph_to_plot$palette <- groups_palette
        }
      }
      graph_to_plot <- igraph::set_vertex_attr(graph_to_plot, "color", 
                                               groups[,1], groups[,2])
      vertex.color <- NULL
    }
    
    igraph::plot.igraph(graph_to_plot,
                        vertex.label.color = vertex.label.color,
                        vertex.label.family = vertex.label.family,
                        vertex.label.cex = vertex.label.cex,
                        vertex.label.dist = vertex.label.dist,
                        edge.color = edge.color,
                        vertex.label = vertex_label,
                        vertex.frame.color = vertex.frame.color,
                        vertex.color = vertex.color,
                        vertex.size = vertex_size,
                        edge.width = edge_width,
                        layout = l * zoom,
                        rescale = ifelse(zoom != 1, FALSE, TRUE),
                        main = title,
                        ...)


    # Adding a legend
    if (legend) {
      # Legend nodes size
      scale_vertex_size <- seq(node_scaling_min, node_scaling_max,
                               length.out = nb_row_legend)
      if (degree_node_scaling) {
        label_legend_node <- seq(
          min(deg), max(deg), length.out = nb_row_legend) %>%
          round %>%
          as.character
      } else {label_legend_node <- as.character(scale_vertex_size)}

      legend_nodes_size <- graphics::legend(
        "bottomright", label_legend_node,
        pt.cex = scale_vertex_size/200, col = 'white', 
        pch = 21, title = "Degree", bty = "n", 
        y.intersp = 0.7, cex = legend_cex)
      x_nodes_size <- (legend_nodes_size$text$x + legend_nodes_size$rect$left)/2
      y_nodes_size <- legend_nodes_size$text$y
      graphics::symbols(x_nodes_size, y_nodes_size, 
                        circles = scale_vertex_size/200, inches = FALSE,
                        add = TRUE, bg = 'gray', fg = 'white')
      
      # Legend nodes colors if groups provided
      if (!is.null(groups)) {
        graphics::legend("topleft", as.character(unique(groups[ ,2])),
                         pt.bg = graph_to_plot$palette, title = names(groups)[2], 
                         pt.cex = 2, pch = 21, col = "white",
                         bty = "n", cex = legend_cex)
      }

      # Legend vertex
      scale_edge_width <- seq(edge_scaling_min, edge_scaling_max,
                              length.out = nb_row_legend) %>% signif
      label_legend_edge <- seq(
        min(edge), max(edge), length.out = nb_row_legend) %>%
        signif %>%
        as.character
      graphics::legend("topright", label_legend_edge, col='gray', title = "Weight",
                       lwd = scale_edge_width, bty = "n", cex = legend_cex)
    }
  # }
  return(l)
}


# Removing errors about dplyr data-variables
utils::globalVariables(c("avg_sil_width", "k", "sub_module"))

#' Detect sub clusters
#'
#' Use a partitioning around medoid (PAM, or k-medoid) clustering method to
#' detect clusters into a provided module using the strength matrix of the
#' network
#'
#' @param network matrix or data.frame, strength of gene co-expression (edge
#' values).
#' @param seq_k vector, sequence of k number of cluster to test
#' @param fit_plot boolean, does the plot with silhouette coefficient depending
#' on the k tested should be plotted.
#' @param ... any other parameter compatible with the
#' \code{\link[cluster]{pam}} function.
#'
#' @importFrom magrittr %>%
#' @importFrom cluster pam
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr filter select
#' @importFrom graphics plot
#'
#' @return data.frame, a two cols table with the gene id in the first one, and
#' the cluster number assignation in the second one.
#'
#' @examples
#' df <- kuehne_expr[1:24, 1:350]
#' net <- build_net(df, n_threads = 1)
#' mods <- detect_modules(df, net$network)
#' net_mod_1 <- net$network[mods$modules$`1`, mods$modules$`1`]
#' get_sub_clusters(net_mod_1)
#'
#' @export

get_sub_clusters <- function(network, seq_k = 1:15, fit_plot = TRUE,
                             ...) {
  # Checking args
  is_network(network)
  if(!is.numeric(seq_k))
    stop("seq_k must be a numeric vector")
  if(length(seq_k) < 2)
    stop("seq_k must at leat have two power values to test")
  if (any(lapply(seq_k, function(x) x < 1 | x %% 1 != 0) %>% unlist))
    stop("seq_k must contain only whole numbers >= 1")
  if(!is.logical(fit_plot))
    stop("fit_plot must be a boolean")

  # Computing the k-medoid
  list_k_tests <- lapply(seq_k, function(k) {
    k_res <- cluster::pam(network, k, diss = TRUE, ...)
    return(k_res)
  }) %>% stats::setNames(paste0("k_", seq_k))

  # Summarizing needed results for plot and returning optimal k
  df_k_tests <- list_k_tests %>%
    lapply(`[[`, "silinfo") %>%
    lapply(`[[`, "avg.width") %>%
    unlist %>%
    data.frame(k = seq_k[-1], avg_sil_width = .)

  # If asked, plotting the silhouette coefficient against k tested
  if (fit_plot) {
    graphics::plot(df_k_tests, type="b", pch = 19,
         frame = FALSE, xlab="Number of clusters K",
         ylab="Average silhouette width")
  }

  # Isolating optimal k regarding the silhouette coefficient
  opti_k <- df_k_tests %>%
    dplyr::filter(avg_sil_width == max(avg_sil_width)) %>%
    dplyr::select(k) %>% as.numeric()

  # Formatting the table to return the gene id association to sub_module
  clusters <- list_k_tests[[opti_k]]$clustering %>%
    data.frame(sub_module = .) %>%
    tibble::rownames_to_column("gene") %>%
    dplyr::mutate(sub_module = as.character(sub_module))

  return(clusters)
}
