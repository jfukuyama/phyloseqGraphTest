

#' Performs graph-based permutation tests
#'
#' Performs graph-based tests for one-way designs.
#'
#' @param physeq A phyloseq object.
#' @param sampletype A string giving the column name of the sample to
#' be tested. This should be a factor with two or more levels.
#' @param grouping Either a string with the name of a sample data
#' column or a factor of length equal to the number of samples in
#' physeq. These are the groups of samples whose labels should be
#' permuted and are used for repeated measures designs. Default is no
#' grouping (each group is of size 1).
#' @param distance A distance, see \code{\link[phyloseq]{distance}} for a
#' list of the possible methods.
#' @param type One of "mst", "knn", "threshold". If "mst", forms the
#' minimum spanning tree of the sample points. If "knn", forms a
#' directed graph with links from each node to its k nearest
#' neighbors. If "threshold", forms a graph with edges between every
#' pair of samples within a certain distance.
#' @param max.dist For type "threshold", the maximum distance between
#' two samples such that we put an edge between them.
#' @param knn For type "knn", the number of nearest neighbors.
#' @param keep.isolates In the returned network, keep the unconnected
#' points?
#' @param nperm The number of permutations to perform.
#' @param nedges If using "threshold.nedges", the number of edges to use.
#'
#'
#' @importFrom igraph graph.adjacency minimum.spanning.tree
#' get.edgelist V<- E<- V E induced_subgraph
#' @import phyloseq
#' @import ggplot2
#' @import ggnetwork
#' 
#' @return A list with the observed number of pure edges, the vector
#' containing the number of pure edges in each permutation, the
#' permutation p-value, the graph used for testing, and a vector with
#' the sample types used for the test.
#' @examples
#' library(phyloseq)
#' data(enterotype)
#' gt = graph_perm_test(enterotype, sampletype = "SeqTech", type = "mst")
#' gt
#' @export
graph_perm_test = function(physeq, sampletype, grouping = 1:nsamples(physeq),
    distance = "jaccard", type = c("mst", "knn", "threshold.value", "threshold.nedges"),
    max.dist = .4, knn = 1, nedges = nsamples(physeq), keep.isolates = TRUE, nperm = 499) {
    type = match.arg(type)
    # make the network
    d = distance(physeq, method = distance, type = "samples")
    if(!validGrouping(sample_data(physeq), sampletype, grouping)) {
        stop("Not a valid grouping, all values of sampletype must
              be the same within each level of grouping")
    }
    switch(type,
           "threshold.value" = {
               neighbors = as.matrix(d) <= max.dist
               diag(neighbors) = 0
               net = graph.adjacency(neighbors, mode = "undirected", add.colnames = "name")   
           },
           "threshold.nedges" = {
               threshold = sort(as.vector(d))[nedges]
               neighbors = as.matrix(d) <= threshold
               diag(neighbors) = 0
               net = graph.adjacency(neighbors, mode = "undirected", add.colnames = "name")
           },
           "knn" = {
               neighbors = t(apply(as.matrix(d),1, function(x) {
                   r = rank(x)
                   nvec = ((r > 1) & (r < (knn + 2))) + 0
               }))
               neighbors = neighbors + t(neighbors)
               net = graph.adjacency(neighbors, mode = "undirected",
                   add.colnames = "name", weighted = TRUE)
           },
           "mst" = {
               gr = graph.adjacency(as.matrix(d), mode = "undirected", weighted = TRUE,
                   add.colnames = "name")
               net = minimum.spanning.tree(gr, algorithm = "prim")
           }           
           )
    el = get.edgelist(net)
    sampledata = data.frame(sample_data(physeq))
    elTypes = el
    elTypes[,1] = sampledata[el[,1], sampletype]
    elTypes[,2] = sampledata[el[,2], sampletype]
    observedPureEdges = apply(elTypes, 1, function(x) x[1] == x[2])
    edgeType = sapply(observedPureEdges, function(x) if(x) "pure" else "mixed")
    # set these attributes for plotting later
    if(is.factor(sampledata[,sampletype]))
        V(net)$sampletype = as.character(sampledata[,sampletype])
    else
        V(net)$sampletype = sampledata[,sampletype]
    E(net)$edgetype = edgeType
    
    # find the number of pure edges for the non-permuted data
    nobserved = sum(observedPureEdges)
    origSampleData = sampledata[,sampletype]
    names(origSampleData) = rownames(sampledata)
    # find the permutation distribution of the number of pure edges
    permvec = numeric(nperm)
    for(i in 1:nperm) {
        sampledata[,sampletype] = permute(sampledata, grouping, sampletype)
        elTypes = el
        elTypes[,1] = sampledata[el[,1], sampletype]
        elTypes[,2] = sampledata[el[,2], sampletype]
        permPureEdges = apply(elTypes, 1, function(x) x[1] == x[2])
        permvec[i] = sum(permPureEdges)
    }
    pval = (sum(permvec >= nobserved) + 1) / (nperm + 1)
    if(!keep.isolates) {
        degrees = igraph::degree(net)
        net = igraph::induced_subgraph(net, which(degrees > 0))
    }
    out = list(observed = nobserved, perm = permvec, pval = pval,
               net = net, sampletype = origSampleData, type = type)
    class(out) = "psgraphtest"
    return(out)
}

#' Print psgraphtest objects
#' @param x \code{psgraphtest} object.
#' @param ... Not used
#' @method print psgraphtest
#' @export
print.psgraphtest <- function(x, ...) {
    cat("Output from graph_perm_test\n")
    cat("---------------------------\n")
    cat(paste("Observed test statistic: ", x$observed, " pure edges", "\n", sep = ""))
    cat(paste(nrow(get.edgelist(x$net)), " total edges in the graph", "\n", sep = ""))
    cat(paste("Permutation p-value: ", x$pval, "\n", sep = ""))
}

#' Permute labels
#'
#' Permutes sample labels, respecting repeated measures.
#'
#' @param sampledata Data frame describing the samples.
#' @param grouping Grouping for repeated measures.
#' @param sampletype The sampletype used for testing (a column of sampledata).
#' @return A permuted set of labels where the permutations are done
#'     over the levels of grouping.
#' @keywords internal
permute = function(sampledata, grouping, sampletype) {
    if(length(grouping) != nrow(sampledata)) {
        grouping = sampledata[,grouping]
    }
    x = as.character(sampledata[,sampletype])
    # gives the original mapping between grouping variables and sampletype
    labels = tapply(x, grouping, function(x) x[1])
    # permute the labels of the groupings
    names(labels) = sample(names(labels))
    return(labels[as.character(grouping)])
}

#' Check for valid grouping
#'
#' Grouping should describe a repeated measures design, so this
#' function tests whether all of the levels of grouping have the same
#' value of sampletype.
#'
#' @param sd Data frame describing the samples.
#' @param sampletype The sampletype used for testing.
#' @param grouping Grouping for repeated measures.
#' @return TRUE or FALSE for valid or invalid grouping.
#' @keywords internal
validGrouping = function(sd, sampletype, grouping) {
    if(!(sampletype %in% colnames(sd))) {
        stop("\'sampletype\' must be a column names of the sample data")
    }
    if(!(all(grouping %in% colnames(sd))) && (length(grouping) != nrow(sd))) {
        stop("\'grouping\' must be either a column name of the sample data
             or a vector with number of elements equal to the number of samples")
    }
    sd = data.frame(sd)
    if(length(grouping) != nrow(sd)) {
        grouping = sd[,grouping]
    }
    valid = all(tapply(sd[,sampletype], grouping, FUN = function(x)
        length(unique(x)) == 1))
    return(valid)
}

#' Plots the graph used for testing
#'
#' When using the graph_perm_test function, a graph is created. This
#' function will plot the graph used for testing with nodes colored by
#' sample type and edges marked as pure or mixed.
#'
#' @param graphtest The output from graph_perm_test.
#' @return A ggplot object created by ggnetwork.
#' @examples
#' library(phyloseq)
#' data(enterotype)
#' gt = graph_perm_test(enterotype, sampletype = "SeqTech")
#' plot_test_network(gt)
#' @export
plot_test_network = function(graphtest) {
    if(graphtest$type == "mst")
        layout = igraph::layout_(graphtest$net, igraph::with_kk())
    else
        layout = igraph::layout_(graphtest$net, igraph::with_fr())
    ggplot(new_fortify.igraph(graphtest$net),
      aes_string(x = "x", y = "y", xend = "xend", yend = "yend"), layout = layout) +
      geom_edges(aes_string(linetype = "edgetype")) +
      geom_nodes(aes_string(color = "sampletype")) +
      scale_linetype_manual(values = c(3,1)) + theme_blank()
}


#' Plots the permutation distribution
#'
#' Plots a histogram of the permutation distribution of the number of
#' pure edges and a mark showing the observed number of pure edges. 
#'
#' @param graphtest The output from graph_perm_test.
#' @param bins The number of bins to use for the histogram.
#' @importFrom utils packageVersion
#' @return A ggplot object.
#' @examples
#' library(phyloseq)
#' data(enterotype)
#' gt = graph_perm_test(enterotype, sampletype = "SeqTech")
#' plot_permutations(gt)
#' @export
plot_permutations = function(graphtest, bins = 30) {
    p = qplot(graphtest$perm, geom = "histogram", bins = bins)
    if(packageVersion("ggplot2") >= "2.2.1.9000") {
        ymax = max(ggplot_build(p)$layout$panel_scales_y[[1]]$get_limits())
    } else {
        ymax = ggplot_build(p)$layout$panel_ranges[[1]][["y.range"]][2]
    }
    p + geom_segment(aes(x = graphtest$observed, y = 0,
                         xend = graphtest$observed, yend = ymax / 10), color = "red") +
        geom_point(aes(x = graphtest$observed, y = ymax / 10), color = "red") +
        xlab("Number of pure edges")                         
}


#' Fortify method for networks of class \code{\link[igraph:igraph-package]{igraph}}
#'
#' This is copied with very slight modification from
#' https://github.com/briatte/ggnetwork/blob/master/R/fortify-igraph.R,
#' as that version is not on CRAN yet.
#'
#' @param model an object of class \code{\link[igraph:igraph-package]{igraph}}.
#' @param data not used by this method.
#' @param layout a function call to an
#'   \code{\link[igraph:igraph-package]{igraph}} layout function, such as
#'   \code{\link[igraph]{layout_nicely}} (the default), or a 2 column matrix
#'   giving the x and y coordinates for the vertices.
#'   See \code{\link[igraph]{layout_}} for details.
#' @inheritParams format_fortify
#' @param ... additional parameters for the \code{\link[igraph]{layout_}} function
#'
#' @return a \code{\link[base]{data.frame}} object.
new_fortify.igraph <- function(
  model,
  data = NULL,
  layout = igraph::nicely(),
  arrow.gap = ifelse(igraph::is.directed(model), 0.025, 0),
  by = NULL,
  scale = TRUE,
  stringsAsFactors = getOption("stringsAsFactors", FALSE),
  ...
) {
  # node placement
  if (inherits(layout, "matrix")) {
    if (nrow(layout) != igraph::gorder(model)) {
      stop("layout matrix dimensions do not match network size")
    }
    nodes <- layout[, 1:2]
  } else {
    nodes <- igraph::layout_(model, layout, ...)
  }

  format_fortify(
    model = model,
    nodes = nodes,
    weights = "none",
    arrow.gap = arrow.gap,
    by = by,
    scale = scale,
    stringsAsFactors = stringsAsFactors,
    .list_vertex_attributes_fun = igraph::list.vertex.attributes,
    .get_vertex_attributes_fun = igraph::get.vertex.attribute,
    .list_edges_attributes_fun = igraph::list.edge.attributes,
    .get_edges_attributes_fun = igraph::get.edge.attribute,
    .as_edges_list_fun = igraph::as_edgelist
  )
}



#' format_fortify
#'
#' a unified function to format \code{\link[network]{network}} or
#' \code{\link[igraph:igraph-package]{igraph}} object. Copied with
#' very slight modification from
#' https://github.com/briatte/ggnetwork/blob/master/R/utilities.R to
#' fix the same CRAN problem as new_fortify.igraph.
#'
#' @param model an object of class \code{\link[network]{network}}
#'   or \code{\link[igraph:igraph-package]{igraph}}.
#' @param nodes a nodes object from a call to fortify.
#' @param weights the name of an edge attribute to use as edge weights when
#'   computing the network layout, if the layout supports such weights (see
#'   'Details').
#'   Defaults to \code{NULL} (no edge weights).
#' @param arrow.gap a parameter that will shorten the network edges in order to
#'   avoid overplotting edge arrows and nodes; defaults to \code{0} when the
#'   network is undirected (no edge shortening), or to \code{0.025} when the
#'   network is directed. Small values near \code{0.025} will generally achieve
#'   good results when the size of the nodes is reasonably small.
#' @param by a character vector that matches an edge attribute, which will be
#' @param by a character vector that matches an edge attribute, which will be
#'   used to generate a data frame that can be plotted with
#'   \code{\link[ggplot2]{facet_wrap}} or \code{\link[ggplot2]{facet_grid}}. The
#'   nodes of the network will appear in all facets, at the same coordinates.
#'   Defaults to \code{NULL} (no faceting).
#' @param scale whether to (re)scale the layout coordinates. Defaults to
#'   \code{TRUE}, but should be set to \code{FALSE} if \code{layout} contains
#'   meaningful spatial coordinates, such as latitude and longitude.
#' @param stringsAsFactors whether vertex and edge attributes should be
#'   converted to factors if they are of class \code{character}. Defaults to
#'   the value of \code{getOption("stringsAsFactors")}, which is \code{FALSE}
#'   by default: see \code{\link[base]{data.frame}}.
#' @param .list_vertex_attributes_fun a "list vertex attributes" function.
#' @param .get_vertex_attributes_fun a "get vertex attributes" function.
#' @param .list_edges_attributes_fun a "get edges attributes" function.
#' @param .get_edges_attributes_fun a "get edges attributes" function.
#' @param .as_edges_list_fun a "as edges list" function.
#'
#' @return a \code{\link[base]{data.frame}} object.
#'
#' @keywords internal
#'
format_fortify <- function(
  model,
  nodes = NULL,
  weights = NULL,
  arrow.gap = 0,
  by = NULL,
  scale = TRUE,
  stringsAsFactors = getOption("stringsAsFactors", FALSE),
  .list_vertex_attributes_fun = NULL,
  .get_vertex_attributes_fun = NULL,
  .list_edges_attributes_fun = NULL,
  .get_edges_attributes_fun = NULL,
  .as_edges_list_fun = NULL
) {
  # store coordinates
  nodes <- data.frame(nodes)
  colnames(nodes) <- c("x", "y")

  # rescale coordinates
  if (scale) {
    nodes$x <- scale_safely(nodes$x)
    nodes$y <- scale_safely(nodes$y)
  }

  # import vertex attributes
  if (length(.list_vertex_attributes_fun(model)) > 0) {
    nodes <- cbind.data.frame(
      nodes,
      sapply(
        X = .list_vertex_attributes_fun(model),
        Y = model,
        FUN = function(X, Y) .get_vertex_attributes_fun(Y, X),
        simplify = FALSE
      ),
      stringsAsFactors = stringsAsFactors
    )
  }

  # edge list
  if (inherits(model, "igraph")) {
    edges <- .as_edges_list_fun(model, names = FALSE)
  } else {
    edges <- .as_edges_list_fun(model, attrname = weights)
  }

  # edge list (if there are duplicated rows)
  if (nrow(edges[, 1:2, drop = FALSE]) > nrow(unique(edges[, 1:2, drop = FALSE]))) {
    warning("duplicated edges detected")
  }

  edges <- data.frame(nodes[edges[, 1], c("x", "y")], nodes[edges[, 2], c("x", "y")])
  colnames(edges) <- c("x", "y", "xend", "yend")

  # arrow gap (thanks to @heike and @ethen8181 for their work on this issue)
  if (arrow.gap > 0) {
    x.length <- edges$xend - edges$x
    y.length <- edges$yend - edges$y
    arrow.gap <- arrow.gap / sqrt(x.length^2 + y.length^2)
    edges$xend <- edges$x + (1 - arrow.gap) * x.length
    edges$yend <- edges$y + (1 - arrow.gap) * y.length
  }

  # import edge attributes
  if (length(.list_edges_attributes_fun(model)) > 0) {
    edges <- cbind.data.frame(
      edges,
      sapply(
        X = .list_edges_attributes_fun(model),
        Y = model,
        FUN = function(X, Y) .get_edges_attributes_fun(Y, X),
        simplify = FALSE
      ),
      stringsAsFactors = stringsAsFactors
    )
  }

  if (nrow(edges) > 0) {
    # drop "na" columns created by 'network' methods
    # this is to ensure consistency with 'igraph' methods
    if ("na" %in% colnames(nodes)) nodes$na <- NULL
    if ("na" %in% colnames(edges)) edges$na <- NULL

    # merge edges and nodes data
    edges <- merge(nodes, edges, by = c("x", "y"), all = TRUE)

    # add missing columns to nodes data
    nodes$xend <- nodes$x
    nodes$yend <- nodes$y
    # names(nodes) <- names(edges)[1:ncol(nodes)] # columns are already named from 'nodes' and 'edges'

    # make nodes data of identical dimensions to edges data
    nodes[, setdiff(names(edges), names(nodes))] <- NA

    # panelize nodes (for temporal networks)
    if (!is.null(by)) {
      nodes <- lapply(sort(unique(edges[, by])), function(x) {
        y <- nodes
        y[, by] <- x
        y
      })
      nodes <- do.call(rbind, nodes)
    }

    # return a data frame with network.size(model) + network.edgecount(model) rows,
    # or length(unique(edges[, by ])) * network.size(model) + network.edgecount(model)
    # rows if the nodes have been panelized

    # [NOTE] `edges` has to be passed first to `rbind` in order for the edge
    # attributes (e.g. factor) to be preserved in the result; this should not
    # affect plotting the result, but differs from `ggnetwork` 0.5.1: `nodes`
    # is now at the end of the result rather at the beginning

    return(unique(rbind(edges[!is.na(edges$xend), ], nodes)))
  } else {
    # add missing columns to nodes data
    nodes$xend <- nodes$x
    nodes$yend <- nodes$y
    return(nodes)
  }
}

#' Rescale x to (0, 1), except if x is constant
#'
#' Copied from https://github.com/briatte/ggnetwork/blob/f3b8b84d28a65620a94f7aecd769c0ea939466e3/R/utilities.R so as to fix a problem with the cran version of ggnetwork.

#' @param x a vector to rescale
#' @param scale the scale on which to rescale the vector
#'
#' @return The rescaled vector, coerced to a vector if necessary.
#'   If the original vector was constant, all of its values are replaced by 0.5.
#'
#' @author Kipp Johnson
scale_safely <- function(x, scale = diff(range(x))) {
  if (!scale) {
    x <- rep(0.5, length.out = length(x))
  } else {
    x <- scale(x, center = min(x), scale = scale)
  }
  as.vector(x)
}
