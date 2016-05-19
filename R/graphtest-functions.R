

#' Performs graph-based permutation tests
#'
#' Performs graph-based tests for one-way designs.
#'
#' @param physeq A phyloseq object.
#' @param sampletype A string giving the column name of the sample to
#' be tested. This should be a factor with two or more levels.
#' @param distance A distance, see \link{phyloseq::distance} for a
#' list of the possible methods.
#' @param type One of "mst", "knn", "threshold". If "mst", forms the
#' minimum spanning tree of the sample points. If "knn", forms a
#' directed graph with links from each node to its k nearest
#' neighbors. If "threshold", forms a graph with edges between every
#' pair of samples within a certain distance.
#' @param max.dist For type "threshold", the maximum distance between
#' two samples such that we put an edge between them.
#' @param knn For type "knn", the number of nearest neighbors.
#' @param nperm The number of permutations to perform.
#'
#'
#' @importFrom igraph graph.adjacency minimum.spanning.tree get.edgelist
#' @import phyloseq
#' @import ggplot2
#' 
#' @return A list with the observed number of pure edges, the vector
#' containing the number of pure edges in each permutation, the
#' permutation p-value, the graph used for testing, and a vector with
#' the sample types used for the test.
#' @export
graph_perm_test = function(physeq, sampletype,
    distance = "jaccard", type = c("mst", "knn", "threshold.value", "threshold.nedges"),
    max.dist = .4, knn = 1, nedges = nsamples(physeq), nperm = 99) {
    type = match.arg(type)
    # make the network
    d = distance(physeq, method = distance, type = "samples")
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
               neighbors = t(apply(d,1, function(x) {
                   r = rank(x)
                   nvec = ((r > 1) & (r < (knn + 2))) + 0
               }))
               net = graph.adjacency(neighbors, mode = "directed", add.colnames = "name")
           },
           "mst" = {
               gr = graph.adjacency(as.matrix(d), mode = "undirected", weighted = TRUE,
                   add.colnames = "name")
               net = minimum.spanning.tree(gr, algorithm = "prim")
           }           
    )
    el = get.edgelist(net)
    sampledata = data.frame(sample_data(physeq))
    # find the number of pure edges for the non-permuted data
    observedPureEdges = apply(el, 1, function(x)
        sampledata[x[1], sampletype] == sampledata[x[2], sampletype])
    nobserved = sum(observedPureEdges)
    origSampleData = sampledata[,sampletype]
    names(origSampleData) = rownames(sampledata)
    # find the permutation distribution of the number of pure edges
    permvec = numeric(nperm)
    for(i in 1:nperm) {
        sampledata[,sampletype] = sampledata[sample(nrow(sampledata)), sampletype]
        permPureEdges = apply(el, 1, function(x)
            sampledata[x[1], sampletype] == sampledata[x[2], sampletype])
        permvec[i] = sum(permPureEdges)
    }
    pval = (sum(permvec >= nobserved) + 1) / (nperm + 1)
    return(list(observed = nobserved, perm = permvec, pval = pval,
                net = net, sampletype = origSampleData))
}


#' Plots the graph used for testing
#'
#' When using the graph_perm_test function, a graph is created. This
#' function will plot the graph used for testing with nodes colored by
#' sample type and edges marked as pure or mixed.
#'
#' @param graphtest The output from graph_perm_test.
#' @return A ggplot object.
#'
#' @importFrom igraph V layout.fruchterman.reingold get.edgelist
#' @export
plot_test_network = function(graphtest) {
    net = graphtest$net
    sampletype = graphtest$sampletype
    nodes = layout.fruchterman.reingold(net)
    rownames(nodes) = V(net)$name
    nodes = data.frame(nodes)
    names(nodes) = c("x", "y")
    nodes$sampletype = sampletype
    el = get.edgelist(net)
    edgeDF = t(apply(el, 1, function(x) {
        unlist(c(nodes[x[1],1:2], nodes[x[2],1:2]))
    }))

    pure = apply(el, 1, function(x) {
        if(sampletype[x[1]] == sampletype[x[2]])
            return("pure")
        return("mixed")
    })
    edgeDF = data.frame(edgeDF, edgetype = pure)
    names(edgeDF)[1:4] = c("x", "y", "xend", "yend")
    p = ggplot(edgeDF) +
        geom_segment(aes(x = x, y = y, xend = xend, yend = yend, linetype = edgetype)) + 
        geom_point(aes(x = x, y = y, color = sampletype), data = nodes) +
        scale_linetype_manual(values = c(3,1))
    p = p + theme(axis.text.y = element_blank(), axis.text.x = element_blank(),
        axis.ticks = element_blank(), axis.title.x = element_blank(),
        axis.title.y = element_blank())
    print(p)
}


#' Plots the permutation distribution
#'
#' Plots a histogram of the permutation distribution of the number of
#' pure edges and a mark showing the observed number of pure edges. 
#'
#' @param graphtest The output from graph_perm_test.
#' @return A ggplot object.
#'
#' @export
plot_permutations = function(graphtest) {
    p = qplot(graphtest$perm, geom = "histogram")
    ymax = ggplot_build(p)$panel$ranges[[1]]$y.range[2]
    p = p + geom_segment(aes(x = graphtest$observed, y = 0,
                         xend = graphtest$observed, yend = ymax / 10), color = "red") +
        geom_point(aes(x = graphtest$observed, y = ymax / 10), color = "red") +
        xlab("Number of pure edges")                         
    print(p)
}
