tree.isbinary = function(phy) {
    if (Ntip(phy) == ((Nnode(phy) + 1) / 2))
        return (TRUE)
    else
        return (FALSE)
}


drop.tip = function(phy, tip) {
    stopifnot(is.tree(phy))
    stopifnot(tree.isbinary(phy))
    stopifnot(class(tip) == "character")
    nodebits = integer(Nnode(phy))
    nodebits[1L:Ntip(phy)][!(tiplabels(phy) %in% tip)] = 1L
    edge = t(.Call(rcl_subgraph, phy, nodebits))-1L # subtract 1 to match C indexing
    singletons = as.integer(names(which(table(edge[, 1L]) == 1L)))  # identify singletons
    remap = integer(length(singletons))
    for (i in 1L:length(singletons)) {
        s = singletons[i]
        ix = match(s, edge[, 1L])
        while (edge[ix, 2L] %in% singletons) {
            s = edge[ix, 2L]
            ix = match(s, edge[, 1L])
        }
        remap[i] = edge[ix, 2L]
    }
    edge = edge[-match(singletons, edge[, 1L]), ]
    edge[match(singletons, edge[, 2L], nomatch=0), 2L] = remap[singletons %in% edge[, 2L]]
    if (any(table(edge[, 1L]) == 1L))
        stop("You've discovered a bug in rcl::drop.tip!")
    brlen = .Call(rcl_subgraph_brlen, phy, edge)
    alist = split(edge, edge[, 1L])
    children = matrix(-1L, Nnode(phy), 2L)
    sapply(alist, function(z) {
        children[z[1L]+1L, ] <<- z[3L:4L]
    })
    newroot = min(edge[, 1L])
    newick = .Call(rcl_subgraph_newick, phy, newroot, children, brlen)
    newphy = .Call(rcl_read_newick, newick)
    class(newphy) = "tree"
    return (newphy)
}


extract.clade = function(phy, node) {
    stopifnot(is.tree(phy))
    storage.mode(node) = "integer"
    stopifnot(node > Ntip(phy))
    stopifnot(node < Nnode(phy))
    clade = .Call(rcl_subtree, phy, node)
    class(clade) = "tree"
    return (clade)
}


mrca = function(phy, node1, node2) {
    stopifnot(is.tree(phy))
    storage.mode(node1) = "integer"
    storage.mode(node2) = "integer"
    stopifnot(all(node1 >= 1))
    stopifnot(all(node1 <= Nnode(phy)))
    stopifnot(all(node2 >= 1))
    stopifnot(all(node2 <= Nnode(phy)))
    .Call(rcl_mrca, phy, node1, node2)
}


tree.traverse = function(node, phy, visit=c("ALL_NODES", "INTERNAL_NODES_ONLY"), order=c("PREORDER", "POSTORDER")) {
    stopifnot(is.tree(phy))
    storage.mode(node) = "integer"
    if (node <= 0 || node > Nnode(phy))
        stop("Invalid node index")


    order = match.arg(order)
    visit = match.arg(visit)

    order = switch(order, PREORDER=0L, POSTORDER=1L)
    visit = switch(visit, ALL_NODES=0L, INTERNAL_NODES_ONLY=1L)

    travers = .Call(rcl_tree_traverse, phy, node, order, visit)

    tree.step = function(skip=NULL) {
        if (!is.null(skip)) {
            if (skip <= 0 || skip > Nnode(phy)) {
                .Call(rcl_tree_reset, travers)
                stop("Invalid node index")
            }
            node = try(.Call(rcl_tree_jump, skip, travers))
            if (inherits(node, "try-error")) {
                .Call(rcl_tree_reset, travers)
                return (NULL)
            }
        } else
            node = .Call(rcl_tree_step, travers)
        if (!node) {
            .Call(rcl_tree_reset, travers)
            return (NULL)
        }
        return (node)
    }

    return (tree.step)
}
