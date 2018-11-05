##' Read a phylogeny in newick format
##'
##' Builds an R object representation of a phylogeny from its newick description
##'
##' @param file
##' The name of the file containing the newick description
read.newick = function(file) {
    stopifnot(file.exists(file))
    newick = scan(file, what=character(), quiet=TRUE)
    tree = .Call(rcl_read_newick, newick)
    class(tree) = "tree"
    return (tree)
}


is.tree = function(x) {
    return (inherits(x, "tree"))
}


print.tree = function(x, ...) {
    str(x)
}


write.newick = function(phy, file="") {
    stopifnot(is.tree(phy))
    cat(.Call(rcl_write_newick, phy), "\n", sep="", file=file)
}


root = function(phy) {
    attr(phy, "root")
}


Ntip = function(phy) {
    attr(phy, "Ntip")
}


Nnode = function(phy) {
    attr(phy, "Nnode")
}


tiplabels = function(phy) {
    if (is.null(tip.label <- attr(phy, "tip.label"))) {
        stopifnot(is.tree(phy))
        tip.label = .Call(rcl_tiplabels, phy)
        attr(phy, "tip.label") = tip.label
    }
    return (tip.label)
}


brlens = function(phy) {
    if (is.null(brlen <- attr(phy, "brlen"))) {
        stopifnot(is.tree(phy))
        brlen = .Call(rcl_node_brlens, phy)
        attr(phy, "brlen") = brlen
    }
    return (brlen)
}


ages = function(phy) {
    if (is.null(age <- attr(phy, "age"))) {
        stopifnot(is.tree(phy))
        age = .Call(rcl_node_ages, phy)
        attr(phy, "age") = age
    }
    return (age)
}


ancestors = function(phy) {
    if (is.null(anc <- attr(phy, "ancestor"))) {
        stopifnot(is.tree(phy))
        anc = lapply(1L:Nnode(phy), function(n) .Call(rcl_ancestors, phy, n)[-1L])
        attr(phy, "ancestor") = anc
    }
    return (anc)
}


children = function(phy) {
    if (is.null(kids <- attr(phy, "children"))) {
        stopifnot(is.tree(phy))
        kids = lapply(1L:Nnode(phy), function(n) .Call(rcl_children, phy, n))
        attr(phy, "children") = kids
    }
    return (kids)
}




