bhc = function(data, alpha=1, beta=1) {
    storage.mode(data) = "integer"
    stopifnot(inherits(data, "matrix"))
    res = .Call(rcl_bhc, as.numeric(alpha), as.numeric(beta), data)
    res$Ntip = ncol(data)
    res
}

list.clusters = function(x) {

    tips = 1L:x$Ntip
    tipsleft = x$Ntip

    parent = x$parent
    cutoff = x$logodds

    res = list()
    g = 1L

    while (tipsleft > 0) {
        node = tips[1L]

        while (node != x$root && cutoff[parent[node]] > log(0.5))
            node = parent[node]

        begin = x$seqpos[node]
        end = x$seqpos[x$lastvisit[node]]
        desc = x$seq[begin:end]
        desc = desc[desc <= x$Ntip]
        tipsleft = tipsleft - length(desc)
        tips = tips[-match(desc, tips)]
        res[[g]] = desc
        g = g + 1L
    }

    res
}
