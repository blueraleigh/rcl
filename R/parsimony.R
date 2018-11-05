.check.states = function(phy, data, NSTATES) {
    storage.mode(NSTATES) = "integer"
    counts = rep(1L, Ntip(phy))
    fdata = factor(data, levels=0L:(NSTATES-1L))
    downpass = xtabs(counts ~ fdata + names(fdata))
    if (sum(downpass) != Ntip(phy))
        stop("Some taxa are polymorphic.")
    if (any(rowSums(downpass) == 0L))
        warning("Not all states are present among the terminal nodes.")
    attr(downpass, "call") = NULL
    names(attr(downpass, "dimnames")) = NULL
    rownames(downpass) = NULL
    downpass = downpass[, tiplabels(phy)]
    downpass = cbind(downpass, matrix(0L, nrow=NSTATES, ncol=Nnode(phy)-Ntip(phy)))
    return (downpass)
}

# Calculate the Fitch parsimony score given an
# r-state character.
pscore.fitch = function(phy, data, NSTATES, intermediates=FALSE) {
    stopifnot(is.tree(phy))
    stopifnot(tree.isbinary(phy))
    storage.mode(data) = "integer"
    stopifnot(!is.null(names(data)))
    downpass = .check.states(phy, data, NSTATES)
    pscore = integer(Nnode(phy))
    .Call(rcl_fitch_pscore, phy, downpass, pscore)
    if (intermediates)
        return (pscore[-(1L:Ntip(phy))])
    return (pscore[root(phy)])
}


# Perform a maximum parsimony reconstruction of ancestral
# states given an r-state character using Fitch's algorithm
mpr.fitch = function(phy, data, NSTATES) {
    stopifnot(is.tree(phy))
    stopifnot(tree.isbinary(phy))
    storage.mode(data) = "integer"
    stopifnot(!is.null(names(data)))
    downpass = .check.states(phy, data, NSTATES)
    uppass = matrix(0L, nrow=NSTATES, ncol=Nnode(phy))
    pscore = integer(Nnode(phy))
    .Call(rcl_fitch_mpr, phy, uppass, downpass, pscore)
    pscore = pscore[root(phy)]
    # Sample a history of character evolution from a
    # maximum parsimony reconstruction of ancestral states
    sample.mpr = function(n) {
        storage.mode(n) = "integer"
        node.state = matrix(0L, nrow=Nnode(phy), ncol=n)
        node.state[1L:Ntip(phy), ] = data[tiplabels(phy)]
        .Call(rcl_fitch_history, phy, node.state, uppass, downpass, n)
        return (node.state)
    }
    return (sample.mpr)
}


# Count the number of MPR reconstructions
count.mpr.fitch = function(phy, data, NSTATES) {
    stopifnot(is.tree(phy))
    stopifnot(tree.isbinary(phy))
    storage.mode(data) = "integer"
    stopifnot(!is.null(names(data)))
    downpass = .check.states(phy, data, NSTATES)
    uppass = matrix(0L, nrow=NSTATES, ncol=Nnode(phy))
    pscore = integer(Nnode(phy))
    .Call(rcl_fitch_mpr, phy, uppass, downpass, pscore)
    .Call(rcl_fitch_count, phy, uppass, downpass)
}


# Calculate the Sankoff parsimony score given an
# r-state character.
pscore.sankoff = function(phy, data, cost, NSTATES, intermediates=FALSE) {
    stopifnot(is.tree(phy))
    stopifnot(tree.isbinary(phy))
    storage.mode(data) = "integer"
    stopifnot(!is.null(names(data)))
    MAXCOST = 10000L
    stopifnot(NROW(cost) == NCOL(cost) && NROW(cost) == NSTATES)
    cost[cost > MAXCOST] = MAXCOST
    diag(cost) = 0L
    storage.mode(cost) = "integer"
    downpass.cost = .check.states(phy, data, NSTATES)
    downpass.cost[downpass.cost == 0L] = MAXCOST
    downpass.cost[downpass.cost == 1L] = 0L
    pscore = integer(Nnode(phy))
    .Call(rcl_sankoff_pscore, phy, cost, downpass.cost, pscore)
    if (intermediates)
        return (pscore[-(1L:Ntip(phy))])
    return (pscore[root(phy)])
}



# Perform a maximum parsimony reconstruction of ancestral
# states given an r-state character using Sankoff's algorithm
mpr.sankoff = function(phy, data, cost, NSTATES) {
    stopifnot(is.tree(phy))
    stopifnot(tree.isbinary(phy))
    storage.mode(data) = "integer"
    stopifnot(!is.null(names(data)))
    MAXCOST = 10000L
    stopifnot(NROW(cost) == NCOL(cost) && NROW(cost) == NSTATES)
    cost[cost > MAXCOST] = MAXCOST
    diag(cost) = 0L
    storage.mode(cost) = "integer"
    downpass.cost = .check.states(phy, data, NSTATES)
    downpass.cost[downpass.cost == 0L] = MAXCOST
    downpass.cost[downpass.cost == 1L] = 0L
    downpass = matrix(0L, nrow=NSTATES, ncol=Nnode(phy))
    uppass = matrix(0L, nrow=NSTATES, ncol=Nnode(phy))
    pscore = integer(Nnode(phy))
    .Call(rcl_sankoff_mpr, phy, cost, uppass, downpass, downpass.cost, pscore)
    pscore = pscore[root(phy)]
    sample.mpr = function(n) {
        storage.mode(n) = "integer"
        node.state = matrix(0L, nrow=Nnode(phy), ncol=n)
        node.state[1L:Ntip(phy), ] = data[tiplabels(phy)]
        .Call(rcl_sankoff_history, phy, cost, node.state, uppass, downpass, downpass.cost, n)
        return (node.state)
    }
    return (sample.mpr)
}


# Count the number of MPR reconstructions
count.mpr.sankoff = function(phy, data, cost, NSTATES) {
    stopifnot(is.tree(phy))
    stopifnot(tree.isbinary(phy))
    storage.mode(data) = "integer"
    stopifnot(!is.null(names(data)))
    MAXCOST = 10000L
    stopifnot(NROW(cost) == NCOL(cost) && NROW(cost) == NSTATES)
    cost[cost > MAXCOST] = MAXCOST
    diag(cost) = 0L
    storage.mode(cost) = "integer"
    downpass.cost = .check.states(phy, data, NSTATES)
    downpass.cost[downpass.cost == 0L] = MAXCOST
    downpass.cost[downpass.cost == 1L] = 0L
    downpass = matrix(0L, nrow=NSTATES, ncol=Nnode(phy))
    uppass = matrix(0L, nrow=NSTATES, ncol=Nnode(phy))
    pscore = integer(Nnode(phy))
    .Call(rcl_sankoff_mpr, phy, cost, uppass, downpass, downpass.cost, pscore)
    .Call(rcl_sankoff_count, phy, cost, uppass, downpass, downpass.cost)
}


sankoff.cost = function(phy, data, cost, NSTATES) {
    stopifnot(is.tree(phy))
    stopifnot(tree.isbinary(phy))
    storage.mode(data) = "integer"
    stopifnot(!is.null(names(data)))
    MAXCOST = 10000L
    stopifnot(NROW(cost) == NCOL(cost) && NROW(cost) == NSTATES)
    cost[cost > MAXCOST] = MAXCOST
    diag(cost) = 0L
    storage.mode(cost) = "integer"
    downpass.cost = .check.states(phy, data, NSTATES)
    downpass.cost[downpass.cost == 0L] = MAXCOST
    downpass.cost[downpass.cost == 1L] = 0L
    uppass.cost = matrix(0L, nrow=NSTATES, ncol=Nnode(phy))
    uppass.cost[, 1L:Ntip(phy)] = downpass.cost[, 1L:Ntip(phy)]
    pscore = integer(Nnode(phy))
    .Call(rcl_sankoff_cost, phy, cost, uppass.cost, downpass.cost, pscore)
    return (uppass.cost)
}
