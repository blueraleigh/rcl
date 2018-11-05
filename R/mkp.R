.check.state.probs = function(phy, data, NSTATES) {
    storage.mode(NSTATES) = "integer"
    storage.mode(data) = "double"
    stopifnot(nrow(data) == NSTATES || ncol(data) == NSTATES)
    if (sum(data) != Ntip(phy))
        stop("Some state probabilities are incorrect.")
    if (ncol(data) == NSTATES)
        data = t(data)
    data = data[, tiplabels(phy)]
    data = cbind(data, matrix(0, nrow=NSTATES, ncol=Nnode(phy)-Ntip(phy)))
    return (data)
}


.check.covarion.states = function(phy, data, NSTATES) {
    storage.mode(NSTATES) = "integer"
    stopifnot(class(data) == "numeric" || class(data) == "integer")
    storage.mode(data) == "integer"
    stopifnot(!is.null(names(data)))
    clk = .check.states(phy, data, NSTATES)
    clk = clk[, tiplabels(phy)]
    clk = cbind(clk, matrix(0L, nrow=NSTATES, ncol=Nnode(phy)-Ntip(phy)))
    clk = rbind(clk, matrix(0L, nrow=NSTATES, ncol=Nnode(phy)))
    for (i in 1L:Nnode(phy)) {
        hidden.state = which(clk[, i] == 1L)
        obs.state = hidden.state + NSTATES
        clk[obs.state, i] = 1L
    }
    storage.mode(clk) = "double"
    return (clk)
}


.check.covarion.state.probs = function(phy, data, NSTATES) {
    storage.mode(NSTATES) = "integer"
    storage.mode(data) = "double"
    stopifnot(nrow(data) == NSTATES || ncol(data) == NSTATES)
    if (sum(data) != Ntip(phy))
        stop("Some state probabilities are incorrect.")
    if (ncol(data) == NSTATES)
        data = t(data)
    data = data[, tiplabels(phy)]
    data = cbind(data, matrix(0, nrow=NSTATES, ncol=Nnode(phy)-Ntip(phy)))
    data = rbind(data, data)
    return (data)
}


# transition rates are numbered starting from
# 1 and proceeding across columns in row-major
# ordering
make.qmat = function(NSTATES) {
    nrate = NSTATES * (NSTATES-1L)

    # not all rates may be free to vary so npar
    # may be different than nrate
    npar = nrate

    # the i-th value returns the index of the i-th rate
    # in the "rate" argument passed to the "qmat" function
    rate.idx = 1L:nrate

    Qmat = matrix(0, NSTATES, NSTATES)
    rownames(Qmat) = colnames(Qmat) = as.character(0L:(NSTATES-1L))
    Qmat.free = rep(TRUE, nrate)    # rates that are free to vary (i.e. not zero)
    Qmat.coo = matrix(0L, nrow=nrate, ncol=2L)

    cntr = 1L
    for (i in 1L:NSTATES) {
        for (j in 1L:NSTATES) {
            if (i == j)
                next
            Qmat.coo[cntr, ] = c(i, j)
            cntr = cntr + 1L
        }
    }
    rm(cntr, i, j)

    qmat = function(rate) {
        diag(Qmat) = 0
        Qmat[Qmat.coo[Qmat.free, , drop=FALSE]] = rate[rate.idx[Qmat.free]]
        Qmat[Qmat.coo[!Qmat.free, , drop=FALSE]] = 0
        diag(Qmat) = -rowSums(Qmat)
        return (Qmat)
    }

    class(qmat) = "qmat"
    return (qmat)
}


is.qmat = function(x) {
    inherits(x, "qmat")
}


print.qmat = function(x, ...) {
    cat("CTMC transition rate matrix.\nNumbers refer to parameter indices.\n")
    npar = environment(x)$npar
    Qmat = x(1L:npar)
    diag(Qmat) = "--"
    class(Qmat) = "table"
    print(Qmat, zero.print=".", na.print="NA")
}


# constrain the transition rate represented by row and
# col. v represents the index of a parameter. if it is 0
# then that rate is fixed at zero and not estimated
constrain.qmat = function(qmat, row, col, v) {
    stopifnot(is.qmat(qmat))

    Qmat.coo = get("Qmat.coo", environment(qmat), inherits=FALSE)
    rate.idx = get("rate.idx", environment(qmat), inherits=FALSE)

    A = as.integer(row) == Qmat.coo[, 1L]
    B = as.integer(col) == Qmat.coo[, 2L]

    idx = which(A & B)

    if (!length(idx) || length(idx) > 1L)
        stop("Invalid constraint")

    rate.idx[idx] = as.integer(v)

    Qmat.free = rate.idx > 0L
    npar = length(unique(rate.idx[Qmat.free]))

    assign("rate.idx", rate.idx, environment(qmat))
    assign("Qmat.free", Qmat.free, environment(qmat))
    assign("npar", npar, environment(qmat))
}

#constrain.qmat = function(qmat, row, col, v) {
#    stopifnot(is.qmat(qmat))
#    rho = as.list(environment(qmat))
#
#    NSTATES = rho$NSTATES
#    nrate = rho$nrate
#    rate.idx = rho$rate.idx
#    Qmat = rho$Qmat
#    Qmat.coo = rho$Qmat.coo
#
#    A = as.integer(row) == Qmat.coo[, 1L]
#    B = as.integer(col) == Qmat.coo[, 2L]
#
#    idx = which(A & B)
#
#    if (!length(idx) || length(idx) > 1L)
#        stop("Invalid constraint")
#
#    rate.idx[idx] = as.integer(v)
#
#    Qmat.free = rate.idx > 0L
#    npar = length(unique(rate.idx[Qmat.free]))
#
#    rm (idx, A, B, row, col, v, rho)
#
#    qmat = function(rate) {
#        if (length(rate) != npar)
#            stop(sprintf("Expected %d parameters, received %d", npar, length(rate)))
#        if (any(rate.idx > npar))
#            stop("Improperly constrained rate matrix")
#        diag(Qmat) = 0
#        Qmat[Qmat.coo[Qmat.free, , drop=FALSE]] = rate[rate.idx[Qmat.free]]
#        Qmat[Qmat.coo[!Qmat.free, , drop=FALSE]] = 0
#        diag(Qmat) = -rowSums(Qmat)
#        return (Qmat)
#    }
#
#    class(qmat) = "qmat"
#    return (qmat)
#}


make.mk = function(phy, data, NSTATES) {
    stopifnot(is.tree(phy))
    if (class(data) == "numeric" || class(data) == "integer") {
        stopifnot(!is.null(names(data)))
        clk = .check.states(phy, data, NSTATES)
        storage.mode(clk) = "double"
    } else if (class(data) == "matrix") {
        stopifnot(!is.null(rownames(data)) || !is.null(colnames(data)))
        clk = .check.state.probs(phy, data, NSTATES)
    } else {
        stop("Unrecognized data input.")
    }
#    qmat = make.qmat(NSTATES)
    mk.lik = function(rate) {
        .Call(rcl_mkp_loglk, phy, qmat(rate), clk)
    }
    class(mk.lik) = "mk"
    return (mk.lik)
}


is.mk = function(x) {
    inherits(x, "mk")
}

make.mk.covarion = function(phy, data, NSTATES) {
    stopifnot(is.tree(phy))
    if (class(data) == "numeric" || class(data) == "integer") {
        clk = .check.covarion.states(phy, data, NSTATES)
    } else if (class(data) == "matrix") {
        stopifnot(!is.null(rownames(data)) || !is.null(colnames(data)))
        clk = .check.covarion.state.probs(phy, data, NSTATES)
    } else {
        stop("Unrecognized data input.")
    }
    #qmat = make.qmat(2L*NSTATES)

    #for (i in 1L:NSTATES) {
        # first NSTATES are the hidden "off" states
        # these can only change to "on" and cannot
        # transition to any other observable state
    #    for (j in 1L:(2L*NSTATES)) {
    #        if (i == j)
    #            next
    #        if (j == (i + NSTATES)) {
    #            qmat = constrain.qmat(qmat, i, j, 1L)
    #            qmat = constrain.qmat(qmat, j, i, 2L)
    #        } else {
    #            qmat = constrain.qmat(qmat, i, j, 0L)
    #            qmat = constrain.qmat(qmat, j, i, 0L)
    #        }
    #    }
    #}
    #cntr = 1L
    #for (i in (NSTATES+1L):(2L*NSTATES)) {
    #    for (j in (NSTATES+1L):(2L*NSTATES)) {
    #        if (i == j)
    #            next
    #        qmat = constrain.qmat(qmat, i, j, 2L+cntr)
    #        cntr = cntr+1L
    #    }
    #}

    #rm (i, j, cntr)

    mk.lik = function(rate) {
        .Call(rcl_mkp_loglk, phy, qmat(rate), clk)
    }
    class(mk.lik) = c("mk.covarion", "mk")
    return (mk.lik)
}

is.mk.covarion = function(x) {
    inherits(x, "mk.covarion")
}


# bind a rate matrix to a likelihood function
bind.qmat = function(qmat, lik) {
    stopifnot(is.mk(lik))
    stopifnot(is.qmat(qmat))
    assign("qmat", qmat, environment(lik))
}


#constrain.mk = function(lik, row, col, v) {
#    clk = environment(lik)$clk
#    phy = environment(lik)$phy
#    qmat = constrain.qmat(environment(lik)$qmat, row, col, v)
#    mk.lik = function(rate) {
#        .Call(rcl_mkp_loglk, phy, qmat(rate), clk)
#    }
#    return (mk.lik)
#}


asr.mk = function(rate, lik, exclude.tips=TRUE) {
    stopifnot(is.mk(lik))
    rho = environment(lik)
    .Call(rcl_mkp_asr, rho$phy, rho$qmat(rate), rho$clk, as.integer(exclude.tips))
}


find.mle = function(lik, init, method=c("Nelder-Mead",
    "BFGS", "CG", "L-BFGS-B", "SANN", "Brent"), ...)
{
    UseMethod("find.mle")
}


find.mle.mk = function(lik, init, method=c("Nelder-Mead",
    "BFGS", "CG", "L-BFGS-B", "SANN", "Brent"), ...)
{
    method = match.arg(method)

    # @param x logged parameter estimates
    # @param v standard deviation of those estimates
    # @return param estimates and confidence interval on unlogged scale
    par.trans = function(x, v) {
        l = x - 1.96*v
        u = x + 1.96*v
        lower = exp(x) * exp(l - x)
        upper = exp(x) * exp(u - x)
        estimate = exp(x)
        return (cbind(estimate, lower, upper))
    }

    qmat = environment(lik)$qmat
    if (is.null(qmat))
        stop("Likelihood function has no rate matrix!")
    npar = get("npar", environment(qmat), inherits=FALSE)

    if (length(init) != npar)
        stop(sprintf("Expected %d parameters, received %d", npar, length(init)))

    stopifnot(all(init > 0))

    fit = optim(log(init), function(par) -lik(exp(par)), method=method, hessian=TRUE, control=list(...))

    fit$par = par.trans(fit$par, sqrt(diag(solve(fit$hessian))))
    fit$loglk = -fit$value
    fit$value = NULL

    return (fit)
}


find.mle.mk.covarion = function(lik, init, method=c("Nelder-Mead",
    "BFGS", "CG", "L-BFGS-B", "SANN", "Brent"), ...)
{
    NextMethod("find.mle")
}


Q = function(lik) {
    stopifnot(is.mk(lik))
    qmat = environment(lik)$qmat
    if (is.null(qmat))
        stop("Likelihood function has no assigned rate matrix!")
    print(qmat)
}


stochastic.map.clade = function(rate, lik) {
    stopifnot(is.mk(lik))
    phy = get("phy", environment(lik), inherits=FALSE)
    qmat = environment(lik)$qmat
    if (is.null(qmat))
        stop("Likelihood function has no rate matrix!")
    qmat = qmat(rate)
    clk = get("clk", environment(lik), inherits=FALSE)
    cache = .Call(rcl_mkp_cache_build, phy, qmat, clk)
    count = matrix(1L, nrow=nrow(qmat), ncol=nrow(qmat))
    mk.smap = function(node, nsim=1L, prior=FALSE) {
        storage.mode(node) = "integer"
        storage.mode(nsim) = "integer"
        stopifnot(node > Ntip(phy) && node <= Nnode(phy))
        .Call(rcl_mkp_smap_clade, phy, node-1L, count, cache, nsim, as.integer(prior))
    }
    class(mk.smap) = c("stochastic.map", "stochastic.map.clade")
    return (mk.smap)
}


stochastic.map.branch = function(rate, lik) {
    stopifnot(is.mk(lik))
    phy = get("phy", environment(lik), inherits=FALSE)
    qmat = environment(lik)$qmat
    if (is.null(qmat))
        stop("Likelihood function has no rate matrix!")
    qmat = qmat(rate)
    clk = get("clk", environment(lik), inherits=FALSE)
    cache = .Call(rcl_mkp_cache_build, phy, qmat, clk)
    count = matrix(1L, nrow=nrow(qmat), ncol=nrow(qmat))
    mk.smap = function(node, nsim=1L, prior=FALSE) {
        storage.mode(node) = "integer"
        storage.mode(nsim) = "integer"
        stopifnot(all(node > 0) && all(node <= Nnode(phy)) && all(node != root(phy)))
        .Call(rcl_mkp_smap_branchset, phy, node-1L, count, cache, nsim, as.integer(prior))
    }
    class(mk.smap) = c("stochastic.map", "stochastic.map.branch")
    return (mk.smap)
}


is.stochastic.map = function(x) {
    inherits(x, "stochastic.map")
}


# By default the stochastic map will count all transitions.
# To constrain it we pass the transitions we want to ignore.
constrain.map = function(map, row, col) {
    stopifnot(is.stochastic.map(map))
    count = get("count", environment(map))
    count[row, col] = 0L
    assign("count", count, environment(map))
}

# Unconstrain a map to count transitions previously ignored
unconstrain.map = function(map, row, col) {
    stopifnot(is.stochastic.map(map))
    count = get("count", environment(map))
    count[row, col] = 1L
    assign("count", count, environment(map))
}


#constrain.map = function(map, row, col) {
#    stopifnot(is.stochastic.map(map))
#    rho = environment(map)
#    stopifnot("mk.smap" %in% ls(rho))
#    phy = rho$phy
#    qmat = rho$qmat
#    clk = rho$clk
#    cache = rho$cache
#    count = rho$count
#    count[row, col] = 0L
#    type = rho$type
#    rm(rho)
#    if (type == "clade") {
#        mk.smap = function(node, nsim=1L, prior=FALSE) {
#            storage.mode(node) = "integer"
#            storage.mode(nsim) = "integer"
#            stopifnot(node > Ntip(phy) && node <= Nnode(phy))
#            .Call(rcl_mkp_smap_clade, phy, node-1L, count, cache, nsim, as.integer(prior))
#        }
#    } else if (type == "branch") {
#        mk.smap = function(node, nsim=1L, prior=FALSE) {
#            storage.mode(node) = "integer"
#            storage.mode(nsim) = "integer"
#            stopifnot(all(node > 0) && all(node <= Nnode(phy)) && all(node != root(phy)))
#            .Call(rcl_mkp_smap_branchset, phy, node-1L, count, cache, nsim, as.integer(prior))
#        }
#    } else {
#        error("Unrecognized type")
#    }
#    return (mk.smap)
#}

