.check.multinomial.data = function(phy, data1, data2) {

    if (is.null(names(data1)))
        stop("Data requires a names component")
    if (NROW(data1) != Ntip(phy) && NCOL(data1) != 1L)
        stop("Invalid argument for data1")
    if (!is.data.frame(data2))
        stop("Invalid argument for data2")
    if (ncol(data2) != 3L)
        stop("Invalid argument for data2")
    if (mode(data2[, 3L]) != "numeric")
        stop("Invalid argument for data2")
    if (mode(data2[, 1L]) != "character")
        stop("Invalid argument for data2")
    if (mode(data2[, 2L]) != "character")
        stop("Invalid argument for data2")

    storage.mode(data2[, 3L]) = "integer"
    storage.mode(data1) = "integer"

    if (!all(names(data1) %in% tiplabels(phy)) || !all(tiplabels(phy) %in% names(data1)))
        stop("Mismatch between tree and data")

    ctab = xtabs(data2[, 3L] ~ data2[, 2L] + data2[, 1L])

    if (!all(colnames(ctab) %in% tiplabels(phy)) || !all(tiplabels(phy) %in% colnames(ctab)))
        stop("Mismatch between tree and data")

    count = data1[tiplabels(phy)]
    ctab = ctab[, tiplabels(phy)]

    if (!all(count > 0))
        stop("Missing data not allowed")

    tipdata = rbind(count, ctab)
    return (tipdata)
}


# data1 is a vector giving total number of observations
# for each species. names correspond to species labels.
# data2 is expected to be a three column data frame
# where the first column are the species labels,
# the second column is the resource category, and
# the third column is the count for that category
make.multinomial.data = function(phy, data1, data2) {
    stopifnot(is.tree(phy))
    tipdata = .check.multinomial.data(phy, data1, data2)
    class(tipdata) = "multinomial.data"
    return (tipdata)
}


# categ is the focal resource category to analyze
make.mk.betabinomial = function(phy, data, categ, NSTATES, model=c("mk", "mk.covarion1", "mk.covarion2")) {
    stopifnot(is.tree(phy))
    stopifnot(inherits(data, "multinomial.data"))
    storage.mode(categ) = "integer"

    NCAT = nrow(data)-1L

    stopifnot(categ <= NCAT && categ > 0L)


    model = match.arg(model)

    clk = if (model == "mk") {
        is.covarion = 0L
        matrix(0, nrow=NSTATES, ncol=Nnode(phy))
    } else {
        is.covarion = 1L
        matrix(0, nrow=2L*NSTATES, ncol=Nnode(phy))
    }

    .Call(rcl_mkp_betabin_init, phy, NSTATES, NCAT, categ, data, clk, is.covarion)

    if (!is.covarion) {
        qmat = matrix(0, NSTATES, NSTATES)
        npar = 1L
    } else {
        qmat = matrix(0, 2L*NSTATES, 2L*NSTATES)
        if (model == "mk.covarion1")
            npar = 2L
        else
            npar = 3L
    }

    # if this is a regular mk model there is
    # a single rate parameter. if this is a covarion model there are
    # either 2 or 3 rate parameters, the last of which is the evolutionary
    # rate of switching between observable states. thus, the first or
    # first two rates are the rates of switching on/off
    mk.lik = function(pars) {
        .Call(rcl_mkp_betabin_qmat, NSTATES, pars, is.covarion, qmat)
        .Call(rcl_mkp_loglk, phy, qmat, clk)
    }

    class(mk.lik) = if (is.covarion) {
        c("mk.betabinomial", "mk.covarion", "mk")
    } else {
        c("mk.betabinomial", "mk")
    }
    return (mk.lik)
}


find.mle.mk.betabinomial = function(lik, init, method=c("Nelder-Mead",
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

    npar = get("npar", environment(lik))

    if (length(init) != npar)
        stop(sprintf("Expected %d parameters, received %d", npar, length(init)))

    stopifnot(all(init > 0))

    fit = optim(log(init), function(par) -lik(exp(par)), method=method, hessian=TRUE, control=list(...))

    fit$par = par.trans(fit$par, sqrt(diag(solve(fit$hessian))))

    fit$loglk = -fit$value
    fit$value = NULL

    return (fit)
}


# categ are the two focal resource categories to analyze
make.mk.dirichletrinomial = function(phy, data, categ, NPART, model=c("mk", "mk.covarion1", "mk.covarion2")) {
    stopifnot(is.tree(phy))
    stopifnot(inherits(data, "multinomial.data"))
    storage.mode(categ) = "integer"

    NCAT = nrow(data)-1L

    stopifnot(length(categ) == 2L)
    stopifnot(all(categ <= NCAT) && all(categ > 0L))

    NSTATES = NPART*NPART

    model = match.arg(model)

    clk = if (model == "mk") {
        is.covarion = 0L
        matrix(0, nrow=NSTATES, ncol=Nnode(phy))
    } else {
        is.covarion = 1L
        matrix(0, nrow=2L*NSTATES, ncol=Nnode(phy))
    }

    bpar.idx = 1L:3L

    # pars is a vector of shape parameters for the beta distribution
    mk.init.lik = function(pars) {
        .Call(rcl_mkp_dirichtri_init, phy, NPART, NCAT, categ, pars, data, clk, is.covarion)
    }

    if (!is.covarion) {
        qmat = matrix(0, NSTATES, NSTATES)
        npar = 4L
    } else {
        qmat = matrix(0, 2L*NSTATES, 2L*NSTATES)
        if (model == "mk.covarion1")
            npar = 5L
        else
            npar = 6L
    }

    # pars is a vector of shape parameters (described above)
    # and evolutionary rate parameters. the rate parameters follow
    # the shape parameters. if this is a regular mk model there is
    # a single rate parameter. if this is a covarion model there are
    # either 2 or 3 rate parameters, the last of which is the evolutionary
    # rate of switching between observable states. thus, the first or
    # first two rates are the rates of switching on/off
    mk.lik = function(pars) {
        mk.init.lik(pars[bpar.idx])
        .Call(rcl_mkp_dirichtri_qmat, NPART, pars, is.covarion, qmat)
        .Call(rcl_mkp_loglk, phy, qmat, clk)
    }

    class(mk.lik) = if (is.covarion) {
        c("mk.dirichletrinomial", "mk.covarion", "mk")
    } else {
        c("mk.dirichletrinomial", "mk")
    }
    return (mk.lik)
}


find.mle.mk.dirichletrinomial = function(lik, init, method=c("Nelder-Mead",
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

    npar = get("npar", environment(lik))

    if (length(init) != npar)
        stop(sprintf("Expected %d parameters, received %d", npar, length(init)))

    stopifnot(all(init > 0))

    fit = optim(log(init), function(par) -lik(exp(par)), method=method, hessian=TRUE, control=list(...))

    pars = par.trans(fit$par, sqrt(diag(solve(fit$hessian))))

    beta.pars = pars[1L:3L, ]
    rate.pars = pars[-(1L:3L), ]

    fit$par = list(
        shape=beta.pars,
        rate=rate.pars
    )

    fit$loglk = -fit$value
    fit$value = NULL

    return (fit)
}
