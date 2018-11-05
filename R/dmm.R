dmm = function(phy, data, ...) {
    stopifnot(is.tree(phy))
    vargs = list(...)
    if (is.null(NGEN <- vargs$niter))
        NGEN = 10000L
    if (is.null(PFREQ <- vargs$print.freq))
        PFREQ = 1000L
    if (is.null(WFREQ <- vargs$write.freq))
        WFREQ = 1000L
    if (is.null(PWEIGHT <- vargs$proposal.weight))
        PWEIGHT = c(10, 1, 1)
    if (is.null(TUNE1 <- vargs$tune.rate))
        TUNE1 = 0.5
    if (is.null(TUNE2 <- vargs$tune.clock))
        TUNE2 = 0.5
    if (is.null(PRIOR1 <- vargs$prior.nclock))
        PRIOR1 = 1
    if (is.null(PRIOR2 <- vargs$prior.rate))
        PRIOR2 = 1
    if (is.null(PRIORONLY <- vargs$prior.only))
        PRIORONLY = FALSE
    if (is.null(NSTATE <- vargs$nstate))
        NSTATE = max(floor(prod(dim(data)) / sum(data > 0L)), 2L)
    if (is.null(burnin <- vargs$burnin))
        burnin = 0.1
    burnin = min(as.numeric(burnin), 0.5)
    BURNIN = burnin * NGEN
    settings = list(
        NGEN=as.integer(NGEN),
        BURNIN=as.integer(BURNIN),
        PFREQ=as.integer(PFREQ),
        WFREQ=as.integer(WFREQ),
        PWEIGHT=as.numeric(PWEIGHT),
        TUNE1=as.numeric(TUNE1),
        TUNE2=as.numeric(TUNE2),
        PRIOR1=as.numeric(PRIOR1),
        PRIOR2=as.numeric(PRIOR2),
        PRIORONLY=as.integer(PRIORONLY),
        NSTATE=as.integer(NSTATE)
    )
    if (is.null(alpha <- vargs$alpha))
        alpha = 0.5

    if (is.null(init.rate <- vargs$init.rate))
        init.rate = runif(1, 0, 0.1)

    nsamples = floor((NGEN - BURNIN) / WFREQ)
    mcmc.out = matrix(0, nsamples, 6L)
    clock.prob = numeric(Nnode(phy))
    branch.rate = numeric(Nnode(phy))
    node.state = matrix(0, NSTATE, Nnode(phy))
    state.prob = matrix(0, NSTATE, nrow(data))

    .Call(rcl_dmm_mcmc, phy, as.numeric(init.rate), as.numeric(alpha), data, state.prob, node.state, clock.prob, branch.rate, mcmc.out, settings)

    return (list(
        mcmc.out=mcmc.out,
        state.prob=state.prob,
        clock.prob=clock.prob,
        branch.rate=branch.rate,
        node.state=node.state
    ))
}
