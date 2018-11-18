bmm = function(phy, data, ...) {
    stopifnot(is.tree(phy))
    stopifnot(!anyNA(data))
    vargs = list(...)
    if (is.null(NGEN <- vargs$niter))
        NGEN = 10000L
    if (is.null(PFREQ <- vargs$print.freq))
        PFREQ = 1000L
    if (is.null(WFREQ <- vargs$write.freq))
        WFREQ = 1000L
    if (is.null(PWEIGHT <- vargs$proposal.weight))
        PWEIGHT = c(1, 2, 1, 2)
    if (is.null(TUNE0 <- vargs$tune.rate))
        TUNE0 = 0.5
    if (is.null(TUNE1 <- vargs$tune.clock))
        TUNE1 = 0.5
    if (is.null(TUNE3 <- vargs$tune.jump))
        TUNE3 = 0.1
    if (is.null(PRIOR0 <- vargs$prior.nclock))
        PRIOR0 = 1
    if (is.null(PRIOR1 <- vargs$prior.rate))
        PRIOR1 = 1
    if (is.null(PRIOR2 <- vargs$prior.njump))
        PRIOR2 = 1
    if (is.null(PRIOR3 <- vargs$prior.jump))
        PRIOR3 = 1
    if (is.null(PRIORONLY <- vargs$prior.only))
        PRIORONLY = FALSE
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
        TUNE0=as.numeric(TUNE0),
        TUNE1=as.numeric(TUNE1),
        TUNE3=as.numeric(TUNE3),
        PRIOR0=as.numeric(PRIOR0),
        PRIOR1=as.numeric(PRIOR1),
        PRIOR2=as.numeric(PRIOR2),
        PRIOR3=as.numeric(PRIOR3),
        PRIORONLY=as.integer(PRIORONLY)
    )

    if (is.null(init.rate <- vargs$init.rate))
        init.rate = runif(1, 0, 0.1)

    nsamples = floor((NGEN - BURNIN) / WFREQ)
    mcmc.out = matrix(0, nsamples, 6L)
    clock.prob = numeric(Nnode(phy))
    jump.prob = numeric(Nnode(phy))
    branch.rate = numeric(Nnode(phy))
    branch.jump = numeric(Nnode(phy))
    node.state = numeric(Nnode(phy))

    .Call(rcl_bmm_mcmc, phy, as.numeric(init.rate), data, node.state, clock.prob, jump.prob, branch.rate, branch.jump, mcmc.out, settings)

    return (list(
        mcmc.out=mcmc.out,
        clock.prob=clock.prob,
        jump.prob=jump.prob,
        branch.rate=branch.rate,
        branch.jump=branch.jump,
        node.state=node.state
    ))
}
