bamm.dmm = function(phy, data, ...) {
    stopifnot(is.tree(phy))
    stopifnot(inherits(data, "bammdata"))
    vargs = list(...)
    TMPDIR = Sys.getenv("R_SESSION_TMPDIR")
    if (is.null(RVAL <- vargs$return))
        RVAL = FALSE
    if (is.null(PWEIGHT <- vargs$proposal.weight))
        PWEIGHT = c(10, 1, 1)
    if (is.null(TUNE2 <- vargs$tune.clockrate))
        TUNE2 = 0.5
    if (is.null(PRIOR2 <- vargs$prior.clockrate))
        PRIOR2 = 1
    if (is.null(PRIOR1 <- vargs$prior.nclock))
        PRIOR1 = 1
    if (is.null(PRIORONLY <- vargs$prior.only))
        PRIORONLY = 0L
    if (is.null(NGEN <- vargs$niter))
        NGEN = 10000L
    if (is.null(PFREQ <- vargs$print.freq))
        PFREQ = 1000L
    if (is.null(WFREQ <- vargs$write.freq))
        WFREQ = 1000L
    if (is.null(MCMCFILE <- vargs$mcmc.file))
        MCMCFILE = sprintf("%s/mcmc.out", TMPDIR)
    if (is.null(CLCKFILE <- vargs$clock.file))
        CLCKFILE = sprintf("%s/clock.out", TMPDIR)
    if (is.null(EVNTFILE <- vargs$event.file))
        EVNTFILE = sprintf("%s/event.out", TMPDIR)
    settings = list(
        PWEIGHT=as.numeric(PWEIGHT),
        TUNE2=as.numeric(TUNE2),
        PRIOR1=as.numeric(PRIOR1),
        PRIOR2=as.numeric(PRIOR2),
        PRIORONLY=as.integer(PRIORONLY),
        NGEN=as.integer(NGEN),
        PFREQ=as.integer(PFREQ),
        WFREQ=as.integer(WFREQ),
        MCMCFILE=as.character(MCMCFILE),
        CLCKFILE=as.character(CLCKFILE),
        EVNTFILE=as.character(EVNTFILE)
    )
    if (is.null(nevent.expected <- vargs$nevent.expected))
        nevent.expected = 1
    .Call(rcl_dmm_mcmc, phy, 1, as.numeric(nevent.expected), data, settings)

    if (RVAL) {
        return (list(
            mcmc.out = read.csv(MCMCFILE, header=FALSE, col.names=c("Iteration", "Nevent", "Nclock", "LogLk", "AcceptFreq")),
            clock.out = read.csv(CLCKFILE, header=FALSE, col.names=c("Iteration", "ClckPos", "ClckRate")),
            event.out = read.csv(EVNTFILE, header=FALSE, col.names=c("Iteration", "EventPos"))
        ))
    }
}



bamm.dmm2 = function(phy, data, ...) {
    stopifnot(is.tree(phy))
    stopifnot(inherits(data, "bammdata"))
    vargs = list(...)
    TMPDIR = Sys.getenv("R_SESSION_TMPDIR")
    if (is.null(RVAL <- vargs$return))
        RVAL = FALSE
    if (is.null(PRIOR0 <- vargs$prior.nevent))
        PRIOR0 = 1
    if (is.null(PRIORONLY <- vargs$prior.only))
        PRIORONLY = 0L
    if (is.null(NGEN <- vargs$niter))
        NGEN = 10000L
    if (is.null(PFREQ <- vargs$print.freq))
        PFREQ = 1000L
    if (is.null(WFREQ <- vargs$write.freq))
        WFREQ = 1000L
    if (is.null(MCMCFILE <- vargs$mcmc.file))
        MCMCFILE = sprintf("%s/mcmc.out", TMPDIR)
    if (is.null(EVNTFILE <- vargs$event.file))
        EVNTFILE = sprintf("%s/event.out", TMPDIR)
    settings = list(
        PRIOR0=as.numeric(PRIOR0),
        PRIORONLY=as.integer(PRIORONLY),
        NGEN=as.integer(NGEN),
        PFREQ=as.integer(PFREQ),
        WFREQ=as.integer(WFREQ),
        MCMCFILE=as.character(MCMCFILE),
        EVNTFILE=as.character(EVNTFILE)
    )
    .Call(rcl_dmm_mcmc2, phy, 1, data, settings)

    if (RVAL) {
        return (list(
            mcmc.out = read.csv(MCMCFILE, header=FALSE, col.names=c("Iteration", "Nevent", "LogLk", "AcceptFreq")),
            event.out = read.csv(EVNTFILE, header=FALSE, col.names=c("Iteration", "EventPos"))
        ))
    }
}


bamm.dmm4 = function(phy, data, ...) {
    stopifnot(is.tree(phy))
    stopifnot(inherits(data, "bammdata"))
    vargs = list(...)
    TMPDIR = Sys.getenv("R_SESSION_TMPDIR")
    if (is.null(RVAL <- vargs$return))
        RVAL = FALSE
    if (is.null(NGEN <- vargs$niter))
        NGEN = 10000L
    if (is.null(PFREQ <- vargs$print.freq))
        PFREQ = 1000L
    if (is.null(WFREQ <- vargs$write.freq))
        WFREQ = 1000L
    if (is.null(MCMCFILE <- vargs$mcmc.file))
        MCMCFILE = sprintf("%s/mcmc.out", TMPDIR)
    if (is.null(EVNTFILE <- vargs$event.file))
        EVNTFILE = sprintf("%s/event.out", TMPDIR)
    if (is.null(PWEIGHT <- vargs$proposal.weight))
        PWEIGHT = c(3, 1, 0)
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
        MCMCFILE=as.character(MCMCFILE),
        EVNTFILE=as.character(EVNTFILE),
        PWEIGHT=as.numeric(PWEIGHT),
        PRIORONLY=as.integer(PRIORONLY)
    )
    if (is.null(nevent.expected <- vargs$nevent.expected))
        nevent.expected = 1
    if (is.null(alpha <- vargs$alpha))
        alpha = 1

    nsamples = floor((NGEN - BURNIN) / WFREQ)
    mcmc.out = matrix(0, nsamples, 6L)
    event.config = matrix(0L, Nnode(phy), nsamples)
    event.prob = numeric(Nnode(phy))

    .Call(rcl_dmm_mcmc4, phy, as.numeric(alpha), as.numeric(nevent.expected), data, event.config, event.prob, mcmc.out, settings)

    return (list(
        mcmc.out = mcmc.out,
        event.config = event.config,
        event.prob = event.prob
    ))
    #if (RVAL) {
    #    return (list(
    #        mcmc.out = read.csv(MCMCFILE, header=FALSE, col.names=c("Iteration", "Nevent", "LogLk", "LogPrior", "Alpha", "AcceptFreq")),
    #        #event.out = read.csv(EVNTFILE, header=FALSE, col.names=c("Iteration", "NodeId"))
    #        event.out = read.csv(EVNTFILE, header=FALSE, col.names=c("NodeId", "Posterior", "Prior"))
    #    ))
    #}
}



bamm.dmm5 = function(phy, data, ...) {
    stopifnot(is.tree(phy))
    stopifnot(inherits(data, "bammdata"))
    vargs = list(...)
    if (is.null(RVAL <- vargs$return))
        RVAL = FALSE
    if (is.null(NGEN <- vargs$niter))
        NGEN = 10000L
    if (is.null(PFREQ <- vargs$print.freq))
        PFREQ = 1000L
    if (is.null(WFREQ <- vargs$write.freq))
        WFREQ = 1000L
    if (is.null(PWEIGHT <- vargs$proposal.weight))
        PWEIGHT = c(3, 1, 1, 1)
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
        PRIORONLY=as.integer(PRIORONLY)
    )
    if (is.null(alpha <- vargs$alpha))
        alpha = 1

    if (is.null(init.rate <- vargs$init.rate))
        init.rate = runif(1, 0, 0.1)

    nsamples = floor((NGEN - BURNIN) / WFREQ)
    mcmc.out = matrix(0, nsamples, 7L)
    event.prob = numeric(Nnode(phy))
    clock.prob = numeric(Nnode(phy))
    branch.rate = numeric(Nnode(phy))

    .Call(rcl_dmm_mcmc5, phy, as.numeric(init.rate), as.numeric(alpha), data, event.prob, clock.prob, branch.rate, mcmc.out, settings)

    return (list(
        mcmc.out = mcmc.out,
        event.prob = event.prob,
        clock.prob = clock.prob,
        branch.rate = branch.rate
    ))
}
