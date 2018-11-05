.check.bammdata = function(phy, categ, data) {

    if (!is.data.frame(data))
        stop("Invalid argument for data")
    if (ncol(data) != 3L)
        stop("Invalid argument for data")
    if (mode(data[, 3L]) != "numeric")
        stop("Invalid argument for data")
    if (mode(data[, 1L]) != "character")
        stop("Invalid argument for data")
    if (mode(data[, 2L]) != "character")
        stop("Invalid argument for data")

    storage.mode(data[, 3L]) = "integer"
    storage.mode(categ) = "integer"

    ctab = xtabs(data[, 3L] ~ data[, 2L] + data[, 1L])

    if (any(categ <= 0) || any(categ > nrow(ctab)))
        stop("Invalid resource categories")

    if (!all(colnames(ctab) %in% tiplabels(phy)) || !all(tiplabels(phy) %in% colnames(ctab)))
        stop("Mismatch between tree and data")

    count = colSums(ctab)
    ctab = ctab[, tiplabels(phy)]

    if (!all(count > 0))
        stop("Missing data not allowed")

    return (rbind(count, ctab[categ, ]))
}


# data is expected to be a three column data frame
# where the first column are the species labels,
# the second column is the resource category, and
# the third column is the count for that category
make.bammdata = function(phy, categ, data) {
    stopifnot(is.tree(phy))
    tipdata = .check.bammdata(phy, categ, data)
    if (nrow(tipdata) == 2L)
        class(tipdata) = c("bammdata.binomial", "bammdata")
    else if (nrow(tipdata) == 3L)
        class(tipdata) = c("bammdata.trinomial", "bammdata")
    else
        stop("Only binomial and trinomial modeling permissible")
    return (tipdata)
}


bamm = function(phy, data, root.p=c("OBS", "FLAT", "ALL"), ...) {
    stopifnot(is.tree(phy))
    stopifnot(inherits(data, "bammdata"))
    vargs = list(...)
    TMPDIR = Sys.getenv("R_SESSION_TMPDIR")
    root.p = match.arg(root.p)
    root.p = switch(root.p, OBS=2L, FLAT=1L, ALL=0L)
    if (is.null(RVAL <- vargs$return))
        RVAL = FALSE
    if (is.null(brks <- vargs$breaks))
        brks = seq(0, 1, 0.1)
    if (is.null(PWEIGHT <- vargs$proposal.weight))
        PWEIGHT = c(1, 1, 1)
    if (is.null(TUNE0 <- vargs$tune.rate))
        TUNE0 = 0.5
    if (is.null(TUNE1 <- vargs$tune.clock))
        TUNE1 = 0.5
    if (is.null(PRIOR1 <- vargs$prior.rate))
        PRIOR1 = 1
    if (is.null(PRIOR2 <- vargs$prior.nclock))
        PRIOR2 = 1
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
    settings = list(
        PWEIGHT=as.numeric(PWEIGHT),
        TUNE0=as.numeric(TUNE0),
        TUNE1=as.numeric(TUNE1),
        PRIOR1=as.numeric(PRIOR1),
        PRIOR2=as.numeric(PRIOR2),
        PRIORONLY=as.integer(PRIORONLY),
        NGEN=as.integer(NGEN),
        PFREQ=as.integer(PFREQ),
        WFREQ=as.integer(WFREQ),
        MCMCFILE=as.character(MCMCFILE),
        CLCKFILE=as.character(CLCKFILE),
        ROOTP=root.p
    )
    if (is.null(init.rate <- vargs$init.rate))
        init.rate = runif(1, 0, 0.1)
    if (!identical(FALSE, is.unsorted(brks)))
        stop("'breaks' must be sorted non-decreasingly and not contain NAs")
    if (max(brks) > 1 || min(brks) < 0)
        stop("'breaks' must be sorted non-decreasingly on (0, 1) and not contain NAs")
    if (inherits(data, "bammdata.binomial"))
        .Call(rcl_bamm_mcmc, phy, as.numeric(init.rate), as.numeric(brks), data, settings, 0L)
    else if (inherits(data, "bammdata.trinomial"))
        .Call(rcl_bamm_mcmc, phy, as.numeric(init.rate), as.numeric(brks), data, settings, 1L)

    if (RVAL) {
        return (list(
            mcmc.out = read.csv(MCMCFILE, header=FALSE, col.names=c("Iteration", "Nclock", "Rate", "Norm", "LogLk", "AcceptFreq")),
            clock.out = read.csv(CLCKFILE, header=FALSE, col.names=c("Iteration", "ClckPos", "ClckRate"))
        ))
    }
}


bamm.bm = function(phy, data, ...) {
    stopifnot(is.tree(phy))
    stopifnot(inherits(data, "bammdata"))
    vargs = list(...)
    TMPDIR = Sys.getenv("R_SESSION_TMPDIR")
    if (is.null(RVAL <- vargs$return))
        RVAL = FALSE
    if (is.null(PWEIGHT <- vargs$proposal.weight))
        PWEIGHT = c(10, 1, 1, 1, 1, 1)
    if (is.null(TUNE0 <- vargs$tune.state))
        TUNE0 = 0.5
    if (is.null(TUNE1 <- vargs$tune.rate))
        TUNE1 = 0.5
    if (is.null(TUNE2 <- vargs$tune.clock))
        TUNE2 = 0.5
    if (is.null(TUNE3 <- vargs$tune.jump))
        TUNE3 = 0.5
    if (is.null(PRIOR2 <- vargs$prior.rate))
        PRIOR2 = 1
    if (is.null(PRIOR3 <- vargs$prior.jump))
        PRIOR3 = 1
    if (is.null(PRIOR4 <- vargs$prior.nclock))
        PRIOR4 = 1
    if (is.null(PRIOR5 <- vargs$prior.njump))
        PRIOR5 = 1
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
    if (is.null(JMPFILE <- vargs$jump.file))
        JMPFILE = sprintf("%s/jump.out", TMPDIR)
    settings = list(
        PWEIGHT=as.numeric(PWEIGHT),
        TUNE0=as.numeric(TUNE0),
        TUNE1=as.numeric(TUNE1),
        TUNE2=as.numeric(TUNE2),
        TUNE3=as.numeric(TUNE3),
        PRIOR2=as.numeric(PRIOR2),
        PRIOR3=as.numeric(PRIOR3),
        PRIOR4=as.numeric(PRIOR4),
        PRIOR5=as.numeric(PRIOR5),
        PRIORONLY=as.integer(PRIORONLY),
        NGEN=as.integer(NGEN),
        PFREQ=as.integer(PFREQ),
        WFREQ=as.integer(WFREQ),
        MCMCFILE=as.character(MCMCFILE),
        CLCKFILE=as.character(CLCKFILE),
        JMPFILE=as.character(JMPFILE)
    )
    if (is.null(init.rate <- vargs$init.rate))
        init.rate = runif(1, 0, 0.1)
    .Call(rcl_bamm_bm_mcmc, phy, as.numeric(init.rate), data, settings)

    if (RVAL) {
        return (list(
            mcmc.out = read.csv(MCMCFILE, header=FALSE, col.names=c("Iteration", "Nclock", "Njump", "Rate", "Norm", "LogLk", "AcceptFreq")),
            clock.out = read.csv(CLCKFILE, header=FALSE, col.names=c("Iteration", "ClckPos", "ClckRate")),
            jump.out = read.csv(JMPFILE, header=FALSE, col.names=c("Iteration", "JumpPos", "JumpSz"))
        ))
    }
}

