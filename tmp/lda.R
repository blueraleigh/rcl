lda = function(data, ...) {

    storage.mode(data) = "integer"
    stopifnot(inherits(data, "matrix"))

    vargs = list(...)

    if (is.null(NSWEEP <- vargs$nsweep))
        NSWEEP = 1000L
    if (is.null(BURNIN <- vargs$burnin))
        BURNIN = 100L
    if (is.null(SAMPLEFREQ <- vargs$sample.freq))
        SAMPLEFREQ = 10L
    if (is.null(ALPHA <- vargs$alpha))
        ALPHA = 1
    if (is.null(BETA <- vargs$beta))
        BETA = 1
    if (is.null(NTOPIC <- vargs$ntopic))
        NTOPIC = 2L

    control = list(
        nsweep = as.integer(NSWEEP),
        burnin = as.integer(BURNIN),
        sample_freq = as.integer(SAMPLEFREQ),
        ntopic = as.integer(NTOPIC),
        alpha = as.numeric(ALPHA),
        beta = as.numeric(BETA)
    )

    L = .Call(rcl_lda, control, data)

    dim(L[[1L]]) = c(NSWEEP, 2)
    dim(L[[2L]]) = c(ncol(data), NTOPIC)
    dim(L[[3L]]) = c(NTOPIC, nrow(data))

    return (L)
}


dmm = function(data, ...) {

    storage.mode(data) = "integer"
    stopifnot(inherits(data, "matrix"))

    vargs = list(...)

    if (is.null(NSWEEP <- vargs$nsweep))
        NSWEEP = 1000L
    if (is.null(BURNIN <- vargs$burnin))
        BURNIN = 100L
    if (is.null(SAMPLEFREQ <- vargs$sample.freq))
        SAMPLEFREQ = 10L
    if (is.null(ALPHA <- vargs$alpha))
        ALPHA = 1
    if (is.null(BETA <- vargs$beta))
        BETA = 1
    if (is.null(NTOPIC <- vargs$ntopic))
        NTOPIC = 2L

    control = list(
        nsweep = as.integer(NSWEEP),
        burnin = as.integer(BURNIN),
        sample_freq = as.integer(SAMPLEFREQ),
        ntopic = as.integer(NTOPIC),
        alpha = as.numeric(ALPHA),
        beta = as.numeric(BETA)
    )

    L = .Call(rcl_gsdmm, control, data)

    dim(L[[1L]]) = c(1L, NTOPIC)
    dim(L[[2L]]) = c(NTOPIC, nrow(data))
    dim(L[[3L]]) = c(ncol(data), NTOPIC)

    return (L)
}
