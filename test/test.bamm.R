library(rcl)

phy = read.newick("snaketree.tre")
dat = readRDS("snakedata.rds")

tipdata = rbind(colSums(dat), colSums(dat[c(13, 16, 18, 20), , drop=FALSE]))
storage.mode(tipdata) = "integer"
class(tipdata) = c("bammdata.binomial", "bammdata")

l = bamm(phy, tipdata, niter=1000000, proposal.weight=c(1, 1, 1), return=TRUE)

tipdata = rbind(colSums(dat), colSums(dat[c(2, 3, 9, 11, 17, 19), ]), colSums(dat[c(13, 16, 18, 20), ]))
storage.mode(tipdata) = "integer"
class(tipdata) = c("bammdata.trinomial", "bammdata")

l = bamm(phy, tipdata, niter=1000000, proposal.weight=c(1, 1, 1), return=TRUE, breaks=c(0, 0.1, 0.6, 1))



median(l$param.out[, 5])

# 5 - .0086645
# 10 - .0096
# 15 - 0.011241
# 20 - .011806
# 25 - .0118355

plot(c(5, 10, 15, 20, 25), c(.0086645, .0096, .011241, .011806, .0118355), type='b')

w = Sys.getenv("R_SESSION_TMPDIR")
read.csv(sprintf("%s/param.out", w), header=FALSE)

curve(dbeta(x, l$param.out[1, 2], l$param.out[1, 3]), 0, 1, col=8, ylim=c(0, 10))
for (i in 2:nrow(l$param.out))
    curve(dbeta(x, l$param.out[i, 2], l$param.out[i, 3]), 0, 1, add=TRUE, col=8)
curve(dbeta(x, mean(l$param.out[, 2]), mean(l$param.out[, 3])), 0, 1, add=TRUE, col=2)


L = plot(phy)
points(L[[1]][331, 1], L[[1]][331, 3])



foo = function(n, a, b, tau) {
    r = matrix(0, n, 2)
    r[, 1] = rgamma(n, a)
    r[, 2] = rgamma(n, b)
    for (i in 1:n) {
        if (runif(1) < .5)
            r[i, 1] = r[i, 1] + rgamma(1, tau[1])
        else
            r[i, 2] = r[i, 2] + rgamma(1, tau[2])
    }
    r = r / rowSums(r)

    return (r[, 1])

}


hist(foo(1000, 11, 4, c(13, 10)), breaks='FD', xlim=c(0, 1))

fullphy = read.newick("/Users/mgrundler/Dropbox/manuscripts/snake_diet/datasets/phylogenies/Tonini_dna.tre")
keep = c(tiplabels(fullphy)[grep("Siphlophis", tiplabels(fullphy))],
tiplabels(fullphy)[grep("Drepanoides", tiplabels(fullphy))],
tiplabels(fullphy)[grep("Mussurana", tiplabels(fullphy))],
tiplabels(fullphy)[grep("Paraphimophis", tiplabels(fullphy))],
tiplabels(fullphy)[grep("Phimophis", tiplabels(fullphy))],
tiplabels(fullphy)[grep("Clelia", tiplabels(fullphy))],
tiplabels(fullphy)[grep("Boiruna", tiplabels(fullphy))],
tiplabels(fullphy)[grep("Rachidelus", tiplabels(fullphy))],
tiplabels(fullphy)[grep("Pseudoboa", tiplabels(fullphy))],
tiplabels(fullphy)[grep("Oxyrhopus", tiplabels(fullphy))])

subphy = drop.tip(phy, tiplabels(fullphy)[!tiplabels(fullphy) %in% keep])
