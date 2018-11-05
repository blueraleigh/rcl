#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include "node.h"
#include "tree.h"

/*
** Incomplete beta function
**
** Evaluate the beta from over (0, x) rather than (0, 1).
** This is just the cumulative distribution of the beta
** distribution multiplied by the beta function.
*/
static double ibeta(double x, double a, double b)
{
    return exp(pbeta(x, a, b, 1, 1) + lbeta(a, b));
}


/*
** Beta binomial probability
**
** Evaluate the probability of k successes in n trials
** given that k is drawn from a binomial distribution
** and p, the parameter of the binomial distribution, is
** drawn from a beta distribution on the interval (l, u)
**
** n     -- number of trials
** k     -- number of successes
** lcoef -- log binomial coefficient
** a,b   -- shape parameters of the beta distribution
** l,u   -- lower and upper bound over which the binomial distribution
**          compounded with the beta distribution is evaluated
*/
static double pbetabin0(int n, int k, double l, double u, double lcoef, double a, double b)
{
    return exp(lcoef + log(ibeta(u, k+a, n-k+b) -
        ibeta(l, k+a, n-k+b)) - lbeta(a, b));
}


/*
** Compute the beta binomial probability over
** the interval (state*step, state*step + step)
*/
static double pbetabin(int state, int n, int k, double step, double lcoef, double a, double b)
{
    double l = state * step;
    return pbetabin0(n, k, l, l+step, lcoef, a, b);
}

/*
** Compute average value of beta distribution over the interval
** (state*step, state*step + step)
*/
static double bavg(int state, double step, double a, double b)
{
    double l = state * step;
    double p = pbeta(l+step, a, b, 1, 0) - pbeta(l, a, b, 1, 0);
    return (ibeta(l+step, a+1, b) - ibeta(l, a+1, b)) / (p*beta(a, b));
}


/*
** data -- Data for the tips. (NCAT+1) by NTIP matrix. The first row contains the
**         number of observations and subsequent rows contain the number of observations
**         in each resource category.
** clk -- Downpass conditional likelihoods. NSTATES by NNODE matrix.
*/
SEXP rcl_mkp_betabin_init(SEXP rtree, SEXP NSTATES, SEXP NCAT, SEXP CATEG, SEXP data, SEXP clk, SEXP iscovarion)
{
    int i;
    int j;
    int n;
    int k;
    int nstate0 = INTEGER(NSTATES)[0];          // number of bins for the beta distribution
    int ncat = INTEGER(NCAT)[0];                // number of resource categories in the data
    int *categ = INTEGER(CATEG);                // resource categories being analyzed
    int is_covarion = INTEGER(iscovarion)[0];
    int *tipdata = INTEGER(data);
    int nstate;

    double step = 1 / (double)nstate0;
    double *dclk = REAL(clk);
    double cdf[nstate0+1];

    struct tree *tree = (struct tree*)R_ExternalPtrAddr(rtree);

    cdf[0] = 0;
    cdf[nstate0] = 1;

    // if this is a covarion model the conditional likelihood matrix
    // for each resource category actually has twice the number
    // of rows as the "true" number of states
    if (is_covarion)
        nstate = 2 * nstate0;
    else
        nstate = nstate0;

    // access into the conditional likelihood array
    #define ELEM1(m, i, j) m[(i) + (j)*nstate]
    #define COLUMN1(m, j) m + (j)*nstate

    // access into the tipdata array
    #define ELEM2(m, i, j) m[(i) + (j)*(ncat+1)]

    for (i = 0; i < tree->ntip; ++i) {

        n = ELEM2(tipdata, 0, i);
        k = ELEM2(tipdata, categ[0], i);

        for (j = 1; j < LENGTH(CATEG); ++j)
            k += ELEM2(tipdata, categ[j], i);

        for (j = 1; j < nstate0; ++j)
            cdf[j] = pbeta(j*step, k+1, n-k+1, 1, 0);

        for (j = 0; j < nstate0; ++j)
            ELEM1(dclk, j, i) = cdf[j+1] - cdf[j];

        if (is_covarion)
            memcpy(COLUMN1(dclk, i) + nstate0, COLUMN1(dclk, i), nstate0 * sizeof(double));

    }

    return R_NilValue;
}


SEXP rcl_mkp_betabin_qmat(SEXP nstate, SEXP pars, SEXP iscovarion, SEXP QMAT)
{
    int i;
    int j;
    int lda;
    int npar = LENGTH(pars);
    int nstate0 = INTEGER(nstate)[0];
    int is_covarion = INTEGER(iscovarion)[0];
    double f = 1/(double)(nstate0-1);
    double lam;
    double mu;
    double v;
    double *qmat = REAL(QMAT);

    if (is_covarion) {
        lda = 2 * nstate0;
        if (npar == 2) {
            lam = REAL(pars)[0];    // "on" rate
            mu = REAL(pars)[0];     // "off" rate
            v = REAL(pars)[1];
        } else {
            lam = REAL(pars)[0];
            mu = REAL(pars)[1];
            v = REAL(pars)[2];
        }
    } else {
        lda = nstate0;
        v = REAL(pars)[0];
    }

    memset(qmat, 0, lda * lda * sizeof(double));

    for (i = 0; i < (nstate0-1); ++i) {
        if (is_covarion) {
            qmat[i + (i+nstate0) * lda] = lam;
            qmat[(i+nstate0) + i * lda] = mu;
            qmat[i + i * lda] -= lam;
            qmat[(i+nstate0) + (i+nstate0) * lda] -= mu;
        }
        for (j = (i+1); j < nstate0; ++j) {
            if (is_covarion) {
                qmat[(i+nstate0) + (j+nstate0) * lda] = v * f;
                qmat[(j+nstate0) + (i+nstate0) * lda] = v * f;
                qmat[(i+nstate0) + (i+nstate0) * lda] -= qmat[(i+nstate0) + (j+nstate0) * lda];
                qmat[(j+nstate0) + (j+nstate0) * lda] -= qmat[(j+nstate0) + (i+nstate0) * lda];
            } else {
                qmat[i + j * lda] = v * f;
                qmat[j + i * lda] = v * f;
                qmat[i + i * lda] -= qmat[i + j * lda];
                qmat[j + j * lda] -= qmat[j + i * lda];
            }
        }
    }

    if (is_covarion) {
        qmat[(nstate0-1) + (2*nstate0-1) * lda] = lam;
        qmat[(2*nstate0-1) + (nstate0-1) * lda] = mu;
        qmat[(nstate0-1) + (nstate0-1) * lda] -= lam;
        qmat[(2*nstate0-1) + (2*nstate0-1) * lda] -= mu;
    }

    return R_NilValue;
}



/*
** Dirichlet trinomial probability
**
** Evaluate the probability of 3-way count data in n trials
** given that counts are drawn from a trinomial distribution
** and p, the parameters of the trinomial distribution, are
** drawn from a particular region of dirichlet distribution.
**
** We evaluate this probability by factoring the dirichlet-trinomial into
** a product of a marginal and conditional beta-binomial distribution.
**
** n1,n2,n3 -- counts
** a,b,c    -- parameters of the dirichlet distribution
** lnorm    -- log normalization constant of dirichlet
** lcoef    -- log multinomial coefficent
** l1,u1    -- lower and upper bound over which the marginal binomial distribution
**             compounded with the beta distribution is evaluated
** l2,u2    -- lower and upper bound over which the conditional binomial distribution
**             compounded with the beta distribution is evaluated
*/
static double pdirichtri0(int n1, int n2, int n3,
    double l1, double u1, double l2, double u2, double lcoef,
    double lnorm, double a, double b, double c)
{
    int n = n1 + n2 + n3;
    double A = ibeta(u1, n1+a, n-n1+b+c) - ibeta(l1, n1+a, n-n1+b+c);
    double B = ibeta(u2, n2+b, n3+c) - ibeta(l2, n2+b, n3+c);
    return exp(lcoef + lnorm + log(A) + log(B));
}


/*
** npart -- number of partitions for each distribution
** state -- evolutionary state [0, npart*npart-1]
*/
static double pdirichtri(int state, int npart, int n1, int n2, int n3,
    double lcoef, double lnorm, double a, double b, double c)
{
    double l1, u1, l2, u2;
    double step = 1 / (double)npart;

    // unravel the state index to get partition index for each distribution
    int i = state % npart;
    int j = (state-i) / npart;

    l1 = i * step;
    u1 = l1 + step;
    l2 = j * step;
    u2 = l2 + step;

    return pdirichtri0(n1, n2, n3, l1, u1, l2, u2, lcoef, lnorm, a, b, c);
}


/*
** Compute the log normalization constant for a 3-variate
** dirichlet distribution with parameters a, b, and c
*/
static double dirichlet_lnorm(double a, double b, double c)
{
    return lgamma(a+b+c) - (lgamma(a) + lgamma(b) + lgamma(c));
}

/*
** Compute the amount of probability contained in a cell
** of a partition of a 3-variate Dirichlet distribution.
**
** The strategy is to factor the Dirichlet into a product
** of a marginal Beta distribution and a conditional Beta
** distribution and then partition this product space. The
** amount of probability a cell in such a partition contains
** is easily found using the incomplete beta function.
**
** npart   -- the number of equal-width bins over (0, 1) in each direction
** i,j     -- the cell index for the marginal and conditional distribution
** a, b, c -- the dirichlet parameters
** lnorm   -- the log normalization constant for the dirichlet distribution
*/
static double dirichlet_cell_prob(int npart, int i, int j, double a, double b, double c, double lnorm)
{
    double l1, u1, l2, u2;
    double step = 1 / (double)npart;

    l1 = i * step, u1 = l1 + step;
    l2 = j * step, u2 = l2 + step;

    double A = ibeta(u1, a, b+c) - ibeta(l1, a, b+c);
    double B = ibeta(u2, b, c) - ibeta(l2, b, c);

    return exp(lnorm + log(A) + log(B));
}


static double dmultinom(int state, int npart, int n1, int n2, int n3,
    double a, double b, double c)
{
    double l1, u1, l2, u2;
    double step = 1 / (double)npart;

    // unravel the state index to get partition index for each distribution
    int i = state % npart;
    int j = (state-i) / npart;

    l1 = i * step;
    u1 = l1 + step;
    l2 = j * step;
    u2 = l2 + step;

    // bavg(i, step, a, b+c);  // marginal beta
    // bavg(j, step, b, c);    // conditional beta

    return dbinom(n1, n1+n2+n3, bavg(i, step, a, b+c), 0) *
        dbinom(n2, n2+n3, bavg(j, step, b, c), 0);
}


SEXP rcl_mkp_dirichtri_qmat(SEXP nstate, SEXP pars, SEXP iscovarion, SEXP QMAT)
{
    int i, j, lda, npar = LENGTH(pars);
    int npart = INTEGER(nstate)[0];
    int nstate0 = npart * npart;
    int is_covarion = INTEGER(iscovarion)[0];
    double prob[nstate0];
    double a = REAL(pars)[0];
    double b = REAL(pars)[1];
    double c = REAL(pars)[2];
    double lnorm = dirichlet_lnorm(a, b, c);
    double lam, mu, v;
    double *qmat = REAL(QMAT);

    // Compute the amount of probability in each cell. This is
    // used for drawing the next state when a transition occurs
    for (i = 0; i < npart; ++i) {
        for (j = 0; j < npart; ++j)
            prob[i + j * npart] = dirichlet_cell_prob(npart, i, j, a, b, c, lnorm);
    }

    if (is_covarion) {
        lda = 2 * nstate0;
        if (npar == 4) {
            lam = REAL(pars)[3];    // "on" rate
            mu = REAL(pars)[4];     // "off" rate
            v = REAL(pars)[5];
        } else {
            lam = REAL(pars)[3];
            mu = REAL(pars)[3];
            v = REAL(pars)[4];
        }
    } else {
        lda = nstate0;
        v = REAL(pars)[3];
    }

    memset(qmat, 0, lda * lda * sizeof(double));

    for (i = 0; i < (nstate0-1); ++i) {
        if (is_covarion) {
            qmat[i + (i+nstate0) * lda] = lam;
            qmat[(i+nstate0) + i * lda] = mu;
            qmat[i + i * lda] -= lam;
            qmat[(i+nstate0) + (i+nstate0) * lda] -= mu;
        }
        for (j = (i+1); j < nstate0; ++j) {
            if (is_covarion) {
                qmat[(i+nstate0) + (j+nstate0) * lda] = v * (prob[j]) / (1-prob[i]);
                qmat[(j+nstate0) + (i+nstate0) * lda] = v * (prob[i]) / (1-prob[j]);
                qmat[(i+nstate0) + (i+nstate0) * lda] -= qmat[(i+nstate0) + (j+nstate0) * lda];
                qmat[(j+nstate0) + (j+nstate0) * lda] -= qmat[(j+nstate0) + (i+nstate0) * lda];
            } else {
                qmat[i + j * lda] = v * (prob[j]) / (1-prob[i]);
                qmat[j + i * lda] = v * (prob[i]) / (1-prob[j]);
                qmat[i + i * lda] -= qmat[i + j * lda];
                qmat[j + j * lda] -= qmat[j + i * lda];
            }
        }
    }

    if (is_covarion) {
        qmat[(nstate0-1) + (2*nstate0-1) * lda] = lam;
        qmat[(2*nstate0-1) + (nstate0-1) * lda] = mu;
        qmat[(nstate0-1) + (nstate0-1) * lda] -= lam;
        qmat[(2*nstate0-1) + (2*nstate0-1) * lda] -= mu;
    }

    return R_NilValue;
}


/*
** pars -- shape parameters for the dirichlet distribution
** data -- Data for the tips. 3 by NTIP matrix. Each row contains the
**         number of counts in each resource category.
** clk -- Downpass conditional likelihoods. NSTATES by NNODE matrix.
*/
SEXP rcl_mkp_dirichtri_init(SEXP rtree, SEXP NPART, SEXP NCAT, SEXP CATEG, SEXP pars, SEXP data, SEXP clk, SEXP iscovarion)
{
    int i, j, n, n1, n2,
        ncat = INTEGER(NCAT)[0],
        npart = INTEGER(NPART)[0],              // number of bins for each beta distribution
        nstate0 = npart * npart,
        is_covarion = INTEGER(iscovarion)[0],
        categ1 = INTEGER(CATEG)[0],
        categ2 = INTEGER(CATEG)[1],
        *tipdata = INTEGER(data),
        nstate;

    double *dpars = REAL(pars),
           *dclk = REAL(clk), lcoef,
           a = dpars[0], b = dpars[1], c = dpars[2],
           lnorm = dirichlet_lnorm(a, b, c), prob[nstate0];

    struct tree *tree = (struct tree*)R_ExternalPtrAddr(rtree);

    for (i = 0; i < npart; ++i) {
        for (j = 0; j < npart; ++j)
            prob[i + j * npart] = dirichlet_cell_prob(npart, i, j, a, b, c, lnorm);
    }

    // if this is a covarion model the conditional likelihood matrix
    // actually has twice the number of rows as the "true" number of states
    if (is_covarion)
        nstate = 2 * nstate0;
    else
        nstate = nstate0;

    #undef ELEM1
    #undef COLUMN1
    #undef ELEM2
    // access into the conditional likelihood array
    #define ELEM1(m, i, j) m[(i) + (j)*nstate]
    #define COLUMN1(m, j) m + (j)*nstate

    // access into the tipdata array
    #define ELEM2(m, i, j) m[(i) + (j)*(ncat+1)]

    for (i = 0; i < tree->ntip; ++i) {

        n = ELEM2(tipdata, 0, i);
        n1 = ELEM2(tipdata, categ1, i);
        n2 = ELEM2(tipdata, categ2, i);
        // log multinomial coefficient
        lcoef = lgamma(n+1) - (lgamma(n1+1) + lgamma(n2+1) + lgamma(n-n1-n2+1));

        for (j = 0; j < nstate0; ++j) {
            ELEM1(dclk, j, i) = pdirichtri(j, npart, n1, n2, n-n1-n2, lcoef, lnorm, a, b, c) / prob[j];
            //ELEM1(dclk, j, i) = dmultinom(j, npart, n1, n2, n-n1-n2, a, b, c);
        }

        if (is_covarion)
            memcpy(COLUMN1(dclk, i) + nstate0, COLUMN1(dclk, i), nstate0 * sizeof(double));

    }

    return R_NilValue;
}
