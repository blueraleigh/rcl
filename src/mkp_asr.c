#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include "node.h"
#include "tree.h"
#include "matrix.h"

/*
** Cache for marginal ancestral state reconstruction.
**
** Each cache is simply an array of doubles in a specified
** order that allows retrieval of the cached item by node
** index.
*/

#define COLUMN(m, j) m + nstate * (j)
#define ELEM(m, i, j) m[(i) + nstate * (j)]

struct mkp_cache {

    /* Cache of transition probability matrices for
    ** each node's branch.
    */
    double *prob;

    /* Cache of the downpass conditional likelihoods
    ** for each node. These give the likelihood of the
    ** data contained in a node's subtree conditional
    ** on the node being in a particular state.
    */
    double *dclk;

    /* Conditional likelihoods are scaled to avoid underflow
    ** errors and we keep track of the logarithms of these
    ** scaling factors in a separate cache.
    */
    double *dclk_lq;

    /* Sum of dclk_lq */
    double ld;

    /* Log likelihood */
    double loglk;
};


static struct mkp_cache *mkp_cache_alloc(int nstate, struct tree *tree)
{
    int n = nstate * nstate * tree->nnode;
    int m = nstate * tree->nnode;
    struct mkp_cache *cache = Calloc(1, struct mkp_cache);
    cache->prob = Calloc(n, double);
    cache->dclk = Calloc(m, double);
    cache->dclk_lq = Calloc(tree->nnode, double);
    cache->ld = 0;
    cache->loglk = 0;
    return cache;
}


static void mkp_cache_free(struct mkp_cache *cache)
{
    Free(cache->prob);
    Free(cache->dclk);
    Free(cache->dclk_lq);
    Free(cache);
}


static void mkp_cache_init(int nstate, double *clk, struct mkp_cache *cache, struct tree *tree)
{
    memcpy(cache->dclk, clk, nstate * tree->nnode * sizeof(double));
}


/*
** Normalize conditional log likelihoods to state probability
** distribution.
*/
static void mkp_normalize(int nstate, double *clk)
{
    int j;
    double maxp = R_NegInf, tot = 0;

    for (j = 0; j < nstate; ++j) {
        if (clk[j] > maxp)
            maxp = clk[j];
    }

    for (j = 0; j < nstate; ++j) {
        if (ISNAN(clk[j]))
            clk[j] = 0.0;
        else
            tot += clk[j] = exp(clk[j] - maxp);
    }

    for (j = 0; j < nstate; ++j)
        clk[j] /= tot;
}


static double mkp_branch_prob(int nstate, double *R, struct mkp_cache *cache, struct node *node, double *out)
{
    int ione = 1;
    double sf, isf;

    double *Pr = cache->prob + nstate * nstate * node->index;

    double *init = COLUMN(cache->dclk, node->index);

    // compute transition probability matrix
    matrix_exponential(R, node->brlen, nstate, Pr);

    // multiply that matrix with the initial conditional likelihoods
    // in the init vector, storing the new conditional likelihoods
    // in the out vector
    matrix_vector_product(Pr, nstate, nstate, init, out);

    // compute the scale factor, which is just the sum of the new
    // conditional likelihoods
    sf = F77_CALL(dasum)(&nstate, out, &ione);

    // rescale the conditional likelihoods by this factor
    isf = 1 / sf;
    F77_CALL(dscal)(&nstate, &isf, out, &ione);

    return sf;
}


static void mkp_node_lk(int nstate, double *R, struct mkp_cache *cache, struct node *node)
{
    int j;
    double sf, lsf, rsf;
    double lf[nstate], rt[nstate];

    // do the left branch
    lsf = mkp_branch_prob(nstate, R, cache, node->children[0], lf);

    // do the right branch
    rsf = mkp_branch_prob(nstate, R, cache, node->children[1], rt);

    // to form the conditional likelihoods for the ancestral node
    // we multiply the contributions from each daughter
    for (j = 0; j < nstate; ++j)
        ELEM(cache->dclk, j, node->index) = lf[j] * rt[j];

    sf = log(lsf) + log(rsf);
    cache->dclk_lq[node->index] = sf;
    cache->ld += sf;
}


/* compute the likelihood and populate the cache in the process */
static void mkp_loglk(int nstate, double *R, struct mkp_cache *cache, struct tree *tree)
{
    int j, ione = 1;
    double p = 0, tot = 0;
    struct tree_traversal t;
    struct node *node;

    t = tree_traverse(POSTORDER, INTERNAL_NODES_ONLY, tree->root, tree);
    node = tree_step(&t);

    while (node != NULL) {
        mkp_node_lk(nstate, R, cache, node);
        node = tree_step(&t);
    }

    // sum the (scaled) conditional likelihoods at the root of the tree, weighting
    // each term by its proportional contribution to the total sum
    tot = F77_CALL(dasum)(&nstate, COLUMN(cache->dclk, tree->root->index), &ione);

    for (j = 0; j < nstate; ++j)
        p += (ELEM(cache->dclk, j, tree->root->index) / tot) * ELEM(cache->dclk, j, tree->root->index);

    cache->loglk = log(p);

    // finally, add in the log of the total scale factor to get the unscaled log likelihood
    cache->loglk += cache->ld;
}



static double mkp_branch_prob_from_cache(int nstate, struct mkp_cache *cache, struct node *node, double *out)
{
    int ione = 1;
    double sf, isf;

    // get the pre-computed transition probability matrix
    double *Pr = cache->prob + nstate * nstate * node->index;

    double *init = COLUMN(cache->dclk, node->index);

    // multiply that matrix with the initial conditional likelihoods
    // in the init vector, storing the new conditional likelihoods
    // in the out vector
    matrix_vector_product(Pr, nstate, nstate, init, out);

    // compute the scale factor, which is just the sum of the new
    // conditional likelihoods
    sf = F77_CALL(dasum)(&nstate, out, &ione);

    // rescale the conditional likelihoods by this factor
    isf = 1 / sf;
    F77_CALL(dscal)(&nstate, &isf, out, &ione);

    return sf;
}


static double mkp_node_lk_from_cache(int nstate, struct mkp_cache *cache, struct node *node)
{
    int j;
    double lsf, rsf;
    double lf[nstate], rt[nstate];

    // do the left branch
    lsf = mkp_branch_prob_from_cache(nstate, cache, node->children[0], lf);

    // do the right branch
    rsf = mkp_branch_prob_from_cache(nstate, cache, node->children[1], rt);

    // to form the conditional likelihoods for the ancestral node
    // we multiply the contributions from each daughter
    for (j = 0; j < nstate; ++j)
        ELEM(cache->dclk, j, node->index) = lf[j] * rt[j];

    return log(lsf) + log(rsf);
}



static double mkp_asr_compute(
    int nstate,
    struct mkp_cache *cache,
    struct node *node,
    struct tree *tree,
    int state)
{
    int j;
    int ione = 1;
    double v;
    double p = 0;
    double tot = 0;
    double loglk = 0;
    double oldsum = 0;
    double newsum = 0;

    // to compute the marginal probability that a node is in a given state
    // we set all conditional likelihoods to zero except the conditional for
    // the state of interest and then recalculate the likehood for the tree
    v = ELEM(cache->dclk, state, node->index);
    memset(COLUMN(cache->dclk, node->index), 0, nstate * sizeof(double));
    ELEM(cache->dclk, state, node->index) = v;

    // we can do this by recalculating the conditional likelihoods
    // for all nodes on the path back to root from current node
    while (node != tree->root) {
        node = node->parent;
        oldsum += cache->dclk_lq[node->index];
        newsum += mkp_node_lk_from_cache(nstate, cache, node);
    }

    tot = F77_CALL(dasum)(&nstate, COLUMN(cache->dclk, tree->root->index), &ione);

    for (j = 0; j < nstate; ++j)
        p += (ELEM(cache->dclk, j, tree->root->index) / tot) * ELEM(cache->dclk, j, tree->root->index);

    loglk = log(p);

    return loglk + (cache->ld - oldsum + newsum);
}


static void mkp_asr(
    int nstate,
    int exclude_tips,
    double *asr,
    struct mkp_cache *cache,
    struct tree *tree)
{
    int j;
    double cpy_dclk[nstate * tree->nnode];
    memcpy(cpy_dclk, cache->dclk, nstate * tree->nnode * sizeof(double));

    struct tree_traversal t;
    struct node *node;

    if (exclude_tips)
        t = tree_traverse(POSTORDER, INTERNAL_NODES_ONLY, tree->root, tree);
    else
        t = tree_traverse(POSTORDER, ALL_NODES, tree->root, tree);
    node = tree_step(&t);

    while (node != NULL) {
        for (j = 0; j < nstate; ++j) {
            ELEM(asr, j, node->index) = mkp_asr_compute(
                nstate,
                cache,
                node,
                tree,
                j
            );
            // copy the original conditional likelihoods back in as these get
            // modified by mkp_asr_compute
            memcpy(cache->dclk, cpy_dclk, nstate * tree->nnode * sizeof(double));
        }
        mkp_normalize(nstate, COLUMN(asr, node->index));
        node = tree_step(&t);
    }
}


SEXP rcl_mkp_asr(SEXP rtree, SEXP R, SEXP clk, SEXP exclude_tips)
{
    struct tree *tree = (struct tree*)R_ExternalPtrAddr(rtree);
    int nstate = INTEGER(getAttrib(R, R_DimSymbol))[0];
    struct mkp_cache *cache = mkp_cache_alloc(nstate, tree);
    SEXP asr = PROTECT(allocVector(REALSXP, nstate * tree->nnode));
    SEXP dim = PROTECT(allocVector(INTSXP, 2));
    INTEGER(dim)[0] = nstate;
    INTEGER(dim)[1] = tree->nnode;
    setAttrib(asr, R_DimSymbol, dim);
    mkp_cache_init(nstate, REAL(clk), cache, tree);
    mkp_loglk(nstate, REAL(R), cache, tree);
    mkp_asr(nstate, INTEGER(exclude_tips)[0], REAL(asr), cache, tree);
    mkp_cache_free(cache);
    UNPROTECT(2);
    return asr;
}
