#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include "node.h"
#include "tree.h"
#include "matrix.h"
#include "atblspc.h"

/*
** Cache for stochastic character mapping.
**
** Each cache is simply an array of doubles in a specified
** order that allows retrieval of the cached item by node
** index.
*/

#define COLUMN(m, j) m + nstate * (j)
#define ELEM(m, i, j) m[(i) + nstate * (j)]


// nstate by nstate matrices for O(1) weighted random
// sampling using the alias method. each entry gives
// the probability of sampling a state for a descendant node
// given an ancestral node state corresponding to the column
// index
struct nodeprob {
    int *alias;
    double *prob;
};


struct rootprob {
    int *alias;
    double *prob;
};


struct stateprob {
    int *alias;
    double *prob;
};


struct mkp_cache {

    /* CTMC rate matrix */
    double *R;

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

    /* Conditional probability distribution on node states
    ** given a sampled ancestral node state and downpass
    ** likelihoods (i.e. data)
    */
    struct nodeprob *nodep;

    /* Conditional probability distribution on node states
    ** given a sampled ancestral node state. These are only
    ** made conditional on the fitted transition probability
    ** matrix, not the observed data at the tips of the tree,
    ** and are referred to as "priors" for that reason.
    */
    struct nodeprob *nodep_prior;

    /* Probability distribution of states at the root
    ** implied by the conditional likelihoods.
    */
    struct rootprob rootp;

    /* Probability distribution on states given that
    ** a transition out of a state has just occurred
    */
    struct stateprob statep;

    /* Sum of dclk_lq */
    double ld;

    /* Log likelihood */
    double loglk;

    int nstate;
    int nnode;
};


static struct mkp_cache *mkp_cache_alloc(int nstate, struct tree *tree)
{
    int i;
    int n = nstate * nstate;
    int m = nstate * tree->nnode;
    struct mkp_cache *cache = Calloc(1, struct mkp_cache);
    cache->R = Calloc(n, double);
    cache->dclk = Calloc(m, double);
    cache->dclk_lq = Calloc(tree->nnode, double);
    cache->nodep = Calloc(tree->nnode, struct nodeprob);
    cache->nodep_prior = Calloc(tree->nnode, struct nodeprob);
    for (i = 0; i < tree->nnode; ++i) {
        cache->nodep[i].alias = Calloc(n, int);
        cache->nodep[i].prob = Calloc(n, double);
        cache->nodep_prior[i].alias = Calloc(n, int);
        cache->nodep_prior[i].prob = Calloc(n, double);
    }
    cache->rootp.alias = Calloc(nstate, int);
    cache->rootp.prob = Calloc(nstate, double);
    cache->statep.alias = Calloc(n, int);
    cache->statep.prob = Calloc(n, double);
    cache->ld = 0;
    cache->loglk = 0;
    cache->nstate = nstate;
    cache->nnode = tree->nnode;
    return cache;
}


static void mkp_cache_free(struct mkp_cache *cache)
{
    int i;
    Free(cache->R);
    Free(cache->dclk);
    Free(cache->dclk_lq);
    for (i = 0; i < cache->nnode; ++i) {
        Free(cache->nodep[i].alias);
        Free(cache->nodep[i].prob);
        Free(cache->nodep_prior[i].alias);
        Free(cache->nodep_prior[i].prob);
    }
    Free(cache->nodep);
    Free(cache->nodep_prior);
    Free(cache->rootp.alias);
    Free(cache->rootp.prob);
    Free(cache->statep.alias);
    Free(cache->statep.prob);
    Free(cache);
}


static void mkp_cache_init(int nstate, double *R, double *clk, struct mkp_cache *cache, struct tree *tree)
{
    int i, j;
    double tmp[nstate];
    memcpy(cache->dclk, clk, nstate * tree->nnode * sizeof(double));
    memcpy(cache->R, R, nstate * nstate * sizeof(double));;
    for (i = 0; i < nstate; ++i) {
        for (j = 0; j < nstate; ++j)
            tmp[j] = ELEM(R, i, j);
        tmp[i] = 0;
        atblspc_set(nstate, tmp,
            COLUMN(cache->statep.alias, i), COLUMN(cache->statep.prob, i));
    }
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
    int i, j, ione = 1;
    double sf, isf;
    double tmp[nstate];
    double Pr[nstate * nstate];

    double *init = COLUMN(cache->dclk, node->index);

    // compute transition probability matrix
    matrix_exponential(R, node->brlen, nstate, Pr);

    // multiply that matrix with the initial conditional likelihoods
    // in the init vector, storing the new conditional likelihoods
    // in the out vector
    matrix_vector_product(Pr, nstate, nstate, init, out);

    // now compute the probability of sampling a state at this node
    // conditional on (data and) a sampled ancestral state
    for (i = 0; i < nstate; ++i) {
        for (j = 0; j < nstate; ++j)
            tmp[j] = (ELEM(Pr, i, j) * init[j]) / out[i];
        // compute alias structure
        atblspc_set(nstate, tmp,
            COLUMN(cache->nodep[node->index].alias, i), COLUMN(cache->nodep[node->index].prob, i));
        for (j = 0; j < nstate; ++j)
            tmp[j] = ELEM(Pr, i, j);
        atblspc_set(nstate, tmp,
            COLUMN(cache->nodep_prior[node->index].alias, i), COLUMN(cache->nodep_prior[node->index].prob, i));
    }

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
    double p = 0, tot = 0, rootp[nstate];
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
    for (j = 0; j < nstate; ++j) {
        p += (ELEM(cache->dclk, j, tree->root->index) / tot) * ELEM(cache->dclk, j, tree->root->index);
        rootp[j] = log(ELEM(cache->dclk, j, tree->root->index)) + cache->ld;
    }
    cache->loglk = log(p);

    mkp_normalize(nstate, rootp);
    atblspc_set(nstate, rootp, cache->rootp.alias, cache->rootp.prob);

    // finally, add in the log of the total scale factor to get the unscaled log likelihood
    cache->loglk += cache->ld;
}


static void mkp_sampleroot(int nstate, int *state, struct mkp_cache *cache, struct node *root)
{
    int j = (int)(unif_rand() * nstate);
    j = unif_rand() < cache->rootp.prob[j] ? j : cache->rootp.alias[j];
    state[root->index] = j;
}


static void mkp_samplenode(int nstate, int *state, struct mkp_cache *cache, struct node *node)
{
    int i = state[node->parent->index];
    int j = (int)(unif_rand() * nstate);

    int *alias = COLUMN(cache->nodep[node->index].alias, i);
    double *prob = COLUMN(cache->nodep[node->index].prob, i);

    j = unif_rand() < prob[j] ? j : alias[j];

    state[node->index] = j;
}


static void mkp_samplenode_prior(int nstate, int *state, struct mkp_cache *cache, struct node *node)
{
    int i = state[node->parent->index];
    int j = (int)(unif_rand() * nstate);

    int *alias = COLUMN(cache->nodep_prior[node->index].alias, i);
    double *prob = COLUMN(cache->nodep_prior[node->index].prob, i);

    j = unif_rand() < prob[j] ? j : alias[j];

    state[node->index] = j;
}


static int mkp_samplestate(int nstate, int state, struct mkp_cache *cache)
{
    int j = (int)(unif_rand() * nstate);
    int *alias = COLUMN(cache->statep.alias, state);  // state is the parent state
    double *prob = COLUMN(cache->statep.prob, state);
    return unif_rand() < prob[j] ? j : alias[j];
}


static int mkp_draw_history(
    int nstate,             // nbr of character states
    int begin,
    int end,
    double brlen,
    int *count,             // count[i, j] whether or not i->j transitions are being counted
    struct mkp_cache *cache)
{
    int tmp, begin0, incr, redo = 1;
    double t, rate;

    if (begin != end) {
        // draw history conditional on at least one change
        while (redo) {
            incr = 0;
            begin0 = begin;
            rate = -ELEM(cache->R, begin0, begin0);
            // draw first change such that it is guaranteed to occur within
            // time available (see Nielsen 2001)
            t = -log(1 - unif_rand() * (1-exp(-rate * brlen))) / rate;
            while (t < brlen) {
                tmp = mkp_samplestate(nstate, begin0, cache); // choose state;
                if (ELEM(count, begin0, tmp))
                    incr += 1;
                begin0 = tmp;
                rate = -ELEM(cache->R, begin0, begin0);
                t += rexp(1/rate); // rexp is parameterized by scale (= 1/rate)
            }
            redo = (begin0 == end) ? 0 : 1;
        }
    } else {
        while (redo) {
            incr = 0;
            begin0 = begin;
            rate = -ELEM(cache->R, begin0, begin0);
            t = rexp(1/rate);
            while (t < brlen) {
                tmp = mkp_samplestate(nstate, begin0, cache);
                if (ELEM(count, begin0, tmp))
                    incr += 1;
                begin0 = tmp;
                rate = -ELEM(cache->R, begin0, begin0);
                t += rexp(1/rate);
            }
            redo = (begin0 == end) ? 0 : 1;
        }
    }
    return incr;
}


/*
** Apply mapping to single branch. This assumes that states have already been
** sampled for nodes
*/
static int mkp_map_branch(
    int nstate,
    int *state,             // state[node->index] returns the character state of node
    int *count,             // count[i, j] returns whether or not i->j transitions are counted
    struct mkp_cache *cache,
    struct node *node)      // branch that the map is being applied to
{
    // assert(node->parent);

    return mkp_draw_history(nstate, state[node->parent->index],
        state[node->index], node->brlen, count, cache);
}


static int mkp_map_branchset(
    int nstate,
    int nbranch,
    int useprior,
    int *branches,      /* array of node indices to map */
    int *count,
    struct mkp_cache *cache,
    struct tree *tree)
{
    int i, incr = 0, state[tree->nnode];
    struct tree_traversal t;
    struct node *node;

    t = tree_traverse(PREORDER, ALL_NODES, tree->root, tree);
    node = tree_step(&t);

    // sample state for root
    mkp_sampleroot(nstate, state, cache, tree->root);

    node = tree_step(&t);

    if (!useprior) {
        // sample state for each node conditional on parent state
        while (node) {
            mkp_samplenode(nstate, state, cache, node);
            node = tree_step(&t);
        }
    } else {
        while (node) {
            mkp_samplenode_prior(nstate, state, cache, node);
            node = tree_step(&t);
        }
    }

    for (i = 0; i < nbranch; ++i)
        incr += mkp_map_branch(nstate, state, count, cache, tree->node[branches[i]]);
    return incr;
}


static int mkp_map_clade(
    int nstate,
    int useprior,
    int *count,             // count[i, j] returns whether or not i->j transitions are counted
    struct mkp_cache *cache,
    struct node *root,      // clade that the map is being applied to
    struct tree *tree)
{
    int incr = 0, state[tree->nnode];
    struct tree_traversal t;
    struct node *node;

    t = tree_traverse(PREORDER, ALL_NODES, tree->root, tree);
    node = tree_step(&t);

    // sample state for root
    mkp_sampleroot(nstate, state, cache, tree->root);

    node = tree_step(&t);

    if (!useprior) {
        // sample state for each node conditional on parent state
        while (node) {
            mkp_samplenode(nstate, state, cache, node);
            node = tree_step(&t);
        }
    } else {
        while (node) {
            mkp_samplenode_prior(nstate, state, cache, node);
            node = tree_step(&t);
        }
    }

    // given the sampled node states generate a character transformation history
    // along each branch

    t = tree_traverse(PREORDER, INTERNAL_NODES_ONLY, root, tree);
    node = tree_step(&t);

    while (node) {

        // do left branch
        incr += mkp_draw_history(nstate, state[node->index],
            state[node->children[0]->index], node->children[0]->brlen, count, cache);

        // do right branch
        incr += mkp_draw_history(nstate, state[node->index],
            state[node->children[1]->index], node->children[1]->brlen, count, cache);

        node = tree_step(&t);
    }
    return incr;
}


void rcl_mkp_cache_free(SEXP c)
{
    struct mkp_cache *cache = (struct mkp_cache*)R_ExternalPtrAddr(c);
    mkp_cache_free(cache);
    R_ClearExternalPtr(c);
}


SEXP rcl_mkp_cache_build(SEXP rtree, SEXP R, SEXP clk)
{
    struct tree *tree = (struct tree*)R_ExternalPtrAddr(rtree);
    int nstate = INTEGER(getAttrib(R, R_DimSymbol))[0];
    struct mkp_cache *cache = mkp_cache_alloc(nstate, tree);
    mkp_cache_init(nstate, REAL(R), REAL(clk), cache, tree);
    mkp_loglk(nstate, REAL(R), cache, tree);
    SEXP c = PROTECT(R_MakeExternalPtr(cache, R_NilValue, R_NilValue));
    R_RegisterCFinalizer(c, &rcl_mkp_cache_free);
    UNPROTECT(1);
    return c;
}


SEXP rcl_mkp_smap_clade(SEXP rtree, SEXP clade, SEXP count, SEXP rcache, SEXP nsim, SEXP fromprior)
{
    struct tree *tree = (struct tree*)R_ExternalPtrAddr(rtree);
    struct node *root = tree->node[INTEGER(clade)[0]];
    struct mkp_cache *cache = (struct mkp_cache*)R_ExternalPtrAddr(rcache);
    int i, nstate = cache->nstate, n = INTEGER(nsim)[0], useprior = INTEGER(fromprior)[0];
    int *cp, *vp;
    SEXP v;
    v = PROTECT(allocVector(INTSXP, n));
    vp = INTEGER(v);
    cp = INTEGER(count);
    GetRNGstate();
    for (i = 0; i < n; ++i)
        vp[i] = mkp_map_clade(nstate, useprior, cp, cache, root, tree);
    PutRNGstate();
    UNPROTECT(1);
    return v;
}


SEXP rcl_mkp_smap_branchset(SEXP rtree, SEXP branches, SEXP count, SEXP rcache, SEXP nsim, SEXP fromprior)
{
    struct tree *tree = (struct tree*)R_ExternalPtrAddr(rtree);
    struct mkp_cache *cache = (struct mkp_cache*)R_ExternalPtrAddr(rcache);
    int i, nstate = cache->nstate, useprior = INTEGER(fromprior)[0];
    SEXP v = PROTECT(allocVector(INTSXP, INTEGER(nsim)[0]));
    GetRNGstate();
    for (i = 0; i < INTEGER(nsim)[0]; ++i)
        INTEGER(v)[i] = mkp_map_branchset(nstate, LENGTH(branches), useprior,
            INTEGER(branches), INTEGER(count), cache, tree);
    PutRNGstate();
    UNPROTECT(1);
    return v;
}
