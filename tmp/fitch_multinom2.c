#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include "tree.h"


struct state {
    /* number of tips in state */
    int n;

    /* indices of tips in state */
    int *id;

    /* amalgamated count data from tips in state */
    int *x;

    /* sum of logs of multinomial coefficients for
    ** each tip in state */
    double lmultinomcoef;

    /* log likelihood of data in state */
    double loglk;

    /* pointer to next state in the set */
    struct state *next;
};


static struct state *alloc_state(int ncat)
{
    struct state *state;
    state = Calloc(1, struct state);
    state->n = 0;
    state->id = NULL;
    state->x = Calloc(ncat, int);
    state->next = NULL;
    return state;
}


static void free_state(struct state *state)
{
    if (state->id)
        Free(state->id);
    Free(state->x);
    Free(state);
}


/* Compute the logarithm of the multinomial coefficient
** for the count data in x */
static double lmultinomfn(int size, int *x)
{
    int i;
    int n;
    double lcoef;

    n = 0;
    lcoef = 0;

    for (i = 0; i < size; ++i) {
        n += x[i];
        lcoef -= lgammafn(x[i] + 1);
    }

    lcoef += lgammafn(n+1);

    return lcoef;
}


static double merged_loglk(int size, struct state *lfstate, struct state *rtstate)
{
    int i;
    int x;
    double loglk;

    x = 0;
    loglk = lfstate->lmultinomcoef + rtstate->lmultinomcoef;

    for (i = 0; i < size; ++i) {
        x += lfstate->x[i] + rtstate->x[i];
        loglk += lgammafn(lfstate->x[i] + rtstate->x[i] + 0.5);
    }

    loglk += lgammafn(0.5 * size) - lgammafn(0.5 * size + x) - size * lgammafn(0.5);

    return loglk;
}


static double state_loglk(int size, struct state *state)
{
    int i;
    int x;
    double loglk;

    x = 0;
    loglk = state->lmultinomcoef;

    for (i = 0; i < size; ++i) {
        x += state->x[i];
        loglk += lgammafn(state->x[i] + 0.5);
    }

    loglk += lgammafn(0.5 * size) - lgammafn(0.5 * size + x) - size * lgammafn(0.5);

    return loglk;
}


static void run_downpass(int ncat, int *data, struct tree *tree, int *nevent, double *lnl, int *uppass)
{
    /* Number of tip groups in each node's subtree */
    int ngroup[tree->nnode];

    /* Log likelihood of data in each node's subtree */
    double loglk[tree->nnode];

    int i;
    int j;
    int n = 0;
    double _loglk;
    double aic1, aic2, w1, w2;
    struct tree_traversal t;
    struct node *node;
    struct node *lfdesc;
    struct node *rtdesc;
    struct state *lfstate;
    struct state *rtstate;
    struct state **nodestate;

    /* bitmask for the uppass stage. during the uppass,
    ** if the indicator variable is 1 that node inherits
    ** the state of its ancestor */
    int merge[tree->nnode];

    memset(merge, 0, tree->nnode * sizeof(int));

    nodestate = Calloc(tree->nnode, struct state*);

    for (i = 0; i < tree->ntip; ++i) {
        nodestate[i] = alloc_state(ncat);
        memcpy(nodestate[i]->x, data + i * ncat, ncat * sizeof(int));
        ngroup[i] = 1;
        nodestate[i]->n = 1;
        nodestate[i]->lmultinomcoef = lmultinomfn(ncat, nodestate[i]->x);
        nodestate[i]->loglk = state_loglk(ncat, nodestate[i]);
        loglk[i] = nodestate[i]->loglk;
    }

    t = tree_traverse(POSTORDER, INTERNAL_NODES_ONLY, tree->root, tree);

    node = tree_step(&t);

    while (node) {
        nodestate[node->index] = alloc_state(ncat);

        lfdesc = node->children[0];
        rtdesc = node->children[1];
        lfstate = nodestate[lfdesc->index];
        rtstate = nodestate[rtdesc->index];

        _loglk = loglk[lfdesc->index] + loglk[rtdesc->index];
        aic2 = -2 * _loglk + 2 * (ngroup[lfdesc->index] + ngroup[rtdesc->index]);
        _loglk -= lfstate->loglk;
        _loglk -= rtstate->loglk;
        _loglk += merged_loglk(ncat, lfstate, rtstate);
        aic1 = -2 * _loglk + 2 * (ngroup[lfdesc->index] + ngroup[rtdesc->index] - 1);
        w1 = exp(-0.5 * (aic1 - fmin(aic1, aic2)));
        w2 = exp(-0.5 * (aic2 - fmin(aic1, aic2)));

        if (unif_rand() < (w1 / (w1 + w2))) {
            ngroup[node->index] = ngroup[lfdesc->index] + ngroup[rtdesc->index] - 1;
            nodestate[node->index]->n = lfstate->n + rtstate->n;
            nodestate[node->index]->lmultinomcoef = lfstate->lmultinomcoef + rtstate->lmultinomcoef;
            for (i = 0; i < ncat; ++i)
                nodestate[node->index]->x[i] = lfstate->x[i] + rtstate->x[i];
            nodestate[node->index]->loglk = state_loglk(ncat, nodestate[node->index]);
            loglk[node->index] = _loglk;
            merge[lfdesc->index] = 1;
            merge[rtdesc->index] = 1;
        } else {
            ++n;
            ngroup[node->index] = ngroup[lfdesc->index] + ngroup[rtdesc->index];
            loglk[node->index] = loglk[lfdesc->index] + loglk[rtdesc->index];
            if (unif_rand() < 0.5) {
                nodestate[node->index]->n = lfstate->n;
                nodestate[node->index]->lmultinomcoef = lfstate->lmultinomcoef;
                memcpy(nodestate[node->index]->x, lfstate->x, ncat * sizeof(int));
                nodestate[node->index]->loglk = lfstate->loglk;
            } else {
                nodestate[node->index]->n = rtstate->n;
                nodestate[node->index]->lmultinomcoef = rtstate->lmultinomcoef;
                memcpy(nodestate[node->index]->x, rtstate->x, ncat * sizeof(int));
                nodestate[node->index]->loglk = rtstate->loglk;
            }
        }

        node = tree_step(&t);
    }

    /* begin uppass */
    t = tree_traverse(PREORDER, ALL_NODES, tree->root, tree);
    node = tree_step(&t);
    node = tree_step(&t);
    memcpy(uppass + ncat * tree->root->index, nodestate[tree->root->index]->x, ncat * sizeof(int));
    while (node) {
        if (merge[node->index]) {
            memcpy(nodestate[node->index]->x, nodestate[node->parent->index]->x, ncat * sizeof(int));
            memcpy(uppass + ncat * node->index, nodestate[node->parent->index]->x, ncat * sizeof(int));
        } else {
            memcpy(uppass + ncat * node->index, nodestate[node->index]->x, ncat * sizeof(int));
        }
        node = tree_step(&t);
    }

    /* end uppass */

    for (i = 0; i < tree->nnode; ++i)
        free_state(nodestate[i]);
    Free(nodestate);

    *lnl = loglk[tree->root->index];
    *nevent = n;
}


SEXP rcl_fitch_multinom2(SEXP rtree, SEXP data)
{
    int ncat;
    struct tree *tree;

    int *nevent;
    int *nodestate;
    double *loglk;

    SEXP dim;
    SEXP result;
    result = PROTECT(allocVector(VECSXP, 3));
    dim = PROTECT(allocVector(INTSXP, 2));
    ncat = INTEGER(getAttrib(data, R_DimSymbol))[0];
    tree = (struct tree*)R_ExternalPtrAddr(rtree);
    SET_VECTOR_ELT(result, 0, allocVector(REALSXP, 1));
    SET_VECTOR_ELT(result, 1, allocVector(INTSXP, 1));
    SET_VECTOR_ELT(result, 2, allocVector(INTSXP, ncat * tree->nnode));

    INTEGER(dim)[0] = ncat;
    INTEGER(dim)[1] = tree->nnode;

    loglk = REAL(VECTOR_ELT(result, 0));
    nevent = INTEGER(VECTOR_ELT(result, 1));
    nodestate = INTEGER(VECTOR_ELT(result, 2));

    GetRNGstate();
    run_downpass(ncat, INTEGER(data), tree, nevent, loglk, nodestate);
    PutRNGstate();

    setAttrib(VECTOR_ELT(result, 2), R_DimSymbol, dim);

    UNPROTECT(2);
    return result;
}
