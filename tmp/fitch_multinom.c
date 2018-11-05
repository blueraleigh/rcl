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
        loglk += lgammafn(lfstate->x[i] + rtstate->x[i] + 1);
    }

    loglk += lgammafn(1 * size) - lgammafn(1 * size + x) - size * lgammafn(1);

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
        loglk += lgammafn(state->x[i] + 1);
    }

    loglk += lgammafn(1 * size) - lgammafn(1 * size + x) - size * lgammafn(1);

    return loglk;
}


static double run_downpass(int ncat, int *data, struct tree *tree)
{
    /* Number of tip groups in each node's subtree */
    int ngroup[tree->nnode];

    /* Log likelihood of data in each node's subtree */
    double loglk[tree->nnode];

    int i;
    int j;
    int n = 0;
    double _loglk;
    double aic;
    double aicmin;
    struct tree_traversal t;
    struct node *node;
    struct node *lfdesc;
    struct node *rtdesc;
    struct state *lfstate;
    struct state *rtstate;
    struct state *k = NULL;
    struct state *p = NULL;

    int needsfree[tree->nnode];
    memset(needsfree, 0, tree->nnode * sizeof(int));

    struct state **nodestate;

    nodestate = Calloc(tree->nnode, struct state*);

    for (i = 0; i < tree->ntip; ++i) {
        nodestate[i] = alloc_state(ncat);
        memcpy(nodestate[i]->x, data + i * ncat, ncat * sizeof(int));
        ngroup[i] = 1;
        nodestate[i]->n = 1;
        nodestate[i]->lmultinomcoef = lmultinomfn(ncat, nodestate[i]->x);
        nodestate[i]->loglk = state_loglk(ncat, nodestate[i]);
        needsfree[i] = 1;
        loglk[i] = nodestate[i]->loglk;
    }

    t = tree_traverse(POSTORDER, INTERNAL_NODES_ONLY, tree->root, tree);

    node = tree_step(&t);

    while (node) {
        lfdesc = node->children[0];
        rtdesc = node->children[1];
        lfstate = nodestate[lfdesc->index];
        rtstate = nodestate[rtdesc->index];

        _loglk = loglk[lfdesc->index] + loglk[rtdesc->index];
        //Rprintf("%f   ", _loglk);
        aic = -2 * _loglk + 2 * (ngroup[lfdesc->index] + ngroup[rtdesc->index]);
        aicmin = aic;

        while (lfstate) {

            _loglk -= lfstate->loglk;

            while (rtstate) {

                _loglk -= rtstate->loglk;

                // merge states, calc likelihood of merged state, add to loglk, check for improvement
                _loglk += merged_loglk(ncat, lfstate, rtstate);
                //Rprintf("%f\n", _loglk);
                // check for improvement and record if so
                aic = -2 * _loglk + 2 * (ngroup[lfdesc->index] + ngroup[rtdesc->index] - 1);

                if (aic < aicmin || (aic - aicmin) < 3) {
                    k = lfstate;
                    p = rtstate;
                    aicmin = aic;
                }

                _loglk += rtstate->loglk;

                rtstate = rtstate->next;
            }

            _loglk += lfstate->loglk;

            lfstate = lfstate->next;
            rtstate = nodestate[rtdesc->index];
        }

        lfstate = nodestate[lfdesc->index];

        if (k && p) {
            ngroup[node->index] = ngroup[lfdesc->index] + ngroup[rtdesc->index] - 1;
            nodestate[node->index] = alloc_state(ncat);
            nodestate[node->index]->n = k->n + p->n;
            nodestate[node->index]->lmultinomcoef = k->lmultinomcoef + p->lmultinomcoef;
            for (i = 0; i < ncat; ++i)
                nodestate[node->index]->x[i] = k->x[i] + p->x[i];
            nodestate[node->index]->loglk = state_loglk(ncat, nodestate[node->index]);
            needsfree[node->index] = 1;
        } else {
            ++n;
            ngroup[node->index] = ngroup[lfdesc->index] + ngroup[rtdesc->index];
            while (lfstate) {
                if (!lfstate->next) {
                    lfstate->next = nodestate[rtdesc->index];
                    break;
                }
                lfstate = lfstate->next;
            }
            nodestate[node->index] = nodestate[lfdesc->index];
        }

        loglk[node->index] = -0.5 * (aicmin - 2 * (ngroup[node->index]));
        k = NULL;
        p = NULL;
        node = tree_step(&t);
    }
    //
    lfstate = nodestate[tree->root->index];

    while (lfstate) {
        for (i = 0; i < ncat; ++i)
            Rprintf("%d  ", lfstate->x[i]);
        Rprintf("\n");
        lfstate = lfstate->next;
    }
    Rprintf("\n\n");
    lfstate = nodestate[tree->root->children[1]->index];

    //while (lfstate) {
    //    for (i = 0; i < ncat; ++i)
    //        Rprintf("%d  ", lfstate->x[i]);
    //    Rprintf("\n");
    //    lfstate = lfstate->next;
    //}
    //Rprintf("\n\n");

    //
    for (i = 0; i < tree->nnode; ++i) {
        if (needsfree[i])
            free_state(nodestate[i]);
    }
    Free(nodestate);

    return loglk[tree->root->index];
}


SEXP rcl_fitch_multinom(SEXP rtree, SEXP data)
{
    int ncat;
    struct tree *tree;
    GetRNGstate();
    ncat = INTEGER(getAttrib(data, R_DimSymbol))[0];
    tree = (struct tree*)R_ExternalPtrAddr(rtree);
    PutRNGstate();
    return ScalarReal(run_downpass(ncat, INTEGER(data), tree));
}
