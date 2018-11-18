#ifndef BMM_MODEL_H
#define BMM_MODEL_H

/* Random local clock Brownian motion with jumps model of continuous trait evolution */

#include <R.h>
#include <Rmath.h>
#include "tree.h"


static int PRIORONLY;

struct model {
    int nclock;             // Number of local clocks on the tree.
    int njump;              // Number of jumps on the tree.
    int nitem;              // Number of items on the nextclockid stack.
    int accept;             // Number accepted moves.
    int reject;             // Number rejected moves.
    double rate;            // Overall rate of state change.
    double loglk;           // Sum of the previous two log likelihoods.

    /* nnode length array of continuous node states */
    double *x;

    /* expected displacement across each branch */
    double *mu;

    /* expected variance across each branch */
    double *var;

    /* Bit vector to indicate whether phenotypic jump occurs across branch */
    int *jump;

    /*
    ** Bit vector to indicate whether a branch has a rate
    ** of evolution that differs from its parent (i.e. initiates
    ** a new random local clock)
    */
    int *clock;

    /* nnode length vector with the clock id for each node */
    int *clockid;

    /* next clock id stack */
    int *nextclockid;

    /*
    ** Storage array that holds the node indices where each local
    ** clock originates. */
    int *clockpos;

    /*
    ** Rate multipliers for each branch. These are only different from
    ** 1 for the nodes where random local clocks originate.
    */
    double *ratemultipl;

    /* Relative branch rates. Always 1 for the root. The
    ** relative rate for a branch is equal to the relative
    ** rate of its parent branch multiplied by a scalar, which
    ** is 1 for every branch that does not initiate a new
    ** random local clock. */
    double *rltvrate;

    /*
    ** Normalization constant for relative branch rates is
    ** given by treelen / treescal
    **
    ** The actual rate on a branch is given by the product
    ** of the overall rate of evolution, the normalization
    ** constant, and the relative rate of the branch. The
    ** normalization constant ensures that the sum of the
    ** product of each branch length and its normalized relative
    ** rate equals the sum of all branch lengths. This allows
    ** the rate of evolution to vary over phylogeny while
    ** preserving the total absolute time in which evolution has
    ** had to occur.
    */
    double treescal;

    /* Sum of all branch lengths */
    double treelen;

    struct tree *tree;

    struct tree_traversal t;
};


static double model_nodeloglk(struct node *node, struct model *model)
{
    double u;       // contrast
    double mu;      // contrast expectation
    double v1;
    double v2;
    double var;     // contrast variance
    double delta;
    double norm;
    struct node *lfdesc;
    struct node *rtdesc;

    norm = model->treelen / model->treescal;

    lfdesc = node->children[0];
    rtdesc = node->children[1];

    if (!lfdesc->degree)
        model->var[lfdesc->index] = (norm * model->rate * model->rltvrate[lfdesc->index] * lfdesc->brlen);

    if (!rtdesc->degree)
        model->var[rtdesc->index] = (norm * model->rate * model->rltvrate[rtdesc->index] * rtdesc->brlen);

    v1 = model->var[lfdesc->index];

    v2 = model->var[rtdesc->index];

    var = v1 + v2;

    mu = model->mu[lfdesc->index] - model->mu[rtdesc->index];

    u = model->x[lfdesc->index] - model->x[rtdesc->index];

    delta = u - mu;

    model->var[node->index] =
        (norm * model->rate * model->rltvrate[node->index] * node->brlen) + (v1*v2) / (v1 + v2);

    model->x[node->index] =
        (v2 * model->x[lfdesc->index] + v1 * model->x[rtdesc->index]) / (v1 + v2);

    return -0.5 * (M_LN_2PI + log(var)) - 0.5 * (delta * delta / var);
}


static double model_loglk(struct model *model)
{
    if (PRIORONLY)
        return 0;

    double loglk = 0;
    struct node *node;

    node = tree_step(&(model->t));

    while (node) {

        loglk += model_nodeloglk(node, model);

        node = tree_step(&(model->t));
    }

    tree_reset(&(model->t));

    return loglk;
}


static struct model *model_init(struct tree *tree)
{
    struct model *model;

    model = Calloc(1, struct model);
    model->nclock = 0;
    model->njump = 0;
    model->nitem = tree->nnode-1;
    model->accept = 0;
    model->reject = 0;
    model->rate = 0;
    model->loglk = 0;
    model->treescal = 0;
    model->treelen = 0;

    model->x = Calloc(tree->nnode, double);
    memset(model->x, 0, tree->nnode * sizeof(double));

    model->mu = Calloc(tree->nnode, double);
    memset(model->mu, 0, tree->nnode * sizeof(double));

    model->var = Calloc(tree->nnode, double);
    memset(model->var, 0, tree->nnode * sizeof(double));

    model->jump = Calloc(tree->nnode, int);
    memset(model->jump, 0, tree->nnode * sizeof(int));

    model->rltvrate = Calloc(tree->nnode, double);
    model->ratemultipl = Calloc(tree->nnode, double);
    model->clock = Calloc(tree->nnode, int);
    model->clockid = Calloc(tree->nnode, int);
    model->clockpos = Calloc(tree->nnode, int);
    model->nextclockid = Calloc(tree->nnode-1, int);
    memset(model->nextclockid, 0, (tree->nnode-1) * sizeof(int));

    model->tree = tree;

    model->t = tree_traverse(POSTORDER, INTERNAL_NODES_ONLY, tree->root, tree);

    int i, j;
    for (i = 0; i < tree->nnode; ++i) {
        model->treelen += tree->node[i]->brlen;
        model->ratemultipl[i] = 1;
        model->rltvrate[i] = 1;
        model->clock[i] = 0;
        model->clockid[i] = 0;
        model->clockpos[i] = 0;
    }
    model->treescal = model->treelen;

    int clockid = tree->nnode - 1;
    for (i = 0; i < (tree->nnode-1); ++i)
        model->nextclockid[i] = clockid--;

    return model;
}


static void model_set(double rate, double *tipdata, struct model *model)
{
    memcpy(model->x, tipdata, model->tree->ntip * sizeof(double));
    model->rate = rate;

    if (PRIORONLY)
        return;

    model->loglk = model_loglk(model);
}


static void model_free(struct model *model)
{
    Free(model->x);
    Free(model->mu);
    Free(model->var);
    Free(model->jump);
    Free(model->rltvrate);
    Free(model->ratemultipl);
    Free(model->clock);
    Free(model->clockpos);
    Free(model->clockid);
    Free(model->nextclockid);
    Free(model);
}


#endif
