#ifndef DMM_MODEL_H
#define DMM_MODEL_H

/* Random local clock Dirichlet-multinomial model of resource use evolution */

#include <R.h>
#include <Rmath.h>
#include "tree.h"


static int PRIORONLY;

struct model {
    int ncat;               // Number of resource use categories.
    int nstate;             // Number of states in the model.
    int nclock;             // Number of local clocks on the tree.
    int nitem;              // Number of items on the nextclockid stack.
    int accept;             // Number accepted moves.
    int reject;             // Number rejected moves.
    double rate;            // Overall rate of state change.
    double dataloglk;       // Log likelihood for sampled count data given state of each tip.
    double treeloglk;       // Log likelihood for history of state change.
    double loglk;           // Sum of the previous two log likelihoods.

    /* ncat by ntip matrix of sampled resource use counts */
    int *data;

    /* ncat by nstate matrix of cumulative resource use counts
    ** for each state. Column index corresponds to state id. */
    int *count;

    /* nnode array that records state id for each node. */
    int *stateid;

    /* Sum of log multinomial coefficient for count data in each tip */
    double multinomlcoef;

    /* Log likelihood of the data contributed by
    ** cumulative counts in each state */
    double *stateloglk;

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

    /* Hyperparameter of Dirichlet prior on multinomial parameters
    ** of proportional resource use vector accompanying each event of
    ** character change
    */
    double alpha;

    struct tree *tree;

    struct tree_traversal t;
};


static double model_branchloglk(struct node *node, struct model *model)
{
    double rate;
    double norm;
    double D;
    double F;
    double pii;
    double pij;

    norm = model->treelen / model->treescal;
    rate = norm * model->rate * model->rltvrate[node->index];

    D = exp(-rate*model->nstate*node->brlen);
    F = 1 - D;

    pii = log(D + F / model->nstate);
    pij = log(F / model->nstate);

    return model->stateid[node->index] == model->stateid[node->parent->index] ? pii : pij;
}


static double model_treeloglk(struct model *model)
{
    double loglk = 0;
    struct node *node;

    node = tree_step(&(model->t));

    while (node) {

        loglk += model_branchloglk(node->children[0], model);
        loglk += model_branchloglk(node->children[1], model);

        node = tree_step(&(model->t));
    }

    tree_reset(&(model->t));

    return loglk - log(model->nstate);
}


/* Log likelihood of count data associated with an state */
static double model_stateloglk(int stateid, struct model *model)
{
    int i;
    int n;
    double loglk;

    n = 0;
    loglk = 0;
    for (i = 0; i < model->ncat; ++i) {
        n += model->count[i + stateid * model->ncat];
        loglk += lgammafn(
                model->count[i + stateid * model->ncat] + model->alpha) - lgammafn(model->alpha);
    }
    loglk += lgammafn(model->ncat * model->alpha) - lgammafn(n + model->ncat * model->alpha);

    return loglk;
}


static double model_dataloglk(struct model *model)
{
    int i;
    double loglk;
    loglk = 0;
    for (i = 0; i < model->nstate; ++i)
        loglk += model_stateloglk(i, model);
    return loglk;
}


static struct model *model_init(int ncat, int nstate, struct tree *tree)
{
    struct model *model;

    model = Calloc(1, struct model);
    model->ncat = ncat;
    model->nstate = nstate;
    model->nclock = 0;
    model->nitem = tree->nnode-1;
    model->accept = 0;
    model->reject = 0;
    model->rate = 0;
    model->dataloglk = 0;
    model->treeloglk = 0;
    model->loglk = 0;
    model->treescal = 0;
    model->treelen = 0;
    model->alpha = 0;
    model->multinomlcoef = 0;

    model->data = Calloc(ncat * tree->ntip, int);
    memset(model->data, 0, ncat * tree->ntip * sizeof(int));

    model->count = Calloc(ncat * nstate, int);
    memset(model->count, 0, ncat * nstate * sizeof(int));

    model->stateid = Calloc(tree->nnode, int);
    memset(model->stateid, 0, tree->nnode * sizeof(int));

    model->stateloglk = Calloc(nstate, double);
    memset(model->stateloglk, 0, nstate * sizeof(double));

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


static void model_set(double rate, double alpha, int *tipdata, struct model *model)
{
    memcpy(model->data, tipdata, model->ncat * model->tree->ntip * sizeof(int));
    model->rate = rate;
    model->alpha = alpha;

    if (PRIORONLY)
        return;

    int i, j, n;

    for (i = 0; i < model->tree->ntip; ++i) {
        n = 0;
        for (j = 0; j < model->ncat; ++j) {
            n += model->data[j + i * model->ncat];
            model->multinomlcoef -= lgammafn(model->data[j + i * model->ncat] + 1);
        }
        model->multinomlcoef += lgammafn(n + 1);
    }

    for (i = 0; i < model->tree->ntip; ++i) {
        for (j = 0; j < model->ncat; ++j)
            model->count[j + 0 * model->ncat] += model->data[j + i * model->ncat];
    }

    n = 0;
    for (i = 0; i < model->ncat; ++i) {
        n += model->count[i + 0 * model->ncat];
        model->stateloglk[0] += lgammafn(
                model->count[i + 0 * model->ncat] + model->alpha) - lgammafn(model->alpha);
    }
    model->stateloglk[0] += lgammafn(model->ncat * model->alpha) - lgammafn(n + model->ncat * model->alpha);

    model->dataloglk = model->multinomlcoef + model->stateloglk[0];
    model->treeloglk = model_treeloglk(model);
    model->loglk = model->dataloglk + model->treeloglk;
}


static void model_free(struct model *model)
{
    Free(model->data);
    Free(model->count);
    Free(model->stateid);
    Free(model->stateloglk);
    Free(model->rltvrate);
    Free(model->ratemultipl);
    Free(model->clock);
    Free(model->clockpos);
    Free(model->clockid);
    Free(model->nextclockid);
    Free(model);
}


#endif
