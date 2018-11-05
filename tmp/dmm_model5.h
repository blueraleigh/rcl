#ifndef DMM_MODEL5_H
#define DMM_MODEL5_H

/* Random local clock Dirichlet-multinomial model of resource use evolution */

#include <R.h>
#include <Rmath.h>
#include "tree.h"


#define ELEM(m, i, j, lda) (m)[(i) + (j) * (lda)]
#define COLUMN(m, j, lda) (m) + (j) * (lda)

static int PRIORONLY;

struct model {
    int ncat;               // Number of resource use categories.
    int nevent;             // Number of events of character change on the tree.
    int nclock;             // Number of local clocks on the tree.
    int nitem1;             // Number of items on the nexteventid stack.
    int nitem2;             // Number of items on the nextclockid stack.
    int accept;             // Number accepted moves.
    int reject;             // Number rejected moves.
    double rate;            // Overall rate of state change.
    double dataloglk;       // Log likelihood for sampled count data given state of each tip.
    double treeloglk;       // Log likelihood for history of state change.
    double loglk;           // Sum of the previous two log likelihoods.

    /* ncat by ntip matrix of sampled resource use counts */
    int *data;

    /* ncat by nnode matrix of cumulative resource use counts
    ** for each event. Column index corresponds to event id. */
    int *count;

    /* nnode length bit vector to indicate whether at least one event
    ** of character change occurs across each branch. */
    int *event;

    /* nnode length vector with the event id for each node. */
    int *eventid;

    /* next event id stack */
    int *nexteventid;

    /* Sum of log multinomial coefficient for count data in each tip */
    double multinomlcoef;

    /* Log likelihood of the data in each event contributed by
    ** cumulative counts in each event */
    double *eventloglk;

    /* Running mean of the indicator variable on each branch,
    ** which is also the marginal probability that a change
    ** occurs on each branch */
    double *margprobevnt;

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

    /* Parameter of Dirichlet hyperprior on multinomial prior
    ** on proportional resource use vector accompanying each event of
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

    norm = model->treelen / model->treescal;
    rate = norm * model->rate * model->rltvrate[node->index];

    return model->event[node->index] ? log1p(-exp(-rate * node->brlen)) : -rate * node->brlen;
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

    return loglk;
}


/* Log likelihood of count data associated with an event */
static double model_eventloglk(int eventid, struct model *model)
{
    int i;
    int n;
    double loglk;

    loglk = 0;
    n = 0;
    for (i = 0; i < model->ncat; ++i) {
        n += model->count[i + eventid * model->ncat];
        loglk += lgammafn(
                model->count[i + eventid * model->ncat] + model->alpha) - lgammafn(model->alpha);
    }
    loglk += lgammafn(model->ncat * model->alpha) - lgammafn(n + model->ncat * model->alpha);

    return loglk;
}


static struct model *model_init(int ncat, struct tree *tree)
{
    struct model *model;

    model = Calloc(1, struct model);
    model->ncat = ncat;
    model->nevent = 0;
    model->nclock = 0;
    model->nitem1 = tree->nnode-1;
    model->nitem2 = tree->nnode-1;
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

    model->count = Calloc(ncat * tree->nnode, int);
    memset(model->count, 0, ncat * tree->nnode * sizeof(int));

    model->event = Calloc(tree->nnode, int);
    memset(model->event, 0, tree->nnode * sizeof(int));

    model->eventid = Calloc(tree->nnode, int);
    memset(model->eventid, 0, tree->nnode * sizeof(int));

    model->nexteventid = Calloc(tree->nnode-1, int);
    memset(model->nexteventid, 0, (tree->nnode-1) * sizeof(int));

    model->eventloglk = Calloc(tree->nnode, double);
    memset(model->eventloglk, 0, tree->nnode * sizeof(double));

    model->rltvrate = Calloc(tree->nnode, double);
    model->ratemultipl = Calloc(tree->nnode, double);
    model->clock = Calloc(tree->nnode, int);
    model->clockid = Calloc(tree->nnode, int);
    model->clockpos = Calloc(tree->nnode, int);
    model->nextclockid = Calloc(tree->nnode-1, int);
    memset(model->nextclockid, 0, (tree->nnode-1) * sizeof(int));

    model->tree = tree;

    model->t = tree_traverse(POSTORDER, INTERNAL_NODES_ONLY, tree->root, tree);

    int i;
    for (i = 0; i < tree->nnode; ++i) {
        model->treelen += tree->node[i]->brlen;
        model->ratemultipl[i] = 1;
        model->rltvrate[i] = 1;
        model->clock[i] = 0;
        model->clockid[i] = 0;
        model->clockpos[i] = 0;
    }
    model->treescal = model->treelen;

    int eventid = tree->nnode - 1;
    int clockid = tree->nnode - 1;
    for (i = 0; i < (tree->nnode-1); ++i) {
        model->nexteventid[i] = eventid--;
        model->nextclockid[i] = clockid--;
    }

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
        model->eventloglk[0] += lgammafn(
                model->count[i + 0 * model->ncat] + model->alpha) - lgammafn(model->alpha);
    }
    model->eventloglk[0] += lgammafn(model->ncat * model->alpha) - lgammafn(n + model->ncat * model->alpha);

    model->dataloglk = model->multinomlcoef + model->eventloglk[0];
    model->treeloglk = model_treeloglk(model);
    model->loglk = model->dataloglk + model->treeloglk;
}


static void model_free(struct model *model)
{
    Free(model->data);
    Free(model->count);
    Free(model->event);
    Free(model->eventid);
    Free(model->nexteventid);
    Free(model->eventloglk);
    Free(model->rltvrate);
    Free(model->ratemultipl);
    Free(model->clock);
    Free(model->clockpos);
    Free(model->clockid);
    Free(model->nextclockid);
    Free(model);
}


#endif
