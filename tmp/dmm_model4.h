#ifndef DMM_MODEL4_H
#define DMM_MODEL4_H

/* Random local clock Dirichlet-multinomial model of resource use evolution */

#include <R.h>
#include <Rmath.h>
#include "tree.h"
#include "atblspc.h"


#define ELEM(m, i, j, lda) (m)[(i) + (j) * (lda)]
#define COLUMN(m, j, lda) (m) + (j) * (lda)

static int PRIORONLY;

struct model {
    int ncat;           // Number of resource use categories.
    int nevent;         // Number of events of character change on the tree.
    int nitem;          // Number of ids on the next event id stack
    int accept;         // Number accepted moves.
    int reject;         // Number rejected moves.
    double loglk;       // Model loglikelihood.
    double logprior;    // Model log prior on event number.

    /* ncat by ntip matrix of sampled resource use counts */
    int *data;

    /* ncat by nnode matrix of cumulative resource use counts
    ** for each event. Column index corresponds to event id. */
    int *count;

    /* nnode length bit vector indicating whether an event occurs
    ** on each branch */
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
    ** which is also the marginal posterior probability
    ** that an event occurs on each branch */
    double *margprobevnt;

    /* Sum of all branch lengths */
    double treelen;

    /* Parameter of Dirichlet hyperprior on multinomial prior
    ** on proportional resource use vector accompanying each event of
    ** character change
    */
    double alpha;

    /* Shape parameter of gamma hyperprior on rate of Poisson process prior.
    ** The rate parameter of gamma hyperprior is set equal to 1. Under this
    ** formulation the prior distribution on the number of events is a
    ** negative binomial with r = beta and p = 0.5, implying that the
    ** expected number of events equals beta.
    */
    double beta;

    struct tree *tree;

    struct tree_traversal t;

    /* Structure for constant time sampling of branches with probability
    ** equal to relative branch length */
    struct atblspc *atblspc;
};


static struct model *model_init(int ncat, struct tree *tree)
{
    struct model *model;

    model = Calloc(1, struct model);
    model->ncat = ncat;
    model->nevent = 0;
    model->nitem = tree->nnode-1;
    model->accept = 0;
    model->reject = 0;
    model->loglk = 0;
    model->logprior = 0;
    model->treelen = 0;
    model->alpha = 0;
    model->beta = 0;
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

    model->margprobevnt = Calloc(tree->nnode, double);
    memset(model->margprobevnt, 0, tree->nnode * sizeof(double));

    model->tree = tree;

    model->t = tree_traverse(POSTORDER, INTERNAL_NODES_ONLY, tree->root, tree);

    model->atblspc = atblspc_alloc(tree->nnode);

    int i;
    int eventid;
    double w[tree->nnode];
    for (i = 0; i < tree->nnode; ++i)
        model->treelen += tree->node[i]->brlen;

    for (i = 0; i < tree->nnode; ++i)
        w[i] = tree->node[i]->brlen / model->treelen;

    eventid = tree->nnode - 1;
    for (i = 0; i < (tree->nnode-1); ++i)
        model->nexteventid[i] = eventid--;

    atblspc_init(w, model->atblspc);

    return model;
}


static void model_set(double alpha, double beta, int *tipdata, struct model *model)
{
    memcpy(model->data, tipdata, model->ncat * model->tree->ntip * sizeof(int));
    model->alpha = alpha;
    model->beta = beta;

    int i;
    int j;
    int n;

    model->logprior = dnbinom(0, model->beta, 0.5, 1);

    if (PRIORONLY)
        return;

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

    model->loglk = model->multinomlcoef + model->eventloglk[0];
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


static void model_updatealpha(double alpha, struct model *model)
{
    if (PRIORONLY)
        return;

    int i;
    model->alpha = alpha;
    for (i = 0; i < model->tree->nnode; ++i) {
        model->loglk -= model->eventloglk[i];
        model->eventloglk[i] = model_eventloglk(i, model);
        model->loglk += model->eventloglk[i];
    }
}


static void model_free(struct model *model)
{
    atblspc_free(model->atblspc);
    Free(model->data);
    Free(model->count);
    Free(model->eventid);
    Free(model->event);
    Free(model->nexteventid);
    Free(model->eventloglk);
    Free(model->margprobevnt);
    Free(model);
}


#endif
