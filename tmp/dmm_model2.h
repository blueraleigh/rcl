#ifndef DMM_MODEL2_H
#define DMM_MODEL2_H

/* Dirichlet-multinomial model of resource use evolution */

#include <R.h>
#include <Rmath.h>
#include "tree.h"


#define ELEM(m, i, j, lda) (m)[(i) + (j) * (lda)]
#define COLUMN(m, j, lda) (m) + (j) * (lda)


struct model {
    int ncat;           // Number of resource use categories.
    int nevent;         // Number of events of character change on the tree.
    int accept;         // Number accepted moves.
    int reject;         // Number rejected moves.
    double loglk;       // Model loglikelihood.
    double logprior;    // Model log prior on event configuration.

    /* ncat by ntip matrix of sampled resource use counts */
    int *data;

    /* nnode length bit vector to indicate whether an event
    ** of character change occurs on each branch
    */
    int *event;

    /*
    ** Storage array that holds the node indices where each
    ** event originates. Only the first nevent positions are
    ** in use but it can hold up to nnode events.
    */
    int *eventpos;

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


static struct model *model_init(int ncat, struct tree *tree)
{
    struct model *model;

    model = Calloc(1, struct model);
    model->ncat = ncat;
    model->nevent = 0;
    model->accept = 0;
    model->reject = 0;
    model->loglk = 0;
    model->logprior = 0;
    model->treelen = 0;
    model->alpha = 0;

    model->data = Calloc(ncat * tree->ntip, int);
    memset(model->data, 0, ncat * tree->ntip * sizeof(int));

    model->event = Calloc(tree->nnode, int);
    memset(model->event, 0, tree->nnode * sizeof(int));

    model->eventpos = Calloc(tree->nnode, int);
    memset(model->eventpos, 0, tree->nnode * sizeof(int));

    model->tree = tree;

    model->t = tree_traverse(POSTORDER, INTERNAL_NODES_ONLY, tree->root, tree);

    int i;
    for (i = 0; i < tree->nnode; ++i)
        model->treelen += tree->node[i]->brlen;

    return model;
}


static void model_set(double alpha, int *tipdata, struct model *model)
{
    memcpy(model->data, tipdata, model->ncat * model->tree->ntip * sizeof(int));
    model->alpha = alpha;
}


static void model_free(struct model *model)
{
    Free(model->data);
    Free(model->event);
    Free(model->eventpos);
    Free(model);
}


#endif
