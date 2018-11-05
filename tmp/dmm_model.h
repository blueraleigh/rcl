#ifndef DMM_MODEL_H
#define DMM_MODEL_H

/* Random local clock Dirichlet-multinomial model of resource use evolution */

#include <R.h>
#include <Rmath.h>
#include "tree.h"


#define ELEM(m, i, j, lda) (m)[(i) + (j) * (lda)]
#define COLUMN(m, j, lda) (m) + (j) * (lda)


struct model {
    int ncat;           // Number of resource use categories.
    int nevent;         // Number of events of character change on the tree.
    int nclck;          // Number of local clocks on the tree.
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


    /* Relative branch rates. Always 1 for the root. */
    double *rltvrate;


    /*
    ** Rate multipliers for each branch. The relative
    ** rate of a branch is equal to the relative rate
    ** of its parent branch multiplied by a scalar.
    */
    double *ratemultipl;


    /*
    ** Bit vector to indicate whether a branch has a rate
    ** of evolution that differs from its parent.
    */
    int *rlclck;


    /*
    ** Storage array that holds the rate multipliers associated
    ** with each local clock. Only the first nclck positions are
    ** in use but it can hold up to nnode clocks.
    */
    double *clckrate;

    /*
    ** Storage array that holds the node indices where each local
    ** clock originates. Only the first nclck positions are
    ** in use but it can hold up to nnode clocks.
    */
    int *clckpos;

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

    /* Shape parameter of gamma hyperprior on rate of Poisson process prior.
    ** The rate parameter of gamma hyperprior is set equal to 1. Under this
    ** formulation the prior distribution on the number of events is a
    ** negative binomial with r = beta and p = 0.5, implying that the
    ** expected number of events equals beta.
    */
    double beta;

    struct tree *tree;

    struct tree_traversal t;
};


static struct model *model_init(int ncat, struct tree *tree)
{
    struct model *model;

    model = Calloc(1, struct model);
    model->ncat = ncat;
    model->nevent = 0;
    model->nclck = 0;
    model->accept = 0;
    model->reject = 0;
    model->loglk = 0;
    model->logprior = 0;
    model->treescal = 0;
    model->treelen = 0;
    model->alpha = 0;
    model->beta = 0;

    model->data = Calloc(ncat * tree->ntip, int);
    memset(model->data, 0, ncat * tree->ntip * sizeof(int));

    model->event = Calloc(tree->nnode, int);
    memset(model->event, 0, tree->nnode * sizeof(int));

    model->eventpos = Calloc(tree->nnode, int);
    memset(model->eventpos, 0, tree->nnode * sizeof(int));

    model->rltvrate = Calloc(tree->nnode, double);
    model->ratemultipl = Calloc(tree->nnode, double);
    model->rlclck = Calloc(tree->nnode, int);
    model->clckrate = Calloc(tree->nnode, double);
    model->clckpos = Calloc(tree->nnode, int);

    model->tree = tree;

    model->t = tree_traverse(POSTORDER, INTERNAL_NODES_ONLY, tree->root, tree);

    int i;
    for (i = 0; i < tree->nnode; ++i) {
        model->treelen += tree->node[i]->brlen;
        model->ratemultipl[i] = 1;
        model->rltvrate[i] = 1;
        model->rlclck[i] = 0;
        model->clckrate[i] = 0;
        model->clckpos[i] = 0;
    }
    model->treescal = model->treelen;

    return model;
}


static void model_set(double alpha, double beta, int *tipdata, struct model *model)
{
    memcpy(model->data, tipdata, model->ncat * model->tree->ntip * sizeof(int));
    model->alpha = alpha;
    model->beta = beta;
}


static void model_free(struct model *model)
{
    Free(model->data);
    Free(model->event);
    Free(model->eventpos);
    Free(model->rltvrate);
    Free(model->ratemultipl);
    Free(model->rlclck);
    Free(model->clckrate);
    Free(model->clckpos);
    Free(model);
}


#endif
