#ifndef MODEL_H
#define MODEL_H

/* Random local clock Poisson process model of discrete trait evolution */

#include <R.h>
#include <Rmath.h>
#include "tree.h"


#define ELEM(m, i, j, lda) (m)[(i) + (j) * (lda)]
#define COLUMN(m, j, lda) (m) + (j) * (lda)


struct model {
    int K;              // Number of bins (=states) for the beta distribution.
    int nclck;          // Number of local clocks on the tree.
    int accept;         // Number accepted moves.
    int reject;         // Number rejected moves.
    double loglk;       // Model loglikelihood.
    double step;        // Width of each bin.
    double *clk;        // Conditional likelihood matrix. Dimension K by Nnode.


    /* Overall rate of evolution */
    double rate;


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


    struct tree *tree;

    struct tree_traversal t;
};


static struct model *model_init(int k, struct tree *tree)
{
    struct model *model;

    model = Calloc(1, struct model);
    model->K = k;
    model->nclck = 0;
    model->accept = 0;
    model->reject = 0;
    model->loglk = 0;
    model->rate = 0;
    model->treescal = 0;
    model->treelen = 0;
    model->step = 1 / (double)k;

    model->clk = Calloc(k * tree->nnode, double);
    memset(model->clk, 0, k * tree->nnode * sizeof(double));

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


static void model_set(double rate, double *brks, int *tipdata, struct model *model)
{
    int i;
    int j;
    int n;
    int k;
    double cdf[model->K+1];

    model->rate = rate;

    cdf[0] = 0;
    cdf[model->K] = 1;

    for (i = 0; i < model->tree->ntip; ++i) {
        n = ELEM(tipdata, 0, i, 2);
        k = ELEM(tipdata, 1, i, 2);

        for (j = 1; j < model->K; ++j)
            cdf[j] = pbeta(brks[j], k+1, n-k+1, 1, 0);

        for (j = 0; j < model->K; ++j)
            ELEM(model->clk, j, i, model->K) = cdf[j+1] - cdf[j];

    }
}


static void model_set2(double rate, double *brks, int *tipdata, struct model *model)
{
    int i;
    int j;
    int k;
    int n;
    int n1;
    int n2;
    int nbin;

    // TODO: find a better way
    nbin = (int)sqrt(model->K);

    double cdf1[nbin+1];
    double cdf2[nbin+1];

    model->rate = rate;

    cdf1[0] = cdf2[0] = 0;
    cdf1[nbin] = cdf2[nbin] = 1;

    for (i = 0; i < model->tree->ntip; ++i) {
        n = ELEM(tipdata, 0, i, 3);
        n1 = ELEM(tipdata, 1, i, 3);
        n2 = ELEM(tipdata, 2, i, 3);

        for (j = 1; j < nbin; ++j) {
            cdf1[j] = pbeta(brks[j], n1+1, n-n1+1, 1, 0);
            cdf2[j] = pbeta(brks[j], n2+1, n-n1-n2+1, 1, 0);
        }

        for (j = 0; j < nbin; ++j) {
            for (k = 0; k < nbin; ++k)
                ELEM(model->clk, j+k*nbin, i, model->K) = (cdf1[j+1] - cdf1[j]) * (cdf2[k+1] - cdf2[k]);
        }
    }
}


static void model_free(struct model *model)
{
    Free(model->clk);
    Free(model->rltvrate);
    Free(model->ratemultipl);
    Free(model->rlclck);
    Free(model->clckrate);
    Free(model->clckpos);
    Free(model);
}


#endif
