#ifndef BAMM_BM_MODEL_H
#define BAMM_BM_MODEL_H

/* Random local clock Brownian motion with jumps model of continuous trait evolution */

#include <R.h>
#include <Rmath.h>
#include "tree.h"


#define ELEM(m, i, j, lda) (m)[(i) + (j) * (lda)]
#define COLUMN(m, j, lda) (m) + (j) * (lda)


struct model {
    int nclck;          // Number of local clocks on the tree.
    int njmp;           // Number of jumps on the tree.
    int accept;         // Number accepted moves.
    int reject;         // Number rejected moves.
    double loglk;       // Model loglikelihood.

    /*
    ** 2 by ntip matrix of counts. First row
    ** is the total number of counts, second
    ** row is the number of counts in the focal
    ** category.
    */
    int *data;

    double *x;          // length nnode array of log probabilities.


    /*
    ** Expected displacement from ancestral state.
    ** Only non-zero at positions where *jump is 1.
    ** Prior is normal distribution that is
    ** centered on zero with arbitrary variance.
    ** length nnode array.
    */
    double *displ;


    /*
    ** Bit vector to indicate whether a phenotypic jump
    ** occurs on a branch.
    */
    int *jump;


    /*
    ** Storage array that holds the node indices where each
    ** jump occurs. Only the first njmp positions are
    ** in use but it can hold up to nnode jumps.
    */
    int *jmppos;


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


static struct model *model_init(struct tree *tree)
{
    struct model *model;

    model = Calloc(1, struct model);
    model->nclck = 0;
    model->accept = 0;
    model->reject = 0;
    model->loglk = 0;
    model->treescal = 0;
    model->treelen = 0;
    model->rate = 0;

    model->data = Calloc(2 * tree->ntip, int);
    memset(model->data, 0, 2 * tree->ntip * sizeof(int));

    model->x = Calloc(tree->nnode, double);
    memset(model->x, 0, tree->nnode * sizeof(double));

    model->rltvrate = Calloc(tree->nnode, double);
    model->ratemultipl = Calloc(tree->nnode, double);
    model->rlclck = Calloc(tree->nnode, int);
    model->clckrate = Calloc(tree->nnode, double);
    model->clckpos = Calloc(tree->nnode, int);

    model->jmppos = Calloc(tree->nnode, int);
    model->displ = Calloc(tree->nnode, double);
    model->jump = Calloc(tree->nnode, int);
    memset(model->displ, 0, tree->nnode * sizeof(double));
    memset(model->jump, 0, tree->nnode * sizeof(int));
    memset(model->jmppos, 0, tree->nnode * sizeof(int));

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


static void model_set(double rate, int *tipdata, struct model *model)
{
    int i;
    double p;
    double v1;
    double v2;
    double x1;
    double x2;
    struct node *node;

    memcpy(model->data, tipdata, 2 * model->tree->ntip * sizeof(int));
    model->rate = rate;

    for (i = 0; i < model->tree->ntip; ++i) {
        p = (double)(ELEM(tipdata, 1, i, 2)+1) / (double)(ELEM(tipdata, 0, i, 2) + 2);
        model->x[i] = log(p / (1-p));
    }

    node = tree_step(&(model->t));

    while (node) {
        v1 = node->children[0]->brlen;
        v2 = node->children[1]->brlen;
        x1 = model->x[node->children[0]->index];
        x2 = model->x[node->children[1]->index];
        model->x[node->index] = (v2*x1 + v1*x2) / (v1 + v2);
        node = tree_step(&(model->t));
    }

    tree_reset(&(model->t));
}


static void model_free(struct model *model)
{
    Free(model->data);
    Free(model->x);
    Free(model->displ);
    Free(model->jump);
    Free(model->jmppos);
    Free(model->rltvrate);
    Free(model->ratemultipl);
    Free(model->rlclck);
    Free(model->clckrate);
    Free(model->clckpos);
    Free(model);
}


#endif
