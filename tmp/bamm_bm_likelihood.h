#ifndef BAMM_BM_LIKELIHOOD_H
#define BAMM_BM_LIKELIHOOD_H

#include <R.h>
#include <Rmath.h>
#include "node.h"
#include "tree.h"
#include "matrix.h"
#include "bamm_bm_model.h"


static int PRIORONLY;


/* Log likelihood contributed by one branch */
static double branch_loglk(double norm, struct node *node, struct model *model)
{
    if (!node)
        return 0;

    double x;
    double y;
    double var;
    double mu;

    y = model->x[node->parent->index];

    x = model->x[node->index];

    mu = model->displ[node->index];
    var = norm * model->rate * model->rltvrate[node->index] * node->brlen;

    if (node->degree)
        return dnorm(x-y, mu, sqrt(var), 1);

    double p;
    double n;
    double k;

    p = exp(x);
    p /= 1+p;

    n = (double)ELEM(model->data, 0, node->index, 2);
    k = (double)ELEM(model->data, 1, node->index, 2);

    return dnorm(x-y, mu, sqrt(var), 1) + dbinom(k, n, p, 1);
}


static double model_loglk(struct model *model)
{
    if (PRIORONLY)
        return 0;

    double loglk = 0;
    double norm;
    struct node *node;

    norm = model->treelen / model->treescal;

    // model->t assumed to have been initialized like so
    // t = tree_traverse(POSTORDER, INTERNAL_NODES_ONLY, tree->root, tree);
    node = tree_step(&(model->t));

    while (node) {
        loglk += branch_loglk(norm, node->children[0], model);
        loglk += branch_loglk(norm, node->children[1], model);
        node = tree_step(&(model->t));
    }

    tree_reset(&(model->t));

    return loglk;
}

#endif
