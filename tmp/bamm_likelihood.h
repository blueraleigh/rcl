#ifndef LIKELIHOOD_H
#define LIKELIHOOD_H

#include <R.h>
#include <Rmath.h>
#include <R_ext/BLAS.h>
#include "node.h"
#include "tree.h"
#include "matrix.h"
#include "bamm_model.h"

#define ALLROOT 0
#define FLATROOT 1
#define OBSROOT 2


static int PRIORONLY;
static int ROOTP;


static double rootp(double ls, struct model *model)
{
    int i;
    int size = model->K;
    double *sclk = COLUMN(model->clk, model->tree->root->index, model->K);
    double loglk = 0, p = 0, tot = 0;

    switch (ROOTP) {
        case ALLROOT:
            for (i = 0; i < size; ++i)
                p += sclk[i];
            loglk = log(p) + ls;
            break;
        case FLATROOT:
            for (i = 0; i < size; ++i)
                p += sclk[i] / size;
            loglk = log(p) + ls;
            break;
        case OBSROOT:
            for (i = 0; i < size; ++i)
                tot += sclk[i];
            for (i = 0; i < size; ++i)
                p += (sclk[i] / tot) * sclk[i];
            loglk = log(p) + ls;
            break;
    }

    return loglk;
}


static double mkp_branch_prob(int nstate, double brlen, double rate, double *init, double *out)
{
    int i;
    int j;
    int ione = 1;
    double sf;
    double isf;
    double Pr[nstate * nstate];

    double D = exp(-rate*brlen);
    double F = 1 - D;
    double prob = 1/(double)(nstate);

    for (i = 0; i < nstate; ++i)
        Pr[i + i * nstate] = D + F * prob;
    for (i = 0; i < (nstate-1); ++i) {
        for (j = (i+1); j < nstate; ++j) {
            Pr[i + j * nstate] = F * prob;
            Pr[j + i * nstate] = F * prob;
        }
    }

    // multiply that matrix with the initial conditional likelihoods
    // in the init vector, storing the new conditional likelihoods
    // in the out vector
    matrix_vector_product(Pr, nstate, nstate, init, out);

    // compute the scale factor, which is just the sum of the new
    // conditional likelihoods
    sf = F77_CALL(dasum)(&nstate, out, &ione);

    // rescale the conditional likelihoods by this factor
    isf = 1 / sf;
    F77_CALL(dscal)(&nstate, &isf, out, &ione);

    return sf;
}


static double mkp_node_lk(double norm, struct node *node, struct model *model)
{
    int i;
    double *init;
    double rate;
    double lsf;
    double rsf;
    double lf[model->K];
    double rt[model->K];

    // do the left branch
    init = COLUMN(model->clk, node->children[0]->index, model->K);
    rate = norm * model->rate * model->rltvrate[node->children[0]->index];
    lsf = mkp_branch_prob(model->K, node->children[0]->brlen, rate,
        init, lf);

    // do the right branch
    init = COLUMN(model->clk, node->children[1]->index, model->K);
    rate = norm * model->rate * model->rltvrate[node->children[1]->index];
    rsf = mkp_branch_prob(model->K, node->children[1]->brlen, rate,
        init, rt);

    // to form the conditional likelihoods for the ancestral node
    // we multiply the contributions from each daughter
    for (i = 0; i < model->K; ++i)
        ELEM(model->clk, i, node->index, model->K) = lf[i] * rt[i];

    // add up the logarithms of the scale factors as we will
    // use this for computing the final (unscaled) log likelihood
    return log(lsf) + log(rsf);
}


static double model_loglk(struct model *model)
{
    if (PRIORONLY)
        return 0;

    double ls = 0;
    double norm;
    struct node *node;

    norm = model->treelen / model->treescal;

    // model->t assumed to have been initialized like so
    // t = tree_traverse(POSTORDER, INTERNAL_NODES_ONLY, tree->root, tree);
    node = tree_step(&(model->t));

    while (node) {
        ls += mkp_node_lk(norm, node, model);
        node = tree_step(&(model->t));
    }

    tree_reset(&(model->t));

    return rootp(ls, model);
}

#endif
