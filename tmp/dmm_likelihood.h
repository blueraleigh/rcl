#ifndef DMM_LIKELIHOOD_H
#define DMM_LIKELIHOOD_H

#include <R.h>
#include <Rmath.h>
#include "node.h"
#include "tree.h"
#include "dmm_model.h"

static int PRIORONLY;

static void accumulate_counts(double *loglk, int *count, struct node *node, struct model *model)
{
    if (node->degree) {
        if (!model->event[node->children[0]->index])
            accumulate_counts(loglk, count, node->children[0], model);
        if (!model->event[node->children[1]->index])
            accumulate_counts(loglk, count, node->children[1], model);
    } else {
        int i;
        int k;
        int n = 0;
        for (i = 0; i < model->ncat; ++i) {
            k = model->data[i + node->index * model->ncat];
            count[i] += k;
            n += k;
            *loglk -= lgammafn(k+1);
        }
        *loglk += lgammafn(n+1);
    }
}



/* Compute the log likelihood of the terminal nodes that inherit
** an event of character change beginning at *node under the
** Dirichlet-multinomial model
*/
static double event_loglk(struct node *node, struct model *model)
{
    int i;
    int n;
    int count[model->ncat];
    double loglk;

    memset(count, 0, model->ncat * sizeof(int));
    loglk = 0;

    accumulate_counts(&loglk, count, node, model);

    n = 0;
    for (i = 0; i < model->ncat; ++i) {
        n += count[i];
        loglk += lgammafn(count[i] + model->alpha) - lgammafn(model->alpha);
    }

    loglk += lgammafn(model->ncat * model->alpha) - lgammafn(n + model->ncat * model->alpha);

    return loglk;
}


static double model_loglk(struct model *model)
{
    if (PRIORONLY)
        return 0;

    int i;
    double loglk = 0;
    loglk += event_loglk(model->tree->root, model);
    for (i = 0; i < model->nevent; ++i)
        loglk += event_loglk(model->tree->node[model->eventpos[i]], model);
    return loglk;
}


static double model_logprior(struct model *model)
{
    int i;
    int k;
    double norm;
    double brlen;
    double p;
    double logprior = 0;
    struct node *node;

    norm = model->treelen / model->treescal;

    node = tree_step(&(model->t));

    while (node) {
        brlen = norm * model->rltvrate[node->children[0]->index] * node->children[0]->brlen;
        p = brlen / model->treelen;
        k = model->event[node->children[0]->index];
        logprior += k * log(p) - lgammafn(k+1);
        brlen = norm * model->rltvrate[node->children[1]->index] * node->children[1]->brlen;
        p = brlen / model->treelen;
        k = model->event[node->children[1]->index];
        logprior += k * log(p) - lgammafn(k+1);
        node = tree_step(&(model->t));
    }
    logprior += lgammafn(model->nevent + 1);

    tree_reset(&(model->t));
    /*
    for (i = 0; i < model->tree->nnode; ++i) {
        node = model->tree->node[i];
        brlen = norm * model->rltvrate[node->index] * node->brlen;
        p = brlen / model->treelen;
        k = model->event[node->index];
        logprior += dbinom(k, model->nevent, p, 1);
    }
    */
    logprior += dnbinom(model->nevent, model->beta, 0.5, 1);
    //Rprintf("%f\n", logprior);
    return logprior;
}

#endif
