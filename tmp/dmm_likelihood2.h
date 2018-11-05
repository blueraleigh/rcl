#ifndef DMM_LIKELIHOOD2_H
#define DMM_LIKELIHOOD2_H

#include <R.h>
#include <Rmath.h>
#include "node.h"
#include "tree.h"
#include "dmm_model2.h"

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

#endif
