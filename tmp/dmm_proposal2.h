#ifndef DMM_PROPOSAL2_H
#define DMM_PROPOSAL2_H

#include <R.h>
#include <Rmath.h>
#include "node.h"
#include "tree.h"
#include "dmm_model2.h"
#include "dmm_likelihood2.h"

/*
** Move types
**     add an event of character change
**     remove an event of character change
*/


static double PRIOR0;   // Poisson prior on number of events


// Add/remove an event of character change
static void propose0(struct model *model)
{
    int i;
    int j;
    int k;
    int evnt;
    int eventpos[model->nevent];
    double q;
    double p;
    double loglk;
    struct node *node;

    i = (int)((model->tree->nnode - 1) * unif_rand());

    if (i >= model->tree->root->index)
        i += 1;

    node = model->tree->node[i];


    if (model->event[i]) {
        // Remove the event. Proposal ratio is (N-K+1) / K,
        // where N is nbr branches and K is nbr events

        for (j = 0; j < model->nevent; ++j)
            if (model->eventpos[j] == i)
                evnt = j;

        eventpos[evnt] = model->eventpos[evnt];

        for (k = evnt; k < (model->nevent - 1); ++k) {
            eventpos[k+1] = model->eventpos[k+1];
            model->eventpos[k] = model->eventpos[k+1];
        }

        q = log((double)(model->tree->nnode - model->nevent) / (double)model->nevent);

        // Prior ratio
        p = log(model->nevent / PRIOR0);

        model->nevent -= 1;
        model->event[i] = 0;

        loglk = model_loglk(model);


        if (unif_rand() < exp(loglk - model->loglk + p + q)) {
            // accept
            model->loglk = loglk;
            model->accept += 1;
        } else {
            //reject
            model->nevent += 1;
            model->event[i] = 1;
            for (k = evnt; k < model->nevent; ++k)
                model->eventpos[k] = eventpos[k];
            model->reject += 1;
        }

    } else {
        // Add an event. Proposal ratio is (K+1) / (N - K), where K = nevent and N = nnode-1
        q = log((double)(model->nevent + 1) / (double)(model->tree->nnode - model->nevent - 1));

        p = log(PRIOR0 / (model->nevent + 1));

        model->eventpos[model->nevent] = i;
        model->nevent += 1;
        model->event[i] = 1;

        loglk = model_loglk(model);


        if (unif_rand() < exp(loglk - model->loglk + p + q)) {
            // accept
            model->loglk = loglk;
            model->accept += 1;
        } else {
            //reject
            model->nevent -= 1;
            model->event[i] = 0;
            model->reject += 1;
        }
    }
}


static void propose(struct model *model)
{
    propose0(model);
}

#endif
