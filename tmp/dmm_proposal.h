#ifndef DMM_PROPOSAL_H
#define DMM_PROPOSAL_H

#include <R.h>
#include <Rmath.h>
#include "node.h"
#include "tree.h"
#include "dmm_model.h"
#include "dmm_likelihood.h"

/*
** Move types
**   within dimensional moves
**     shift multipl  -- adjust rate multiplier for a local clock
**   trans dimensional moves
**     add a local random clock
**     remove a local random clock
**     add an event of character change
**     remove an event of character change
*/


#define SHFTMLT 2
#define RCLCK 1
#define REVNT 0

static int PROPOSAL;
static double PWEIGHT[3];
static double TUNE2;    // Tuning param for adjusting rate multipliers
static double PRIOR2;   // Gamma prior of rate multipliers
static double PRIOR1;   // Poisson prior on number of local clocks


static int choose_proposal()
{
    int i;
    double w = unif_rand();
    for (i = 0; i < 3; ++i) {
        if (w < PWEIGHT[i])
            return i;
    }
    return 2;
}


/*
** Set the rate multiplier for a branch and adjust relative branch
** rates accordingly. Return the difference in the sum of the product
** of branch lengths and relative rates for the affected subtree under
** the new and old configuration of relative rates to update the
** normalization constant.
*/
static double set_ratemultipl(double r, struct node *node, struct model *model)
{
    double opsum;
    double npsum;
    struct tree_traversal t;
    struct node *nd;

    t = tree_traverse(ALL_NODES, PREORDER, node, model->tree);
    nd = tree_step(&t);

    opsum = nd->brlen * model->rltvrate[nd->index];

    model->ratemultipl[nd->index] = r;
    model->rltvrate[nd->index] = r * model->rltvrate[nd->parent->index];

    npsum = nd->brlen * model->rltvrate[nd->index];

    nd = tree_step(&t);

    while (nd) {
        opsum += nd->brlen * model->rltvrate[nd->index];
        model->rltvrate[nd->index] = model->ratemultipl[nd->index] * model->rltvrate[nd->parent->index];
        npsum += nd->brlen * model->rltvrate[nd->index];
        nd = tree_step(&t);
    }

    return npsum - opsum;
}


// Adjust the rate of a random local clock
static void propose2(struct model *model)
{
    int i;
    int j;
    double r;
    double adjr;
    double q;       // proposal ratio
    double p;       // prior ratio
    double U;
    double oscal;
    double logprior;
    struct node *node;

    if (model->nclck) {
        i = (int)(model->nclck * unif_rand());

        r = model->clckrate[i];
        j = model->clckpos[i];
        node = model->tree->node[j];

        U = runif(TUNE2, 1/TUNE2);
        q = -log(U);
        adjr = r * U;

        oscal = model->treescal;
        model->treescal = oscal + set_ratemultipl(adjr, node, model);
        logprior = model_logprior(model);

        // log prior ratio on event config plus log prior ratio on rate multiplier
        p = logprior - model->logprior + (PRIOR2-1)*log(adjr / r) - (PRIOR2) * (adjr - r);

        if (unif_rand() < exp(p + q)) {
            // accept
            model->logprior = logprior;
            model->clckrate[i] = adjr;
            model->accept += 1;
        } else {
            set_ratemultipl(r, node, model);
            model->treescal = oscal;
            model->reject += 1;
        }
    }
}


// Add/remove a random local clock
static void propose1(struct model *model)
{
    int i;
    int j;
    int k;
    int clck;
    int clckpos[model->nclck];
    double q;
    double p;
    double r;
    double oscal;
    double logprior;
    double clckrate[model->nclck];
    struct node *node;

    i = (int)((model->tree->nnode - 1) * unif_rand());

    if (i >= model->tree->root->index)
        i += 1;

    node = model->tree->node[i];

    oscal = model->treescal;

    if (model->rlclck[i]) {
        // Remove the clock. Proposal ratio is (N-K+1) / K,
        // where N is nbr branches and K is nbr clocks

        for (j = 0; j < model->nclck; ++j)
            if (model->clckpos[j] == i)
                clck = j;

        clckrate[clck] = model->clckrate[clck];
        clckpos[clck] = model->clckpos[clck];

        for (k = clck; k < (model->nclck - 1); ++k) {
            clckrate[k+1] = model->clckrate[k+1];
            clckpos[k+1] = model->clckpos[k+1];
            model->clckrate[k] = model->clckrate[k+1];
            model->clckpos[k] = model->clckpos[k+1];
        }


        q = log((double)(model->tree->nnode - model->nclck) / (double)model->nclck);

        r = model->ratemultipl[i];

        model->treescal = oscal + set_ratemultipl(1, node, model);
        logprior = model_logprior(model);

        // log prior ratio on event config plus log prior ratio on clock count
        p = logprior - model->logprior + log(model->nclck / PRIOR1);

        if (unif_rand() < exp(p + q)) {
            // accept
            model->logprior = logprior;
            model->rlclck[i] = 0;
            model->nclck -= 1;
            model->accept += 1;
        } else {
            //reject
            set_ratemultipl(r, node, model);
            model->treescal = oscal;
            for (k = clck; k < model->nclck; ++k) {
                model->clckrate[k] = clckrate[k];
                model->clckpos[k] = clckpos[k];
            }
            model->reject += 1;
        }

    } else {
        // Add a clock. Proposal ratio is (K+1) / (N - K), where K = nclck and N = nnode-1
        q = log((double)(model->nclck + 1) / (double)(model->tree->nnode - model->nclck - 1));

        // Draw new rate multiplier from prior gamma distribution with mean 1.
        r = rgamma(PRIOR2, 1/PRIOR2);

        model->treescal = oscal + set_ratemultipl(r, node, model);
        logprior = model_logprior(model);

        // log prior ratio on event config plus log prior ratio on clock count
        p = logprior - model->logprior + log(PRIOR1 / (model->nclck + 1));

        if (unif_rand() < exp(p + q)) {
            // accept
            model->logprior = logprior;
            model->rlclck[i] = 1;
            model->clckrate[model->nclck] = r;
            model->clckpos[model->nclck] = i;
            model->nclck += 1;
            model->accept += 1;
        } else {
            //reject
            set_ratemultipl(1, node, model);
            model->treescal = oscal;
            model->reject += 1;
        }
    }
}


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
    double logprior;
    double brlen;   // new
    double norm;    // new
    struct node *node;

    i = (int)((model->tree->nnode - 1) * unif_rand());

    if (i >= model->tree->root->index)
        i += 1;

    node = model->tree->node[i];
    norm = model->treelen / model->treescal;


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

        // new block
        logprior = model->logprior;
        logprior -= dnbinom(model->nevent, model->beta, 0.5, 1);
        logprior -= lgammafn(model->nevent + 1);
        brlen = norm * model->rltvrate[node->index] * node->brlen;
        p = brlen / model->treelen;
        k = model->event[node->index];
        logprior -= k * log(p) - lgammafn(k+1);
        // end

        model->nevent -= 1;
        model->event[i] = 0;

        loglk = model_loglk(model);
        //logprior = model_logprior(model);

        // new block
        logprior += dnbinom(model->nevent, model->beta, 0.5, 1);
        logprior += lgammafn(model->nevent + 1);
        k = model->event[node->index];
        logprior += k * log(p) - lgammafn(k+1);
        //end


        // log prior ratio on event config
        p = logprior - model->logprior;

        if (unif_rand() < exp(loglk - model->loglk + p + q)) {
            // accept
            model->loglk = loglk;
            model->logprior = logprior;
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

        // new block
        logprior = model->logprior;
        logprior -= dnbinom(model->nevent, model->beta, 0.5, 1);
        logprior -= lgammafn(model->nevent + 1);
        brlen = norm * model->rltvrate[node->index] * node->brlen;
        p = brlen / model->treelen;
        k = model->event[node->index];
        logprior -= k * log(p) - lgammafn(k+1);
        // end

        model->eventpos[model->nevent] = i;
        model->nevent += 1;
        model->event[i] = 1;

        loglk = model_loglk(model);
        //logprior = model_logprior(model);

        // new block
        logprior += dnbinom(model->nevent, model->beta, 0.5, 1);
        logprior += lgammafn(model->nevent + 1);
        k = model->event[node->index];
        logprior += k * log(p) - lgammafn(k+1);
        //end

        p = logprior - model->logprior;

        if (unif_rand() < exp(loglk - model->loglk + p + q)) {
            // accept
            model->loglk = loglk;
            model->logprior = logprior;
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
    PROPOSAL = choose_proposal();
    switch (PROPOSAL) {
        case REVNT:
            propose0(model);
            break;
        case RCLCK:
            propose1(model);
            break;
        case SHFTMLT:
            propose2(model);
            break;
    }
}

#endif
