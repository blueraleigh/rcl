#ifndef PROPOSAL_H
#define PROPOSAL_H

#include <R.h>
#include <Rmath.h>
#include "node.h"
#include "tree.h"
#include "bamm_model.h"
#include "bamm_likelihood.h"

/*
** Move types
**   within dimensional moves
**     shift shape    -- adjust the shape parameters of the beta distribution
**     shift rate     -- adjust the overall rate of evolution by a small amount
**     shift multipl  -- adjust rate multiplier for a local clock
**   trans dimensional moves
**     add a local random clock
**     remove a local random clock
*/

#define SHFTRT 0
#define SHFTMLT 1
#define RCLCK 2

static int PROPOSAL;
static double PWEIGHT[3];
static double TUNE0;    // Tuning param for adjusting global rate of evolution
static double TUNE1;    // Tuning param for adjusting rate multipliers
static double PRIOR1;   // Gamma prior of rate multipliers
static double PRIOR2;   // Poisson prior on number of local clocks



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


// Perturb the overall rate of evolution by a small amount.
static void propose0(struct model *model)
{
    double q;
    double p;
    double U;
    double rate;
    double loglk;

    U = runif(TUNE0, 1/TUNE0);
    rate = model->rate;
    model->rate = rate * U;
    loglk = model_loglk(model);
    q = -log(U);
    p = 0;

    if (unif_rand() < exp(loglk - model->loglk + p + q)) {
        model->loglk = loglk;
        model->accept += 1;
    } else {
        model->rate = rate;
        model->reject += 1;
    }

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
static void propose1(struct model *model)
{
    int i;
    int j;
    double r;
    double adjr;
    double q;       // proposal ratio
    double p;       // prior ratio
    double U;
    double oscal;
    double loglk;
    struct node *node;

    if (model->nclck) {
        i = (int)(model->nclck * unif_rand());

        r = model->clckrate[i];
        j = model->clckpos[i];
        node = model->tree->node[j];

        U = runif(TUNE1, 1/TUNE1);
        q = -log(U);
        adjr = r * U;
        p = (PRIOR1-1)*log(adjr / r) - (PRIOR1) * (adjr - r);

        oscal = model->treescal;
        model->treescal = oscal + set_ratemultipl(adjr, node, model);
        loglk = model_loglk(model);

        if (unif_rand() < exp(loglk - model->loglk + p + q)) {
            // accept
            model->loglk = loglk;
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
static void propose2(struct model *model)
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
    double loglk;
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
        // Prior ratio
        p = log(model->nclck / PRIOR2);

        r = model->ratemultipl[i];

        model->treescal = oscal + set_ratemultipl(1, node, model);
        loglk = model_loglk(model);

        if (unif_rand() < exp(loglk - model->loglk + p + q)) {
            // accept
            model->loglk = loglk;
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
        // Prior ratio
        p = log(PRIOR2 / (model->nclck + 1));
        // Draw new rate multiplier from prior gamma distribution with mean 1.
        r = rgamma(PRIOR1, 1/PRIOR1);

        model->treescal = oscal + set_ratemultipl(r, node, model);
        loglk = model_loglk(model);

        if (unif_rand() < exp(loglk - model->loglk + p + q)) {
            // accept
            model->loglk = loglk;
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


static void propose(struct model *model)
{
    PROPOSAL = choose_proposal();
    switch (PROPOSAL) {
        case SHFTRT:
            propose0(model);
            break;
        case SHFTMLT:
            propose1(model);
            break;
        case RCLCK:
            propose2(model);
            break;
    }
}

#endif
