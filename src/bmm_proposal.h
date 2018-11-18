#ifndef BMM_PROPOSAL_H
#define BMM_PROPOSAL_H

#include <R.h>
#include <Rmath.h>
#include "node.h"
#include "tree.h"
#include "bmm_model.h"

/*
** Move types
**   within dimensional moves
**     adjust the rate of a local clock or the global rate
**     change size of jump
**   trans dimensional moves
**     add a local random clock
**     remove a local random clock
**     add a jump
**     remove a jump
*/


#define CLOCKNBR 0
#define EVORATE 1
#define BRNCHJMP 2
#define TUNEJMP 3

static int PROPOSAL;
static double PWEIGHT[4];
static double TUNE0;    // Tuning param for adjusting global rate
static double TUNE1;    // Tuning param for adjusting rate multipliers
static double PRIOR0;   // Poisson prior on number of local clocks
static double PRIOR1;   // Gamma prior of rate multipliers
static double PRIOR2;   // Poisson prior on number of jumps
static double PRIOR3;   // Normal prior on variance of jumps
static double TUNE3;    // Tuning param for adjusting jump sizes


static int choose_proposal()
{
    int i;
    double w = unif_rand();
    for (i = 0; i < 4; ++i) {
        if (w < PWEIGHT[i])
            return i;
    }
    return 3;
}


/*
** Set the rate multiplier for the random local clock that initiates at node
** and adjust relative branch rates accordingly. Return the difference in
** the sum of the product of branch lengths and relative rates for the affected
** subtree under the new and old configuration of relative rates to update the
** normalization constant.
*/
static double set_ratemultipl(int clockid, double r, struct node *node, struct model *model)
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
    model->clockid[nd->index] = clockid;

    npsum = nd->brlen * model->rltvrate[nd->index];

    nd = tree_step(&t);

    while (nd) {
        model->clockid[nd->index] = clockid;
        opsum += nd->brlen * model->rltvrate[nd->index];
        model->rltvrate[nd->index] = model->ratemultipl[nd->index] * model->rltvrate[nd->parent->index];
        npsum += nd->brlen * model->rltvrate[nd->index];
        nd = tree_step(&t);
    }

    return npsum - opsum;
}


static void add_clock(double r, struct node *node, struct model *model)
{
    int newid;

    newid = model->nextclockid[--(model->nitem)];

    model->clockpos[newid] = node->index;
    model->nclock += 1;
    model->treescal += set_ratemultipl(newid, r, node, model);
    model->clock[node->index] = 1;
}


static void delete_clock(struct node *node, struct model *model)
{
    int oldid;
    int newid;

    oldid = model->clockid[node->index];
    newid = model->clockid[node->parent->index];

    model->nextclockid[(model->nitem)++] = oldid;
    model->nclock -= 1;
    model->treescal += set_ratemultipl(newid, 1, node, model);
    model->clock[node->index] = 0;
}


static void update_clock(double r, struct node *node, struct model *model)
{
    model->treescal += set_ratemultipl(model->clockid[node->index], r, node, model);
}


static void set_jump(double r, struct node *node, struct model *model)
{
    struct tree_traversal t;
    struct node *nd;

    t = tree_traverse(ALL_NODES, PREORDER, node, model->tree);

    nd = tree_step(&t);

    while (nd) {
        model->mu[nd->index] = r;
        nd = tree_step(&t);
        if (nd && model->jump[nd->index])
            nd = tree_jump(nd, &t);
    }
}


static void add_jump(double r, struct node *node, struct model *model)
{
    model->njump += 1;
    set_jump(r, node, model);
    model->jump[node->index] = 1;
}


static void delete_jump(struct node *node, struct model *model)
{
    model->njump -= 1;
    set_jump(0, node, model);
    model->jump[node->index] = 0;
}


static void update_jump(double r, struct node *node, struct model *model)
{
    set_jump(r, node, model);
}


// Add/remove a jump
static void propose2(struct model *model)
{
    int i;
    int j;
    int k;
    double q;
    double p;
    double r;
    double loglk;
    struct node *node;

    i = (int)((model->tree->nnode - 1) * unif_rand());

    if (i >= model->tree->root->index)
        i += 1;

    node = model->tree->node[i];


    if (model->jump[i]) {
        // Remove the jump. Proposal ratio is (N-K+1) / K,
        // where N is nbr branches and K is nbr jumps

        // log prior ratio jump count
        p = log(model->njump / PRIOR2);

        r = model->mu[node->index];

        q = log((double)(model->tree->nnode - model->njump) / (double)model->njump);

        delete_jump(node, model);

        // compute likelihood
        loglk = model_loglk(model);

        if (unif_rand() < exp(loglk - model->loglk + p + q)) {
            // accept
            model->loglk = loglk;
            model->accept += 1;
        } else {
            //reject
            add_jump(r, node, model);
            model->reject += 1;
        }

    } else {
        // Add a jump. Proposal ratio is (K+1) / (N - K), where K = njump and N = nnode-1
        q = log((double)(model->njump + 1) / (double)(model->tree->nnode - model->njump - 1));

        // log prior ratio on jump count
        p = log(PRIOR2 / (model->njump + 1));

        // Draw new jump from prior normal distribution with mean 0.
        r = rnorm(0, sqrt(PRIOR3));

        add_jump(r, node, model);
        // compute likelihood
        loglk = model_loglk(model);

        if (unif_rand() < exp(loglk - model->loglk + p + q)) {
            // accept
            model->loglk = loglk;
            model->accept += 1;
        } else {
            //reject
            delete_jump(node, model);
            model->reject += 1;
        }
    }
}


// Adjust the size of a jump
static void propose3(struct model *model)
{
    int i;
    double r;
    double p;       // prior ratio
    double U;
    double loglk;
    struct node *node;

    i = (int)((model->tree->nnode - 1) * unif_rand());

    if (i >= model->tree->root->index)
        i += 1;

    if (model->jump[i]) {
        node = model->tree->node[i];

        r = model->mu[node->index];

        U = rnorm(0, TUNE3);

        // log prior ratio on jump size
        p = -0.5 * PRIOR3 * (U * (2 * r + U));

        update_jump(r + U, node, model);

        loglk = model_loglk(model);

        if (unif_rand() < exp(loglk - model->loglk + p)) {
            // accept
            model->loglk = loglk;
            model->accept += 1;
        } else {
            update_jump(r, node, model);
            model->reject += 1;
        }
    }
}

// Adjust the overall rate or the rate of a random local clock
static void propose1(struct model *model)
{
    int i;
    int clockid;
    double r;
    double q;       // proposal ratio
    double p;       // prior ratio
    double U;
    double loglk;
    struct node *node;

    i = (int)((model->tree->nnode - 1) * unif_rand());

    if (i >= model->tree->root->index)
        i += 1;

    clockid = model->clockid[i];

    // clockid = 0 is reserved for clock at root, in which case we update
    // the global rate of evolution
    if (clockid) {
        node = model->tree->node[model->clockpos[clockid]];

        r = model->ratemultipl[node->index];

        U = runif(TUNE1, 1/TUNE1);
        q = -log(U);

        update_clock(r * U, node, model);

        // log prior ratio on rate multiplier
        p = (PRIOR1-1)*log(U) - (PRIOR1) * (r * (U - 1));

        // compute likelihood
        loglk = model_loglk(model);

        if (unif_rand() < exp(loglk - model->loglk + p + q)) {
            // accept
            model->loglk = loglk;
            model->accept += 1;
        } else {
            update_clock(r, node, model);
            model->reject += 1;
        }
    } else {
        r = model->rate;
        U = runif(TUNE1, 1/TUNE1);
        q = -log(U);

        model->rate *= U;

        loglk = model_loglk(model);

        if (unif_rand() < exp(loglk - model->loglk + q)) {
            // accept
            model->loglk = loglk;
            model->accept += 1;
        } else {
            model->rate = r;
            model->reject += 1;
        }
    }
}


// Add/remove a random local clock
static void propose0(struct model *model)
{
    int i;
    int j;
    int k;
    double q;
    double p;
    double r;
    double loglk;
    struct node *node;

    i = (int)((model->tree->nnode - 1) * unif_rand());

    if (i >= model->tree->root->index)
        i += 1;

    node = model->tree->node[i];


    if (model->clock[i]) {
        // Remove the clock. Proposal ratio is (N-K+1) / K,
        // where N is nbr branches and K is nbr clocks

        // log prior ratio clock count
        p = log(model->nclock / PRIOR0);

        r = model->ratemultipl[node->index];

        delete_clock(node, model);

        q = log((double)(model->tree->nnode - model->nclock) / (double)model->nclock);

        // compute likelihood
        loglk = model_loglk(model);

        if (unif_rand() < exp(loglk - model->loglk + p + q)) {
            // accept
            model->loglk = loglk;
            model->accept += 1;
        } else {
            //reject
            add_clock(r, node, model);
            model->reject += 1;
        }

    } else {
        // Add a clock. Proposal ratio is (K+1) / (N - K), where K = nclock and N = nnode-1
        q = log((double)(model->nclock + 1) / (double)(model->tree->nnode - model->nclock - 1));

        // log prior ratio on clock count
        p = log(PRIOR0 / (model->nclock + 1));

        // Draw new rate multiplier from prior gamma distribution with mean 1.
        r = rgamma(PRIOR1, 1/PRIOR1);

        add_clock(r, node, model);

        // compute likelihood
        loglk = model_loglk(model);

        if (unif_rand() < exp(loglk - model->loglk + p + q)) {
            // accept
            model->loglk = loglk;
            model->accept += 1;
        } else {
            //reject
            delete_clock(node, model);
            model->reject += 1;
        }
    }
}


static void propose(struct model *model)
{
    PROPOSAL = choose_proposal();
    switch (PROPOSAL) {
        case CLOCKNBR:
            propose0(model);
            break;
        case EVORATE:
            propose1(model);
            break;
        case BRNCHJMP:
            propose2(model);
            break;
        case TUNEJMP:
            propose3(model);
            break;
    }
}

#endif
