#ifndef DMM_PROPOSAL_H
#define DMM_PROPOSAL_H

#include <R.h>
#include <Rmath.h>
#include "node.h"
#include "tree.h"
#include "dmm_model.h"

/*
** Move types
**   within dimensional moves
**     adjust the rate of a local clock or the global rate
**     change the state for a node
**   trans dimensional moves
**     add a local random clock
**     remove a local random clock
*/


#define CHNGSTATE 0
#define CLOCKNBR 1
#define EVORATE 2

static int PROPOSAL;
static double PWEIGHT[3];
static double TUNE1;    // Tuning param for adjusting global rate
static double TUNE2;    // Tuning param for adjusting rate multipliers
static double PRIOR1;   // Poisson prior on number of local clocks
static double PRIOR2;   // Gamma prior of rate multipliers


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


// Adjust the overall rate or the rate of a random local clock
static void propose2(struct model *model)
{
    int i;
    int clockid;
    double r;
    double q;       // proposal ratio
    double p;       // prior ratio
    double U;
    double loglk;
    double treeloglk;
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

        U = runif(TUNE2, 1/TUNE2);
        q = -log(U);

        update_clock(r * U, node, model);

        // log prior ratio on rate multiplier
        p = (PRIOR2-1)*log(U) - (PRIOR2) * (r * (U - 1));

        // compute likelihood
        treeloglk = model_treeloglk(model);
        loglk = model->loglk - model->treeloglk + treeloglk;

        if (unif_rand() < exp(loglk - model->loglk + p + q)) {
            // accept
            model->treeloglk = treeloglk;
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

        treeloglk = model_treeloglk(model);
        loglk = model->loglk - model->treeloglk + treeloglk;

        if (unif_rand() < exp(loglk - model->loglk + q)) {
            // accept
            model->treeloglk = treeloglk;
            model->loglk = loglk;
            model->accept += 1;
        } else {
            model->rate = r;
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
    double q;
    double p;
    double r;
    double loglk;
    double treeloglk;
    struct node *node;

    i = (int)((model->tree->nnode - 1) * unif_rand());

    if (i >= model->tree->root->index)
        i += 1;

    node = model->tree->node[i];


    if (model->clock[i]) {
        // Remove the clock. Proposal ratio is (N-K+1) / K,
        // where N is nbr branches and K is nbr clocks

        // log prior ratio clock count
        p = log(model->nclock / PRIOR1);

        r = model->ratemultipl[node->index];

        delete_clock(node, model);

        q = log((double)(model->tree->nnode - model->nclock) / (double)model->nclock);

        // compute likelihood
        treeloglk = model_treeloglk(model);
        loglk = model->loglk - model->treeloglk + treeloglk;

        if (unif_rand() < exp(loglk - model->loglk + p + q)) {
            // accept
            model->treeloglk = treeloglk;
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
        p = log(PRIOR1 / (model->nclock + 1));

        // Draw new rate multiplier from prior gamma distribution with mean 1.
        r = rgamma(PRIOR2, 1/PRIOR2);

        add_clock(r, node, model);

        // compute likelihood
        treeloglk = model_treeloglk(model);
        loglk = model->loglk - model->treeloglk + treeloglk;

        if (unif_rand() < exp(loglk - model->loglk + p + q)) {
            // accept
            model->treeloglk = treeloglk;
            model->loglk = loglk;
            model->accept += 1;
        } else {
            //reject
            delete_clock(node, model);
            model->reject += 1;
        }
    }
}


static void change_state(int old, int new, struct node *node, struct model *model)
{
    int i;
    int n = 0;
    int m = 0;

    if (node->parent)
        model->treeloglk -= model_branchloglk(node, model);

    if (node->degree) {
        model->treeloglk -= model_branchloglk(node->children[0], model);
        model->treeloglk -= model_branchloglk(node->children[1], model);
    }

    if (!node->degree) {

        model->dataloglk -= model->stateloglk[old];
        model->dataloglk -= model->stateloglk[new];

        model->stateloglk[old] = 0;
        model->stateloglk[new] = 0;

        for (i = 0; i < model->ncat; ++i) {
            model->count[i + old * model->ncat] -= model->data[i + node->index * model->ncat];
            model->count[i + new * model->ncat] += model->data[i + node->index * model->ncat];

            n += model->count[i + old * model->ncat];
            model->stateloglk[old] += lgammafn(
                    model->count[i + old * model->ncat] + model->alpha) - lgammafn(model->alpha);

            m += model->count[i + new * model->ncat];
            model->stateloglk[new] += lgammafn(
                    model->count[i + new * model->ncat] + model->alpha) - lgammafn(model->alpha);
        }

        model->stateloglk[old] += lgammafn(model->ncat * model->alpha) - lgammafn(n + model->ncat * model->alpha);
        model->stateloglk[new] += lgammafn(model->ncat * model->alpha) - lgammafn(m + model->ncat * model->alpha);

        model->dataloglk += model->stateloglk[old];
        model->dataloglk += model->stateloglk[new];

    }

    model->stateid[node->index] = new;

    if (node->parent)
        model->treeloglk += model_branchloglk(node, model);

    if (node->degree) {
        model->treeloglk += model_branchloglk(node->children[0], model);
        model->treeloglk += model_branchloglk(node->children[1], model);
    }

    model->loglk = model->dataloglk + model->treeloglk;
}


// Change node state
static void propose0(struct model *model)
{
    int i;
    int old;
    int new;
    double loglk;
    struct node *node;

    i = (int)(model->tree->nnode * unif_rand());

    node = model->tree->node[i];

    old = model->stateid[node->index];

    new = (int)((model->nstate-1) * unif_rand());

    if (new >= old)
        ++new;

    loglk = model->loglk;

    change_state(old, new, node, model);

    if (unif_rand() < exp(model->loglk - loglk)) {
        model->accept += 1;
    } else {
        change_state(new, old, node, model);
        model->reject += 1;
    }
}


static void propose(struct model *model)
{
    PROPOSAL = choose_proposal();
    switch (PROPOSAL) {
        case CHNGSTATE:
            propose0(model);
            break;
        case CLOCKNBR:
            propose1(model);
            break;
        case EVORATE:
            propose2(model);
            break;
    }
}

#endif
