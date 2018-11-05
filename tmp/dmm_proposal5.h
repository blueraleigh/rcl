#ifndef DMM_PROPOSAL5_H
#define DMM_PROPOSAL5_H

#include <R.h>
#include <Rmath.h>
#include "node.h"
#include "tree.h"
#include "dmm_model5.h"

/*
** Move types
**   within dimensional moves
**     adjust the rate of a local clock or the global rate
**     add or delete an event of evolutionary change
**   trans dimensional moves
**     add a local random clock
**     remove a local random clock
*/


#define EVENTNBR 0
#define EVENTMOVE 1
#define CLOCKNBR 2
#define EVORATE 3

static int PROPOSAL;
static double PWEIGHT[4];
static double TUNE1;    // Tuning param for adjusting global rate
static double TUNE2;    // Tuning param for adjusting rate multipliers
static double PRIOR1;   // Poisson prior on number of local clocks
static double PRIOR2;   // Gamma prior of rate multipliers


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

    newid = model->nextclockid[--(model->nitem2)];

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

    model->nextclockid[(model->nitem2)++] = oldid;
    model->nclock -= 1;
    model->treescal += set_ratemultipl(newid, 1, node, model);
    model->clock[node->index] = 0;
}


static void update_clock(double r, struct node *node, struct model *model)
{
    model->treescal += set_ratemultipl(model->clockid[node->index], r, node, model);
}


// Adjust the overall rate or the rate of a random local clock
static void propose3(struct model *model)
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
static void propose2(struct model *model)
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


static void add_event(struct node *node, struct model *model)
{
    int i;
    int n1;
    int n2;
    int oldid;
    int newid;
    struct tree_traversal t;

    model->nevent += 1;
    model->event[node->index] = 1;

    if (PRIORONLY)
        return;

    oldid = model->eventid[node->index];
    model->dataloglk -= model->eventloglk[oldid];

    newid = model->nexteventid[--(model->nitem1)];

    t = tree_traverse(ALL_NODES, PREORDER, node, model->tree);

    node = tree_step(&t);

    while (node) {

        if (!node->degree) {

            for (i = 0; i < model->ncat; ++i) {
                model->count[i + newid * model->ncat] += model->data[i + node->index * model->ncat];
                model->count[i + oldid * model->ncat] -= model->data[i + node->index * model->ncat];
            }

        }

        model->eventid[node->index] = newid;

        node = tree_step(&t);

        while (node && model->event[node->index])
            node = tree_jump(node, &t);
    }

    model->eventloglk[newid] = 0;
    model->eventloglk[oldid] = 0;

    n1 = 0, n2 = 0;
    for (i = 0; i < model->ncat; ++i) {
        n1 += model->count[i + newid * model->ncat];
        n2 += model->count[i + oldid * model->ncat];
        model->eventloglk[newid] += lgammafn(
            model->count[i + newid * model->ncat] + model->alpha) - lgammafn(model->alpha);
        model->eventloglk[oldid] += lgammafn(
            model->count[i + oldid * model->ncat] + model->alpha) - lgammafn(model->alpha);
    }

    model->eventloglk[newid] += lgammafn(model->ncat * model->alpha) - lgammafn(n1 + model->ncat * model->alpha);
    model->eventloglk[oldid] += lgammafn(model->ncat * model->alpha) - lgammafn(n2 + model->ncat * model->alpha);

    model->dataloglk += model->eventloglk[newid];
    model->dataloglk += model->eventloglk[oldid];
}


static void delete_event(struct node *node, struct model *model)
{
    int i;
    int n;
    int oldid;
    int newid;
    struct tree_traversal t;

    model->nevent -= 1;
    model->event[node->index] = 0;

    if (PRIORONLY)
        return;

    oldid = model->eventid[node->index];
    newid = model->eventid[node->parent->index];

    model->nexteventid[(model->nitem1)++] = oldid;

    t = tree_traverse(ALL_NODES, PREORDER, node, model->tree);

    node = tree_step(&t);

    model->dataloglk -= model->eventloglk[oldid];
    model->dataloglk -= model->eventloglk[newid];

    while (node) {

        if (!node->degree) {

            for (i = 0; i < model->ncat; ++i)
                model->count[i + newid * model->ncat] += model->data[i + node->index * model->ncat];

        }

        model->eventid[node->index] = newid;

        node = tree_step(&t);

        while (node && model->event[node->index])
            node = tree_jump(node, &t);

    }

    model->eventloglk[newid] = 0;
    model->eventloglk[oldid] = 0;

    n = 0;
    for (i = 0; i < model->ncat; ++i) {
        n += model->count[i + newid * model->ncat];
        model->count[i + oldid * model->ncat] = 0;
        model->eventloglk[newid] += lgammafn(
            model->count[i + newid * model->ncat] + model->alpha) - lgammafn(model->alpha);
    }

    model->eventloglk[newid] += lgammafn(model->ncat * model->alpha) - lgammafn(n + model->ncat * model->alpha);

    model->dataloglk += model->eventloglk[newid];
}


// Add/remove an event of character change
static void propose0(struct model *model)
{
    int i;
    double loglk;
    double branchloglk;
    double treeloglk;
    struct node *node;

    i = (int)((model->tree->nnode - 1) * unif_rand());

    if (i >= model->tree->root->index)
        ++i;

    node = model->tree->node[i];

    if (model->event[i]) {
        // Delete the event

        branchloglk = model_branchloglk(node, model);
        delete_event(node, model);

        treeloglk = model->treeloglk - branchloglk + model_branchloglk(node, model);
        loglk = model->dataloglk + model->treeloglk;

        if (unif_rand() < exp(loglk - model->loglk)) {
            model->treeloglk = treeloglk;
            model->loglk = loglk;
            model->accept += 1;
        } else {
            add_event(node, model);
            model->reject += 1;
        }

    } else {
        // Add an event

        branchloglk = model_branchloglk(node, model);
        add_event(node, model);

        treeloglk = model->treeloglk - branchloglk + model_branchloglk(node, model);
        loglk = model->dataloglk + model->treeloglk;

        if (unif_rand() < exp(loglk - model->loglk)) {
            model->treeloglk = treeloglk;
            model->loglk = loglk;
            model->accept += 1;
        } else {
            delete_event(node, model);
            model->reject += 1;
        }

    }
}


// Move (global) an event of character change
static void propose1(struct model *model)
{
    int i;
    int j;
    double branchloglk;
    double treeloglk;
    double loglk;
    struct node *node1;
    struct node *node2;

    i = (int)((model->tree->nnode - 1) * unif_rand());
    j = (int)((model->tree->nnode - 1) * unif_rand());

    if (i >= model->tree->root->index)
        ++i;
    if (j >= model->tree->root->index)
        ++j;

    node1 = model->tree->node[i];
    node2 = model->tree->node[j];

    if (model->event[i] && !model->event[j]) {
        branchloglk = model_branchloglk(node1, model) + model_branchloglk(node2, model);

        delete_event(node1, model);
        add_event(node2, model);

        treeloglk = model->treeloglk - branchloglk + model_branchloglk(node1, model) + model_branchloglk(node2, model);
        loglk = model->dataloglk + model->treeloglk;

        if (unif_rand() < exp(loglk - model->loglk)) {
            model->treeloglk = treeloglk;
            model->loglk = loglk;
            model->accept += 1;
        } else {
            delete_event(node2, model);
            add_event(node1, model);
            model->reject += 1;
        }

    } else if (model->event[j] && !model->event[i]) {
        branchloglk = model_branchloglk(node1, model) + model_branchloglk(node2, model);

        delete_event(node2, model);
        add_event(node1, model);

        treeloglk = model->treeloglk - branchloglk + model_branchloglk(node1, model) + model_branchloglk(node2, model);
        loglk = model->dataloglk + model->treeloglk;

        if (unif_rand() < exp(loglk - model->loglk)) {
            model->treeloglk = treeloglk;
            model->loglk = loglk;
            model->accept += 1;
        } else {
            delete_event(node1, model);
            add_event(node2, model);
            model->reject += 1;
        }
    }
}


static void propose(struct model *model)
{
    PROPOSAL = choose_proposal();
    switch (PROPOSAL) {
        case EVENTNBR:
            propose0(model);
            break;
        case EVENTMOVE:
            propose1(model);
            break;
        case CLOCKNBR:
            propose2(model);
            break;
        case EVORATE:
            propose3(model);
            break;
    }
}

#endif
