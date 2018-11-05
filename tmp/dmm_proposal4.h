#ifndef DMM_PROPOSAL4_H
#define DMM_PROPOSAL4_H

#include <R.h>
#include <Rmath.h>
#include "node.h"
#include "tree.h"
#include "dmm_model4.h"

static double PWEIGHT[3];

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

    model->logprior = dnbinom(model->nevent, model->beta, 0.5, 1);

    if (PRIORONLY)
        return;

    oldid = model->eventid[node->index];
    model->loglk -= model->eventloglk[oldid];

    if (!model->nitem)
        error("Max events reached!");
    else
        newid = model->nexteventid[--(model->nitem)];

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

    model->loglk += model->eventloglk[newid];
    model->loglk += model->eventloglk[oldid];
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

    model->logprior = dnbinom(model->nevent, model->beta, 0.5, 1);

    if (PRIORONLY)
        return;

    oldid = model->eventid[node->index];
    newid = model->eventid[node->parent->index];

    model->nexteventid[(model->nitem)++] = oldid;

    t = tree_traverse(ALL_NODES, PREORDER, node, model->tree);

    node = tree_step(&t);

    model->loglk -= model->eventloglk[oldid];
    model->loglk -= model->eventloglk[newid];

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

    model->loglk += model->eventloglk[newid];
}


static void propose_bd(struct model *model)
{
    int i;
    double q;
    double nbranch;
    double nevent;
    double loglk;
    double logprior;
    struct node *node;

    loglk = model->loglk;
    logprior = model->logprior;

    nevent = (double)model->nevent;
    nbranch = (double)model->tree->nnode - 1;

    i = (int)(nbranch * unif_rand());

    if (i >= model->tree->root->index)
        ++i;

    node = model->tree->node[i];

    if (model->event[i]) {
        // Delete the event
        delete_event(node, model);

        q = log((nbranch-nevent+1) / nevent);

        if (unif_rand() < exp(model->loglk - loglk + model->logprior - logprior + q)) {
            model->accept += 1;
        } else {
            add_event(node, model);
            model->reject += 1;
        }

    } else {
        // Add an event
        add_event(node, model);

        q = log((nevent+1) / (nbranch-nevent));

        if (unif_rand() < exp(model->loglk - loglk + model->logprior - logprior + q)) {
            model->accept += 1;
        } else {
            delete_event(node, model);
            model->reject += 1;
        }

    }
}


static void propose_swap(struct model *model)
{
    int i;
    int j;
    int nbranch;
    double loglk;
    struct node *node1;
    struct node *node2;

    loglk = model->loglk;

    nbranch = model->tree->nnode - 1;

    i = (int)(nbranch * unif_rand());
    j = (int)(nbranch * unif_rand());

    if (i >= model->tree->root->index)
        ++i;
    if (j >= model->tree->root->index)
        ++j;

    node1 = model->tree->node[i];
    node2 = model->tree->node[j];

    if (model->event[i] && !model->event[j]) {
        delete_event(node1, model);
        add_event(node2, model);

        if (unif_rand() < exp(model->loglk - loglk)) {
            model->accept += 1;
        } else {
            delete_event(node2, model);
            add_event(node1, model);
            model->reject += 1;
        }

    } else if (model->event[j] && !model->event[i]) {
        delete_event(node2, model);
        add_event(node1, model);

        if (unif_rand() < exp(model->loglk - loglk)) {
            model->accept += 1;
        } else {
            delete_event(node1, model);
            add_event(node2, model);
            model->reject += 1;
        }
    }
}


static void propose_alpha(struct model *model)
{
    double q;
    double U;
    double loglk;
    double alpha;

    U = runif(0.5, 2);
    q = -log(U);
    alpha = model->alpha;
    loglk = model->loglk;

    model_updatealpha(alpha*U, model);

    if (unif_rand() < exp(model->loglk - loglk + q)) {
        model->accept += 1;
    } else {
        model_updatealpha(alpha, model);
        model->reject += 1;
    }
}


static void propose(struct model *model)
{
    double U = unif_rand();
    if (U < PWEIGHT[0])
        propose_bd(model);
    else if (U < PWEIGHT[1])
        propose_swap(model);
    else
        propose_alpha(model);
}

#endif
