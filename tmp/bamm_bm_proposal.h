#ifndef BAMM_BM_PROPOSAL_H
#define BAMM_BM_PROPOSAL_H

#include <R.h>
#include <Rmath.h>
#include "node.h"
#include "tree.h"
#include "bamm_bm_model.h"
#include "bamm_bm_likelihood.h"

/*
** Move types
**   within dimensional moves
**     shift trait    -- adjust a node state by a small amount
**     shift rate     -- adjust the overall rate of evolution by a small amount
**     shift multipl  -- adjust rate multiplier for a local clock
**     shift jump     -- adjust displacement for jump
**   trans dimensional moves
**     add a local random clock
**     remove a local random clock
**     add a jump
**     remove a jump
*/


#define SHFTST 0
#define SHFTRT 1
#define SHFTMLT 2
#define SHFTJMP 3
#define RCLCK 4
#define RJMP 5

static int PROPOSAL;
static double PWEIGHT[6];
static double TUNE0;    // Tuning param for adjusting node state (sd of normal distribution)
static double TUNE1;    // Tuning param for adjusting global rate of evolution
static double TUNE2;    // Tuning param for adjusting rate multipliers
static double TUNE3;    // Tuning param for adjusting jump sizes
static double PRIOR2;   // Gamma prior of rate multipliers
static double PRIOR3;   // Variance of normal prior on jump size
static double PRIOR4;   // Poisson prior on number of local clocks
static double PRIOR5;   // Poisson prior on number of jumps



static int choose_proposal()
{
    int i;
    double w = unif_rand();
    for (i = 0; i < 6; ++i) {
        if (w < PWEIGHT[i])
            return i;
    }
    return 5;
}


// Adjust a node state
static void propose0(struct model *model)
{
    int i;
    double x;
    double q;
    double p;
    double norm;
    double loglk;
    struct node *node;

    norm = model->treelen / model->treescal;

    // choose random node
    i = (int)(model->tree->nnode * unif_rand());

    node = model->tree->node[i];
    x = model->x[i];
    loglk = model->loglk;
    if (node->degree) {
        loglk -= branch_loglk(norm, node->children[0], model);
        loglk -= branch_loglk(norm, node->children[1], model);
        if (node != model->tree->root)
            loglk -= branch_loglk(norm, node, model);
        model->x[i] += rnorm(0, TUNE0);
        loglk += branch_loglk(norm, node->children[0], model);
        loglk += branch_loglk(norm, node->children[1], model);
        if (node != model->tree->root)
            loglk += branch_loglk(norm, node, model);
    } else {
        loglk -= branch_loglk(norm, node, model);
        model->x[i] += rnorm(0, TUNE0);
        loglk += branch_loglk(norm, node, model);
    }

    p = q = 0;
    if (unif_rand() < exp(loglk - model->loglk + p + q)) {
        model->loglk = loglk;
        model->accept += 1;
    } else {
        model->x[i] = x;
        model->reject += 1;
    }
}


// Perturb the overall rate of evolution by a small amount.
static void propose1(struct model *model)
{
    double q;
    double p;
    double U;
    double rate;
    double loglk;

    U = runif(TUNE1, 1/TUNE1);
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
    double loglk;
    struct node *node;

    if (model->nclck) {
        i = (int)(model->nclck * unif_rand());

        r = model->clckrate[i];
        j = model->clckpos[i];
        node = model->tree->node[j];

        U = runif(TUNE2, 1/TUNE2);
        q = -log(U);
        adjr = r * U;
        p = (PRIOR2-1)*log(adjr / r) - (PRIOR2) * (adjr - r);

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


// Adjust size of jump
static void propose3(struct model *model)
{
    int i;
    int k;
    double r;
    double adjr;
    double q;       // proposal ratio
    double p;       // prior ratio
    double norm;
    double loglk;
    struct node *node;

    if (model->njmp) {

        norm = model->treelen / model->treescal;

        loglk = model->loglk;

        i = (int)(model->njmp * unif_rand());

        k = model->jmppos[i];
        r = model->displ[k];
        adjr = r + rnorm(0, TUNE3);
        node = model->tree->node[k];

        loglk -= branch_loglk(norm, node, model);
        model->displ[k] = adjr;
        loglk += branch_loglk(norm, node, model);

        p = (1 / (2 * PRIOR3)) * (r - adjr) * (r - adjr);
        q = 0;

        if (unif_rand() < exp(loglk - model->loglk + p + q)) {
            model->loglk = loglk;
            model->accept += 1;
        } else {
            model->displ[k] = r;
            model->reject += 1;
        }
    }
}


// Add/remove a random local clock
static void propose4(struct model *model)
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
        p = log(model->nclck / PRIOR4);

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
        p = log(PRIOR4 / (model->nclck + 1));
        // Draw new rate multiplier from prior gamma distribution with mean 1.
        r = rgamma(PRIOR2, 1/PRIOR2);

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


// Add/remove a jump
static void propose5(struct model *model)
{
    int i;
    int j;
    int k;
    int jmp;
    int jmppos[model->njmp];
    double q;
    double p;
    double norm;
    double loglk;
    double jmpsz;
    struct node *node;

    i = (int)((model->tree->nnode - 1) * unif_rand());

    if (i >= model->tree->root->index)
        i += 1;

    node = model->tree->node[i];
    norm = model->treelen / model->treescal;

    if (model->jump[i]) {
        // Remove the jump. Proposal ratio is (N-K+1) / K,
        // where N is nbr branches and K is nbr jumps
        for (j = 0; j < model->njmp; ++j)
            if (model->jmppos[j] == i)
                jmp = j;

        jmpsz = model->displ[i];

        jmppos[jmp] = model->jmppos[jmp];

        for (k = jmp; k < (model->njmp - 1); ++k) {
            jmppos[k+1] = model->jmppos[k+1];
            model->jmppos[k] = model->jmppos[k+1];
        }

        loglk = model->loglk;
        loglk -= branch_loglk(norm, node, model);
        model->displ[i] = 0;
        loglk += branch_loglk(norm, node, model);

        q = log((double)(model->tree->nnode - model->njmp) / (double)model->njmp);
        // Prior ratio
        p = log(model->njmp / PRIOR5);

        if (unif_rand() < exp(loglk - model->loglk + p + q)) {
            // accept
            model->loglk = loglk;
            model->jump[i] = 0;
            model->njmp -= 1;
            model->accept += 1;
        } else {
            //reject
            model->displ[i] = jmpsz;
            for (k = jmp; k < model->njmp; ++k)
                model->jmppos[k] = jmppos[k];
            model->reject += 1;
        }
    } else {
        // Add a clock. Proposal ratio is (K+1) / (N - K), where K = njmp and N = nnode-1
        q = log((double)(model->njmp + 1) / (double)(model->tree->nnode - model->njmp - 1));
        // Prior ratio
        p = log(PRIOR5 / (model->njmp + 1));
        // Draw new jump from prior normal distribution with mean 0.

        loglk = model->loglk;
        loglk -= branch_loglk(norm, node, model);
        model->displ[i] = rnorm(0, sqrt(PRIOR3));
        loglk += branch_loglk(norm, node, model);

        if (unif_rand() < exp(loglk - model->loglk + p + q)) {
            // accept
            model->loglk = loglk;
            model->jmppos[model->njmp] = i;
            model->jump[i] = 1;
            model->njmp += 1;
            model->accept += 1;
        } else {
            //reject
            model->displ[i] = 0;
            model->reject += 1;
        }
    }

}


static void propose(struct model *model)
{
    PROPOSAL = choose_proposal();
    switch (PROPOSAL) {
        case SHFTST:
            propose0(model);
            break;
        case SHFTRT:
            propose1(model);
            break;
        case SHFTMLT:
            propose2(model);
            break;
        case SHFTJMP:
            propose3(model);
            break;
        case RCLCK:
            propose4(model);
            break;
        case RJMP:
            propose5(model);
            break;
    }
}

#endif
