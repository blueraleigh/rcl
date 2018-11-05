#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include "node.h"
#include "tree.h"

struct model {

    /* Number of resource categories */
    int ncat;

    int *data;

    /* Log probability of data in each node's subtree under
    ** hypothesis that each observation drawn i.i.d. from
    ** same distribution */
    double *h0;

    /* Log likelihood of data in each node's subtree. */
    double *loglk;

    /* Log likelihood of full data */
    double LnLk;

    /* Sum of all branch lengths descended from each node */
    double *treelen;

    struct tree *tree;

    struct tree_traversal t;

    /* Dirichlet prior parameter */
    double alpha;

    /* Rate of Poisson process. */
    double rate;
};


static struct model *model_alloc(int ncat, int ntip, int nnode)
{
    struct model *model = Calloc(1, struct model);
    model->ncat = ncat;
    model->alpha = 0;
    model->rate = 0;
    model->data = Calloc(ncat * ntip, int);
    model->loglk = Calloc(nnode, double);
    model->h0 = Calloc(nnode, double);
    model->treelen = Calloc(nnode, double);
    return model;
}


static void model_free(struct model *model)
{
    Free(model->data);
    Free(model->loglk);
    Free(model->h0);
    Free(model->treelen);
    Free(model);
}


static void accumulate_counts(double *loglk, int *count, struct node *node, struct model* model)
{
    if (node->degree) {
        accumulate_counts(loglk, count, node->children[0], model);
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


/* Compute the log likelihood of the terminal nodes that descend
** from *node assuming their data is drawn i.i.d. from a single
** Dirichlet-multinomial distribution
*/
static double model_h0(struct node *node, struct model *model)
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


static void model_set(double alpha, double rate, int *data, struct tree *tree, struct model *model)
{
    memcpy(model->data, data, model->ncat * tree->ntip * sizeof(int));
    model->tree = tree;
    model->t = tree_traverse(POSTORDER, INTERNAL_NODES_ONLY, tree->root, tree);
    model->alpha = alpha;
    model->rate = rate;

    int i;
    for (i = 0; i < tree->nnode; ++i)
        model->h0[i] = model_h0(tree->node[i], model);

    for (i = 0; i < tree->ntip; ++i)
        model->loglk[i] = model->h0[i];

    struct node *node;
    memset(model->treelen, 0, tree->nnode * sizeof(double));
    node = tree_step(&(model->t));
    while (node) {
        model->treelen[node->index] += model->treelen[node->children[0]->index];
        model->treelen[node->index] += model->treelen[node->children[1]->index];
        model->treelen[node->index] += node->children[0]->brlen;
        model->treelen[node->index] += node->children[1]->brlen;
        node = tree_step(&(model->t));
    }

    tree_reset(&(model->t));
}


/* Log likelihood. *h0 are precomputed by assumption */
static double model_loglk(struct model *model)
{
    int i;
    int lfdesc;
    int rtdesc;
    double pk;
    double notpk;
    double lnlk;
    double lnli;
    double lnlj;
    struct node *node;

    // model->t assumed to have been initialized like so
    // t = tree_traverse(POSTORDER, INTERNAL_NODES_ONLY, tree->root, tree);
    node = tree_step(&(model->t));

    while (node) {
        lfdesc = node->children[0]->index;
        rtdesc = node->children[1]->index;
        lnli = model->loglk[lfdesc];
        lnlj = model->loglk[rtdesc];
        pk = -model->rate * model->treelen[node->index];
        notpk = log1p(-exp(pk));
        lnlk = pk + model->h0[node->index] + logspace_add(notpk, lnli + lnlj);
        model->loglk[node->index] = lnlk;
        node = tree_step(&(model->t));
    }

    tree_reset(&(model->t));

    return model->loglk[model->tree->root->index];
}


static void run(int niter, struct model *model)
{
    int i;
    double q = 0;
    double U;
    double rate;
    double loglk;

    for (i = 0; i < niter; ++i) {

        rate = model->rate;
        U = runif(0.5, 2);
        q = -log(U);
        model->rate = U * rate;
        loglk = model_loglk(model);

        if (unif_rand() < exp(loglk - model->LnLk + q))
            model->LnLk = loglk;
        else
            model->rate = rate;

        if (i % 1024 == 0) {
            Rprintf("%f    ", model->rate);
            Rprintf("%f\n", model->LnLk);
        }
    }
}


SEXP rcl_dmm_mcmc3(SEXP rtree, SEXP alpha, SEXP rate, SEXP data, SEXP nstep)
{
    int niter;
    int ncat;
    int *tipdata;
    struct tree *tree;
    struct model *model;

    GetRNGstate();

    ncat = INTEGER(getAttrib(data, R_DimSymbol))[0];
    tipdata = INTEGER(data);
    tree = (struct tree*)R_ExternalPtrAddr(rtree);

    niter = INTEGER(nstep)[0];

    model = model_alloc(ncat, tree->ntip, tree->nnode);
    model_set(REAL(alpha)[0], REAL(rate)[0], tipdata, tree, model);

    model->LnLk = model_loglk(model);

    run(niter, model);

    model_free(model);

    PutRNGstate();

    return R_NilValue;
}
