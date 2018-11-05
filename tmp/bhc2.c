#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include "node.h"
#include "tree.h"

void rcl_tree_free(SEXP rtree);

/*
struct node {
    int index;
    int nk;
    double dk;
    double pk;
    double h1;
    double h2;
    double rk
    double loglk;
    struct node *lfdesc;
    struct node *rtdesc;
    struct node *parent;
};


static struct node *alloc_node(int id)
{
    struct node *node;
    node = Calloc(1, struct node);
    node->index = id;
    node->nk = 0;
    node->dk = 0;
    node->pk = 0;
    node->h1 = 0;
    node->h2 = 0;
    node->rk = 0;
    node->loglk = 0;
    node->lfdecs = NULL;
    node->rtdesc = NULL;
    node->parent = NULL;
}

struct tree_traversal {
    int size;
    int order;
    int n;
    int m;
    struct node **in;
    struct node **out;
};

#define PREORDER 0
#define POSTORDER 1
#define ALL_NODES 0
#define INTERNAL_NODES_ONLY 1

static void tree_traverse(int order, int visit, struct node *node, struct tree_traversal *t)
{
    t->order = order;
    t->in[(t->n)++] = node;
    while (t->n) {
        t->out[t->m] = t->in[--(t->n)];
        if (t->out[t->m]->degree) {
            if (visit == ALL_NODES) {
                t->in[(t->n)++] = t->out[t->m]->lfdesc;
                t->in[(t->n)++] = t->out[(t->m)++]->rtdesc;
            } else {
                if (t->out[t->m]->rtdesc->degree)
                    t->in[(t->n)++] = t->out[t->m]->lfdesc;
                if (t->out[t->m]->lfdesc->degree)
                    t->in[(t->n)++] = t->out[t->m]->rtdesc;
                ++(t->m);
            }
        }
    }
}


static struct node *tree_step(struct tree_traversal *t)
{
    switch (t->order) {
        case POSTORDER:
            if (t->m)
                return t->out[--(t->m)];
            return NULL;
        case PREORDER:
            if (t->n != t->m)
                return t->out[(t->n)++];
            return NULL;
    }
    return NULL;
}



static void add_child(struct node *child, struct node *parent)
{
    if (parent->lfdesc)
        parent->rtdesc = child;
    else
        parent->lfdesc = child;
    parent->nk += child->nk
}
*/

/* Compute the logarithm of the multinomial coefficient
** for the count data in x */
static double lmultinomfn(int size, int *x)
{
    int i;
    int n;
    double lcoef;

    n = 0;
    lcoef = 0;

    for (i = 0; i < size; ++i) {
        n += x[i];
        lcoef -= lgammafn(x[i] + 1);
    }

    lcoef += lgammafn(n+1);

    return lcoef;
}

/* dirichlet multinomial log probability */
static double lpmfn(int size, double lcoef, int *count)
{
    int i;
    int x;
    double loglk;

    x = 0;
    loglk = lcoef;

    for (i = 0; i < size; ++i) {
        x += count[i];
        loglk += lgammafn(count[i] + 1);
    }

    loglk += lgammafn(1 * size) - lgammafn(1 * size + x) - size * lgammafn(1);

    return loglk;
}


static double lpmfn2(int size, double lcoef1, double lcoef2, int *count1, int *count2)
{
    int i;
    int x;
    double loglk;

    x = 0;
    loglk = lcoef1 + lcoef2;

    for (i = 0; i < size; ++i) {
        x += count1[i] + count2[i];
        loglk += lgammafn(count1[i] + count2[i] + 1);
    }

    loglk += lgammafn(1 * size) - lgammafn(1 * size + x) - size * lgammafn(1);

    return loglk;
}


//static double cluster(int nobs, int ncat, int *data, int *merge, double *height, double alpha)
static struct node *cluster(int nobs, int ncat, int *data, double alpha, double *lnlk, SEXP colnames)
{
    int i;
    int j;
    int k;

    int nnode = 2*nobs - 1;
    //int count[ncat * nnode];

    //int nk[nnode];
    //int skip[nnode];

    /* merge workspace
    ** wsp[0] = log(alpha)
    ** wsp[2] = dk
    ** wsp[3] = pk
    ** wsp[4] = h1
    ** wsp[5] = h2
    */
    double wsp[6];
    wsp[0] = log(alpha);

    double best;
    //double lmultinomcoef[nnode];
    //double dk[nnode];
    //double pk[nnode];
    //double h1[nnode];
    //double h2[nnode];
    //double loglk[nnode];
    //double rk[nnode][nnode];

    int *count = Calloc(ncat * nnode, int);
    int *nk = Calloc(nnode, int);
    int *skip = Calloc(nnode, int);
    double *lmultinomcoef = Calloc(nnode, double);
    double *dk = Calloc(nnode, double);
    double *pk = Calloc(nnode, double);
    double *h1 = Calloc(nnode, double);
    double *h2 = Calloc(nnode, double);
    double *loglk = Calloc(nnode, double);
    double *rk = Calloc(nnode * nnode, double);




    int lfdesc;
    int rtdesc;
    struct node *nodearr[nnode];

    memset(skip, 0, nnode * sizeof(int));
    memset(count, 0, ncat * nnode * sizeof(int));
    memcpy(count, data, ncat * nobs * sizeof(int));

    for (i = 0; i < nobs; ++i) {
        nodearr[i] = node_alloc();
        nk[i] = 1;
        dk[i] = wsp[0];
        pk[i] = 0;
        lmultinomcoef[i] = lmultinomfn(ncat, data + i * ncat);
        h1[i] = lpmfn(ncat, lmultinomcoef[i], data + i * ncat);
        h2[i] = 0;
        loglk[i] = h1[i];
        //rk[i][i] = 0;
        rk[i + i*nnode] = 0;
        nodearr[i]->label = Calloc(strlen(CHAR(STRING_ELT(colnames, i)))+1, char);
        strcpy(nodearr[i]->label, CHAR(STRING_ELT(colnames, i)));
        nodearr[i]->brlen = 0;
    }

    for (i = 1; i < nobs; ++i) {
        for (j = 0; j < i; ++j) {
            wsp[1] = wsp[0] + lgammafn(nk[i] + nk[j]);
            wsp[2] = logspace_add(wsp[1], dk[i] + dk[j]);
            wsp[3] = wsp[1] - wsp[2];
            wsp[4] = wsp[3] + lpmfn2(ncat, lmultinomcoef[i], lmultinomcoef[j], count + i*ncat, count + j*ncat);
            wsp[5] = log1p(-exp(wsp[3])) + loglk[i] + loglk[j];
            //rk[i][j] = rk[j][i] = wsp[4] - logspace_add(wsp[4], wsp[5]);
            rk[i + j*nnode] = rk[j + i*nnode] = wsp[4] - logspace_add(wsp[4], wsp[5]);
        }
    }

    for (i = nobs; i < nnode; ++i) {
        best = R_NegInf;
        nodearr[i] = node_alloc();

        for (j = 1; j < i; ++j) {
            for (k = 0; k < j; ++k) {
                if (!skip[j] && !skip[k] && rk[j + k*nnode] > best) {
                    lfdesc = j;
                    rtdesc = k;
                    //best = rk[j][k];
                    best = rk[j + k*nnode];
                }
            }
        }

        node_add_child(nodearr[i], nodearr[lfdesc]);
        node_add_child(nodearr[i], nodearr[rtdesc]);

        skip[lfdesc] = skip[rtdesc] = 1;

        //merge[(i - nobs) + 0 * (nobs-1)] = lfdesc < nobs ? -(lfdesc+1) : (lfdesc - nobs + 1);
        //merge[(i - nobs) + 1 * (nobs-1)] = rtdesc < nobs ? -(rtdesc+1) : (rtdesc - nobs + 1);

        memcpy(count + i*ncat, count + lfdesc*ncat, ncat * sizeof(int));
        for (j = 0; j < ncat; ++j)
            count[j + i*ncat] += count[j + rtdesc*ncat];

        lmultinomcoef[i] = lmultinomcoef[lfdesc] + lmultinomcoef[rtdesc];
        nk[i] = nk[lfdesc] + nk[rtdesc];
        wsp[1] = wsp[0] + lgammafn(nk[i]);
        dk[i] = logspace_add(wsp[1], dk[lfdesc] + dk[rtdesc]);
        pk[i] = wsp[1] - dk[i];
        h1[i] = pk[i] + lpmfn(ncat, lmultinomcoef[i], count + i*ncat);
        h2[i] = log1p(-exp(pk[i])) + loglk[lfdesc] + loglk[rtdesc];
        loglk[i] = logspace_add(h1[i], h2[i]);
        rk[i + i*nnode] = h1[i] - loglk[i];

        //height[i - nobs] = rk[i][i];

        nodearr[i]->brlen = -rk[i + i*nnode];

        for (j = 0; j < i; ++j) {
            if (skip[j])
                continue;
            wsp[1] = wsp[0] + lgammafn(nk[i] + nk[j]);
            wsp[2] = logspace_add(wsp[1], dk[i] + dk[j]);
            wsp[3] = wsp[1] - wsp[2];
            wsp[4] = wsp[3] + lpmfn2(ncat, lmultinomcoef[i], lmultinomcoef[j], count + i*ncat, count + j*ncat);
            wsp[5] = log1p(-exp(wsp[3])) + loglk[i] + loglk[j];
            rk[i + j*nnode] = rk[j + i*nnode] = wsp[4] - logspace_add(wsp[4], wsp[5]);
        }

    }

    *lnlk = loglk[nnode-1];
    Free(count);
    Free(nk);
    Free(skip);
    Free(lmultinomcoef);
    Free(dk);
    Free(pk);
    Free(h1);
    Free(h2);
    Free(loglk);
    Free(rk);
    return nodearr[nnode-1];
}


SEXP rcl_bhc2(SEXP alpha, SEXP data)
{
    int ncat;
    int nobs;
    double loglk;

    //SEXP res;
    //SEXP merge;
    //SEXP height;
    //SEXP dim;

    ncat = INTEGER(getAttrib(data, R_DimSymbol))[0];
    nobs = INTEGER(getAttrib(data, R_DimSymbol))[1];

    //dim = PROTECT(allocVector(INTSXP, 2));
    //INTEGER(dim)[0] = nobs-1;
    //INTEGER(dim)[1] = 2;
    //merge = PROTECT(allocVector(INTSXP, 2 * (nobs-1)));
    //height = PROTECT(allocVector(REALSXP, nobs-1));
    //cluster(nobs, ncat, INTEGER(data), INTEGER(merge), REAL(height), REAL(alpha)[0]);

    //setAttrib(merge, R_DimSymbol, dim);

    //res = PROTECT(allocVector(VECSXP, 2));

    //SET_VECTOR_ELT(res, 0, merge);
    //SET_VECTOR_ELT(res, 1, height);

    //UNPROTECT(4);
    //return res;

    struct node *root;
    struct tree *tree;
    SEXP colnames;

    colnames = VECTOR_ELT(getAttrib(data, R_DimNamesSymbol), 1);

    root = cluster(nobs, ncat, INTEGER(data), REAL(alpha)[0], &loglk, colnames);

    tree = tree_build_tree(root);

    SEXP rtree = PROTECT(R_MakeExternalPtr(tree, R_NilValue, R_NilValue));

    R_RegisterCFinalizer(rtree, &rcl_tree_free);
    setAttrib(rtree, install("root"), ScalarInteger(tree->root->index+1));
    setAttrib(rtree, install("Ntip"), ScalarInteger(tree->ntip));
    setAttrib(rtree, install("Nnode"), ScalarInteger(tree->nnode));
    setAttrib(rtree, install("LnLk"), ScalarReal(loglk));
    UNPROTECT(1);

    return rtree;
}
