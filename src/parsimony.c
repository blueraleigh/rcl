#include <R.h>
#include <Rinternals.h>
#include "node.h"
#include "tree.h"


/* (inefficient) Parsimony routines for a single r-state characer */


/*
** Fitch parsimony
*/

/*
** Perform a union and intersection operation on
** two r-length binary vectors
**
** @param r number of character states
** @param a left binary vector
** @param b right binary vector
** @param x vector to hold intersection of a and b
** @param u vector to hold union of a and b
** @return 1 if we can take an intersection, 0 otherwise
*/
static int xoru(int r, int *a, int *b, int *x, int *u)
{
    int xtype, n = 0;
    for (int i = 0; i < r; ++i) {
        xtype = a[i] + b[i];
        switch (xtype) {
            case 0:
                x[i] = u[i] = 0;
                break;
            case 1:
                x[i] = 0;  // intersection
                u[i] = 1;  // union
                break;
            case 2:
                x[i] = u[i] = 1;
                n = 1;
                break;
        }
    }
    return n;
}


/*
** Determine if two r-length binary vectors are identical,
** disjoint or somewhere in between
**
** @return 1 if identical, 0 if disjoint, -1 if in between
*/
static int iddj(int r, int *a, int *b)
{
    int size_a = 0, size_b = 0, id = 0;
    for (int i = 0; i < r; ++i) {
        if (a[i] && b[i])
            id += 1;
        if (a[i])
            size_a += 1;
        if (b[i])
            size_b += 1;
    }
    if (id == 0)
        return 0;
    if ((size_a == size_b) && (id == size_a))
        return 1;
    return -1;
}


/*
** Uniformly sample a 1-valued index from an r-length
** binary vector (reservoir sampling)
*/
static int choose_state(int r, int *a)
{
    int j;
    double k = 1;
    for (int i = 0; i < r; ++i) {
        if (a[i]) {
            if (unif_rand() < 1/k)
                j = i;
            k += 1;
        }
    }
    return j;
}


/*
** Perform the Fitch parsimony downpass algorithm
**
** @param ddata binary vector to hold the downpass state sets
**  at each internal node. The i-th position will correspond
**  to the downpass state set for the node with i-th index.
**  i.e. This is a r x tree->nnode matrix stored in column major
**  as a flat array. This same convention applies to all other
**  variables in this file with names like udata, dcost, ucost.
** @param pscores vector of parsimony scores. Each value is
**  the number of parsimony changes occurring in the subclade
**  rooted at the node with the corresponding index number.
** @param r the number of character states
*/
static void fitch_downpass(struct tree *tree, int *ddata, int *pscores, int r)
{
    int parent, lfchild, rtchild, pscore;
    int x[r];
    int u[r];
    memset(pscores, 0, tree->nnode * sizeof(int));
    struct tree_traversal t = tree_traverse(POSTORDER, INTERNAL_NODES_ONLY, tree->root, tree);
    struct node *node = tree_step(&t);
    while (node != NULL) {
        pscore = 0;
        parent = node->index;
        lfchild = node->children[0]->index;
        rtchild = node->children[1]->index;
        if (xoru(r, ddata + lfchild*r, ddata + rtchild*r, x, u)) {
            memcpy(ddata + parent*r, x, r * sizeof(int));
        } else {
            pscore += 1;
            memcpy(ddata + parent*r, u, r * sizeof(int));
        }
        pscores[parent] += pscores[lfchild];
        pscores[parent] += pscores[rtchild];
        pscores[parent] += pscore;
        node = tree_step(&t);
    }
}


/*
** Perform the Fitch uppass algorithm
**
** @param udata binary vector to hold the uppass state sets
**  at each internal node. The i-th position will correspond
**  to the uppass state set for the node with i-th index
** @param ddata the downpass state sets resulting from the
**  Fitch downpass algorithm
*/
static void fitch_uppass(struct tree *tree, int *udata, int *ddata, int r)
{
    int focal, parent, lfchild, rtchild, diff;
    int x[r];
    int u[r];
    memcpy(udata, ddata, (tree->nnode * r) * sizeof(int));
    struct tree_traversal t = tree_traverse(PREORDER, INTERNAL_NODES_ONLY, tree->root, tree);
    struct node *node = tree_step(&t);
    node = tree_step(&t);   // step past the root because its final state set is already complete
    while (node != NULL) {
        focal = node->index;
        parent = node->parent->index;
        lfchild = node->children[0]->index;
        rtchild = node->children[1]->index;
        xoru(r, ddata + focal*r, udata + parent*r, x, u);
        if (iddj(r, udata + parent*r, x) == 1) {
            memcpy(udata + focal*r, x, r * sizeof(int));
        } else if (iddj(r, ddata + lfchild*r, ddata + rtchild*r) == 0) {
            memcpy(udata + focal*r, u, r * sizeof(int));
        } else {
            xoru(r, ddata + lfchild*r, ddata + rtchild*r, x, u);
            xoru(r, udata + parent*r, u, x, u);
            xoru(r, ddata + focal*r, x, x, u);
            memcpy(udata + focal*r, u, r * sizeof(int));
        }
        node = tree_step(&t);
    }
}


SEXP rcl_fitch_pscore(SEXP rtree, SEXP ddata, SEXP pscore)
{
    int r = INTEGER(getAttrib(ddata, R_DimSymbol))[0];
    struct tree *tree = (struct tree*)R_ExternalPtrAddr(rtree);
    fitch_downpass(tree, INTEGER(ddata), INTEGER(pscore), r);
    return R_NilValue;
}


SEXP rcl_fitch_mpr(SEXP rtree, SEXP up, SEXP down, SEXP pscore)
{
    int r = INTEGER(getAttrib(down, R_DimSymbol))[0];
    struct tree *tree = (struct tree*)R_ExternalPtrAddr(rtree);
    int *ddata = INTEGER(down);
    int *udata = INTEGER(up);
    fitch_downpass(tree, ddata, INTEGER(pscore), r);
    fitch_uppass(tree, udata, ddata, r);
    return R_NilValue;
}


/*
** Count the number of MPR reconstructions
**
** Algorithm modified from Mesquite project
** https://github.com/MesquiteProject/MesquiteCore/blob/master/Source/mesquite/parsimony/lib/MPRProcessor.java
*/

#define ELEM(m, i, j) m[(i) + tree->nnode * (j)]

#define DO_CHILD(d) do {                                                    \
    child = node->children[(d)]->index;                                     \
    for (int j = 0; j < r; ++j) {                                           \
        if (udata[j + parent * r]) {                                        \
            np = 0;                                                         \
            if (ddata[j + child * r]) {                                     \
                np += ELEM(nh, child, j);                                   \
            } else {                                                        \
                for (int k = 0; k < r; ++k) {                               \
                    if (ddata[k + child * r] || udata[k + child * r])       \
                        np += ELEM(nh, child, k);                           \
                }                                                           \
            }                                                               \
            ELEM(nh, parent, j) *= np;                                      \
        }                                                                   \
    }                                                                       \
} while (0)

static double fitch_count(struct tree *tree, int *udata, int *ddata, int r)
{
    int parent, child;
    double np, nbr = 0;
    double nh[tree->nnode * r];
    memset(nh, 0, (tree->nnode * r) * sizeof(double));
    struct tree_traversal t;
    struct node *node;

    for (int i = 0; i < tree->ntip; ++i) {
        for (int j = 0; j < r; ++j) {
            if (ddata[j + i * r])
                ELEM(nh, i, j) = 1;
        }
    }

    t = tree_traverse(POSTORDER, INTERNAL_NODES_ONLY, tree->root, tree);
    node = tree_step(&t);

    while (node != NULL) {
        parent = node->index;
        for (int j = 0; j < r; ++j) {
            if (udata[j + parent * r])
                ELEM(nh, parent, j) = 1;
        }
        DO_CHILD(0);
        DO_CHILD(1);
        node = tree_step(&t);
    }

    for (int j = 0; j < r; ++j)
        nbr += ELEM(nh, tree->root->index, j);

    return nbr;
}

#undef DO_CHILD


/*
** Sample histories of character evolution from a maximum parsimony reconstruction
** of ancestral states.
**
** @param nodestate a matrix to store the sampled states at each node in each sample
** @param nsample the number of samples to take
*/
static void fitch_history(struct tree *tree, int *nodestate, int *udata, int *ddata, int r, int nsample)
{
    int n = tree->nnode - tree->ntip;
    int focal, parent, cstate, pstate;

    struct tree_traversal t;
    struct node *node;

    for (int k = 0; k < nsample; ++k) {
        nodestate[tree->root->index + k * tree->nnode] = choose_state(r, udata + r * tree->root->index);

        t = tree_traverse(PREORDER, INTERNAL_NODES_ONLY, tree->root, tree);
        node = tree_step(&t);
        node = tree_step(&t);

        while (node != NULL) {
            focal = node->index;
            parent = node->parent->index;
            pstate = nodestate[parent + k * tree->nnode];
            if (ddata[pstate + focal*r]) {
                cstate = pstate;
            } else {
                if (udata[pstate + focal*r]) {
                    // temporarily modify the downpass state set
                    // to include the parent state, which is a
                    // valid MPR state-to-state transition in this case.
                    ddata[pstate + r*focal] = 1;
                    cstate = choose_state(r, ddata + r*focal);
                    ddata[pstate + r*focal] = 0;
                } else {
                    cstate = choose_state(r, ddata + r*focal);
                }
            }
            nodestate[focal + k*tree->nnode] = cstate;
            node = tree_step(&t);
        }
    }
}


SEXP rcl_fitch_history(SEXP rtree, SEXP nodestate, SEXP up, SEXP down, SEXP nsample)
{
    GetRNGstate();
    int r = INTEGER(getAttrib(down, R_DimSymbol))[0];
    struct tree *tree = (struct tree*)R_ExternalPtrAddr(rtree);
    int *ddata = INTEGER(down);
    int *udata = INTEGER(up);
    int *state = INTEGER(nodestate);
    fitch_history(tree, state, udata, ddata, r, INTEGER(nsample)[0]);
    PutRNGstate();
    return R_NilValue;
}


SEXP rcl_fitch_count(SEXP rtree, SEXP up, SEXP down)
{
    int r = INTEGER(getAttrib(down, R_DimSymbol))[0];
    struct tree *tree = (struct tree*)R_ExternalPtrAddr(rtree);
    int *ddata = INTEGER(down);
    int *udata = INTEGER(up);
    SEXP nbr = PROTECT(allocVector(REALSXP, 1));
    REAL(nbr)[0] = fitch_count(tree, udata, ddata, r);
    UNPROTECT(1);
    return nbr;
}


/*
** Sankoff parsimony
*/

#define MAX_SANKOFF_COST 10000  // needs to be the same as max on R side

/*
** Determine the cost associated with each state assignment to an
** internal node
**
** @param cost matrix of state-to-state transition costs
** @param pstate cost of each state assignment to parent
** @param lfstate cost of each state assignment to left descendant
** @param rtstate cost of each state assignment to right descendant
*/
static int sankoff_cost(int r, int *cost, int *pstate, int *lfstate, int *rtstate)
{
    int ts, min_lf, min_rt, mint = MAX_SANKOFF_COST;
    for (int s = 0; s < r; ++s) {
        min_lf = MAX_SANKOFF_COST;
        for (int t = 0; t < r; ++t) {
            ts = cost[s + t * r] + lfstate[t];
            if (ts < min_lf)
                min_lf = ts;
        }
        min_rt = MAX_SANKOFF_COST;
        for (int t = 0; t < r; ++t) {
            ts = cost[s + t * r] + rtstate[t];
            if (ts < min_rt)
                min_rt = ts;
        }
        pstate[s] = min_lf + min_rt;
        if (pstate[s] < mint)
            mint = pstate[s];
    }
    return mint;
}


/*
** Convert Sankoff downpass scores to 1s and 0s
**
** Once the sankoff downpass algorithm is completed the set of preliminary
** node state sets includes all the states at a node that have minimum cost.
** This function steps through the tree and identifies those states and
** assigns them membership in the state set by giving them a 1. Non-membership
** is indicated by a 0.
*/
static void sankoff_bitmask(struct tree *tree, int r, int *ddata, int *dcost, int *pscores)
{
    int mint, *cost, *state;
    struct tree_traversal t;
    struct node *node;
    t = tree_traverse(PREORDER, ALL_NODES, tree->root, tree);
    node = tree_step(&t);
    while (node != NULL) {
        cost = dcost + r * node->index;
        state = ddata + r * node->index;
        mint = pscores[node->index];
        for (int s = 0; s < r; ++s)
            state[s] = (cost[s] > mint) ? 0 : 1;
        node = tree_step(&t);
    }
}


/*
** Sankoff downpass algorithm
*/
static void sankoff_downpass(struct tree *tree, int *cost, int *dcost, int *pscores, int r)
{
    int parent, lfchild, rtchild;
    memset(pscores, 0, tree->nnode * sizeof(int));

    struct tree_traversal t;
    struct node *node;

    t = tree_traverse(POSTORDER, INTERNAL_NODES_ONLY, tree->root, tree);
    node = tree_step(&t);

    while (node != NULL) {
        parent = node->index;
        lfchild = node->children[0]->index;
        rtchild = node->children[1]->index;
        pscores[parent] = sankoff_cost(r, cost, dcost + parent*r, dcost + lfchild*r, dcost + rtchild*r);
        node = tree_step(&t);
    }
}


static void sankoff_uppass(struct tree *tree, int *cost, int *udata, int *ddata, int *dcost, int *pscores, int r)
{
    sankoff_bitmask(tree, r, ddata, dcost, pscores);
    memcpy(udata, ddata, (tree->nnode * r) * sizeof(int));

    struct tree_traversal t;
    struct node *node;
    int *state, *pstate;
    int tmpstate[r];
    int mint;

    t = tree_traverse(PREORDER, INTERNAL_NODES_ONLY, tree->root, tree);
    node = tree_step(&t);
    node = tree_step(&t);

    while (node != NULL) {
        state = dcost + r * node->index;
        pstate = udata + r * node->parent->index;
        memset(tmpstate, 0, r * sizeof(int));
        for (int s = 0; s < r; ++s) {
            mint = MAX_SANKOFF_COST;
            if (pstate[s]) {
                for (int q = 0; q < r; ++q) {
                    if ((cost[s + q * r] + state[q]) < mint)
                        mint = cost[s + q * r] + state[q];
                }
                for (int q = 0; q < r; ++q) {
                    if ((cost[s + q * r] + state[q]) == mint)
                        tmpstate[q] = 1;
                }
            }
        }
        memcpy(udata + r * node->index, tmpstate, r * sizeof(int));
        node = tree_step(&t);
    }
}


SEXP rcl_sankoff_pscore(SEXP rtree, SEXP costmatrix, SEXP dcost, SEXP pscore)
{
    int r = INTEGER(getAttrib(dcost, R_DimSymbol))[0];
    struct tree *tree = (struct tree*)R_ExternalPtrAddr(rtree);
    sankoff_downpass(tree, INTEGER(costmatrix), INTEGER(dcost), INTEGER(pscore), r);
    return R_NilValue;
}


SEXP rcl_sankoff_mpr(SEXP rtree, SEXP costmatrix, SEXP up, SEXP down, SEXP downcost, SEXP pscore)
{
    int r = INTEGER(getAttrib(down, R_DimSymbol))[0];
    struct tree *tree = (struct tree*)R_ExternalPtrAddr(rtree);
    int *ddata = INTEGER(down);
    int *dcost = INTEGER(downcost);
    int *udata = INTEGER(up);
    sankoff_downpass(tree, INTEGER(costmatrix), dcost, INTEGER(pscore), r);
    sankoff_uppass(tree, INTEGER(costmatrix), udata, ddata, dcost, INTEGER(pscore), r);
    return R_NilValue;
}


static void sankoff_history(struct tree *tree, int *nodestate, int *cost, int *udata, int *ddata, int *dcost, int r, int nsample)
{
    int n = tree->nnode - tree->ntip;
    int focal, parent, cstate, pstate, mint;

    int *state;
    int tmpstate[r];

    struct tree_traversal t;
    struct node *node;

    for (int k = 0; k < nsample; ++k) {
        nodestate[tree->root->index + k * tree->nnode] = choose_state(r, udata + r * tree->root->index);

        t = tree_traverse(PREORDER, INTERNAL_NODES_ONLY, tree->root, tree);
        node = tree_step(&t);
        node = tree_step(&t);

        while (node != NULL) {
            memset(tmpstate, 0, r * sizeof(int));
            mint = MAX_SANKOFF_COST;
            focal = node->index;
            parent = node->parent->index;
            pstate = nodestate[parent + k * tree->nnode];

            state = dcost + r * focal;

            for (int q = 0; q < r; ++q) {
                if ((cost[pstate + q * r] + state[q]) < mint)
                    mint = cost[pstate + q * r] + state[q];
            }
            for (int q = 0; q < r; ++q) {
                if ((cost[pstate + q * r] + state[q]) == mint)
                    tmpstate[q] = 1;
            }

            cstate = choose_state(r, tmpstate);

            nodestate[focal + k*tree->nnode] = cstate;
            node = tree_step(&t);
        }
    }
}


SEXP rcl_sankoff_history(SEXP rtree, SEXP costmatrix, SEXP nodestate, SEXP up, SEXP down, SEXP dcost, SEXP nsample)
{
    GetRNGstate();
    int r = INTEGER(getAttrib(down, R_DimSymbol))[0];
    struct tree *tree = (struct tree*)R_ExternalPtrAddr(rtree);
    int *ddata = INTEGER(down);
    int *udata = INTEGER(up);
    int *state = INTEGER(nodestate);
    sankoff_history(tree, state, INTEGER(costmatrix), udata, ddata, INTEGER(dcost), r, INTEGER(nsample)[0]);
    PutRNGstate();
    return R_NilValue;
}


#define DO_CHILD(d) do {                                                    \
    child = node->children[(d)]->index;                                     \
    state = dcost + child * r;                                              \
    for (int s = 0; s < r; ++s) {                                           \
        mint = MAX_SANKOFF_COST;                                            \
        if (udata[s + parent * r]) {                                        \
            np = 0;                                                         \
            for (int q = 0; q < r; ++q) {                                   \
                if ((cost[s + q * r] + state[q]) < mint)                    \
                    mint = cost[s + q * r] + state[q];                      \
            }                                                               \
            for (int q = 0; q < r; ++q) {                                   \
                if ((cost[s + q * r] + state[q]) == mint)                   \
                    np += ELEM(nh, child, q);                               \
            }                                                               \
            ELEM(nh, parent, s) *= np;                                      \
        }                                                                   \
    }                                                                       \
} while (0)

static double sankoff_count(struct tree *tree, int *cost, int *udata, int *ddata, int *dcost, int r)
{
    int parent, child;
    double np, nbr = 0;
    double nh[tree->nnode * r];
    memset(nh, 0, (tree->nnode * r) * sizeof(double));
    struct tree_traversal t;
    struct node *node;

    for (int i = 0; i < tree->ntip; ++i) {
        for (int j = 0; j < r; ++j) {
            if (ddata[j + i * r])
                ELEM(nh, i, j) = 1;
        }
    }

    int *state;
    int mint;

    t = tree_traverse(POSTORDER, INTERNAL_NODES_ONLY, tree->root, tree);
    node = tree_step(&t);

    while (node != NULL) {
        parent = node->index;
        for (int j = 0; j < r; ++j) {
            if (udata[j + parent * r])
                ELEM(nh, parent, j) = 1;
        }
        DO_CHILD(0);
        DO_CHILD(1);
        node = tree_step(&t);
    }

    for (int j = 0; j < r; ++j)
        nbr += ELEM(nh, tree->root->index, j);

    return nbr;
}

#undef DO_CHILD


SEXP rcl_sankoff_count(SEXP rtree, SEXP costmatrix, SEXP up, SEXP down, SEXP dcost)
{
    int r = INTEGER(getAttrib(down, R_DimSymbol))[0];
    struct tree *tree = (struct tree*)R_ExternalPtrAddr(rtree);
    int *ddata = INTEGER(down);
    int *udata = INTEGER(up);
    SEXP nbr = PROTECT(allocVector(REALSXP, 1));
    REAL(nbr)[0] = sankoff_count(tree, INTEGER(costmatrix), udata, ddata, INTEGER(dcost), r);
    UNPROTECT(1);
    return nbr;
}


static void sankoff_uppass2(struct tree *tree, int *cost, int *ucost, int *dcost, int r)
{
    memcpy(ucost, dcost, (tree->nnode * r) * sizeof(int));

    struct tree_traversal t;
    struct node *node;
    int *state, *pstate;
    int tmpcost[r];
    int mint, ts;

    t = tree_traverse(PREORDER, INTERNAL_NODES_ONLY, tree->root, tree);
    node = tree_step(&t);
    node = tree_step(&t);

    while (node != NULL) {
        state = dcost + r * node->index;
        pstate = ucost + r * node->parent->index;
        memset(tmpcost, 0, r * sizeof(int));

        for (int s = 0; s < r; ++s) {
            mint = MAX_SANKOFF_COST;
            for (int t = 0; t < r; ++t) {
                ts = cost[s + t * r] + state[t];
                if (ts < mint)
                    mint = ts;
            }
            tmpcost[s] = mint;
        }

        for (int s = 0; s < r; ++s) {
            mint = MAX_SANKOFF_COST;
            for (int t = 0; t < r; ++t) {
                ts = (pstate[t]-tmpcost[t]) + cost[t + s * r] + state[s];
                if (ts < mint)
                    mint = ts;
            }
            ucost[s + r * node->index] = mint;
        }

        node = tree_step(&t);
    }
}


SEXP rcl_sankoff_cost(SEXP rtree, SEXP costmatrix, SEXP upcost, SEXP downcost, SEXP pscores)
{
    int r = INTEGER(getAttrib(downcost, R_DimSymbol))[0];
    struct tree *tree = (struct tree*)R_ExternalPtrAddr(rtree);
    sankoff_downpass(tree, INTEGER(costmatrix), INTEGER(downcost), INTEGER(pscores), r);
    sankoff_uppass2(tree, INTEGER(costmatrix), INTEGER(upcost), INTEGER(downcost), r);
    return R_NilValue;
}
