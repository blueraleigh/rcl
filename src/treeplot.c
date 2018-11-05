#include <R.h>
#include <Rinternals.h>

#include "tree.h"


SEXP rcl_plot_ctree(SEXP rtree, SEXP ages, SEXP direction)
{
    struct tree *tree = (struct tree*)R_ExternalPtrAddr(rtree);
    int nnode = tree->nnode;
    int ntip = tree->ntip;
    int d = INTEGER(direction)[0];
    SEXP coord = PROTECT(allocVector(REALSXP, nnode * 4));
    SEXP bar = PROTECT(allocVector(REALSXP, (nnode-ntip) * 4));
    double *segs = REAL(coord);
    double *bars = REAL(bar);
    double *age = REAL(ages);
    double a, b;
    struct tree_traversal t;
    struct node *node;

    double maxage = 0;
    for (int i = 0; i < tree->nnode; ++i) {
        if (age[tree->node[i]->index] > maxage)
            maxage = age[tree->node[i]->index];
    }

    for (int i = 0; i < ntip; ++i) {
        node = tree->node[i];
        switch (d) {
            case 0: // right
                segs[node->index + 0 * nnode] = age[node->index];                   // x0
                segs[node->index + 1 * nnode] = age[node->index] - node->brlen;     // x1
                segs[node->index + 2 * nnode] = (double)i+1;                        // y0
                segs[node->index + 3 * nnode] = (double)i+1;                        // y1
                break;
            case 1:  // left
                segs[node->index + 0 * nnode] = maxage - age[node->index];
                segs[node->index + 1 * nnode] = maxage - age[node->index] + node->brlen;
                segs[node->index + 2 * nnode] = (double)i+1;
                segs[node->index + 3 * nnode] = (double)i+1;
                break;
            case 2:  // up
                segs[node->index + 0 * nnode] = (double)i+1;
                segs[node->index + 1 * nnode] = (double)i+1;
                segs[node->index + 2 * nnode] = age[node->index];
                segs[node->index + 3 * nnode] = age[node->index] - node->brlen;
                break;
            case 3:  // down
                segs[node->index + 0 * nnode] = (double)i+1;
                segs[node->index + 1 * nnode] = (double)i+1;
                segs[node->index + 2 * nnode] = maxage - age[node->index];
                segs[node->index + 3 * nnode] = maxage - age[node->index] + node->brlen;
                break;
        }
    }

    t = tree_traverse(POSTORDER, INTERNAL_NODES_ONLY, tree->root, tree);
    node = tree_step(&t);

    while (node != NULL) {
        switch (d) {
            case 0:
                segs[node->index + 0 * nnode] = age[node->index];
                segs[node->index + 1 * nnode] = age[node->index] - node->brlen;
                a = segs[node->children[0]->index + 2 * nnode];
                b = segs[node->children[node->degree-1]->index + 2 * nnode];
                segs[node->index + 2 * nnode] = (a + b) / 2;
                segs[node->index + 3 * nnode] = (a + b) / 2;
                bars[(node->index-ntip) + 0 * (nnode-ntip)] = age[node->index];
                bars[(node->index-ntip) + 1 * (nnode-ntip)] = age[node->index];
                bars[(node->index-ntip) + 2 * (nnode-ntip)] = a;
                bars[(node->index-ntip) + 3 * (nnode-ntip)] = b;
                break;
            case 1:
                segs[node->index + 0 * nnode] = maxage - age[node->index];
                segs[node->index + 1 * nnode] = maxage - age[node->index] + node->brlen;
                a = segs[node->children[0]->index + 2 * nnode];
                b = segs[node->children[node->degree-1]->index + 2 * nnode];
                segs[node->index + 2 * nnode] = (a + b) / 2;
                segs[node->index + 3 * nnode] = (a + b) / 2;
                bars[(node->index-ntip) + 0 * (nnode-ntip)] = maxage - age[node->index];
                bars[(node->index-ntip) + 1 * (nnode-ntip)] = maxage - age[node->index];
                bars[(node->index-ntip) + 2 * (nnode-ntip)] = a;
                bars[(node->index-ntip) + 3 * (nnode-ntip)] = b;
                break;
            case 2:
                segs[node->index + 2 * nnode] = age[node->index];
                segs[node->index + 3 * nnode] = age[node->index] - node->brlen;
                a = segs[node->children[0]->index + 0 * nnode];
                b = segs[node->children[node->degree-1]->index + 0 * nnode];
                segs[node->index + 0 * nnode] = (a + b) / 2;
                segs[node->index + 1 * nnode] = (a + b) / 2;
                bars[(node->index-ntip) + 2 * (nnode-ntip)] = age[node->index];
                bars[(node->index-ntip) + 3 * (nnode-ntip)] = age[node->index];
                bars[(node->index-ntip) + 0 * (nnode-ntip)] = a;
                bars[(node->index-ntip) + 1 * (nnode-ntip)] = b;
                break;
            case 3:
                segs[node->index + 2 * nnode] = maxage - age[node->index];
                segs[node->index + 3 * nnode] = maxage - age[node->index] + node->brlen;
                a = segs[node->children[0]->index + 0 * nnode];
                b = segs[node->children[node->degree-1]->index + 0 * nnode];
                segs[node->index + 0 * nnode] = (a + b) / 2;
                segs[node->index + 1 * nnode] = (a + b) / 2;
                bars[(node->index-ntip) + 2 * (nnode-ntip)] = maxage - age[node->index];
                bars[(node->index-ntip) + 3 * (nnode-ntip)] = maxage - age[node->index];
                bars[(node->index-ntip) + 0 * (nnode-ntip)] = a;
                bars[(node->index-ntip) + 1 * (nnode-ntip)] = b;
                break;
        }
        node = tree_step(&t);
    }
    SEXP dim1 = PROTECT(allocVector(INTSXP, 2));
    SEXP dim2 = PROTECT(allocVector(INTSXP, 2));
    INTEGER(dim1)[0] = nnode;
    INTEGER(dim2)[0] = nnode - ntip;
    INTEGER(dim1)[1] = 4;
    INTEGER(dim2)[1] = 4;
    setAttrib(coord, R_DimSymbol, dim1);
    setAttrib(bar, R_DimSymbol, dim2);
    SEXP ret = PROTECT(allocVector(VECSXP, 2));
    SET_VECTOR_ELT(ret, 0, coord);
    SET_VECTOR_ELT(ret, 1, bar);
    UNPROTECT(5);
    return ret;
}


SEXP rcl_plot_ptree(SEXP rtree, SEXP step)
{
    struct tree *tree = (struct tree*)R_ExternalPtrAddr(rtree);
    int nnode = tree->nnode;
    int ntip = tree->ntip;
    SEXP theta = PROTECT(allocVector(REALSXP, nnode * 3));
    double *th = REAL(theta);
    double a, b, vstep = REAL(step)[0];
    struct tree_traversal t;
    struct node *node;
    for (int i = 0; i < ntip; ++i) {
        node = tree->node[i];
        th[node->index + 0 * nnode] = i * vstep;
        th[node->index + 1 * nnode] = 0;
        th[node->index + 2 * nnode] = 0;
    }
    t = tree_traverse(POSTORDER, INTERNAL_NODES_ONLY, tree->root, tree);
    node = tree_step(&t);
    while (node != NULL) {
        a = th[node->children[0]->index];
        b = th[node->children[node->degree-1]->index];
        th[node->index + 0 * nnode] = (a + b) / 2;
        th[node->index + 1 * nnode] = a;
        th[node->index + 2 * nnode] = b;
        node = tree_step(&t);
    }
    SEXP dim = PROTECT(allocVector(INTSXP, 2));
    INTEGER(dim)[0] = nnode;
    INTEGER(dim)[1] = 3;
    setAttrib(theta, R_DimSymbol, dim);
    UNPROTECT(2);
    return theta;
}
