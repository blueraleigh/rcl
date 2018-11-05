#include <R.h>
#include <Rinternals.h>
#include "node.h"
#include "tree.h"


void rcl_tree_free(SEXP rtree);
SEXP rcl_read_newick(SEXP newick);
SEXP rcl_build_tree(SEXP rtree);
SEXP rcl_write_newick(SEXP rtree);
SEXP rcl_subtree(SEXP rtree, SEXP node);


void rcl_tree_free(SEXP rtree)
{
    if (TYPEOF(rtree) == NILSXP)
        return;
    struct tree *tree = (struct tree*)R_ExternalPtrAddr(rtree);
    tree_free(tree);
    R_ClearExternalPtr(rtree);
}


SEXP rcl_read_newick(SEXP newick)
{
    struct tree *tree = tree_read(CHAR(STRING_ELT(newick, 0)));
    SEXP rtree = PROTECT(R_MakeExternalPtr(tree, R_NilValue, R_NilValue));
    R_RegisterCFinalizer(rtree, &rcl_tree_free);
    setAttrib(rtree, install("root"), ScalarInteger(tree->root->index+1));
    setAttrib(rtree, install("Ntip"), ScalarInteger(tree->ntip));
    setAttrib(rtree, install("Nnode"), ScalarInteger(tree->nnode));
    UNPROTECT(1);
    return rtree;
}


SEXP rcl_build_tree(SEXP rtree)
{
    struct tree *tree = (struct tree*)R_ExternalPtrAddr(rtree);
    struct node *node;
    struct tree_traversal t;
    int ntip = tree->ntip;
    int nnode = tree->nnode;
    int nodeidx;
    SEXP children;
    if (tree_isbinary(tree)) {
        children = PROTECT(allocVector(INTSXP, 2 * nnode));
        memset(INTEGER(children), 0, sizeof(int) * (2 * nnode));
    }
    else {
        children = PROTECT(allocVector(VECSXP, nnode));
    }
    SEXP parent = PROTECT(allocVector(INTSXP, nnode));
    SEXP brlen = PROTECT(allocVector(REALSXP, nnode));
    SEXP age = PROTECT(allocVector(REALSXP, nnode));
    SEXP seq_int = PROTECT(allocVector(INTSXP, nnode-ntip));
    SEXP seq_all = PROTECT(allocVector(INTSXP, nnode));
    SEXP seq_pos = PROTECT(allocVector(INTSXP, nnode));
    SEXP lastvisit_int = PROTECT(allocVector(INTSXP, nnode-ntip));
    SEXP lastvisit_all = PROTECT(allocVector(INTSXP, nnode));
    SEXP Ntip = PROTECT(allocVector(INTSXP, 1));
    SEXP Nnode = PROTECT(allocVector(INTSXP, 1));
    SEXP root = PROTECT(allocVector(INTSXP, 1));
    SEXP tiplabel = PROTECT(allocVector(STRSXP, ntip));
    INTEGER(Ntip)[0] = ntip;
    INTEGER(Nnode)[0] = nnode;
    INTEGER(root)[0] = ntip+1;
    memset(INTEGER(parent), 0, nnode * sizeof(int));
    memset(REAL(brlen), 0, nnode * sizeof(double));

    memcpy(INTEGER(seq_pos), tree->seq_pos, nnode * sizeof(int));
    memcpy(INTEGER(lastvisit_int), tree->lastvisit_int, (nnode - ntip) * sizeof(int));
    memcpy(INTEGER(lastvisit_all), tree->lastvisit_all, nnode * sizeof(int));

    t = tree_traverse(PREORDER, ALL_NODES, tree->root, tree);
    node = tree_step(&t);
    if (tree_isbinary(tree)) {
        while (node != NULL) {
            nodeidx = node->index;
            REAL(brlen)[nodeidx] = node->brlen;
            INTEGER(seq_all)[INTEGER(seq_pos)[nodeidx]] = nodeidx + 1;
            if (node->parent != NULL)
                INTEGER(parent)[nodeidx] = node->parent->index + 1;
            if (nodeidx >= ntip) {
                // store children in column-major matrix
                INTEGER(children)[nodeidx + 0 * nnode] = node->children[0]->index + 1;
                INTEGER(children)[nodeidx + 1 * nnode] = node->children[1]->index + 1;
                INTEGER(seq_int)[nodeidx - ntip] = nodeidx + 1;
            } else {
                SET_STRING_ELT(tiplabel, nodeidx, mkChar(node->label));
            }
            node = tree_step(&t);
        }
        SEXP dims = PROTECT(allocVector(INTSXP, 2));
        INTEGER(dims)[0] = nnode;
        INTEGER(dims)[1] = 2;
        setAttrib(children, R_DimSymbol, dims);
    } else {
        while (node != NULL) {
            nodeidx = node->index;
            REAL(brlen)[nodeidx] = node->brlen;
            INTEGER(seq_all)[INTEGER(seq_pos)[nodeidx]] = nodeidx + 1;
            if (node->parent != NULL)
                INTEGER(parent)[nodeidx] = node->parent->index + 1;
            if (nodeidx >= ntip) {
                SET_VECTOR_ELT(children, nodeidx - tree->ntip, allocVector(INTSXP, node->degree));
                INTEGER(seq_int)[nodeidx - ntip] = nodeidx + 1;
                for (int j = 0; j < node->degree; ++j)
                    INTEGER(VECTOR_ELT(children, nodeidx - tree->ntip))[j] = node->children[j]->index + 1;
            } else {
                SET_STRING_ELT(tiplabel, nodeidx, mkChar(node->label));
            }
            node = tree_step(&t);
        }
    }
    for (int i = 0; i < nnode; ++i) {
        node = tree->node[i];
        REAL(age)[i] = 0.0;
        while (node != NULL) {
            REAL(age)[i] += node->brlen;
            node = node->parent;
        }
        INTEGER(seq_pos)[i] += 1;
        INTEGER(lastvisit_all)[i] += 1;
        if (i >= ntip)
            INTEGER(lastvisit_int)[i - ntip] += 1;
    }
    SEXP phy = PROTECT(allocVector(VECSXP, 13));
    SEXP names = PROTECT(allocVector(STRSXP, 13));
    SET_STRING_ELT(names, 0, mkChar("children"));
    SET_STRING_ELT(names, 1, mkChar("parent"));
    SET_STRING_ELT(names, 2, mkChar("brlen"));
    SET_STRING_ELT(names, 3, mkChar("tip.label"));
    SET_STRING_ELT(names, 4, mkChar("Ntip"));
    SET_STRING_ELT(names, 5, mkChar("Nnode"));
    SET_STRING_ELT(names, 6, mkChar("root"));
    SET_STRING_ELT(names, 7, mkChar("preorder_int"));
    SET_STRING_ELT(names, 8, mkChar("preorder_all"));
    SET_STRING_ELT(names, 9, mkChar("preorder_pos"));
    SET_STRING_ELT(names, 10, mkChar("lastvisit_int"));
    SET_STRING_ELT(names, 11, mkChar("lastvisit_all"));
    SET_STRING_ELT(names, 12, mkChar("age"));
    SET_VECTOR_ELT(phy, 0, children);
    SET_VECTOR_ELT(phy, 1, parent);
    SET_VECTOR_ELT(phy, 2, brlen);
    SET_VECTOR_ELT(phy, 3, tiplabel);
    SET_VECTOR_ELT(phy, 4, Ntip);
    SET_VECTOR_ELT(phy, 5, Nnode);
    SET_VECTOR_ELT(phy, 6, root);
    SET_VECTOR_ELT(phy, 7, seq_int);
    SET_VECTOR_ELT(phy, 8, seq_all);
    SET_VECTOR_ELT(phy, 9, seq_pos);
    SET_VECTOR_ELT(phy, 10, lastvisit_int);
    SET_VECTOR_ELT(phy, 11, lastvisit_all);
    SET_VECTOR_ELT(phy, 12, age);
    setAttrib(phy, R_NamesSymbol, names);
    if (tree_isbinary(tree))
        UNPROTECT(16);
    else
        UNPROTECT(15);
    return phy;
}


SEXP rcl_write_newick(SEXP rtree)
{
    struct tree *tree = (struct tree*)R_ExternalPtrAddr(rtree);
    char *newick = node_print_newick(tree->root);
    SEXP ret = PROTECT(allocVector(STRSXP, 1));
    SET_STRING_ELT(ret, 0, mkChar(newick));
    Free(newick);
    UNPROTECT(1);
    return ret;
}


SEXP rcl_subtree(SEXP rtree, SEXP node)
{
    struct tree *tree = (struct tree*)R_ExternalPtrAddr(rtree);
    char *newick = node_print_newick(tree->node[INTEGER(node)[0]-1]);
    struct tree *subtree = tree_read(newick);
    subtree->root->brlen = 0;
    SEXP subrtree = PROTECT(R_MakeExternalPtr(subtree, R_NilValue, R_NilValue));
    R_RegisterCFinalizer(subrtree, &rcl_tree_free);
    setAttrib(subrtree, install("root"), ScalarInteger(subtree->root->index+1));
    setAttrib(subrtree, install("Ntip"), ScalarInteger(subtree->ntip));
    setAttrib(subrtree, install("Nnode"), ScalarInteger(subtree->nnode));
    Free(newick);
    UNPROTECT(1);
    return subrtree;
}


SEXP rcl_tiplabels(SEXP rtree)
{
    int i;
    SEXP tiplabel;
    struct tree *tree = (struct tree*)R_ExternalPtrAddr(rtree);

    tiplabel = PROTECT(allocVector(STRSXP, tree->ntip));

    for (i = 0; i < tree->ntip; ++i)
        SET_STRING_ELT(tiplabel, i, mkChar(tree->node[i]->label));

    UNPROTECT(1);
    return tiplabel;
}


SEXP rcl_node_brlens(SEXP rtree)
{
    int i;
    SEXP brlen;
    struct tree *tree = (struct tree*)R_ExternalPtrAddr(rtree);

    brlen = PROTECT(allocVector(REALSXP, tree->nnode));

    for (i = 0; i < tree->nnode; ++i)
        REAL(brlen)[i] = tree->node[i]->brlen;

    UNPROTECT(1);
    return brlen;
}


SEXP rcl_node_ages(SEXP rtree)
{
    int i;
    double *node_age;
    SEXP age;
    struct tree *tree = (struct tree*)R_ExternalPtrAddr(rtree);
    struct node *node;

    age = PROTECT(allocVector(REALSXP, tree->nnode));
    node_age = REAL(age);

    for (i = 0; i < tree->nnode; ++i) {
        node = tree->node[i];
        node_age[i] = 0.0;
        while (node != NULL) {
            node_age[i] += node->brlen;
            node = node->parent;
        }
    }

    UNPROTECT(1);
    return age;
}
