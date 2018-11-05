#include <R.h>
#include <Rinternals.h>
#include "node.h"
#include "tree.h"
#include "sb.h"

SEXP rcl_ancestors(SEXP, SEXP);
SEXP rcl_mrca(SEXP, SEXP, SEXP);
SEXP rcl_subgraph(SEXP, SEXP);
SEXP rcl_subgraph_brlen(SEXP, SEXP);
SEXP rcl_subgraph_newick(SEXP, SEXP, SEXP, SEXP);


SEXP rcl_ancestors(SEXP rtree, SEXP node)
{
    struct tree *tree = (struct tree*)R_ExternalPtrAddr(rtree);
    int *res = NULL;
    struct node *p = tree->node[INTEGER(node)[0]-1];
    while (p != NULL) {
        sb_push(res, p->index + 1);
        p = p->parent;
    }
    SEXP ret = PROTECT(allocVector(INTSXP, sb_count(res)));
    memcpy(INTEGER(ret), res, sb_count(res) * sizeof(int));
    sb_free(res);
    UNPROTECT(1);
    return ret;
}


SEXP rcl_children(SEXP rtree, SEXP node)
{
    int i;
    struct tree *tree = (struct tree*)R_ExternalPtrAddr(rtree);
    struct node *p = tree->node[INTEGER(node)[0]-1];
    SEXP ret = PROTECT(allocVector(INTSXP, p->degree));
    for (i = 0; i < p->degree; ++i)
        INTEGER(ret)[i] = p->children[i]->index + 1;
    UNPROTECT(1);
    return ret;
}


// Given a rooted tree return the subgraph
// (as an adjaceny list) connecting all nodes
// in nodes argument
SEXP rcl_subgraph(SEXP rtree, SEXP nodes)
{
    struct tree *tree = (struct tree*)R_ExternalPtrAddr(rtree);
    int *nodebits = INTEGER(nodes);    // bitmask indicating nodes to connect
    int *pathbits = Calloc(tree->nnode, int);
    memset(pathbits, 0, tree->nnode * sizeof(int));
    int i, n = tree->nnode - tree->ntip;

    struct tree_traversal t;
    struct node *node;

    // ancestor descendant, ancestor descendant, ... , etc
    int *alist = NULL;

    t = tree_traverse(POSTORDER, ALL_NODES, tree->root, tree);
    node = tree_step(&t);

    while (node) {
        if (nodebits[node->index]) {
            // need to include components back to root
            while (node->parent != NULL) {
                // bail out if we already got this path
                if (pathbits[node->index])
                    break;
                sb_push(alist, node->parent->index + 1);
                sb_push(alist, node->index + 1);
                pathbits[node->index] = 1;
                node = node->parent;
            }
        }
        node = tree_step(&t);
    }

    Free(pathbits);
    SEXP ret = PROTECT(allocVector(INTSXP, sb_count(alist)));
    SEXP dim = PROTECT(allocVector(INTSXP, 2));
    INTEGER(dim)[0] = 2;
    INTEGER(dim)[1] = sb_count(alist) / 2;
    memcpy(INTEGER(ret), alist, sb_count(alist) * sizeof(int));
    setAttrib(ret, R_DimSymbol, dim);
    UNPROTECT(2);
    return ret;
}


SEXP rcl_subgraph_brlen(SEXP rtree, SEXP edge)
{
    struct tree *tree = (struct tree*)R_ExternalPtrAddr(rtree);
    int nedge = INTEGER(getAttrib(edge, R_DimSymbol))[0];
    int *alist = INTEGER(edge);
    int parent, child;
    double brlen;
    struct node *node;
    SEXP brlenv = PROTECT(allocVector(REALSXP, tree->nnode));
    memset(REAL(brlenv), 0, tree->nnode * sizeof(double));
    for (int i = 0; i < nedge; ++i) {
        brlen = 0;
        child = alist[i + 1 * nedge];
        parent = alist[i + 0 * nedge];
        node = tree->node[child];
        while (node->index != parent) {
            brlen += node->brlen;
            node = node->parent;
        }
        REAL(brlenv)[child] = brlen;
    }
    UNPROTECT(1);
    return brlenv;
}


static void subgraph_write_newick(struct tree *tree, struct node *node, char **buffer, int *children, double *brlenv)
{
    char lab[50], brlen[50], *buf = *buffer;
    int child;
    memset(brlen, '\0', sizeof(char) * 50);
    memset(lab, '\0', sizeof(char) * 50);
    if (children[node->index + 0 * tree->nnode] < 0) {
        sprintf(brlen, "%f", brlenv[node->index]);
        sb_add(buf, strlen(brlen)+1);
        if (node->degree == 0) {
            sb_add(buf, strlen(node->label)+1);
            sb_add(buf, 2);
            strcat(buf, node->label);
        }
        else {
            sprintf(lab, "%d", node->index+1);
            sb_add(buf, strlen(lab)+1);
            sb_add(buf, 2);
            strcat(buf, lab);
        }
        strcat(buf, ":");
        strcat(buf, brlen);
    } else {
        sb_add(buf, 2);
        strcat(buf, "(");
        *buffer = buf;
        for (int i = 0; i < 2; ++i) {
            child = children[node->index + i * tree->nnode];
            subgraph_write_newick(tree, tree->node[child], buffer, children, brlenv);
            buf = *buffer;
            sb_add(buf, 2);
            if (i < 1)
                strcat(buf, ",");
            *buffer = buf;
        }
        sb_add(buf, 3);
        strcat(buf, "):");
        sprintf(brlen, "%f", brlenv[node->index]);
        sb_add(buf, strlen(brlen)+1);
        strcat(buf, brlen);
    }
    *buffer = buf;
}


SEXP rcl_subgraph_newick(SEXP rtree, SEXP newroot, SEXP children, SEXP brlen)
{
    struct tree *tree = (struct tree*)R_ExternalPtrAddr(rtree);
    char *buffer = NULL;
    sb_push(buffer, '\0');
    subgraph_write_newick(tree, tree->node[INTEGER(newroot)[0]], &buffer, INTEGER(children), REAL(brlen));
    sb_add(buffer, 2);
    strcat(buffer, ";");
    char *newick = Calloc(sb_count(buffer)+1, char);
    strcpy(newick, buffer);
    sb_free(buffer);
    SEXP ret = PROTECT(allocVector(STRSXP, 1));
    SET_STRING_ELT(ret, 0, mkChar(newick));
    Free(newick);
    UNPROTECT(1);
    return ret;
}


SEXP rcl_mrca(SEXP rtree, SEXP a, SEXP b)
{
    int i, n, la, lb;
    struct node *node1, *node2, *mrca;
    struct tree *tree;
    SEXP res;

    tree = (struct tree*)R_ExternalPtrAddr(rtree);
    la = LENGTH(a);
    lb = LENGTH(b);
    n = la < lb ? la : lb;

    res = PROTECT(allocVector(INTSXP, n));

    for (i = 0; i < n; ++i) {
        node1 = tree->node[INTEGER(a)[i]-1];
        node2 = tree->node[INTEGER(b)[i]-1];
        mrca = tree_mrca(node1, node2);
        INTEGER(res)[i] = mrca->index + 1;
    }

    UNPROTECT(1);
    return res;
}


void rcl_free_tree_traverse(SEXP rt)
{
    struct tree_traversal *t = R_ExternalPtrAddr(rt);
    Free(t);
}


SEXP rcl_tree_traverse(SEXP rtree, SEXP node, SEXP order, SEXP visit)
{
    struct tree_traversal *t;
    struct tree *tree;

    SEXP rt;

    t = Calloc(1, struct tree_traversal);
    tree = (struct tree*)R_ExternalPtrAddr(rtree);

    *t = tree_traverse(INTEGER(order)[0], INTEGER(visit)[0], tree->node[INTEGER(node)[0]-1], tree);

    rt = PROTECT(R_MakeExternalPtr(t, R_NilValue, R_NilValue));
    R_RegisterCFinalizer(rt, &rcl_free_tree_traverse);
    UNPROTECT(1);
    return rt;
}


SEXP rcl_tree_reset(SEXP rt)
{
    struct tree_traversal *t = R_ExternalPtrAddr(rt);
    tree_reset(t);
    return R_NilValue;
}


SEXP rcl_tree_step(SEXP rt)
{
    struct tree_traversal *t = R_ExternalPtrAddr(rt);
    struct node *node;
    node = tree_step(t);
    if (node)
        return ScalarInteger(node->index+1);
    return ScalarInteger(0);
}


SEXP rcl_tree_jump(SEXP nodeid, SEXP rt)
{
    struct tree_traversal *t = R_ExternalPtrAddr(rt);
    struct node *node;
    node = tree_jump(t->tree->node[INTEGER(nodeid)[0]-1], t);
    if (node)
        return ScalarInteger(node->index+1);
    return ScalarInteger(0);
}


