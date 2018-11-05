#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>


/* Bayesian hierarchical clustering of multinomial data */

struct hclust {

    /* Number of data points */
    int n;

    /* Number of categories in the multinomial */
    int k;

    /* Current number of clusters*/
    int nclust;

    /*
    ** Input data. k by n matrix
    ** that holds counts in each category.
    */
    int *data;

    /* Root of the cluster tree */
    struct node *root;

    /* Workspace for building the cluster tree */
    struct node **node;

    /* Dirichlet process concentration parameter */
    double alpha;

    /* log(alpha) */
    double lalpha;

    /* Symmetric Dirichlet prior on multinomial counts */
    double beta;

    /* lgamma(beta) */
    double lgbeta;

    /* lgamma(beta*k) */
    double lgbetak;
};


struct node {

    int index;

    /* Number of data points in node's subtree */
    int nk;

    /* Indices of data points in node's subtree */
    int *dataid;

    /* Marginal log likelihood assuming data i.i.d. from single cluster */
    double h1;

    /* Marginal log likelihood assuming data from multiple clusters */
    double h2;

    /* Scalar used to compute pk */
    double dk;

    /* Prior log likelihood that node's data i.i.d. from one cluster */
    double pk;

    /* Posterior log likelihood that node's data i.i.d. from one cluster */
    double rk;

    double loglk;

    /* Left descendant of node */
    struct node *lfdesc;

    /* Right descendant of node */
    struct node *rtdesc;

    /* Ancestor of node */
    struct node *parent;
};


static struct node *node_alloc(int index)
{
    struct node *node;

    node = Calloc(1, struct node);

    node->index = index;
    node->nk = 0;
    node->h1 = 0;
    node->h2 = 0;
    node->dk = 0;
    node->pk = 0;
    node->rk = 0;
    node->dataid = NULL;
    node->rtdesc = NULL;
    node->lfdesc = NULL;
    node->parent = NULL;

    return node;
}


static void node_free(struct node *node)
{
    if (node) {
        node_free(node->lfdesc);
        node_free(node->rtdesc);
        Free(node->dataid);
        Free(node);
    }
}


static void node_setdata(struct node *node)
{
    node->dataid = Calloc(node->nk, int);
    memcpy(node->dataid, node->lfdesc->dataid, node->lfdesc->nk * sizeof(int));
    memcpy(node->dataid + node->lfdesc->nk, node->rtdesc->dataid, node->rtdesc->nk * sizeof(int));
}


static void node_addchild(struct node *child, struct node *parent)
{
    if (parent->lfdesc && parent->rtdesc)
        error("Attempt to add third descendant");

    if (parent->lfdesc)
        parent->rtdesc = child;
    else
        parent->lfdesc = child;

    child->parent = parent;

    parent->nk += child->nk;
}


static void node_rmchild(struct node *child, struct node *parent)
{
    if (!parent->lfdesc && !parent->rtdesc)
        error("Attempt to remove null descendant");

    if (parent->lfdesc == child)
        parent->lfdesc = NULL;
    else
        parent->rtdesc = NULL;

    child->parent = NULL;

    parent->nk -= child->nk;
}


static struct hclust *hclust_alloc(int n, int k, int *data, double alpha, double beta)
{
    int i;
    int j;
    int size;
    int *count;
    double loglk;
    struct node *node;
    struct hclust *hclust;

    hclust = Calloc(1, struct hclust);

    hclust->n = n;
    hclust->k = k;
    hclust->nclust = n;
    hclust->alpha = alpha;
    hclust->lalpha = log(alpha);
    hclust->beta = beta;
    hclust->lgbeta = lgamma(beta);
    hclust->lgbetak = lgamma(k * beta);
    hclust->data = data;
    hclust->root = NULL;
    hclust->node = Calloc(n, struct node*);

    for (i = 0; i < n; ++i) {
        node = node_alloc(i);
        node->nk = 1;
        node->dataid = Calloc(1, int);
        node->dataid[0] = i;
        count = data + i*k;
        loglk = 0;
        size = 0;
        for (j = 0; j < k; ++j) {
            size += count[j];
            loglk += lgamma(count[j] + hclust->beta) - hclust->lgbeta - lgamma(count[j]+1);
        }
        loglk += hclust->lgbetak + lgamma(size+1) - lgamma(size + hclust->k * hclust->beta);
        node->h1 = node->loglk = loglk;
        hclust->node[i] = node;
    }

    return hclust;
}


static void hclust_free(struct hclust *hclust)
{
    node_free(hclust->root);
    Free(hclust->node);
    Free(hclust);
}



/* Log dirichlet multinomial probability */
static double cluster_loglk(struct node *node, struct hclust *hclust)
{
    int i;
    int j;
    int n;
    int count[hclust->k];
    double loglk;

    memset(count, 0, hclust->k * sizeof(int));

    loglk = 0;
    for (i = 0; i < node->lfdesc->nk; ++i) {
        n = 0;
        for (j = 0; j < hclust->k; ++j) {
            count[j] += hclust->data[j + node->lfdesc->dataid[i] * hclust->k];
            n += hclust->data[j + node->lfdesc->dataid[i] * hclust->k];
            loglk -= lgamma(hclust->data[j + node->lfdesc->dataid[i] * hclust->k]+1);
        }
        loglk += lgamma(n+1);
    }
    for (i = 0; i < node->rtdesc->nk; ++i) {
        n = 0;
        for (j = 0; j < hclust->k; ++j) {
            count[j] += hclust->data[j + node->rtdesc->dataid[i] * hclust->k];
            n += hclust->data[j + node->rtdesc->dataid[i] * hclust->k];
            loglk -= lgamma(hclust->data[j + node->rtdesc->dataid[i] * hclust->k]+1);
        }
        loglk += lgamma(n+1);
    }

    n = 0;
    for (j = 0; j < hclust->k; ++j) {
        n += count[j];
        loglk += lgamma(count[j] + hclust->beta) - hclust->lgbeta;
    }

    loglk += hclust->lgbetak - lgamma(n + hclust->k * hclust->beta);
    return loglk;
}


/* Compute node likelihoods */
static void node_compute(struct node *node, struct hclust *hclust)
{
    double x;
    double clnl;

    clnl = cluster_loglk(node, hclust);

    x = hclust->lalpha + lgamma(node->nk);

    node->dk = logspace_add(x, node->lfdesc->dk + node->rtdesc->dk);
    node->pk = x - node->dk;

    node->h1 = node->pk + clnl;
    node->h2 = log1p(-exp(node->pk)) + node->lfdesc->loglk + node->rtdesc->loglk;
    node->loglk = logspace_add(node->h1, node->h2);

    node->rk = node->h1 - node->loglk;
}


static void cluster(struct hclust *hclust)
{
    int i;
    int j;
    int besti;
    int bestj;
    int index = hclust->nclust;
    double mx;

    struct node *node;
    struct node *left;
    struct node *right;

    while (hclust->nclust > 1) {
        mx = R_NegInf;
        node = node_alloc(index++);
        for (i = 0; i < (hclust->nclust-1); ++i) {
            node_addchild(hclust->node[i], node);
            for (j = (i+1); j < hclust->nclust; ++j) {
                node_addchild(hclust->node[j], node);
                node_compute(node, hclust);
                if (node->rk > mx) {
                    besti = i, bestj = j;
                    left = node->lfdesc;
                    right = node->rtdesc;
                    mx = node->rk;
                }
                node_rmchild(hclust->node[j], node);
            }
            node_rmchild(hclust->node[i], node);
        }
        node_addchild(left, node);
        node_addchild(right, node);
        node_setdata(node);
        node_compute(node, hclust);
        hclust->nclust -= 1;
        if (besti < bestj) {
            hclust->node[besti] = node;
            memmove(hclust->node + bestj, hclust->node + bestj + 1, (hclust->n - bestj - 1) * sizeof(struct node*));
        } else {
            hclust->node[bestj] = node;
            memmove(hclust->node + besti, hclust->node + besti + 1, (hclust->n - besti - 1) * sizeof(struct node*));
        }
    }
    hclust->root = node;
}


static void traverse1(int *idx, int *seq, int *seqpos, int *lastvisit, struct node *root)
{
    seq[*idx] = root->index + 1;
    seqpos[root->index] = *idx + 1;
    (*idx)++;
    if (root->lfdesc)
        traverse1(idx, seq, seqpos, lastvisit, root->lfdesc);
    if (root->rtdesc)
        traverse1(idx, seq, seqpos, lastvisit, root->rtdesc);
    lastvisit[root->index] = seq[*idx - 1];
}


static void traverse2(double *loglk, double *logodds, struct node *root)
{
    loglk[root->index] = root->loglk;
    logodds[root->index] = root->rk;
    if (root->lfdesc)
        traverse2(loglk, logodds, root->lfdesc);
    if (root->rtdesc)
        traverse2(loglk, logodds, root->rtdesc);
}


static void traverse3(int nnode, int *children, int *parent, struct node *root)
{
    if (root->parent)
        parent[root->index] = root->parent->index + 1;
    if (root->lfdesc) {
        children[root->index + 0 * nnode] = root->lfdesc->index + 1;
        traverse3(nnode, children, parent, root->lfdesc);
    }
    if (root->rtdesc) {
        children[root->index + 1 * nnode] = root->rtdesc->index + 1;
        traverse3(nnode, children, parent, root->rtdesc);
    }
}



SEXP rcl_bhc(SEXP alpha, SEXP beta, SEXP data)
{
    int n;
    int k;
    int idx = 0;
    struct hclust *hclust;

    SEXP result;
    SEXP names;
    SEXP dim;

    int *seq;
    int *seqpos;
    int *lastvisit;
    int *children;
    int *parent;
    double *loglk;
    double *logodds;

    k = INTEGER(getAttrib(data, R_DimSymbol))[0];
    n = INTEGER(getAttrib(data, R_DimSymbol))[1];
    hclust = hclust_alloc(n, k, INTEGER(data), REAL(alpha)[0], REAL(beta)[0]);
    cluster(hclust);

    result = PROTECT(allocVector(VECSXP, 8));

    SET_VECTOR_ELT(result, 0, allocVector(INTSXP, 2 * hclust->n - 1));
    SET_VECTOR_ELT(result, 1, allocVector(INTSXP, 2 * hclust->n - 1));
    SET_VECTOR_ELT(result, 2, allocVector(INTSXP, 2 * hclust->n - 1));
    SET_VECTOR_ELT(result, 3, allocVector(INTSXP, 2 * (2 * hclust->n - 1)));
    SET_VECTOR_ELT(result, 4, allocVector(INTSXP, 2 * hclust->n - 1));
    SET_VECTOR_ELT(result, 5, allocVector(REALSXP, 2 * hclust->n - 1));
    SET_VECTOR_ELT(result, 6, allocVector(REALSXP, 2 * hclust->n - 1));
    SET_VECTOR_ELT(result, 7, ScalarInteger(hclust->root->index + 1));

    seq = INTEGER(VECTOR_ELT(result, 0));
    seqpos = INTEGER(VECTOR_ELT(result, 1));
    lastvisit = INTEGER(VECTOR_ELT(result, 2));
    children = INTEGER(VECTOR_ELT(result, 3));
    parent = INTEGER(VECTOR_ELT(result, 4));
    loglk = REAL(VECTOR_ELT(result, 5));
    logodds = REAL(VECTOR_ELT(result, 6));

    memset(children, 0, 2 * (2*hclust->n - 1) * sizeof(int));
    memset(parent, 0, (2*hclust->n - 1) * sizeof(int));

    traverse1(&idx, seq, seqpos, lastvisit, hclust->root);
    traverse2(loglk, logodds, hclust->root);
    traverse3(2*hclust->n-1, children, parent, hclust->root);

    hclust_free(hclust);

    names = PROTECT(allocVector(STRSXP, 8));

    SET_STRING_ELT(names, 0, mkChar("seq"));
    SET_STRING_ELT(names, 1, mkChar("seqpos"));
    SET_STRING_ELT(names, 2, mkChar("lastvisit"));
    SET_STRING_ELT(names, 3, mkChar("children"));
    SET_STRING_ELT(names, 4, mkChar("parent"));
    SET_STRING_ELT(names, 5, mkChar("loglk"));
    SET_STRING_ELT(names, 6, mkChar("logodds"));
    SET_STRING_ELT(names, 7, mkChar("root"));

    dim = PROTECT(allocVector(INTSXP, 2));
    INTEGER(dim)[0] = 2 * hclust->n - 1;
    INTEGER(dim)[1] = 2;

    setAttrib(VECTOR_ELT(result, 3), R_DimSymbol, dim);
    setAttrib(result, R_NamesSymbol, names);

    UNPROTECT(3);
    return result;
}

