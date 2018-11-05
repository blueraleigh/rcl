#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include "node.h"
#include "tree.h"
#include "matrix.h"


/**********************************************************************
**
** Markov process model of discrete character macroevolution likelihood
** functions.
**
**********************************************************************/


// macros for column and element acess into a column-major
// matrix (with number of rows equal to lda) stored as
// a flattened array. these are only used in function contexts
// where the local variable lda is defined.
#define COLUMN(m, j) (m) + lda * (j)
#define ELEM(m, i, j) (m)[(i) + lda * (j)]


static void mkp_branch_prob(int nstate, int nsite, double brlen, double *R, double *clk, double *out, double *sf)
{
    int i, ione = 1;

    double Pr[nstate * nstate], isf;

    // compute transition probability matrix
    matrix_exponential(R, brlen, nstate, Pr);

    // multiply that matrix with the initial conditional likelihoods
    // in the clk vector, storing the new conditional likelihoods
    // in the out vector and rescaling
    for (i = 0; i < nsite; ++i) {
        matrix_vector_product(Pr, nstate, nstate, clk + i*nstate, out + i*nstate);
        sf[i] = F77_CALL(dasum)(&nstate, out + i*nstate, &ione);
        isf = 1 / sf[i];
        F77_CALL(dscal)(&nstate, &isf, out + i*nstate, &ione);
    }
}


static void mkp_node_lk(int nstate, int nsite, double *R, double *clk, double *ls, struct node *node)
{
    int i, j, lda = nstate * nsite;
    double *init;
    double lsf[nsite], rsf[nsite], lf[nstate*nsite], rt[nstate*nsite];

    // do the left branch
    init = COLUMN(clk, node->children[0]->index);
    mkp_branch_prob(nstate, nsite, node->children[0]->brlen,
        R, init, lf, lsf);

    // do the right branch
    init = COLUMN(clk, node->children[1]->index);
    mkp_branch_prob(nstate, nsite, node->children[1]->brlen,
        R, init, rt, rsf);

    // to form the conditional likelihoods for the ancestral node
    // we multiply the contributions from each daughter
    for (i = 0; i < nsite; ++i) {
        for (j = 0; j < nstate; ++j)
            ELEM(clk+i*nstate, j, node->index) = lf[j+i*nstate] * rt[j+i*nstate];
        ls[i] += log(lsf[i]) + log(rsf[i]);
    }
}


static double mkp_loglk(int nstate, int nsite, double *R, double *clk, struct tree *tree)
{
    int i, j, ione = 1, lda = nstate * nsite;
    double p, tot, tmp, loglk = 0, ls[nsite];
    struct tree_traversal t;
    struct node *node;

    memset(ls, 0, nsite * sizeof(double));

    t = tree_traverse(POSTORDER, INTERNAL_NODES_ONLY, tree->root, tree);
    node = tree_step(&t);

    while (node != NULL) {
        mkp_node_lk(nstate, nsite, R, clk, ls, node);
        node = tree_step(&t);
    }

    // sum the (scaled) conditional likelihoods at the root of the tree, weighting
    // each term by its proportional contribution to the total sum
    for (i = 0; i < nsite; ++i) {
        p = 0;
        tot = F77_CALL(dasum)(&nstate, COLUMN(clk, tree->root->index) + i*nstate, &ione);
        for (j = 0; j < nstate; ++j) {
            tmp = ELEM(clk+i*nstate, j, tree->root->index);
            p += (tmp / tot) * tmp;
        }
        loglk += log(p) + ls[i];
    }

    return loglk;
}


SEXP rcl_mkpmulti_loglk(SEXP rtree, SEXP NSITE, SEXP R, SEXP clk)
{
    struct tree *tree = (struct tree*)R_ExternalPtrAddr(rtree);
    int nstate = INTEGER(getAttrib(R, R_DimSymbol))[0];
    int nsite = INTEGER(NSITE)[0];
    return ScalarReal(mkp_loglk(nstate, nsite, REAL(R), REAL(clk), tree));
}
