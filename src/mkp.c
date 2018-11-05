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
// matrix (with number of rows equal to nstate) stored as
// a flattened array. these are only used in function contexts
// where the local variable nstate is defined.
#define COLUMN(m, j) m + nstate * (j)
#define ELEM(m, i, j) m[(i) + nstate * (j)]


static double mkp_branch_prob(int nstate, double brlen, double *R, double *init, double *out)
{
    int ione = 1;
    double Pr[nstate * nstate], sf, isf;

    // compute transition probability matrix
    matrix_exponential(R, brlen, nstate, Pr);

    // multiply that matrix with the initial conditional likelihoods
    // in the init vector, storing the new conditional likelihoods
    // in the out vector
    matrix_vector_product(Pr, nstate, nstate, init, out);

    // compute the scale factor, which is just the sum of the new
    // conditional likelihoods
    sf = F77_CALL(dasum)(&nstate, out, &ione);

    // rescale the conditional likelihoods by this factor
    isf = 1 / sf;
    F77_CALL(dscal)(&nstate, &isf, out, &ione);

    return sf;
}


static double mkp_node_lk(int nstate, double *R, double *clk, struct node *node)
{
    int j;
    int repeat;
    int remainder;
    double *init;
    double lsf;
    double rsf;
    double lf[nstate];
    double rt[nstate];

    // do the left branch
    init = COLUMN(clk, node->children[0]->index);
    lsf = mkp_branch_prob(nstate, node->children[0]->brlen,
        R, init, lf);

    // do the right branch
    init = COLUMN(clk, node->children[1]->index);
    rsf = mkp_branch_prob(nstate, node->children[1]->brlen,
        R, init, rt);

    j = 0;
    repeat = nstate / 5;
    remainder = nstate % 5;

    // to form the conditional likelihoods for the ancestral node
    // we multiply the contributions from each daughter
    while (j < repeat) {
        ELEM(clk, j+0, node->index) = lf[j+0] * rt[j+0];
        ELEM(clk, j+1, node->index) = lf[j+1] * rt[j+1];
        ELEM(clk, j+2, node->index) = lf[j+2] * rt[j+2];
        ELEM(clk, j+3, node->index) = lf[j+3] * rt[j+3];
        ELEM(clk, j+4, node->index) = lf[j+4] * rt[j+4];
        j += 5;
    }

    switch (remainder) {
        case 4:
            ELEM(clk, j+3, node->index) = lf[j+3] * rt[j+3];
        case 3:
            ELEM(clk, j+2, node->index) = lf[j+2] * rt[j+2];
        case 2:
            ELEM(clk, j+1, node->index) = lf[j+1] * rt[j+1];
        case 1:
            ELEM(clk, j+0, node->index) = lf[j+0] * rt[j+0];
        case 0: ;
    }

    // add up the logarithms of the scale factors as we will
    // use this for computing the final (unscaled) log likelihood
    return log(lsf) + log(rsf);
}


static double mkp_loglk(int nstate, double *R, double *clk, struct tree *tree)
{
    int j;
    int ione = 1;
    double p = 0;
    double tot = 0;
    double loglk = 0;
    double ls = 0;
    struct tree_traversal t;
    struct node *node;

    t = tree_traverse(POSTORDER, INTERNAL_NODES_ONLY, tree->root, tree);
    node = tree_step(&t);

    while (node != NULL) {
        ls += mkp_node_lk(nstate, R, clk, node);
        node = tree_step(&t);
    }

    // sum the (scaled) conditional likelihoods at the root of the tree, weighting
    // each term by its proportional contribution to the total sum
    tot = F77_CALL(dasum)(&nstate, COLUMN(clk, tree->root->index), &ione);
    for (j = 0; j < nstate; ++j)
        p += (ELEM(clk, j, tree->root->index) / tot) * ELEM(clk, j, tree->root->index);
    loglk = log(p);

    // finally, add in the log of the total scale factor to get the unscaled log likelihood
    loglk += ls;
    return loglk;
}


SEXP rcl_mkp_loglk(SEXP rtree, SEXP R, SEXP clk)
{
    struct tree *tree = (struct tree*)R_ExternalPtrAddr(rtree);
    int nstate = INTEGER(getAttrib(R, R_DimSymbol))[0];
    return ScalarReal(mkp_loglk(nstate, REAL(R), REAL(clk), tree));
}
