#ifndef RCL_ALIAS_TBLSPC_H_
#define RCL_ALIAS_TBLSPC_H_

#include <R.h>
#include <Rinternals.h>

/*
** Alias method for O(1) sampling from a discrete probability distribution
*/

struct atblspc;
struct atblspc *atblspc_alloc(int);
void atblspc_init(double*, struct atblspc*);
void atblspc_free(struct atblspc*);
void atblspc_set(int, double*, int*, double*);
int atblspc_choose(struct atblspc*);
SEXP rcl_atblspc_build(SEXP, SEXP);
void rcl_atblspc_free(SEXP);
SEXP rcl_atblspc_choose(SEXP);

struct atblspc {
    int n;              /* size of the sample space */
    int *alias;
    double *prob;
};

#endif
