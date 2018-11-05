#include "atblspc.h"


struct atblspc *atblspc_alloc(int n)
{
    struct atblspc *a = Calloc(1, struct atblspc);
    a->n = n;
    a->alias = Calloc(n, int);
    a->prob = Calloc(n, double);
    return a;
}


void atblspc_init(double *w, struct atblspc *a)
{
    atblspc_set(a->n, w, a->alias, a->prob);
}


void atblspc_free(struct atblspc *a)
{
    if (a) {
        Free(a->prob);
        Free(a->alias);
        Free(a);
    }
}


int atblspc_choose(struct atblspc *a)
{
    int j = (int)(a->n * unif_rand());
    return (unif_rand() < a->prob[j]) ? j : a->alias[j];
}


SEXP rcl_atblspc_build(SEXP n, SEXP w)
{
    struct atblspc *a = atblspc_alloc(INTEGER(n)[0]);
    atblspc_init(REAL(w), a);
    SEXP r = PROTECT(R_MakeExternalPtr(a, R_NilValue, R_NilValue));
    R_RegisterCFinalizer(r, &rcl_atblspc_free);
    UNPROTECT(1);
    return r;
}


void rcl_atblspc_free(SEXP r)
{
    if (TYPEOF(r) == NILSXP)
        return;
    struct atblspc *a = (struct atblspc*)R_ExternalPtrAddr(r);
    atblspc_free(a);
    R_ClearExternalPtr(r);
}


SEXP rcl_atblspc_choose(SEXP r)
{
    struct atblspc *a = (struct atblspc*)R_ExternalPtrAddr(r);
    return ScalarInteger(atblspc_choose(a));
}


/*
** Trivially modified from original source code here: apps.jcns.fz-juelich.de/ransampl
**   (FreeBSD Copyright (c) 2013 Joachim Wuttke, Forschungszentrum Juelich GmbH)
** with detailed explanation of the method
** here: http://www.keithschwarz.com/darts-dice-coins
*/
void atblspc_set(int n, double *w, int *alias, double *prob)
{
    int i, a, g, S[n], L[n], nS = 0, nL = 0;
    double tot = 0, P[n];

    for (i = 0; i < n; ++i)
        tot += w[i];

    if (!tot)
        error("No non-zero sample weights");

    for (i = 0; i < n; ++i)
        P[i] = w[i] * n / tot;

    for (i = n-1; i >= 0; --i) {
        if (P[i] < 1)
            S[nS++] = i;
        else
            L[nL++] = i;
    }

    while (nS && nL) {
        a = S[--nS];
        g = L[--nL];
        prob[a] = P[a];
        alias[a] = g;
        P[g] = P[g] + P[a] - 1;
        if (P[g] < 1)
            S[nS++] = g;
        else
            L[nL++] = g;
    }

    while (nL)
        prob[L[--nL]] = 1;

    while (nS)
        prob[S[--nS]] = 1;
}
