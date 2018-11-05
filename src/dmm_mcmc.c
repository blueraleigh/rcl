#include <R.h>
#include <Rinternals.h>
#include "dmm_model.h"
#include "dmm_proposal.h"


static int NGEN;
static int BURNIN;
static int PFREQ;
static int WFREQ;
static int NSAMPLE;

static double *R_clockprob;     // marginal probability of new local clock for each branch
static double *R_branchrate;    // mean relative rate for each branch
static double *R_nodestate;     // marginal state probabilities for each node
static double *R_stateprob;     // posterior estimates of distribution over resource use categories for each state
static double *R_mcmcout;


static void mcmc_write(int n, int gen, struct model *model)
{
    R_mcmcout[n + 0 * NSAMPLE] = gen;
    R_mcmcout[n + 1 * NSAMPLE] = model->nclock;
    R_mcmcout[n + 2 * NSAMPLE] = model->rate;
    R_mcmcout[n + 3 * NSAMPLE] = model->loglk;
    R_mcmcout[n + 4 * NSAMPLE] = model->alpha;
    R_mcmcout[n + 5 * NSAMPLE] = (double)model->accept / (double)model->reject;
}


static void mcmc_header(struct model *model)
{
    Rprintf("%-20s", "Iteration");
    Rprintf("%-20s", "Nclock");
    Rprintf("%-20s", "Rate");
    Rprintf("%-20s", "LogLk");
    Rprintf("%-20s", "Alpha");
    Rprintf("%-20s\n", "AcceptFreq");
    Rprintf("%-20d", 0);
    Rprintf("%-20d", model->nclock);
    Rprintf("%-20.6f", model->rate);
    Rprintf("%-20.6f", model->loglk);
    Rprintf("%-20.6f", model->alpha);
    Rprintf("%-20s\n", "NA");
}


static void mcmc_print(int gen, struct model *model)
{
    Rprintf("%-20d", gen);
    Rprintf("%-20d", model->nclock);
    Rprintf("%-20.6f", model->rate);
    Rprintf("%-20.6f", model->loglk);
    Rprintf("%-20.6f", model->alpha);
    Rprintf("%-20.6f\n", (double)model->accept / (double)(model->accept + model->reject));
}


static void mcmc_run(struct model *model)
{
    int i;
    int j;
    int k;
    int n = 0;

    double tot;
    double p[model->ncat];

    mcmc_header(model);

    for (i = 0; i < NGEN; ++i) {
        propose(model);

        if (i % 1024 == 0) {
            model->accept = 0;
            model->reject = 0;
        }

        if (PFREQ && (i+1) % PFREQ == 0)
            mcmc_print(i+1, model);

        if (WFREQ && (i+1) > BURNIN && (i+1) % WFREQ == 0) {
            mcmc_write(n, i+1, model);
            for (j = 0; j < model->tree->nnode; ++j) {
                R_clockprob[j] += (model->clock[j] - R_clockprob[j]) / (n + 1);
                R_branchrate[j] += (model->rltvrate[j] - R_branchrate[j]) / (n + 1);
                for (k = 0; k < model->nstate; ++k) {
                    if (model->stateid[j] == k)
                        R_nodestate[k + j * model->nstate] +=
                            (1 - R_nodestate[k + j * model->nstate]) / (n + 1);
                    else
                        R_nodestate[k + j * model->nstate] +=
                            (0 - R_nodestate[k + j * model->nstate]) / (n + 1);
                }
            }
            for (j = 0; j < model->nstate; ++j) {
                tot = 0;
                for (k = 0; k < model->ncat; ++k)
                    tot += p[k] = model->count[k + j * model->ncat] + model->alpha;
                for (k = 0; k < model->ncat; ++k)
                    R_stateprob[j + k * model->nstate] +=
                        (p[k]/tot - R_stateprob[j + k * model->nstate]) / (n + 1);
            }
            ++n;
        }
    }

}


static SEXP mcmc_settings_get(const char *name, SEXP settings)
{
    int i;
    int len;
    SEXP names;
    names = getAttrib(settings, R_NamesSymbol);
    len = LENGTH(names);
    for (i = 0; i < len; ++i) {
        if (strcmp(name, CHAR(STRING_ELT(names, i))) == 0) {
            return VECTOR_ELT(settings, i);
        }
    }
    return R_NilValue;
}


// Set all the global variables associated with the MCMC
static void mcmc_settings(SEXP settings)
{
    int i;
    double tot = 0;
    memcpy(PWEIGHT, REAL(mcmc_settings_get("PWEIGHT", settings)), 3 * sizeof(double));

    for (i = 0; i < 3; ++i)
        tot += PWEIGHT[i];
    for (i = 0; i < 3; ++i)
        PWEIGHT[i] /= tot;
    for (i = 1; i < 3; ++i)
        PWEIGHT[i] += PWEIGHT[i-1];

    NGEN = INTEGER(mcmc_settings_get("NGEN", settings))[0];
    BURNIN = INTEGER(mcmc_settings_get("BURNIN", settings))[0];
    PFREQ = INTEGER(mcmc_settings_get("PFREQ", settings))[0];
    WFREQ = INTEGER(mcmc_settings_get("WFREQ", settings))[0];
    NSAMPLE = (int)((NGEN - BURNIN) / WFREQ);
    PRIORONLY = INTEGER(mcmc_settings_get("PRIORONLY", settings))[0];
    TUNE1 = REAL(mcmc_settings_get("TUNE1", settings))[0];
    TUNE2 = REAL(mcmc_settings_get("TUNE2", settings))[0];
    PRIOR1 = REAL(mcmc_settings_get("PRIOR1", settings))[0];
    PRIOR2 = REAL(mcmc_settings_get("PRIOR2", settings))[0];
}


SEXP rcl_dmm_mcmc(SEXP rtree, SEXP rate, SEXP alpha, SEXP data, SEXP stateprob, SEXP nodestate, SEXP clockprob, SEXP branchrate, SEXP mcmcout, SEXP settings)
{
    int ncat;
    int nstate;
    int *tipdata;
    struct tree *tree;
    struct model *model;

    GetRNGstate();

    ncat = INTEGER(getAttrib(data, R_DimSymbol))[0];
    tipdata = INTEGER(data);
    tree = (struct tree*)R_ExternalPtrAddr(rtree);

    mcmc_settings(settings);
    nstate = INTEGER(mcmc_settings_get("NSTATE", settings))[0];

    model = model_init(ncat, nstate, tree);
    model_set(REAL(rate)[0], REAL(alpha)[0], tipdata, model);

    R_mcmcout = REAL(mcmcout);
    R_clockprob = REAL(clockprob);
    R_branchrate = REAL(branchrate);
    R_stateprob = REAL(stateprob);
    R_nodestate = REAL(nodestate);

    mcmc_run(model);

    model_free(model);

    PutRNGstate();

    return R_NilValue;
}
