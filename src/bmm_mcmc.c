#include <R.h>
#include <Rinternals.h>
#include "bmm_model.h"
#include "bmm_proposal.h"


static int NGEN;
static int BURNIN;
static int PFREQ;
static int WFREQ;
static int NSAMPLE;

static double *R_clockprob;     // marginal probability of new local clock for each branch
static double *R_jumpprob;      // marginal probability of jump for each branch
static double *R_branchrate;    // mean relative rate for each branch
static double *R_branchjump;    // mean jump size for each branch
static double *R_nodestate;     // mean trait value for each node
static double *R_mcmcout;


static void mcmc_write(int n, int gen, struct model *model)
{
    R_mcmcout[n + 0 * NSAMPLE] = gen;
    R_mcmcout[n + 1 * NSAMPLE] = model->nclock;
    R_mcmcout[n + 2 * NSAMPLE] = model->njump;
    R_mcmcout[n + 3 * NSAMPLE] = model->rate;
    R_mcmcout[n + 4 * NSAMPLE] = model->loglk;
    R_mcmcout[n + 5 * NSAMPLE] = (double)model->accept / (double)(model->accept + model->reject);
}


static void mcmc_header(struct model *model)
{
    Rprintf("%-20s", "Iteration");
    Rprintf("%-20s", "Nclock");
    Rprintf("%-20s", "Njump");
    Rprintf("%-20s", "Rate");
    Rprintf("%-20s", "LogLk");
    Rprintf("%-20s\n", "AcceptFreq");
    Rprintf("%-20d", 0);
    Rprintf("%-20d", model->nclock);
    Rprintf("%-20d", model->njump);
    Rprintf("%-20.6f", model->rate);
    Rprintf("%-20.6f", model->loglk);
    Rprintf("%-20s\n", "NA");
}


static void mcmc_print(int gen, struct model *model)
{
    Rprintf("%-20d", gen);
    Rprintf("%-20d", model->nclock);
    Rprintf("%-20d", model->njump);
    Rprintf("%-20.6f", model->rate);
    Rprintf("%-20.6f", model->loglk);
    Rprintf("%-20.6f\n", (double)model->accept / (double)(model->accept + model->reject));
}


static void mcmc_run(struct model *model)
{
    int i;
    int j;
    int k;
    int n = 0;

    double mu;

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
                R_jumpprob[j] += (model->jump[j] - R_jumpprob[j]) / (n + 1);
                R_branchrate[j] += (model->rltvrate[j] - R_branchrate[j]) / (n + 1);
                R_branchjump[j] += (model->mu[j] - R_branchjump[j]) / (n + 1);
                R_nodestate[j] += (model->x[j] - R_nodestate[j]) / (n + 1);
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
    memcpy(PWEIGHT, REAL(mcmc_settings_get("PWEIGHT", settings)), 4 * sizeof(double));

    for (i = 0; i < 4; ++i)
        tot += PWEIGHT[i];
    for (i = 0; i < 4; ++i)
        PWEIGHT[i] /= tot;
    for (i = 1; i < 4; ++i)
        PWEIGHT[i] += PWEIGHT[i-1];

    NGEN = INTEGER(mcmc_settings_get("NGEN", settings))[0];
    BURNIN = INTEGER(mcmc_settings_get("BURNIN", settings))[0];
    PFREQ = INTEGER(mcmc_settings_get("PFREQ", settings))[0];
    WFREQ = INTEGER(mcmc_settings_get("WFREQ", settings))[0];
    NSAMPLE = (int)((NGEN - BURNIN) / WFREQ);
    PRIORONLY = INTEGER(mcmc_settings_get("PRIORONLY", settings))[0];
    TUNE0 = REAL(mcmc_settings_get("TUNE0", settings))[0];
    TUNE1 = REAL(mcmc_settings_get("TUNE1", settings))[0];
    TUNE3 = REAL(mcmc_settings_get("TUNE3", settings))[0];
    PRIOR0 = REAL(mcmc_settings_get("PRIOR0", settings))[0];
    PRIOR1 = REAL(mcmc_settings_get("PRIOR1", settings))[0];
    PRIOR2 = REAL(mcmc_settings_get("PRIOR2", settings))[0];
    PRIOR3 = REAL(mcmc_settings_get("PRIOR3", settings))[0];
}


SEXP rcl_bmm_mcmc(SEXP rtree, SEXP rate, SEXP data, SEXP nodestate, SEXP clockprob, SEXP jumpprob, SEXP branchrate, SEXP branchjump, SEXP mcmcout, SEXP settings)
{
    double *tipdata;
    struct tree *tree;
    struct model *model;

    GetRNGstate();

    tipdata = REAL(data);
    tree = (struct tree*)R_ExternalPtrAddr(rtree);

    mcmc_settings(settings);

    model = model_init(tree);
    model_set(REAL(rate)[0], tipdata, model);

    R_mcmcout = REAL(mcmcout);
    R_clockprob = REAL(clockprob);
    R_jumpprob = REAL(jumpprob);
    R_branchrate = REAL(branchrate);
    R_branchjump = REAL(branchjump);
    R_nodestate = REAL(nodestate);

    mcmc_run(model);

    model_free(model);

    PutRNGstate();

    return R_NilValue;
}
