#include <R.h>
#include <Rinternals.h>
#include "dmm_model4.h"
#include "dmm_proposal4.h"

static FILE *EVNTFILE;
static FILE *MCMCFILE;

static int NGEN;
static int BURNIN;
static int PFREQ;
static int WFREQ;
static int NSAMPLE;

static int *R_eventconfig;
static double *R_eventprob;
static double *R_mcmcout;


static void mcmc_write(int n, int gen, struct model *model)
{
    int i;

    memcpy(R_eventconfig + n * model->tree->nnode, model->event, model->tree->nnode * sizeof(int));

    R_mcmcout[n + 0 * NSAMPLE] = gen;
    R_mcmcout[n + 1 * NSAMPLE] = model->nevent;
    R_mcmcout[n + 2 * NSAMPLE] = model->loglk;
    R_mcmcout[n + 3 * NSAMPLE] = model->logprior;
    R_mcmcout[n + 4 * NSAMPLE] = model->alpha;
    R_mcmcout[n + 5 * NSAMPLE] = (double)model->accept / (double)model->reject;
}


static void mcmc_header(struct model *model)
{
    Rprintf("%-20s", "Iteration");
    Rprintf("%-20s", "Nevent");
    Rprintf("%-20s", "LogLk");
    Rprintf("%-20s", "LogPr");
    Rprintf("%-20s", "Alpha");
    Rprintf("%-20s\n", "AcceptFreq");
    Rprintf("%-20d", 0);
    Rprintf("%-20d", model->nevent);
    Rprintf("%-20.6f", model->loglk);
    Rprintf("%-20.6f", model->logprior);
    Rprintf("%-20.6f", model->alpha);
    Rprintf("%-20s\n", "NA");
}


static void mcmc_write1(int gen, struct model *model)
{
    char buffer[100];
    sprintf(buffer, "%d,%d,%.6f,%.6f,%.6f,%.6f\n", gen, model->nevent, model->loglk, model->logprior,
         model->alpha, (double)model->accept / (double)(model->accept + model->reject));
    fputs(buffer, MCMCFILE);
}



static void mcmc_write2(struct model *model)
{
    int i;
    char buffer[100];

    for (i = 0; i < model->tree->nnode; ++i) {
        sprintf(buffer, "%d,%f,%f\n", i+1, model->margprobevnt[i], 0.5);
        fputs(buffer, EVNTFILE);
    }
}


static void mcmc_write3(int gen, struct model *model)
{
    int i;
    char buffer[100];

    for (i = 0; i < model->tree->nnode; ++i) {
        if (model->event[i]) {
            sprintf(buffer, "%d,%d\n", gen, i+1);
            fputs(buffer, EVNTFILE);
        }
    }
}


static void mcmc_print(int gen, struct model *model)
{
    Rprintf("%-20d", gen);
    Rprintf("%-20d", model->nevent);
    Rprintf("%-20.6f", model->loglk);
    Rprintf("%-20.6f", model->logprior);
    Rprintf("%-20.6f", model->alpha);
    Rprintf("%-20.6f\n", (double)model->accept / (double)(model->accept + model->reject));
}


static void mcmc_run(struct model *model)
{
    int i;
    int j;
    int n = 0;

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
            //mcmc_write1(i+1, model);
            //mcmc_write3(i+1, model);
            mcmc_write(n, i+1, model);
            for (j = 0; j < model->tree->nnode; ++j)
                model->margprobevnt[j] += (model->event[j] - model->margprobevnt[j]) / (n + 1);
            ++n;
        }
    }

    memcpy(R_eventprob, model->margprobevnt, model->tree->nnode * sizeof(double));
    //mcmc_write2(model);

    fclose(MCMCFILE);
    fclose(EVNTFILE);
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
    MCMCFILE = fopen(CHAR(STRING_ELT(mcmc_settings_get("MCMCFILE", settings), 0)), "w");
    EVNTFILE = fopen(CHAR(STRING_ELT(mcmc_settings_get("EVNTFILE", settings), 0)), "w");

    if (!MCMCFILE)
        error("Cannot open MCMCFILE");
    if (!EVNTFILE)
        error("Cannot open EVNTFILE");
}


SEXP rcl_dmm_mcmc4(SEXP rtree, SEXP alpha, SEXP beta, SEXP data, SEXP eventconfig, SEXP eventprob, SEXP mcmcout, SEXP settings)
{
    int ncat;
    int *tipdata;
    struct tree *tree;
    struct model *model;

    GetRNGstate();

    ncat = INTEGER(getAttrib(data, R_DimSymbol))[0];
    tipdata = INTEGER(data);
    tree = (struct tree*)R_ExternalPtrAddr(rtree);

    mcmc_settings(settings);

    model = model_init(ncat, tree);
    model_set(REAL(alpha)[0], REAL(beta)[0], tipdata, model);

    R_mcmcout = REAL(mcmcout);
    R_eventconfig = INTEGER(eventconfig);
    R_eventprob = REAL(eventprob);

    mcmc_run(model);

    model_free(model);

    PutRNGstate();

    return R_NilValue;
}
