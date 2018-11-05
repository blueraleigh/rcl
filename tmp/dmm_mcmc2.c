#include <R.h>
#include <Rinternals.h>
#include "dmm_model2.h"
#include "dmm_proposal2.h"

static FILE *EVNTFILE;
static FILE *MCMCFILE;

static int NGEN;
static int PFREQ;
static int WFREQ;


static void mcmc_header(struct model *model)
{
    Rprintf("%-20s", "Iteration");
    Rprintf("%-20s", "Nevent");
    Rprintf("%-20s", "LogLk");
    Rprintf("%-20s\n", "AcceptFreq");
    Rprintf("%-20d", 0);
    Rprintf("%-20d", model->nevent);
    Rprintf("%-20.6f", model->loglk);
    Rprintf("%-20s\n", "NA");
}


static void mcmc_write1(int gen, struct model *model)
{
    char buffer[100];
    sprintf(buffer, "%d,%d,%.6f,%.6f\n", gen, model->nevent, model->loglk,
         (double)model->accept / (double)(model->accept + model->reject));
    fputs(buffer, MCMCFILE);
}


static void mcmc_write2(int gen, struct model *model)
{
    int i;
    char buffer[100];

    for (i = 0; i < model->nevent; ++i) {
        sprintf(buffer, "%d,%d\n", gen, model->eventpos[i]+1);
        fputs(buffer, EVNTFILE);
    }
}


static void mcmc_write(int gen, struct model *model)
{
    mcmc_write1(gen, model);
    mcmc_write2(gen, model);
}


static void mcmc_print(int gen, struct model *model)
{
    Rprintf("%-20d", gen);
    Rprintf("%-20d", model->nevent);
    Rprintf("%-20.6f", model->loglk);
    Rprintf("%-20.6f\n", (double)model->accept / (double)(model->accept + model->reject));
}


static void mcmc_run(struct model *model)
{
    int i;

    model->loglk = model_loglk(model);
    mcmc_header(model);

    for (i = 0; i < NGEN; ++i) {
        propose(model);

        if (i % 1024 == 0) {
            model->accept = 0;
            model->reject = 0;
        }

        if ((i+1) % PFREQ == 0)
            mcmc_print(i+1, model);

        if ((i+1) % WFREQ == 0)
            mcmc_write(i+1, model);
    }

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

    PRIOR0 = REAL(mcmc_settings_get("PRIOR0", settings))[0];

    PRIORONLY = INTEGER(mcmc_settings_get("PRIORONLY", settings))[0];

    NGEN = INTEGER(mcmc_settings_get("NGEN", settings))[0];
    PFREQ = INTEGER(mcmc_settings_get("PFREQ", settings))[0];
    WFREQ = INTEGER(mcmc_settings_get("WFREQ", settings))[0];
    MCMCFILE = fopen(CHAR(STRING_ELT(mcmc_settings_get("MCMCFILE", settings), 0)), "w");
    EVNTFILE = fopen(CHAR(STRING_ELT(mcmc_settings_get("EVNTFILE", settings), 0)), "w");

    if (!MCMCFILE)
        error("Cannot open MCMCFILE");
    if (!EVNTFILE)
        error("Cannot open EVNTFILE");
}


SEXP rcl_dmm_mcmc2(SEXP rtree, SEXP alpha, SEXP data, SEXP settings)
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
    model_set(REAL(alpha)[0], tipdata, model);

    mcmc_run(model);

    model_free(model);

    PutRNGstate();

    return R_NilValue;
}
