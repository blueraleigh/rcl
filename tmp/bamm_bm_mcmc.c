#ifndef BAMM_BM_MCMC_H
#define BAMM_BM_MCMC_H

#include <R.h>
#include <Rinternals.h>
#include "bamm_bm_model.h"
#include "bamm_bm_proposal.h"

static FILE *CLCKFILE;
static FILE *JMPFILE;
static FILE *MCMCFILE;


static int NGEN;
static int PFREQ;
static int WFREQ;


static void mcmc_header(struct model *model)
{
    Rprintf("%-20s", "Iteration");
    Rprintf("%-20s", "Nclock");
    Rprintf("%-20s", "Njump");
    Rprintf("%-20s", "Rate");
    Rprintf("%-20s", "Norm");
    Rprintf("%-20s", "LogLk");
    Rprintf("%-20s\n", "AcceptFreq");
    Rprintf("%-20d", 0);
    Rprintf("%-20d", model->nclck);
    Rprintf("%-20d", model->njmp);
    Rprintf("%-20.6f", model->rate);
    Rprintf("%-20.6f", 1.0);
    Rprintf("%-20.6f", model->loglk);
    Rprintf("%-20s\n", "NA");
}


static void mcmc_write1(int gen, struct model *model)
{
    char buffer[100];

    sprintf(buffer, "%d,%d,%d,%.6f,%.6f,%.6f,%.6f\n", gen, model->nclck, model->njmp, model->rate,
         model->treelen / model->treescal, model->loglk,
         (double)model->accept / (double)(model->accept + model->reject));

    fputs(buffer, MCMCFILE);
}



static void mcmc_write2(int gen, struct model *model)
{
    int i;
    char buffer[100];

    for (i = 0; i < model->nclck; ++i) {
        sprintf(buffer, "%d,%d,%.6f\n", gen, model->clckpos[i]+1, model->clckrate[i]);
        fputs(buffer, CLCKFILE);
    }
}


static void mcmc_write3(int gen, struct model *model)
{
    int i;
    char buffer[100];

    for (i = 0; i < model->njmp; ++i) {
        sprintf(buffer, "%d,%d,%.6f\n", gen, model->jmppos[i]+1, model->displ[model->jmppos[i]]);
        fputs(buffer, JMPFILE);
    }
}


static void mcmc_write(int gen, struct model *model)
{
    mcmc_write1(gen, model);
    mcmc_write2(gen, model);
    mcmc_write3(gen, model);
}


static void mcmc_print(int gen, struct model *model)
{
    Rprintf("%-20d", gen);
    Rprintf("%-20d", model->nclck);
    Rprintf("%-20d", model->njmp);
    Rprintf("%-20.6f", model->rate);
    Rprintf("%-20.6f", model->treelen / model->treescal);
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
    fclose(CLCKFILE);
    fclose(JMPFILE);
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
    memcpy(PWEIGHT, REAL(mcmc_settings_get("PWEIGHT", settings)), 6 * sizeof(double));

    for (i = 0; i < 6; ++i)
        tot += PWEIGHT[i];
    for (i = 0; i < 6; ++i)
        PWEIGHT[i] /= tot;
    for (i = 1; i < 6; ++i)
        PWEIGHT[i] += PWEIGHT[i-1];

    TUNE0 = REAL(mcmc_settings_get("TUNE0", settings))[0];
    TUNE1 = REAL(mcmc_settings_get("TUNE1", settings))[0];
    TUNE2 = REAL(mcmc_settings_get("TUNE2", settings))[0];
    TUNE3 = REAL(mcmc_settings_get("TUNE3", settings))[0];
    PRIOR2 = REAL(mcmc_settings_get("PRIOR2", settings))[0];
    PRIOR3 = REAL(mcmc_settings_get("PRIOR3", settings))[0];
    PRIOR4 = REAL(mcmc_settings_get("PRIOR4", settings))[0];
    PRIOR5 = REAL(mcmc_settings_get("PRIOR5", settings))[0];

    PRIORONLY = INTEGER(mcmc_settings_get("PRIORONLY", settings))[0];

    NGEN = INTEGER(mcmc_settings_get("NGEN", settings))[0];
    PFREQ = INTEGER(mcmc_settings_get("PFREQ", settings))[0];
    WFREQ = INTEGER(mcmc_settings_get("WFREQ", settings))[0];
    MCMCFILE = fopen(CHAR(STRING_ELT(mcmc_settings_get("MCMCFILE", settings), 0)), "w");
    CLCKFILE = fopen(CHAR(STRING_ELT(mcmc_settings_get("CLCKFILE", settings), 0)), "w");
    JMPFILE = fopen(CHAR(STRING_ELT(mcmc_settings_get("JMPFILE", settings), 0)), "w");

    if (!MCMCFILE)
        error("Cannot open MCMCFILE");
    if (!CLCKFILE)
        error("Cannot open CLCKFILE");
    if (!JMPFILE)
        error("Cannot open JMPFILE");
}


SEXP rcl_bamm_bm_mcmc(SEXP rtree, SEXP rate, SEXP data, SEXP settings)
{
    int *tipdata;
    double init_rate;
    struct tree *tree;
    struct model *model;

    GetRNGstate();

    tipdata = INTEGER(data);
    init_rate = REAL(rate)[0];
    tree = (struct tree*)R_ExternalPtrAddr(rtree);

    mcmc_settings(settings);

    model = model_init(tree);
    model_set(init_rate, tipdata, model);

    mcmc_run(model);

    model_free(model);

    PutRNGstate();

    return R_NilValue;
}


#endif
