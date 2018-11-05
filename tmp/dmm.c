#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

/*
** Gibbs sampler for Dirichlet multinomial mixture model
*/


struct lda {

    /* Number of documents (=species) */
    int ndoc;

    /* Number of terms (=resource categories) */
    int nterm;

    /* Number of topics (=guilds) */
    int ntop;

    /* Number of Gibbs sweeps */
    int nsweep;

    /*
    ** Number of words in each document (document size).
    **
    ** Sum of frequencies of terms appearing in each doc.
    */
    int *docsz;

    /* Number of unique words in each document (document degree) */
    int *docdg;

    /*
    ** Offset of each document in the data arrays
    **
    ** That is, docdat + docofs[i] returns a
    ** pointer to the first term in document i
    */
    int *docofs;

    /*
    ** Term index associated with each term in each
    ** document. For example,
    **
    **      docdat[docofs[i] + j]
    **
    ** returns the term index associated with the j-th
    ** unique word in document i.
    */
    int *docdat;

    /*
    ** Frequency of each term in each document
    */
    int *doctf;

    /*
    ** Topic assigment for each document.
    */
    int *doctop;


    /*
    ** Each entry (i, j) is number of words with term index
    ** j that are assigned to topic i (across all documents).
    ** ntop by nterm matrix.
    */
    int *nzw;


    /*
    ** Each entry i is the number of documents assigned
    ** to the i-th topic. ntop by 1 matrix.
    */
    int *mz;


    /*
    ** Each entry i is the number of words assigned
    ** to the i-th topic. ntop by 1 matrix.
    */
    int *wz;


    /*
    ** Symmetric dirichlet prior on the topic
    ** distribution over documents
    */
    double alpha;


    /*
    ** Symmetric dirichlet prior on the topic
    ** distribution over words
    */
    double beta;

    /*
    ** ntop by 1 matrix that records the average
    ** distribution over topics.
    */
    double *theta;

    /*
    ** ntop by nterm matrix that records the average
    ** distribution over words for each topic.
    */
    double *phi;


    /*
    ** ndoc by ntop matrix that records the average
    ** distribution over topics for each doc.
    */
    double *docphi;
};


/*
** Draw a random integer from [0, size)
** using the probabilities in *prob.
*/
static int choose_one(int size, double *prob)
{
    int i;
    double p_tot = 1;
    for (i = 0; i < size-1; ++i) {
        if (unif_rand() < prob[i] / p_tot)
            return i;
        else
            p_tot -= prob[i];
    }
    return size-1;
}


static void normalize(int size, double *x)
{
    int j;
    double maxp = R_NegInf, tot = 0;

    for (j = 0; j < size; ++j) {
        if (x[j] > maxp)
            maxp = x[j];
    }

    for (j = 0; j < size; ++j) {
        if (ISNAN(x[j]))
            x[j] = 0.0;
        else
            tot += x[j] = exp(x[j] - maxp);
    }

    for (j = 0; j < size; ++j)
        x[j] /= tot;
}


/*
** Draw a new topic assignment for the i-th document
*/
static void draw(int i, struct lda *lda)
{
    int j;
    int k;
    int z;
    int w;
    double tmp1;
    double tmp2;
    double tot = 0;
    double prob[lda->ntop];

    // current topic assignment for this doc
    z = lda->doctop[i];

    lda->mz[z] -= 1;
    lda->wz[z] -= lda->docsz[i];
    for (j = 0; j < lda->docdg[i]; ++j) {
        w = lda->docdat[lda->docofs[i] + j];
        lda->nzw[z + w * lda->ntop] -= lda->doctf[lda->docofs[i] + j];
    }

    for (k = 0; k < lda->ntop; ++k) {
        tmp1 = tmp2 = 0;

        for (j = 0; j < lda->docdg[i]; ++j) {
            w = lda->docdat[lda->docofs[i] + j];
            tmp1 += log(lda->nzw[k + w * lda->ntop] + lda->beta + j - 1);
        }

        for (j = 0; j < lda->docsz[i]; ++j)
            tmp2 += log(lda->wz[k] + lda->nterm * lda->beta + j - 1);

        prob[k] = log(lda->mz[k] + lda->alpha) + tmp1 - tmp2 - log(lda->ndoc - 1 + lda->ntop * lda->alpha);
        tot += prob[k];
    }

    normalize(lda->ntop, prob);
    z = choose_one(lda->ntop, prob);

    lda->mz[z] += 1;
    lda->wz[z] += lda->docsz[i];
    for (j = 0; j < lda->docdg[i]; ++j) {
        w = lda->docdat[lda->docofs[i] + j];
        lda->nzw[z + w * lda->ntop] += lda->doctf[lda->docofs[i] + j];
    }

    // store new topic assignment for this doc
    lda->doctop[i] = z;
}

/* Perform one sweep of the LDA Gibb's sampler */
static void sweep(int iter, struct lda *lda)
{
    int i;
    for (i = 0; i < lda->ndoc; ++i)
        draw(i, lda);
}


/*
** Update the partial estimate of the average topic
** distribution for the i-th topic
*/
static void doctheta(double n, struct lda *lda)
{
    int k;
    double denom = lda->ndoc + lda->ntop * lda->alpha;
    for (k = 0; k < lda->ntop; ++k)
        lda->theta[k] += (lda->mz[k] + lda->alpha) / (n*denom);
}


/*
** Update the partial estimate of the average distribution over terms
** for the i-th topic
*/
static void topicphi(int i, double n, struct lda *lda)
{
    int w;
    double denom = lda->nterm * lda->beta + lda->wz[i];
    for (w = 0; w < lda->nterm; ++w)
        lda->phi[i + w * lda->ntop] += (lda->nzw[i + w * lda->ntop] + lda->beta) / (n*denom);
}


static void docphi(int i, double n, struct lda *lda)
{
    int k;
    k = lda->doctop[i];
    lda->docphi[i + k * lda->ndoc] += 1;
}


/*
** Initialize the lda structure.
**
** Data is assumed to be an nterm by ndoc matrix where each element
** (i, j) records the count of the i-th term in the j-th doc
*/
static struct lda *lda_init(double alpha, double beta, int nsweep, int ntop, int ndoc, int nterm, int *data)
{
    int i;
    int j;
    int w;
    int z;
    int nwrd = 0;
    int ntrm = 0;
    struct lda *lda;

    lda = Calloc(1, struct lda);
    lda->ntop = ntop;
    lda->ndoc = ndoc;
    lda->nterm = nterm;
    lda->nsweep = nsweep;
    lda->alpha = alpha;
    lda->beta = beta;
    lda->nzw = Calloc(ntop * nterm, int);
    lda->mz = Calloc(ntop, int);
    lda->wz = Calloc(ntop, int);
    lda->docsz = Calloc(ndoc, int);
    lda->docdg = Calloc(ndoc, int);
    lda->docofs = Calloc(ndoc, int);
    lda->doctop = Calloc(ndoc, int);
    lda->theta = Calloc(ntop, double);
    lda->phi = Calloc(ntop * nterm, double);
    lda->docphi = Calloc(ndoc * ntop, double);

    memset(lda->nzw, 0, ntop * nterm * sizeof(int));
    memset(lda->mz, 0, ntop * sizeof(int));
    memset(lda->wz, 0, ntop * sizeof(int));
    memset(lda->docsz, 0, ndoc * sizeof(int));
    memset(lda->docdg, 0, ndoc * sizeof(int));
    memset(lda->theta, 0, ntop * sizeof(double));
    memset(lda->phi, 0, ntop * nterm * sizeof(double));
    memset(lda->docphi, 0, ndoc * ntop * sizeof(double));

    for (j = 0; j < ndoc; ++j) {
        lda->docofs[j] = ntrm;
        for (i = 0; i < nterm; ++i) {
            if (data[i + j * nterm]) {
                lda->docsz[j] += data[i + j * nterm];
                lda->docdg[j] += 1;
                nwrd += data[i + j * nterm];
                ntrm += 1;
            }
        }
    }

    lda->docdat = Calloc(ntrm, int);
    lda->doctf = Calloc(ntrm, int);

    for (j = 0; j < ndoc; ++j) {
        w = 0;
        z = (int)(unif_rand() * ntop);
        lda->doctop[j] = z;
        lda->mz[z] += 1;
        lda->wz[z] += lda->docsz[j];
        for (i = 0; i < nterm; ++i) {
            if (w < lda->docdg[j] && data[i + j * nterm] > 0) {
                lda->docdat[lda->docofs[j] + w] = i;
                lda->doctf[lda->docofs[j] + w] = data[i + j * nterm];
                lda->nzw[z + i * ntop] += data[i + j * nterm];
                w += 1;
            }
        }
    }

    return lda;
}


static void lda_free(struct lda *lda)
{
    Free(lda->docsz);
    Free(lda->docofs);
    Free(lda->docdg);
    Free(lda->docdat);
    Free(lda->doctf);
    Free(lda->doctop);
    Free(lda->nzw);
    Free(lda->mz);
    Free(lda->wz);
    Free(lda->theta);
    Free(lda->phi);
    Free(lda->docphi);
    Free(lda);
}


static void lda_run(int burnin, int nsweep, int sample_freq, struct lda *lda)
{
    int i;
    int iter;
    int nsample = (nsweep - burnin) / sample_freq;
    for (iter = 0; iter < nsweep; ++iter) {
        sweep(iter, lda);
        if ((iter + 1) > burnin && (iter+1) % sample_freq == 0) {
            doctheta((double)nsample, lda);
            for (i = 0; i < lda->ntop; ++i)
                topicphi(i, (double)nsample, lda);
            for (i = 0; i < lda->ndoc; ++i)
                docphi(i, (double)nsample, lda);
        }
    }
}


static SEXP dmm_settings_get(const char *name, SEXP settings)
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


SEXP rcl_gsdmm(SEXP control, SEXP data)
{
    GetRNGstate();

    int ntop;
    int ndoc;
    int nterm;
    int nsweep;
    int burnin;
    int sample_freq;
    double alpha;
    double beta;
    struct lda *lda;
    SEXP result;

    nterm = INTEGER(getAttrib(data, R_DimSymbol))[0];
    ndoc = INTEGER(getAttrib(data, R_DimSymbol))[1];
    ntop = INTEGER(dmm_settings_get("ntopic", control))[0];
    nsweep = INTEGER(dmm_settings_get("nsweep", control))[0];
    burnin = INTEGER(dmm_settings_get("burnin", control))[0];
    sample_freq = INTEGER(dmm_settings_get("sample_freq", control))[0];
    alpha = REAL(dmm_settings_get("alpha", control))[0];
    beta = REAL(dmm_settings_get("beta", control))[0];

    lda = lda_init(alpha, beta, nsweep, ntop, ndoc, nterm, INTEGER(data));

    lda_run(burnin, nsweep, sample_freq, lda);


    result = PROTECT(allocVector(VECSXP, 3));

    SET_VECTOR_ELT(result, 0, allocVector(REALSXP, ntop));
    SET_VECTOR_ELT(result, 1, allocVector(REALSXP, ntop * nterm));
    SET_VECTOR_ELT(result, 2, allocVector(REALSXP, ndoc * ntop));

    memcpy(REAL(VECTOR_ELT(result, 0)), lda->theta, ntop * sizeof(double));
    memcpy(REAL(VECTOR_ELT(result, 1)), lda->phi, ntop * nterm * sizeof(double));
    memcpy(REAL(VECTOR_ELT(result, 2)), lda->docphi, ndoc * ntop * sizeof(double));

    lda_free(lda);

    PutRNGstate();
    UNPROTECT(1);

    return result;
}


