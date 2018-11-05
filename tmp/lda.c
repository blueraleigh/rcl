#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

/*
** Gibbs sampler for Latent Dirichlet Allocation mixture model
*/


struct lda {

    /* Number of documents (=species) */
    int ndoc;

    /* Number of terms (=resource categories) */
    int nterm;

    /* Sum of all term frequencies */
    int nword;

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

    /*
    ** Offset of each document in the data arrays
    **
    ** That is, docdat + docofs[i] returns a
    ** pointer to the first word in document i
    ** and likewise for doctop.
    */
    int *docofs;

    /*
    ** Document data. Each item is the term index
    ** associated with a given word in a given document.
    ** For example, the term index associated with the
    ** j-th word in the i-th document is
    **
    **      docdat[docofs[i] + j]
    */
    int *docdat;

    /*
    ** Topic assigment of each word in each document.
    ** For example, the topic index associated with the
    ** j-th word in the i-th document is
    **
    **      doctop[docofs[i] + j]
    */
    int *doctop;


    /*
    ** Each entry (i, j) is number of words in document
    ** i assigned to topic j. ndoc by ntop matrix. Each
    ** row sum equals corresponding value in dsz.
    */
    int *nmz;


    /*
    ** Each entry (i, j) is number of words with term index
    ** j that are assigned to topic i (across all documents).
    ** ntop by nterm matrix.
    */
    int *nzw;


    /*
    ** Each entry i is number of times topic
    ** i is assigned to words (across all documents).
    ** ntop by 1 matrix.
    */
    int *nz;


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
    ** ndoc by ntop matrix that records the average
    ** distribution over topics for each document.
    */
    double *theta;

    /*
    ** ntop by nterm matrix that records the average
    ** distribution over words for each topic.
    */
    double *phi;


    double *loglk;

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


/*
** Draw a new topic assignment for the j-th word
** in the i-th document
*/
static void draw(int i, int j, struct lda *lda)
{
    int k;
    int z;
    int w;
    double tot = 0;
    double prob[lda->ntop];

    // term index corresponding to the j-th word
    w = lda->docdat[lda->docofs[i] + j];

    // current topic assignment for this word
    z = lda->doctop[lda->docofs[i] + j];
    lda->nmz[i + z * lda->ndoc] -= 1;
    lda->nzw[z + w * lda->ntop] -= 1;
    lda->nz[z] -= 1;

    for (k = 0; k < lda->ntop; ++k) {
        prob[k] = (lda->alpha + lda->nmz[i + k * lda->ndoc])
                        * (lda->beta + lda->nzw[k + w * lda->ntop]);
        prob[k] /= (lda->nz[k] + lda->beta * lda->nterm);
        tot += prob[k];
    }
    for (k = 0; k < lda->ntop; ++k)
        prob[k] /= tot;


    z = choose_one(lda->ntop, prob);
    lda->nmz[i + z * lda->ndoc] += 1;
    lda->nzw[z + w * lda->ntop] += 1;
    lda->nz[z] += 1;

    // store new topic assignment for this word
    lda->doctop[lda->docofs[i] + j] = z;
}


static void compute_loglk(int iter, struct lda *lda)
{
    int i;
    int j;
    double tmp1;
    double tmp2;
    double loglk = 0;
/*
    int zz = 0;
    for (i = 0; i < lda->ndoc; ++i) {
        for (j = 0; j < lda->ntop; ++j)
            zz += lda->nmz[i + j * lda->ndoc];
    }
    Rprintf("%d   ", zz);
    zz = 0;
    for (i = 0; i < lda->ntop; ++i) {
        for (j = 0; j < lda->nterm; ++j)
            zz += lda->nzw[i + j * lda->ntop];
    }
    Rprintf("%d   ", zz);
    zz = 0;
    for (i = 0; i < lda->ntop; ++i) {
        zz += lda->nz[i];
    }
    Rprintf("%d\n", zz);
*/
    // P(w|z, T, alpha, beta)
    for (i = 0; i < lda->ntop; ++i) {
        tmp1 = tmp2 = 0;
        for (j = 0; j < lda->nterm; ++j) {
            tmp1 += lgamma(lda->beta + lda->nzw[i + j * lda->ntop]);
            tmp2 += lda->beta + lda->nzw[i + j * lda->ntop];
        }
        loglk += tmp1 - lgamma(tmp2);
    }

    loglk += lda->ntop * (lgamma(lda->nterm * lda->beta) - lda->nterm * lgamma(lda->beta));

    lda->loglk[iter + 0 * lda->nsweep] = loglk, loglk = 0;

    // P(z|T, alpha, beta)
    for (i = 0; i < lda->ndoc; ++i) {
        tmp1 = tmp2 = 0;
        for (j = 0; j < lda->ntop; ++j) {
            tmp1 += lgamma(lda->alpha + lda->nmz[i + j * lda->ndoc]);
            tmp2 += lda->alpha + lda->nmz[i + j * lda->ndoc];
        }
        loglk += tmp1 - lgamma(tmp2);
    }

    loglk += lda->ndoc * (lgamma(lda->ntop * lda->alpha) - lda->ntop * lgamma(lda->alpha));

    lda->loglk[iter + 1 * lda->nsweep] = loglk;
}


/* Perform one sweep of the LDA Gibb's sampler */
static void sweep(int iter, struct lda *lda)
{
    int i;
    int j;
    for (i = 0; i < lda->ndoc; ++i) {
        for (j = 0; j < lda->docsz[i]; ++j)
            draw(i, j, lda);
    }
    compute_loglk(iter, lda);
}


/*
** Update the partial estimate of the average topic
** distribution for the i-th document
*/
static void doctheta(int i, double n, struct lda *lda)
{
    int k;
    double denom = lda->ntop * lda->alpha + lda->docsz[i];
    for (k = 0; k < lda->ntop; ++k)
        lda->theta[i + k * lda->ndoc] += (lda->nmz[i + k * lda->ndoc] + lda->alpha) / (n*denom);
}


/*
** Update the partial estimate of the average distribution over terms
** for the i-th topic
*/
static void topicphi(int i, double n, struct lda *lda)
{
    int w;
    double denom = lda->nterm * lda->beta + lda->nz[i];
    for (w = 0; w < lda->nterm; ++w)
        lda->phi[i + w * lda->ntop] += (lda->nzw[i + w * lda->ntop] + lda->beta) / (n*denom);
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
    struct lda *lda;

    lda = Calloc(1, struct lda);
    lda->ntop = ntop;
    lda->ndoc = ndoc;
    lda->nterm = nterm;
    lda->nsweep = nsweep;
    lda->alpha = alpha;
    lda->beta = beta;
    lda->nmz = Calloc(ndoc * ntop, int);
    lda->nzw = Calloc(ntop * nterm, int);
    lda->nz = Calloc(ntop, int);
    lda->docsz = Calloc(ndoc, int);
    lda->docofs = Calloc(ndoc, int);
    lda->theta = Calloc(ndoc * ntop, double);
    lda->phi = Calloc(ntop * nterm, double);
    lda->loglk = Calloc(2 * nsweep, double);

    memset(lda->nmz, 0, ndoc * ntop * sizeof(int));
    memset(lda->nzw, 0, ntop * nterm * sizeof(int));
    memset(lda->nz, 0, ntop * sizeof(int));
    memset(lda->docsz, 0, ndoc * sizeof(int));
    memset(lda->theta, 0, ndoc * ntop * sizeof(double));
    memset(lda->phi, 0, ntop * nterm * sizeof(double));

    for (j = 0; j < ndoc; ++j) {
        lda->docofs[j] = nwrd;
        for (i = 0; i < nterm; ++i) {
            lda->docsz[j] += data[i + j * nterm];
            nwrd += data[i + j * nterm];
        }
    }

    lda->nword = nwrd;
    lda->docdat = Calloc(nwrd, int);
    lda->doctop = Calloc(nwrd, int);

    for (j = 0; j < ndoc; ++j) {
        w = 0;
        for (i = 0; i < nterm; ++i) {
            if (w < lda->docsz[j] && data[i + j * nterm] > 0) {
                nwrd = 0;
                while (nwrd < data[i + j * nterm]) {
                    z = (int)(unif_rand() * ntop);
                    lda->docdat[lda->docofs[j] + w] = i;
                    lda->doctop[lda->docofs[j] + w] = z;
                    lda->nmz[j + z * ndoc] += 1;
                    lda->nzw[z + i * ntop] += 1;
                    lda->nz[z] += 1;
                    nwrd += 1;
                    w += 1;
                }
            }
        }
    }

    return lda;
}


static void lda_free(struct lda *lda)
{
    Free(lda->docsz);
    Free(lda->docofs);
    Free(lda->docdat);
    Free(lda->doctop);
    Free(lda->nmz);
    Free(lda->nzw);
    Free(lda->nz);
    Free(lda->theta);
    Free(lda->phi);
    Free(lda->loglk);
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
            for (i = 0; i < lda->ndoc; ++i)
                doctheta(i, (double)nsample, lda);
            for (i = 0; i < lda->ntop; ++i)
                topicphi(i, (double)nsample, lda);
        }
    }
}


static SEXP lda_settings_get(const char *name, SEXP settings)
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


SEXP rcl_lda(SEXP control, SEXP data)
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
    ntop = INTEGER(lda_settings_get("ntopic", control))[0];
    nsweep = INTEGER(lda_settings_get("nsweep", control))[0];
    burnin = INTEGER(lda_settings_get("burnin", control))[0];
    sample_freq = INTEGER(lda_settings_get("sample_freq", control))[0];
    alpha = REAL(lda_settings_get("alpha", control))[0];
    beta = REAL(lda_settings_get("beta", control))[0];

    lda = lda_init(alpha, beta, nsweep, ntop, ndoc, nterm, INTEGER(data));

    lda_run(burnin, nsweep, sample_freq, lda);


    result = PROTECT(allocVector(VECSXP, 3));

    SET_VECTOR_ELT(result, 0, allocVector(REALSXP, 2 * nsweep));
    SET_VECTOR_ELT(result, 1, allocVector(REALSXP, ndoc * ntop));
    SET_VECTOR_ELT(result, 2, allocVector(REALSXP, ntop * nterm));

    memcpy(REAL(VECTOR_ELT(result, 0)), lda->loglk, 2 * nsweep * sizeof(double));
    memcpy(REAL(VECTOR_ELT(result, 1)), lda->theta, ndoc * ntop * sizeof(double));
    memcpy(REAL(VECTOR_ELT(result, 2)), lda->phi, ntop * nterm * sizeof(double));

    lda_free(lda);

    PutRNGstate();
    UNPROTECT(1);

    return result;
}


