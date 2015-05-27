#include <Python.h>
#include "numpy/arrayobject.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#ifdef _OPENMP
#include <omp.h>
#endif

/* This program mainly comes from http://bioinfadmin.cs.ucl.ac.uk/downloads/PSICOV/psicov21.c
I correct a minor mistake and make the memory more safety. */

#define FALSE 0
#define TRUE 1

#define SQR(x) ((x)*(x))
#define MAX(x,y) ((x)>(y)?(x):(y))
#define MIN(x,y) ((x)<(y)?(x):(y))

#define MAXSEQLEN 5000
#define MINEFSEQS (seqlen)

/* Convert AA letter to numeric code (0-21) */
int aanum(int ch)
{
  const static int aacvs[] ={999, 0, 3, 4, 3, 6, 13, 7, 8, 9, 21, 11, 10,
   12, 2, 21, 14, 5, 1, 15, 16, 21, 19, 17, 21, 18, 6};

  return (isalpha(ch) ? aacvs[ch & 31] : 20);
}

/* Allocate matrix */
void* allocmat(int rows, int columns, int size)
{
  int   i,j;
  void  **p, *rp;

  rp = malloc(rows * sizeof(void *) + sizeof(int));

  if (rp == NULL)
    return NULL;

  *((int *)rp) = rows;

  p = rp + sizeof(int);

  for (i = 0; i < rows; i++)
    if ((p[i] = calloc(columns, size)) == NULL){
      for (j=0;j<i;j++)
        free(p[j]);
      free((void *)(p-sizeof(int)));
      return NULL;
    }

  return p;
}

/* Allocate vector */
void* allocvec(int columns, int size)
{
  void *p;

  p = calloc(columns, size);

  if (p == NULL)
    return NULL;

  return p;
}


/* LASSO*/

#define EPS (1.1e-15)
#define BIG (1e9)

int glassofast(const int n, double **S, double **L, const double thr, const int maxit, int approxflg, int warm, double **X, double **W)
{
/*  This subroutine computes the L1 regularized covariance matrix estimate
  using the algorithm described in the paper:
  J. Friedman, T. Hastie, R. Tibshirani:
  Sparse inverse covariance estimation with the graphical lasso
  Biostatistics, 9(3):432-441, July 2008.
  This code is adapted from the Fortran code described in the following report:
  M. A. Sustik & B. Calderhead:
  GLASSOFAST: An efficient GLASSO implementation
  Technical Report TR-12-29, University of Texas at Austin

  NOTE: that when multiple threads are used, we gain a huge time saving by
  avoiding full thread synchronisation when updating elements of the W (covariance)
  matrix. In multithreaded mode, the order of updates to the W matrix at each iteration
  will depend on the order in which threads complete. In practice, this hardly matters,
  because the algorithm is iterative, and in testing still converges to within 6 d.p.
  of the non-threaded code. If a very small degree of non-deterministic behaviour really
  worries you, then set the maximum number of threads to 1 (or compile without OpenMP).
*/
  int i, j, ii, iter, jj;
  double a, b, c, delta, dlx, dw, shr, sum, thrlasso, tmp, wd[MAXSEQLEN*21], wxj[MAXSEQLEN*21];

  for (shr=ii=0; ii<n; ii++)
    for (jj=0; jj<n; jj++)
      shr += fabs(S[ii][jj]);
  
  for (i=0; i<n; i++)
    shr -= fabs(S[i][i]);
  
  if (shr == 0.0)
  {
    /* S is diagonal. */
    for (ii=0; ii<n; ii++)
      for (jj=0; jj<n; jj++)
        W[ii][jj] = X[ii][jj] = 0.0;
    
    for (i=0; i<n; i++)
      W[i][i] = W[i][i] + L[i][i];
    
    for (ii=0; ii<n; ii++)
      for (jj=0; jj<n; jj++)
        X[ii][jj] = 0.0;
    
    for (i=0; i<n; i++)
      X[i][i] = 1.0 / MAX(W[i][i], EPS);

    return 0;
  }
  
  shr *= thr/(n-1);
  thrlasso = shr/n;
  if (thrlasso < 2*EPS)
    thrlasso = 2*EPS;
  
  if (!warm)
  {
    for (ii=0; ii<n; ii++)
      for (jj=0; jj<n; jj++)
      {
        W[ii][jj] = S[ii][jj];
        X[ii][jj] = 0.0;
      }
  }
  else
  {
    for (i=0; i<n; i++)
    {
      for (ii=0; ii<n; ii++)
        X[i][ii] = -X[i][ii]/X[i][i];
      X[i][i] = 0.0;
    }
  }
  
  for (i=0; i<n; i++)
  {
    wd[i] = S[i][i] + L[i][i];
    W[i][i] = wd[i];
  }
  
  for (iter = 1; iter<=maxit; iter++)
  {
    dw = 0.0;

    #pragma omp parallel for default(shared) private(i,j,ii,wxj,a,b,c,dlx,delta,sum)
    for (j=0; j<n; j++)
    {
      for (ii=0; ii<n; ii++)
        wxj[ii] = 0.0;

      for (i=0; i<n; i++)
        if (X[j][i] != 0.0)
          for (ii=0; ii<n; ii++)
            wxj[ii] += W[i][ii] * X[j][i];

      for (;;)
      {
        dlx = 0.0;
    
        for (i=0; i<n; i++)
        {
          if (i != j && L[j][i] < BIG)
          {
            a = S[j][i] - wxj[i] + wd[i] * X[j][i];
            b = fabs(a) - L[j][i];
            if (b <= 0.0)
              c = 0.0;
            else if (a >= 0.0)
              c = b / wd[i];
            else
              c = -b / wd[i];

            delta = c - X[j][i];
            if (delta != 0.0 && (!approxflg || fabs(delta) > 1e-6))
            {
              X[j][i] = c;
            
              for (ii=0; ii<n; ii++)
                wxj[ii] += W[i][ii] * delta;
              
              if (fabs(delta) > dlx)
                dlx = fabs(delta);
            }
          }
        }
      
        if (dlx < thrlasso)
          break;
      }
          
      wxj[j] = wd[j];
      
      for (sum=ii=0; ii<n; ii++)
        sum += fabs(wxj[ii] - W[j][ii]);

      #pragma omp critical
        if (sum > dw)
          dw = sum;

      for (ii=0; ii<n; ii++)
        W[j][ii] = wxj[ii];
      for (ii=0; ii<n; ii++)
        W[ii][j] = wxj[ii];
    }
  
    if (dw <= shr)
      break;
  }

  for (i=0; i<n; i++)
  {
    for (sum=ii=0; ii<n; ii++)
      sum += X[i][ii] * W[i][ii];
    
    tmp = 1.0 / (wd[i] - sum);
    
    for (ii=0; ii<n; ii++)
      X[i][ii] = -tmp * X[i][ii];
    X[i][i] = tmp;
  }
  
  for (i=0; i<n-1; i++)
  {
    for (ii=i+1; ii<n; ii++)
    {
      X[i][ii] = 0.5 * (X[i][ii] + X[ii][i]);
      X[ii][i] = X[i][ii];
    }
  }
  
  return iter;
}

/* Perform Cholesky decomposition on matrix */
int test_cholesky(double **a, const int n) 
{
  int i, j, k, status=0;
  double sum;
  static double *diag;

  if (diag == NULL)
    diag = allocvec(n, sizeof(double));
  if (diag==NULL)
    return 2;

  for (i=0; i<n; i++)
  {
    if (!status)
      for (j=i; j<n; j++)
      {
        sum = a[i][j];
    
        for (k=i-1; k >= 0; k--)
          sum -= a[i][k]*a[j][k];
    
        if (i == j)
        {
          if (sum <= 0.0)
          status = 1;
          
          diag[i] = sqrt(sum);
        }
        else
          a[j][i] = sum / diag[i];
      }
  }

  return status;
}

struct sc_entry
{
  double sc;
  int i, j;
} *sclist;

static PyObject *msapsicov(PyObject *self, PyObject *args, PyObject *kwargs)
{
  PyArrayObject *msa, *psicov;

  int approxflg, shrinkflg, overrideflg, rawscflg, apcflg,  
      minseqsep, maxthread, pseudoc;
  double rhodefault, targfnzero, thresh, idthresh, maxgapf;
  char *blockfn = NULL;

  int *wtcount;
  double *weight, **pa, **pcmat, *pcsum;
  char **aln;

  int a, b, i, j, k, ndim, maxit=10000, initflg=0, npair, nnzero, ncon;
  double wtsum, smean, lambda, lastfnzero, trialrho, rfact, score, fnzero, pcmean,
         pc, scsum, scsumsq, mean, sd, zscore, ppv;
  FILE *ifp;

  static char *kwlist[] = {"msa", "psicov", "approxflg", "shrinkflg", 
                           "overrideflg", "rawscflg", "apcflg", 
                           "rhodefault", "targfnzero", "thresh", "idthresh", 
                           "pseudoc", "minseqsep", "blockfn", "maxgapf", 
                           "maxthread", NULL};

  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OO|iiiiiddddiisdi", kwlist,
                                   &msa, &psicov,
                                   &approxflg, &shrinkflg, &overrideflg, 
                                   &rawscflg, &apcflg,
                                   &rhodefault,&targfnzero,&thresh,&idthresh,
                                   &pseudoc,&minseqsep,
                                   &blockfn,
                                   &maxgapf,
                                   &maxthread
                                   ))
      return NULL;

  /* make sure to have a contiguous and well-behaved array */
  msa = PyArray_GETCONTIGUOUS(msa);

  /* check dimensions */
  long nseqs = msa->dimensions[0], seqlen = msa->dimensions[1];

  /* get pointers to data */
  char *seq = (char *) PyArray_DATA(msa); /*size: number x length */
  double *mut = (double *) PyArray_DATA(psicov);

  if (maxthread!=-1)
    #ifdef _OPENMP
      omp_set_num_threads(maxthread);
    #endif
  
  aln = allocvec(nseqs, sizeof(char *));
  if (aln==NULL){
    return Py_BuildValue("Oi", Py_None, 1);
  }

  weight = allocvec(nseqs, sizeof(double));
  if (weight==NULL){
    free(aln);
    return Py_BuildValue("Oi", Py_None, 1);
  }

  wtcount = allocvec(nseqs, sizeof(int));
  if (wtcount==NULL){
    free(aln);
    free(weight);
    return Py_BuildValue("Oi", Py_None, 1);
  }

  if (!(aln[0] = malloc(seqlen))){
    free(aln);
    free(weight);
    free(wtcount);
    return Py_BuildValue("Oi", Py_None, 1);
  }

  for (j=0; j<seqlen; j++)
    aln[0][j] = aanum(seq[j]);

  for (i=1; i<nseqs; i++)
  {
    if (!(aln[i] = malloc(seqlen))){
      for (j=0;j<i;j++)
        free(aln[j]);
      free(aln);
      free(weight);
      free(wtcount);
      return Py_BuildValue("Oi", Py_None, 1);
    }
    
    for (j=0; j<seqlen; j++)
      aln[i][j] = aanum(seq[i*seqlen+j]);
  }

  /* Calculate sequence weights (use openMP/pthreads if available) */
  if (idthresh < 0.0)
  {
    double meanfracid = 0.0;
  
    #pragma omp parallel for default(shared) private(j,k) reduction(+:meanfracid)
      for (i=0; i<nseqs; i++)
        for (j=i+1; j<nseqs; j++)
        {
          int nids;
          double fracid;

          for (nids=k=0; k<seqlen; k++)
            nids += (aln[i][k] == aln[j][k]);
          
          fracid = (double)nids / seqlen;
          
          meanfracid += fracid;
        }
    
    meanfracid /= 0.5 * nseqs * (nseqs - 1.0);

    idthresh = MIN(0.6, 0.38 * 0.32 / meanfracid);
  }

  #pragma omp parallel for default(shared) private(j,k)
    for (i=0; i<nseqs; i++)
      for (j=i+1; j<nseqs; j++)
      {
        int nthresh = seqlen * idthresh;

        for (k=0; nthresh > 0 && k<seqlen; k++)
          nthresh -= (aln[i][k] != aln[j][k]);
        
        if (nthresh > 0)
        {
          #pragma omp critical
          {
            wtcount[i]++;
            wtcount[j]++;
          }
        }
      }

  for (wtsum=i=0; i<nseqs; i++)
    wtsum += (weight[i] = 1.0 / (1 + wtcount[i]));

  if (wtsum < seqlen && !overrideflg){
      for (j=0;j<nseqs;j++)
        free(aln[j]);
      free(aln);
      free(weight);
      free(wtcount);
      return Py_BuildValue("Oid", Py_None, 2, wtsum);
  }

  pa = allocmat(seqlen, 21, sizeof(double));
  if (pa==NULL){
      for (j=0;j<nseqs;j++)
        free(aln[j]);
      free(aln);
      free(weight);
      free(wtcount);
      return Py_BuildValue("Oi", Py_None, 1);
  }

  /* Calculate singlet frequencies with pseudocount */
  for (i=0; i<seqlen; i++)
  {
    for (a=0; a<21; a++)
      pa[i][a] = pseudoc;
    
    for (k=0; k<nseqs; k++)
    {
      a = aln[k][i];
      if (a < 21)
      pa[i][a] += weight[k];
    }
    
    for (a=0; a<21; a++)
      pa[i][a] /= pseudoc * 21.0 + wtsum;
  }

  double **cmat, **rho, **ww, **wwi, **tempmat;

  ndim = seqlen * 21;

  cmat = allocmat(ndim, ndim, sizeof(double));
  if (cmat==NULL){
      for (j=0;j<nseqs;j++)
        free(aln[j]);
      free(aln);
      free(weight);
      free(wtcount);
      for (j=0;j<seqlen;j++)
        free(pa[j]);
      free((void *)pa - sizeof(int));
      return Py_BuildValue("Oi", Py_None, 1);
  }

  tempmat = allocmat(ndim, ndim, sizeof(double));
  if (tempmat==NULL){
      for (j=0;j<nseqs;j++)
        free(aln[j]);
      free(aln);
      free(weight);
      free(wtcount);
      for (j=0;j<seqlen;j++)
        free(pa[j]);
      free((void *)pa - sizeof(int));
      for (j=0;j<ndim;j++)
        free(cmat[j]);
      free((void *)cmat - sizeof(int));
      return Py_BuildValue("Oi", Py_None, 1);
  }

  /* Form the covariance matrix */
  #pragma omp parallel for default(shared) private(j,k,a,b)
    for (i=0; i<seqlen; i++)
      for (j=i; j<seqlen; j++)
      {
        double pab[21][21];

        for (a=0; a<21; a++)
          for (b=0; b<21; b++)
            if (i == j)
              pab[a][b] = (a == b) ? pa[i][a] : 0.0;
            else
              pab[a][b] = pseudoc / 21.0;
        
        if (i != j)
        {
          for (k=0; k<nseqs; k++)
          {
            a = aln[k][i];
            b = aln[k][j];
            if (a < 21 && b < 21)
              pab[a][b] += weight[k];
          }
          
          for (a=0; a<21; a++)
            for (b=0; b<21; b++)
              pab[a][b] /= pseudoc * 21.0 + wtsum;
        }
        
        for (a=0; a<21; a++)
          for (b=0; b<21; b++)
            if (i != j || a == b)
              cmat[i*21+a][j*21+b] = cmat[j*21+b][i*21+a] = pab[a][b] - pa[i][a] * pa[j][b];
      }

  /* Shrink sample covariance matrix towards shrinkage target F = Diag(1,1,1,...,1) * smean */
  int checkflag=0;
  if (shrinkflg)
  {
    for (smean=i=0; i<ndim; i++)
      smean += cmat[i][i];
    
    smean /= (double)ndim;
    lambda = 0.2;

    for (;;)
    {
      for (i=0; i<ndim; i++)
        memcpy(tempmat[i], cmat[i], ndim*sizeof(double));
      
      /* Test if positive definite using Cholesky decomposition */

      int testc=test_cholesky(tempmat, ndim);
      if (testc==2){
        checkflag=1;
        break;
      }
      if (!testc)
        break;
      
      #pragma omp parallel for default(shared) private(j,a,b)
        for (i=0; i<seqlen; i++)
        for (j=0; j<seqlen; j++)
          for (a=0; a<21; a++)
          for (b=0; b<21; b++)
            if (i != j)
              cmat[i*21+a][j*21+b] *= 1.0 - lambda;
            else if (a == b)
              cmat[i*21+a][j*21+b] = smean * lambda + (1.0 - lambda) * cmat[i*21+a][j*21+b];
    }
  }
  if (checkflag){
    for (j=0;j<nseqs;j++)
      free(aln[j]);
    free(aln);
    free(weight);
    free(wtcount);
    for (j=0;j<seqlen;j++)
      free(pa[j]);
    free((void *)pa - sizeof(int));
    for (j=0;j<ndim;j++)
      free(cmat[j]);
    free((void *)cmat - sizeof(int));
    for(j=0;j<ndim;j++)
      free(tempmat[j]);
    free((void *)tempmat - sizeof(int));
    return Py_BuildValue("Oi", Py_None, 1);
  }

  rho = allocmat(ndim, ndim, sizeof(double));
  if (rho==NULL){
    for (j=0;j<nseqs;j++)
      free(aln[j]);
    free(aln);
    free(weight);
    free(wtcount);
    for (j=0;j<seqlen;j++)
      free(pa[j]);
    free((void *)pa - sizeof(int));
    for (j=0;j<ndim;j++)
      free(cmat[j]);
    free((void *)cmat - sizeof(int));
    for(j=0;j<ndim;j++)
      free(tempmat[j]);
    free((void *)tempmat - sizeof(int));
    return Py_BuildValue("Oi", Py_None, 1);
  }
  ww = allocmat(ndim, ndim, sizeof(double));
  if (ww==NULL){
    for (j=0;j<nseqs;j++)
      free(aln[j]);
    free(aln);
    free(weight);
    free(wtcount);
    for (j=0;j<seqlen;j++)
      free(pa[j]);
    free((void *)pa - sizeof(int));
    for (j=0;j<ndim;j++)
      free(cmat[j]);
    free((void *)cmat - sizeof(int));
    for(j=0;j<ndim;j++)
      free(tempmat[j]);
    free((void *)tempmat - sizeof(int));
    for(j=0;j<ndim;j++)
      free(rho[j]);
    free((void *)rho - sizeof(int));
    return Py_BuildValue("Oi", Py_None, 1);
  }
  wwi = allocmat(ndim, ndim, sizeof(double));
  if (wwi==NULL){
    for (j=0;j<nseqs;j++)
      free(aln[j]);
    free(aln);
    free(weight);
    free(wtcount);
    for (j=0;j<seqlen;j++)
      free(pa[j]);
    free((void *)pa - sizeof(int));
    for (j=0;j<ndim;j++)
      free(cmat[j]);
    free((void *)cmat - sizeof(int));
    for(j=0;j<ndim;j++)
      free(tempmat[j]);
    free((void *)tempmat - sizeof(int));
    for(j=0;j<ndim;j++)
      free(rho[j]);
    free((void *)rho - sizeof(int));
    for(j=0;j<ndim;j++)
      free(ww[j]);
    free((void *)ww - sizeof(int));
    return Py_BuildValue("Oi", Py_None, 1);
  }

  lastfnzero=0.0;

  /* Guess at a reasonable starting rho value if undefined */
  if (rhodefault < 0.0)
    trialrho = MAX(0.001, 1.0 / wtsum);
  else
    trialrho = rhodefault;

  rfact = 0.0;

  if (blockfn[0]=='\0')
    blockfn=NULL;

  double besttd = BIG, bestrho=trialrho;
  for (;;)
  {
    double targdiff;
    
    if (trialrho <= 0.0 || trialrho >= 1.0)
    {
      /* Give up search - recalculate with best rho found so far and exit */
      trialrho = bestrho;
      targfnzero = 0.0;
    }
      
    for (i=0; i<ndim; i++)
      for (j=0; j<ndim; j++)
        rho[i][j] = trialrho;
    
    for (i=0; i<seqlen; i++)
    for (j=0; j<seqlen; j++)
      for (a=0; a<21; a++)
      for (b=0; b<21; b++)
        if ((a != b && i == j) || pa[i][20] > maxgapf || pa[j][20] > maxgapf)
          rho[i*21+a][j*21+b] = BIG;
    
    /* Mask out regions if block-out list provided */
    if (blockfn != NULL)
    {
      ifp = fopen(blockfn, "r");
      
      for (;;)
      {
      if (fscanf(ifp, "%d %d %lf", &i, &j, &score) != 3)
        break;
      
      for (a=0; a<21; a++)
        for (b=0; b<21; b++)
        {
          rho[(i-1)*21+a][(j-1)*21+b] = score;
          rho[(j-1)*21+b][(i-1)*21+a] = score;
        }
      }
      
      fclose(ifp);
    }

    glassofast(ndim, cmat, rho, thresh, maxit, approxflg, initflg, wwi, ww);

    /* Don't attempt interation if too few sequences */
    if (targfnzero <= 0.0 || wtsum < seqlen)
      break;
    
    for (npair=nnzero=i=0; i<ndim; i++)
      for (j=i+1; j<ndim; j++,npair++)
        if (wwi[i][j] != 0.0)
          nnzero++;

    fnzero = (double) nnzero / npair;

    /* Stop iterating if we have achieved the target sparsity level */
    targdiff = fabs(fnzero - targfnzero)/targfnzero;

    if (targdiff < 0.01)
      break;

    if (targdiff < besttd)
    {
      besttd = targdiff;
      bestrho = trialrho;
    }
    
    if (fnzero == 0.0)
    {
      /* As we have guessed far too high, halve rho and try again */
      trialrho *= 0.5;
      continue;
    }
    
    if (lastfnzero > 0.0 && fnzero != lastfnzero)
        rfact = pow(rfact, log(targfnzero / fnzero) / log(fnzero / lastfnzero));

    lastfnzero = fnzero;

    /* Make a small trial step in the appropriate direction */
    if (rfact == 0.0)
      rfact = (fnzero < targfnzero) ? 0.9 : 1.1;
    
    trialrho *= rfact;
  }

  /* Calculate background corrected scores using average product correction */
  pcmat = allocmat(seqlen, seqlen, sizeof(double));
  if (pcmat==NULL){
    for (j=0;j<nseqs;j++)
      free(aln[j]);
    free(aln);
    free(weight);
    free(wtcount);
    for (j=0;j<seqlen;j++)
      free(pa[j]);
    free((void *)pa - sizeof(int));
    for (j=0;j<ndim;j++)
      free(cmat[j]);
    free((void *)cmat - sizeof(int));
    for(j=0;j<ndim;j++)
      free(tempmat[j]);
    free((void *)tempmat - sizeof(int));
    for(j=0;j<ndim;j++)
      free(rho[j]);
    free((void *)rho - sizeof(int));
    for(j=0;j<ndim;j++)
      free(ww[j]);
    free((void *)ww - sizeof(int));
    for(j=0;j<ndim;j++)
      free(wwi[j]);
    free((void *)wwi - sizeof(int));
    return Py_BuildValue("Oi", Py_None, 1);
  }
  pcsum = allocvec(seqlen, sizeof(double));
  if (pcsum==NULL){
    for (j=0;j<nseqs;j++)
      free(aln[j]);
    free(aln);
    free(weight);
    free(wtcount);
    for (j=0;j<seqlen;j++)
      free(pa[j]);
    free((void *)pa - sizeof(int));
    for (j=0;j<ndim;j++)
      free(cmat[j]);
    free((void *)cmat - sizeof(int));
    for(j=0;j<ndim;j++)
      free(tempmat[j]);
    free((void *)tempmat - sizeof(int));
    for(j=0;j<ndim;j++)
      free(rho[j]);
    free((void *)rho - sizeof(int));
    for(j=0;j<ndim;j++)
      free(ww[j]);
    free((void *)ww - sizeof(int));
    for(j=0;j<ndim;j++)
      free(wwi[j]);
    free((void *)wwi - sizeof(int));
    for(j=0;j<seqlen;j++)
      free(pcmat[j]);
    free((void *)pcmat - sizeof(int));
    return Py_BuildValue("Oi", Py_None, 1);
  }
  
  pcmean = 0.0;
  
  for (i=0; i<seqlen; i++)
    for (j=i+1; j<seqlen; j++)
    { 
      for (pc=a=0; a<20; a++)
      for (b=0; b<20; b++)
        pc += fabs(wwi[i*21+a][j*21+b]);

      pcmat[i][j] = pcmat[j][i] = pc;
      pcsum[i] += pc;
      pcsum[j] += pc;

      pcmean += pc;
    }

  pcmean /= seqlen * (seqlen - 1) * 0.5;

  /* Build final list of predicted contacts */
  for (scsum=scsumsq=ncon=i=0; i<seqlen; i++)
  {
    for (j=i; j<seqlen; j++)
      if (pcmat[i][j] > 0.0)
      {
        /* Calculate APC score */
        if (apcflg)
          pcmat[i][j] = pcmat[j][i] = pcmat[i][j] - pcsum[i] * pcsum[j] / SQR(seqlen - 1.0) / pcmean;
        if (j>=i+minseqsep){
          scsum += pcmat[i][j];
          scsumsq += SQR(pcmat[i][j]);
          ncon++;
        }
      }
      else{
        pcmat[i][j]=pcmat[j][i]=0.;
      }
  }

  mean = scsum / ncon;
  sd = 1.25 * sqrt(scsumsq / ncon - SQR(mean)); /* Corrected for extreme-value bias */

  if (!rawscflg)
    for (i=0; i<seqlen; i++)
      for (j=i+minseqsep; j<seqlen; j++)
      {
        zscore = (pcmat[i][j] - mean) / sd;
        ppv = 0.904 / (1.0 + 16.61 * exp(-0.8105 * zscore));
        mut[i*seqlen+j]=mut[j*seqlen+i]=ppv;
      }
  else
    for (i=0; i<seqlen; i++)
      for (j=i+minseqsep; j<seqlen; j++)
          mut[i*seqlen+j]=mut[j*seqlen+i]=pcmat[i][j];

  for (j=0;j<nseqs;j++)
    free(aln[j]);
  free(aln);
  free(weight);
  free(wtcount);
  for (j=0;j<seqlen;j++)
    free(pa[j]);
  free((void *)pa - sizeof(int));
  for (j=0;j<ndim;j++)
    free(cmat[j]);
  free((void *)cmat - sizeof(int));
  for(j=0;j<ndim;j++)
    free(tempmat[j]);
  free((void *)tempmat - sizeof(int));
  for(j=0;j<ndim;j++)
    free(rho[j]);
  free((void *)rho - sizeof(int));
  for(j=0;j<ndim;j++)
    free(ww[j]);
  free((void *)ww - sizeof(int));
  for(j=0;j<ndim;j++)
    free(wwi[j]);
  free((void *)wwi - sizeof(int));
  for(j=0;j<seqlen;j++)
    free(pcmat[j]);
  free((void *)pcmat - sizeof(int));
  free(pcsum);

  return Py_BuildValue("O", psicov);
}

static PyMethodDef Ccorrelation_p_methods[] =
{
   {"msapsicov", (PyCFunction)msapsicov, 
   METH_VARARGS | METH_KEYWORDS, 
   "Return PSICOV matrix calculated for given character \n"
   "array that contains an MSA."},

   {NULL, NULL, 0, NULL}
};



#if PY_MAJOR_VERSION >= 3

static struct PyModuleDef Ccorrelation_pmodule = {
    PyModuleDef_HEAD_INIT,
    "Ccorrelation_p",
    "MSA correlation tools with parallel.",
    -1,
    Ccorrelation_p_methods,
};
PyMODINIT_FUNC PyInit_Ccorrelation_p(void) {
  import_array();
  return PyModule_Create(&Ccorrelation_pmodule);
}
#else
PyMODINIT_FUNC initCcorrelation_p(void) {

  Py_InitModule3("Ccorrelation_p", Ccorrelation_p_methods,
    "MSA correlation tools with parallel.");

  import_array();
}
#endif
