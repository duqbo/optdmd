/*=========================================================
 * XGEQP3_M.c 
 * 
 * MatLab interface for LAPACK QR routine ZGEQP3 or DGEQP3
 *
 * [AREF,IPIV,TAU] = xgeqp3_m(A)
 *
 * Copyright Travis Askham 2017
 *
 * MIT License
 *
 * Permission is hereby granted, free of charge, to any person 
 * obtaining a copy of this software and associated documentation 
 * files (the "Software"), to deal in the Software without restriction, 
 * including without limitation the rights to use, copy, modify, merge, 
 * publish, distribute, sublicense, and/or sell copies of the Software, 
 * and to permit persons to whom the Software is furnished to do so, 
 * subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included 
 * in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS 
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY 
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, 
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE 
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *
 *=======================================================*/

#if !defined(_WIN32)
#define zgeqp3 zgeqp3_
#define dgeqp3 dgeqp3_
#endif

#include "mex.h"
#include "fort.h"
#include "stddef.h"

extern void zgeqp3( 
    ptrdiff_t *M, 
    ptrdiff_t *N, 
    double *A, 
    ptrdiff_t *LDA, 
    ptrdiff_t *JPVT, 
    double *TAU, 
    double *WORK, 
    ptrdiff_t *LWORK, 
    double *RWORK,
    ptrdiff_t *INFO );

extern void dgeqp3( 
    ptrdiff_t *M, 
    ptrdiff_t *N, 
    double *A, 
    ptrdiff_t *LDA, 
    ptrdiff_t *JPVT, 
    double *TAU, 
    double *WORK, 
    ptrdiff_t *LWORK, 
    ptrdiff_t *INFO );

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{

	/* mex interface to LAPACK functions dgeqp3 and zgeqp3 */

    char msg[101];
    int cplx;
    ptrdiff_t info, lda, m, n, minmn, lwork, i;
    ptrdiff_t *JPVT = NULL;
    double *A, *Ain, *TAU, *work, *work1, *rwork=NULL;
	
    /* check if appropriate # of inputs/outputs */
    
    if ((nlhs != 3) || (nrhs != 1)) {
      mexErrMsgIdAndTxt( "MATLAB:xgeqp3_m:invalidNumInputOutput",
              "Expect 1 input argument and return 3 output arguments");
    }
    
    /* check if appropriate input type */
    
    if (!mxIsDouble(prhs[0])){
        mexErrMsgIdAndTxt( "MATLAB:xgeqp3_m:invalidInputType",
              "Expect 1 array of double precision");
    }

    /* figure out array dimensions */
    
    n = mxGetN(prhs[0]);
    lda = mxGetM(prhs[0]);
    m = lda;
    
    minmn = n;
    if (m < n){
        minmn = m;
    }

    /* allocate memory */
    
    cplx = (mxGetPi(prhs[0]));
    if (cplx) {
      A = mat2fort(prhs[0],lda,n);
      TAU = (double *)mxCalloc(2*minmn,sizeof(double));
    } else {
      Ain = mxGetPr(prhs[0]);
      A = (double *)mxCalloc(m*n,sizeof(double));
      for (i=0; i<m*n; i++){
          A[i] = Ain[i];
      }
      TAU = (double *)mxCalloc(minmn,sizeof(double));      
    }
    
    JPVT = (ptrdiff_t *)mxCalloc(n,sizeof(ptrdiff_t));
    
    for (i=0; i<n; i++){
        JPVT[i] = 0;
    }
    
    lwork = -1;
    info = 0;
    
    if (cplx) {
      work1 = (double *)mxCalloc(2,sizeof(double));
      rwork = (double *)mxCalloc(2*n,sizeof(double));
      /* Query zgeqp3 on the value of lwork */
      zgeqp3(&m, &n, A, &lda, JPVT, TAU, work1, &lwork, rwork, &info);
        if (info < 0) {
          sprintf(msg, "Input %d to zgeqp3 had an illegal value",-info);
          mexErrMsgIdAndTxt( "MATLAB:xgeqp3_m:illegalInputTozgeqp3", msg);
        }
      lwork = (ptrdiff_t)(work1[0]);
      work = (double *)mxCalloc(2*lwork,sizeof(double));
      zgeqp3(&m, &n, A, &lda, JPVT, TAU, work, &lwork, rwork, &info);      
      if (info < 0) {
          sprintf(msg, "Input %d to zgeqp3 had an illegal value",-info);
          mexErrMsgIdAndTxt( "MATLAB:xgeqp3_m:illegalInputTozgeqp3", msg);
        }   
    } else {
        work1 = (double *)mxCalloc(1,sizeof(double));
      /* Query dgeqp3 on the value of lwork */
        dgeqp3(&m, &n, A, &lda, JPVT, TAU, work1, &lwork, &info);
        if (info < 0) {
          sprintf(msg, "Input %d to dgeqp3 (1st call) had an illegal value",-info);
          mexErrMsgIdAndTxt( "MATLAB:xgeqp3_m:illegalInputTodgeqp3_1", msg);
      }
      lwork = (ptrdiff_t)(work1[0]);
      work = (double *)mxCalloc(lwork,sizeof(double));
      dgeqp3(&m, &n, A, &lda, JPVT, TAU, work, &lwork, &info);
      if (info < 0) {
          sprintf(msg, "Input %d to dgeqp3 (2nd call) had an illegal value",-info);
          mexErrMsgIdAndTxt( "MATLAB:xgeqp3_m:illegalInputTodgeqp3_2", msg);
      }
    }

    if (cplx) {
      plhs[0] = fort2mat(A,lda,m,n);
      plhs[1] = mxCreateNumericMatrix(0,0,
               ((sizeof(ptrdiff_t) == 8) ? mxINT64_CLASS : mxINT32_CLASS),
                       mxREAL);
      mxSetM(plhs[1],n);
      mxSetN(plhs[1],1);
      mxSetData(plhs[1],JPVT);
      plhs[2] = fort2mat(TAU,minmn,minmn,1);
      mxFree(rwork);
    } else {
      plhs[0] = mxCreateDoubleMatrix(0,0,0);
      mxSetPr(plhs[0],A);
      mxSetM(plhs[0],m);
      mxSetN(plhs[0],n);
      plhs[1] = mxCreateNumericMatrix(0,0,
               ((sizeof(ptrdiff_t) == 8) ? mxINT64_CLASS : mxINT32_CLASS),
                       mxREAL);
      mxSetM(plhs[1],n);
      mxSetN(plhs[1],1);
      mxSetData(plhs[1],JPVT);
      plhs[2] = mxCreateDoubleMatrix(0,0,0);      
      mxSetPr(plhs[2],TAU);
      mxSetM(plhs[2],minmn);
      mxSetN(plhs[2],1);
    }
     
    mxFree(work1);
    mxFree(work);

}
