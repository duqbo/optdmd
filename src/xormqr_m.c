/*=========================================================
 * XORMQR_M.c 
 * 
 * MatLab interface for LAPACK QR multiply routine ZUNMQR or DORMQR
 *
 * C = xormqr_m(SIDE,TRANS,AREF,tau,B,varargin)
 *
 * AREF and tau are as returned by xgeqp3_m (should be of same type,
 * either real or complex double) These store the necessary information
 * for applying the Householder reflectors.
 *
 * SIDE = 'L' or 'R', multiply on left or right
 * TRANS = 'N' or 'T' or 'C', multiply by Q, Q transpose, or Q conjugate
 * transpose
 * varargin{1} = k = number of reflectors defining Q (typically k = 
 *  min(m,n) where (m,n) are dimensions of AREF)
 *
 * C = Q*B or C = B*Q or C = Q'*B or C = B*Q' (determined by SIDE and
 * TRANS)
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
#define zunmqr zunmqr_
#define dormqr dormqr_
#endif

#include "mex.h"
#include "fort.h"
#include "stddef.h"

extern void zunmqr( 
    char * SIDE,
    char * TRANS,
    ptrdiff_t *M, 
    ptrdiff_t *N, 
    ptrdiff_t *K,         
    double *A, 
    ptrdiff_t *LDA, 
    double *TAU, 
    double *C, 
    ptrdiff_t *LDC, 
    double *WORK,
    ptrdiff_t *LWORK,         
    ptrdiff_t *INFO );

extern void dormqr( 
    char * SIDE,
    char * TRANS,
    ptrdiff_t *M, 
    ptrdiff_t *N, 
    ptrdiff_t *K,         
    double *A, 
    ptrdiff_t *LDA, 
    double *TAU, 
    double *C, 
    ptrdiff_t *LDC, 
    double *WORK,
    ptrdiff_t *LWORK,         
    ptrdiff_t *INFO );

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{

	/* mex interface to LAPACK functions zunmqr and dormqr */

    char msg[101];
    char * cp;
    char SIDE[1], TRANS[1];
    int cplx;
    int8_T * k8;
    int16_T * k16;
    int32_T * k32;
    int64_T * k64;
    uint8_T * ku8;
    uint16_T * ku16;
    uint32_T * ku32;
    uint64_T * ku64;
    double * kd;
    ptrdiff_t info, lda, ldc, ma, na, m, n, minmn, lwork, i, k, ltau;
    double *A, *B, *Bin, *TAU, *work, *work1=NULL;
	
    /* check if appropriate # of inputs/outputs */
    
    if ((nlhs != 1) || (nrhs < 5)) {
      mexErrMsgIdAndTxt( "MATLAB:xormqr_m:invalidNumInputOutput",
              "Expect at least 5 input arguments and return 1 output argument");
    }
    
    /* check if appropriate input type */
    
    if ((!mxIsDouble(prhs[2])) || (!mxIsDouble(prhs[3])) || (!mxIsDouble(prhs[4])) ){
        mexErrMsgIdAndTxt( "MATLAB:xormqr_m:invalidInputType",
              "Expect inputs 3,4,5 to be arrays of double precision");
    }
    
    if ((!mxIsChar(prhs[0])) || (!mxIsChar(prhs[1]))){
        mexErrMsgIdAndTxt( "MATLAB:xormqr_m:invalidInputType",
              "Expect inputs 1,2 to be character type");
    }
    
    /* figure out array dimensions */
    
    n = mxGetN(prhs[4]);
    m = mxGetM(prhs[4]);
    ma = mxGetM(prhs[2]);
    na = mxGetN(prhs[2]);
    ltau = mxGetM(prhs[3]);
    ldc = m;
    lda = ma;
    
    
    /* get strings */
    
    cp = (char *)mxGetData(prhs[0]);
    SIDE[0] = cp[0];
    cp = (char *)mxGetData(prhs[1]);
    TRANS[0] = cp[0];

    /* check strings */
    
    if (! (SIDE[0] == 'L' || SIDE[0] == 'R')){
        mexErrMsgIdAndTxt( "MATLAB:xormqr_m:invalidInput",
              "Expect input 1 to be 'L' or 'R'");
    }
    if (! (TRANS[0] == 'T' || TRANS[0] == 'C' || TRANS[0] == 'N')){
        mexErrMsgIdAndTxt( "MATLAB:xormqr_m:invalidInput",
              "Expect input 2 to be 'T', 'C', or 'N'");
    }
    
    minmn = na;
    if (ma < na){
        minmn = ma;
    }
    
    k = minmn; /* default behavior */
    
    /* get k if provided */
    
    if (nrhs > 5) {
        if(mxGetClassID(prhs[5]) == mxDOUBLE_CLASS){
                kd = (double *) mxGetData(prhs[5]);
                k = (ptrdiff_t) kd[0];
        }else if (mxGetClassID(prhs[5]) == mxINT8_CLASS){
                k8 = (int8_T *) mxGetData(prhs[5]);
                k = (ptrdiff_t) k8[0];
        }else if (mxGetClassID(prhs[5]) == mxINT16_CLASS){
                k16 = (int16_T *) mxGetData(prhs[5]);
                k = (ptrdiff_t) k16[0];
        }else if (mxGetClassID(prhs[5]) == mxINT32_CLASS){
                k32 = (int32_T *) mxGetData(prhs[5]);
                k = (ptrdiff_t) k32[0];
        }else if (mxGetClassID(prhs[5]) == mxINT64_CLASS){
                k64 = mxGetData(prhs[5]);
                k = (ptrdiff_t) k64[0];
        }else if (mxGetClassID(prhs[5]) == mxUINT8_CLASS){
                ku8 = mxGetData(prhs[5]);
                k = (ptrdiff_t) ku8[0];
        }else if (mxGetClassID(prhs[5]) == mxUINT16_CLASS){
                ku16 = mxGetData(prhs[5]);
                k = (ptrdiff_t) ku16[0];
        }else if (mxGetClassID(prhs[5]) == mxUINT32_CLASS){
                ku32 = mxGetData(prhs[5]);
                k = (ptrdiff_t) ku32[0];
        }else if (mxGetClassID(prhs[5]) == mxUINT64_CLASS){
                ku64 = mxGetData(prhs[5]);
                k = (ptrdiff_t) ku64[0];
        }else {
                mexErrMsgIdAndTxt( "MATLAB:xormqr_m:invalidInputType",
              "Expect input 6 to be one of double or (signed or unsigned) int8, int16, int32, int64");
        }
    }
    
    /* check array dimensions */
    
    if (SIDE[0] == 'L'){
        if ( ma != m ){
            mexErrMsgIdAndTxt( "MATLAB:xormqr_m:invalidInput",
              "Inputs 3 and 5 should have the same number of rows for left multiply");
        }
    } else {
        if ( ma != n ){
            mexErrMsgIdAndTxt( "MATLAB:xormqr_m:invalidInput",
              "Number of rows in input 3 should equal number of cols in input 5 for right multiply");
        }
    }
        
    if (nrhs > 5){
        if ( k < 0 || k > minmn ){
            mexErrMsgIdAndTxt( "MATLAB:xormqr_m:invalidInput",
              "Input 6 is invalid");
        }
    }
    
    
        
    /* allocate memory */
    
    cplx = (mxGetPi(prhs[2]) || mxGetPi(prhs[3]) || mxGetPi(prhs[4]) );
    
    if (cplx && TRANS[0] == 'T'){
        TRANS[0] = 'C';
    }
    
    if ( !(cplx) && TRANS[0] == 'C'){
        TRANS[0] = 'T';
    }
    
    if (cplx) {
      A = mat2fort(prhs[2],ma,na);
      TAU = mat2fort(prhs[3],ltau,1);
      B = mat2fort(prhs[4],m,n);
    } else {
      A = mxGetPr(prhs[2]);
      TAU = mxGetPr(prhs[3]);
      Bin = mxGetPr(prhs[4]);
      B = (double *)mxCalloc(m*n,sizeof(double));
      for (i=0; i<m*n; i++){
          B[i] = Bin[i];
      }
    }
    
    lwork = -1;
    info = 0;
    
    if (cplx) {
      work1 = (double *)mxCalloc(2,sizeof(double));
      /* Query zunmqr on the value of lwork */
      zunmqr(SIDE, TRANS, &m, &n, &k, A, &lda, TAU, B, &ldc, work1, &lwork, &info);
        if (info < 0) {
          sprintf(msg, "Input %d to zunmqr had an illegal value",-info);
          mexErrMsgIdAndTxt( "MATLAB:xormqr_m:illegalInputTozunmqr", msg);
        }
      lwork = (ptrdiff_t)(work1[0]);
      work = (double *)mxCalloc(2*lwork,sizeof(double));
      zunmqr(SIDE, TRANS, &m, &n, &k, A, &lda, TAU, B, &ldc, work, &lwork, &info);
      if (info < 0) {
          sprintf(msg, "Input %d to zunmqr had an illegal value",-info);
          mexErrMsgIdAndTxt( "MATLAB:xormqr_m:illegalInputTozunmqr", msg);
        }   
    } else {
        work1 = (double *)mxCalloc(1,sizeof(double));
      /* Query dormqr on the value of lwork */
        dormqr(SIDE, TRANS, &m, &n, &k, A, &lda, TAU, B, &ldc, work1, &lwork, &info);
        if (info < 0) {
          sprintf(msg, "Input %d to dormqr (1st call) had an illegal value",-info);
          mexErrMsgIdAndTxt( "MATLAB:xormqr_m:illegalInputTodormqr_1", msg);
      }
      lwork = (ptrdiff_t)(work1[0]);
      work = (double *)mxCalloc(lwork,sizeof(double));
      dormqr(SIDE, TRANS, &m, &n, &k, A, &lda, TAU, B, &ldc, work, &lwork, &info);
      if (info < 0) {
          sprintf(msg, "Input %d to dormqr (2nd call) had an illegal value",-info);
          mexErrMsgIdAndTxt( "MATLAB:xormqr_m:illegalInputTodormqr_2", msg);
      }
    }

    if (cplx) {
      plhs[0] = fort2mat(B,ldc,m,n);
      mxFree(A);
      mxFree(TAU);
    } else {
      plhs[0] = mxCreateDoubleMatrix(0,0,0);
      mxSetPr(plhs[0],B);
      mxSetM(plhs[0],m);
      mxSetN(plhs[0],n);
    }
     
    mxFree(work1);
    mxFree(work);

}
