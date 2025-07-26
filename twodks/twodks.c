#include "mex.h"

// Fortran subroutine declaration
extern void mttwo27_(
    int *n, int *p, int *w, int *u, int *cw, int *cu,
    int *z, int *x, int *lub, int *jdim, int *maxit,
    int *nsmall, int *maxr, int *maxl, int *levl,
    int *levk, int *jck, int *nnode, int *iscale
);

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]) {
    if (nrhs != 5) {
        mexErrMsgIdAndTxt("TwoConstraintKP:InvalidInput",
                          "Expected 5 inputs: P (int32), W (int32), U (int32), CW, CU");
    }

    if (!mxIsInt32(prhs[0]) || !mxIsInt32(prhs[1]) || !mxIsInt32(prhs[2])) {
        mexErrMsgIdAndTxt("TwoConstraintKP:InputType",
                          "P, W, and U must be int32 arrays.");
    }

    int N = (int) mxGetNumberOfElements(prhs[0]);
    int *P = (int *) mxGetData(prhs[0]);
    int *W = (int *) mxGetData(prhs[1]);
    int *U = (int *) mxGetData(prhs[2]);
    int CW = (int)(*mxGetPr(prhs[3]));
    int CU = (int)(*mxGetPr(prhs[4]));

    int X[10000];
    int Z, LUB = 0, JDIM = 10000, MAXIT = 100000;
    int NSMALL = 0, MAXR = 100, MAXL = 100;
    int LEVL = 0, LEVK = 0, JCK = 1, NNODE, ISCALE = 1;

    // Call Fortran solver
    mttwo27_(&N, P, W, U, &CW, &CU, &Z, X, &LUB,
             &JDIM, &MAXIT, &NSMALL, &MAXR, &MAXL,
             &LEVL, &LEVK, &JCK, &NNODE, &ISCALE);

    // Output: X and Z
    plhs[0] = mxCreateDoubleMatrix(N, 1, mxREAL);
    double *outX = mxGetPr(plhs[0]);
    for (int i = 0; i < N; ++i) {
        outX[i] = (double)X[i];
    }

    if (nlhs > 1) {
        plhs[1] = mxCreateDoubleScalar((double)Z);
    }
}
