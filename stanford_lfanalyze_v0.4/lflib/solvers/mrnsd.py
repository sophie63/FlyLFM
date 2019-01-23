# __BEGIN_LICENSE__
#
# Copyright (C) 2010-2013 Stanford University.
# All rights reserved.
#
# __END_LICENSE__

import numpy as np
import time

from scipy.sparse.linalg.interface import LinearOperator

from lflib.lightfield import LightField
from lflib.imageio import save_image
from lflib.linear_operators import LightFieldOperator, RegularizedNormalEquationLightFieldOperator

# ----------------------------------------------------------------------------------------
#                   Modified Residual Norm Steepest Descent SOLVER
# ----------------------------------------------------------------------------------------

def mrnsd_reconstruction(A, b, Rtol = 1e-6, NE_Rtol = 1e-6, max_iter = 100, x0 = None):
    '''
    Modified Residual Norm Steepest Descent
    Nonnegatively constrained steepest descent method.

    Ported from the RestoreTools MATLAB package available at:
    http://www.mathcs.emory.edu/~nagy/RestoreTools/

    Input: A  -  object defining the coefficient matrix.
           b  -  Right hand side vector.

 
    Optional Intputs:
      options - Structure that can have:
           x0      - initial guess (must be strictly positive); default is x0 = 1
           max_iter - integer specifying maximum number of iterations;
                      default is 100
           Rtol    - stopping tolerance for the relative residual,
                     norm(b - A*x)/norm(b)
                     default is 1e-6
           NE_Rtol - stopping tolerance for the relative residual,
                     norm(A.T*b - A.T*A*x)/norm(A.T*b)
                     default is 1e-6

    Output:
          x  -  solution

    Original MATLAB code by J. Nagy, August, 2011

    References:
    [1]  J. Nagy, Z. Strakos.
         "Enforcing nonnegativity in image reconstruction algorithms"
         in Mathematical Modeling, Estimation, and Imaging, 
         David C. Wilson, et.al., Eds., 4121 (2000), pg. 182--190. 
    [2]  L. Kaufman.
         "Maximum likelihood, least squares and penalized least squares for PET",
         IEEE Trans. Med. Imag. 12 (1993) pp. 200--214.
    '''

    # The A operator represents a large, sparse matrix that has dimensions [ nrays x nvoxels ]
    nrays = A.shape[0]
    nvoxels = A.shape[1]

    # Pre-compute some values for use in stopping criteria below
    b_norm = np.linalg.norm(b)
    trAb = A.rmatvec(b)
    trAb_norm = np.linalg.norm(trAb)

    # Start the optimization from the initial volume of a focal stack.
    if x0 != None:
        x = x0
    else:
        x = np.ones(nvoxels)

    Rnrm = np.zeros(max_iter+1);
    Xnrm = np.zeros(max_iter+1);
    NE_Rnrm = np.zeros(max_iter+1);

    eps = np.spacing(1)
    tau = np.sqrt(eps);
    sigsq = tau;
    minx = x.min()

    # If initial guess has negative values, compensate
    if minx < 0:
        x = x - min(0,minx) + sigsq;

    #  Initialize some values before iterations begin.
    Rnrm = np.zeros((max_iter+1, 1))
    Xnrm = np.zeros((max_iter+1, 1))

    # Initial Iteration
    r = b - A.matvec(x)
    g = -(A.rmatvec(r));
    xg = x * g;
    gamma = np.dot(g.T, xg);

    for i in range(max_iter):
        tic = time.time()
        x_prev = x
        
        # STEP 1: MRNSD Update step
        s = - x * g;
        u = A.matvec(s);
        theta = gamma / np.dot(u.T, u);
        neg_ind = np.nonzero(s < 0)

        zero_ratio = -x[neg_ind] / s[neg_ind]
        if zero_ratio.shape[0] == 0:
            alpha = theta
        else:
            alpha = min( theta, zero_ratio.min() );
        x = x + alpha*s;

        g = g + alpha * A.rmatvec(u);
        xg = x * g;
        gamma = np.dot(g.T, xg);

        # STEP 2: Compute residuals and check stopping criteria
        Rnrm[i] = np.sqrt(gamma) / b_norm
        Xnrm[i] = np.linalg.norm(x - x_prev) / nvoxels

        toc = time.time()
        print '\t--> [ MRNSD Iteration %d   (%0.2f seconds) ] ' % (i, toc-tic)
        print '\t      Residual Norm: %0.4g               (tol = %0.2e)  ' % (Rnrm[i], Rtol)
        print '\t        Update Norm: %0.4g                              ' % (Xnrm[i])

        # stop because residual satisfies ||b-A*x|| / ||b||<= Rtol
        if Rnrm[i] <= Rtol:
            break

    return x.astype(np.float32)


# ----------------------------------------------------------------------------------------
#             Weighted Modified Residual Norm Steepest Descent SOLVER
# ----------------------------------------------------------------------------------------

def wmrnsd_reconstruction(A, b,
                         Rtol = 1e-6, NE_Rtol = 1e-6, max_iter = 100, x0 = None,
                         sigmaSq = 0.0, beta = 0.0):
    '''
    Modified Residual Norm Steepest Descent
    Nonnegatively constrained steepest descent method.

    Ported from the RestoreTools MATLAB package available at:
    http://www.mathcs.emory.edu/~nagy/RestoreTools/

    Input: A  -  object defining the coefficient matrix.
           b  -  Right hand side vector.

 
    Optional Intputs:
      options - Structure that can have:
           x0      - initial guess (must be strictly positive); default is x0 = 1
           sigmaSq - the square of the standard deviation for the
                     white Gaussian read noise (variance)
           beta    - Poisson parameter for background light level  

           max_iter - integer specifying maximum number of iterations;
                      default is 100
           Rtol    - stopping tolerance for the relative residual,
                     norm(b - A*x)/norm(b)
                     default is 1e-6
           NE_Rtol - stopping tolerance for the relative residual,
                     norm(A.T*b - A.T*A*x)/norm(A.T*b)
                     default is 1e-6

    Output:
          x  -  solution

    Original MATLAB code by J. Nagy, August, 2011

    References:
    [1]  J. Nagy, Z. Strakos.
         "Enforcing nonnegativity in image reconstruction algorithms"
         in Mathematical Modeling, Estimation, and Imaging, 
         David C. Wilson, et.al., Eds., 4121 (2000), pg. 182--190. 
    [2]  L. Kaufman.
         "Maximum likelihood, least squares and penalized least squares for PET",
         IEEE Trans. Med. Imag. 12 (1993) pp. 200--214.
    '''

    # The A operator represents a large, sparse matrix that has dimensions [ nrays x nvoxels ]
    nrays = A.shape[0]
    nvoxels = A.shape[1]

    # Pre-compute some values for use in stopping criteria below
    b_norm = np.linalg.norm(b)
    trAb = A.rmatvec(b)
    trAb_norm = np.linalg.norm(trAb)

    # Start the optimization from the initial volume of a focal stack.
    if x0 != None:
        x = x0
    else:
        x = np.ones(nvoxels)

    Rnrm = np.zeros(max_iter+1);
    Xnrm = np.zeros(max_iter+1);
    NE_Rnrm = np.zeros(max_iter+1);

    eps = np.spacing(1)
    tau = np.sqrt(eps);
    sigsq = tau;
    minx = x.min()

    # If initial guess has negative values, compensate
    if minx < 0:
        x = x - min(0,minx) + sigsq;

    #  Initialize some values before iterations begin.
    Rnrm = np.zeros((max_iter+1, 1))
    Xnrm = np.zeros((max_iter+1, 1))

    # Initial Iteration
    c = b + sigmaSq;
    b = b - beta;

    r = b - A.matvec(x);
    trAr = A.rmatvec(r);

    wt = np.sqrt(c);

    for i in range(max_iter):
        tic = time.time()
        x_prev = x
        
        # STEP 1: WMRNSD Update step
        v = A.rmatvec(r/c)
        d = x * v;  
        w = A.matvec(d);
  
        w = w/wt;
  
        tau_uc = np.dot(d.T,v) / np.dot(w.T,w);
        neg_ind = np.nonzero(d < 0)

        zero_ratio = -x[neg_ind] / d[neg_ind]
        if zero_ratio.shape[0] == 0:
            tau = tau_uc;
        else:
            tau_bd = np.min( zero_ratio );
            tau = min(tau_uc, tau_bd);

        x = x + tau*d;
        w = w * wt;
  
        r = r - tau*w;
        trAr = A.rmatvec(r);

        # STEP 2: Compute residuals and check stopping criteria
        Rnrm[i] = np.linalg.norm(r) / b_norm
        Xnrm[i] = np.linalg.norm(x - x_prev) / nvoxels
        NE_Rnrm[i] = np.linalg.norm(trAr) / trAb_norm

        toc = time.time()
        print '\t--> [ MRNSD Iteration %d   (%0.2f seconds) ] ' % (i, toc-tic)
        print '\t      Residual Norm: %0.4g               (tol = %0.2e)  ' % (Rnrm[i], Rtol)
        print '\t         Error Norm: %0.4g               (tol = %0.2e)  ' % (NE_Rnrm[i], NE_Rtol)
        print '\t        Update Norm: %0.4g                              ' % (Xnrm[i])

        # stop because residual satisfies ||b-A*x|| / ||b||<= Rtol
        if Rnrm[i] <= Rtol:
            break

        # stop because normal equations residual satisfies ||A'*b-A'*A*x|| / ||A'b||<= NE_Rtol
        if NE_Rnrm[i] <= NE_Rtol:
            break

    return x.astype(np.float32)



#---------------------------------------------------------------------------------------

if __name__ == "__main__":
    pass
