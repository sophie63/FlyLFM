# __BEGIN_LICENSE__
#
# Copyright (C) 2010-2012 Stanford University.
# All rights reserved.
#
# __END_LICENSE__

# Major libraries
import h5py, time
import numpy as np
import pylab as pl

# Classes and functions
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
import rpy2.robjects.numpy2ri as np2ri; np2ri.activate()
from multiprocessing import Pool
import scipy.linalg as la

# Local functions
from group_lasso import _get_gram_factorization, _run_primal_group, center_scale
from r_functions import _get_grplasso_R
from lflib.volume import SimultaneousIterativeReconstruction
from lflib.lightfield import LightField
from lflib.imageio import save_image

# ----------------------------------------------------------------------------------------
#                            ADMM ITERATIVE DECONVOLUTION
# ----------------------------------------------------------------------------------------

#--------------------------------------------------------------------------#
# CLASS optProblem: a base class for optimization problems
#--------------------------------------------------------------------------#

class iterativeOptProblem( object ):
    """ 
    A base class for iterative optimization problems. 
    """

    def __init__(self, data, parameters):
        """ Initialize problem. """
        self.data = data
        self.parameters = parameters
        self.verbose = True
        self.err = np.infty
        self.errs = [0]
        self.iters = 0
        self.diagnostics = dict()
        self.results = dict()
        self._record_all_errors = True

        if 'tol' in parameters.keys():
            self.tol = parameters['tol']
        else:
            self.tol = 1e-5

        if 'max_iters' in parameters.keys():
            self.max_iters = parameters['max_iters']
        else:
            self.max_iters = 1e4

    def solve( self, verbose = self.verbose ):
        """ Solve the problem. """
        # if verbose, keep track of solve time
        if self.verbose:
            self.tic = time.time()

        # Update until the convergence criterion is met.
        while self._convergence_criterion is False:
            # iteration update of all model parameters
            self._full_update()
            # record a list of errors
            if self._record_all_errors:
                self.errs.extend( [self.error] )
            # if verbose, print things to screen
            if self.verbose:
                self._print_to_screen()
            # stop if max iterations reached
            if self.iters >= self.max_iters:
                print "!!!!!! MAXIMUM ITERATIONS REACHED !!!!!!" 
                break
            self.iters += 1

        if self.verbose:
            self.toc = time.time()
            print "\n"
            print "#------------------ PROBLEM SOLVED! -------------------------#"
            print "#"
            print "# Problem converged in ", self.iters, " iterations,"
            print "# to tolerance: ", self.tol
            print "# and final error: ", self.error
            print "# in ", self.toc - self.tic, " seconds."
            print "#"
            print "#------------------------------------------------------------#"

        self.diagnostics = self._get_diagnostics()
        self._get_results()
        self._clean_up()
        
    def _print_to_screen(self):
        """ Print per iteration diagnostics to screen """
        if self.iters == 0:
            print " ITERATION   ERR   CONV CRIT"
            print "-----------------------------------------------"
            print ( self.iters, "      ",
                    self.error, "      ",
                    self._convergence_criterion )
        else:
            print ( self.iters, "      ",
                    self.error, "      ",
                    self._convergence_criterion )

    # ---- The following are abstract methods to be overridden by subclasses ---- # 
   
    def _get_diagnostics( self ):
        """ Record diagnostics related to the problem solution. """
        pass

    def _full_update(self):
        """ Abstract method, replace in subclass """        
        pass

    def _get_results(self):
        """ Abstract method, replace in subclass """        
        self.results['parameters'] = self.parameters

    def _clean_up(self):
        """ Abstract method, replace in subclass """
        pass

    @property
    def error( self ):
        """ Abstract method, replace in subclass """
        return 1./((self.iters+1)**2)

    @property
    def _convergence_criterion( self ):
        """ Abstract method, replace in subclass """
        return 1./((self.iters+1)**2) < self.tol

#--------------------------------------------------------------------------#
# CLASS ADMMdeconvolution: a subclass of optProblem for solving the
# parallel group lasso problem
#--------------------------------------------------------------------------#

class ADMMdeconvolution( iterativeOptProblem ):

    """
    This subclass of class iterativeOptProblem solves the ADMM "sharing" problem:

    .. math:
    \begin{eqnarray*}
    \underset{x_{i,}z_{i}}{\text{minimize}} &  & l\left(\sum_{i=1}^{N}z_{i}-b\right)+\sum_{i=1}^{N}g_{i}(x_{i})\\
    \text{subject to} &  & A_{i}x_{i}-z_{i}=0,\ i=1,\ldots,N\\
    &  & -x_{i}\preceq0,\ i=1,\ldots,N,
    \end{eqnarray*}

    where 

    -- $l(\cdot)$ is a convex loss function, for example the
    standard least squares loss function $l(\cdot)=\|\cdot\|_{2}^{2}$,

    -- $g_i(\cdot)$ is a penalty (or {}``regularization'') function appropriate for
    the $x_i$ subvector, 

    -- the optimization variables $x_i$ and $z_i$ correspond to 
    the columns of the partitioned data matrix $A_i$, which are columns of the 
    original data matrix $A$, such that $A = [A_0 A_1 \cdots A_N]$

    -- the full matrix $A$ is the light field projection matrix, which projects
    the nonnegative volume $x$ into the light field $b$. 

    Because the matrix $A$ is large (on the order of $1e6 \times 1e6$, 
    we cannot easily solve the above problem using direct methods. This formulation solves 
    a version of the ADMM sharing problem (Boyd et al., Foundations and Trends in 
    Machine Learning 2011 -- see especially subsections 7.3 and 8.3).
    """

    def __init__(self, data, parameters):
        """ Initialize the class """
        # initialize the superclass
        optProblem.__init__(self, data, parameters)

        # get the data and its shape
        self.X = data['X']
        self.y = data['y']
        self.n, self.p  = self.X.shape

        # initialize model parameters
        self.lam = parameters['lam']
        self.rho = parameters['rho']
        if 'nonnegative' in parameters.keys():
            self.nonnegative = parameters['nonnegative']
        else:
            self.nonnegative = False

        # convergence tolerances
        self._abs_tol = self.tol
        self._rel_tol = self.tol

        # initialize fit-related things
        self._fit_init()

    def _fit_init(self):
        """ Used to initialize or re-initialize fitting variables """
        # weights to be used by adaptive PGL 
        self._weights = np.ones((self.n_groups,)) 

        # initialize various other things
        self._est_avg = np.zeros( (self.n,) )
        self._shared = np.zeros( (self.n,) )
        self._dual = np.zeros( (self.n,) )
        self._gram_factors = None
        self._primals = dict()
        self._last_coefficients = np.random.random((self.p,))
        self._primal_resid = np.inf
        self._dual_resid = np.inf
        self._primal_tol = self.tol
        self._dual_tol = self.tol

    def solve_path( self, lambda_min_frac, num = 10, lambdas = None ):
        """
        Find solution along a path of lambda values, beginning with the 
        largest lambda (smoothest model), and then solve for each model as lambda 
        decreases using the previous model as a warm start. If a vector of 
        lambdas is not provided explictly (i.e., lambdas = None), then 
        a sequence of lambdas of length 'num' between a calculated 
        value of lambda_max and lamda_min = lambda_min_frac*lambda_max
        will be used. 
        """
        self.path_results = dict()
        if lambdas is None:
            lam_max = self._get_lambda_max()
            lam_min = lam_max * lambda_min_frac
            lambdas = np.linspace( lam_min, lam_max, num=num )

        L = len( lambdas )
        for l in xrange(L):
            self.lam = lambdas[l]
            self.parameters['lam'] = self.lam
            self.solve()
            self._get_results()
            self.path_results[l] = self.results
            self.iters = 0 
            self.pool = Pool()

    def _full_update(self):
        """ A single iteration through all updates """
        if self.iters >= 1:
            self._last_coefficients = self.coefficients
        self._serial_primal_update()
        #self._parallel_primal_update()
        self._sharing_update()
        self._dual_update()
        
    def _serial_primal_update(self):
        """ Solve the primal ADMM updates serially """
        # Cache initial Gram matrix eigendecompositions
        if self._gram_factors is None:
            # verbosity...
            if self.verbose:
                print "#--------------------------------------------------------#."
                print "Getting eigendecompositions of group Gram matrices..."
                tic = time.clock()
                
            # instantiate pool of workers
            self.pool = Pool()

            # organize data for workers
            self.data_to_factor = [ (self.Xs[ t ], t) for t in xrange( self.n_groups ) ]

            # let the workers loose
            results = self.pool.map( _get_gram_factorization, self.data_to_factor )

            # collect results
            self._gram_idxs = []
            self._Qs = dict()
            self._zetas = dict()
            self._gram_factors = dict()

            for r in results:
                self._gram_idxs.append( r[0] )
                self._Qs[ r[0] ] = r[1]
                self._zetas[ r[0] ] = r[2]
                self._gram_factors[ r[0] ] = r[3]

            # report how long things took (if verbose)
            if self.verbose:
                toc = time.clock()
                print "Caching Gram matrix eigendecompositions took ", toc-tic, " seconds."
                print "#--------------------------------------------------------#."

        # once the Gram matrix decompositions and workers are already around,
        # if verbose...
        if self.verbose:
            tic = time.clock()

        for t in xrange( self.n_groups ):
            self._primals[t] =_run_primal_group( ( self.Xs[ t ],
                                                    self._Qs[ t ],
                                                    self._zetas[ t ],
                                                    self._group_ests[:,t],
                                                    self._est_avg,
                                                    self._shared,
                                                    self._dual,
                                                    self.lam*self._weights[t],
                                                    self.rho,
                                                    self._gram_factors[ t ],
                                                    t ) )[1]

        # report back how long things took (if verbose)
        if self.verbose:
            self._primal_timing = time.clock()-tic
            
        return self._primals

    def _parallel_primal_update(self):
        """ Solve the primal updates for grouped variables in parallel """
        # Cache initial Gram matrix eigendecompositions
        if self._gram_factors is None:
            # verbosity...
            if self.verbose:
                print "#--------------------------------------------------------#."
                print "Getting eigendecompositions of group Gram matrices..."
                tic = time.clock()
                
            # instantiate pool of workers
            self.pool = Pool()

            # organize data for workers
            self.data_to_factor = [ (self.Xs[ t ], t) for t in xrange( self.n_groups ) ]

            # let the workers loose
            results = self.pool.map( _get_gram_factorization, self.data_to_factor )

            # collect results
            self._gram_idxs = []
            self._Qs = dict()
            self._zetas = dict()
            self._gram_factors = dict()

            for r in results:
                self._gram_idxs.append( r[0] )
                self._Qs[ r[0] ] = r[1]
                self._zetas[ r[0] ] = r[2]
                self._gram_factors[ r[0] ] = r[3]

            # report how long things took (if verbose)
            if self.verbose:
                toc = time.clock()
                print "Caching Gram matrix eigendecompositions took ", toc-tic, " seconds."
                print "#--------------------------------------------------------#."

        # once the Gram matrix decompositions and workers are already around,
        # if verbose...
        if self.verbose:
            tic = time.time()

        # organize problems for workers
        self.problems = ( [ ( self.Xs[ t ],
                              self._Qs[ t ],
                              self._zetas[ t ],
                              self._group_ests[:,t],
                              self._est_avg,
                              self._shared,
                              self._dual,
                              self.lam*self._weights[t],
                              self.rho,
                              self._gram_factors[ t ],
                              t ) for t in xrange( self.n_groups ) ] )

        # release the workers!
        results = self.pool.map( _run_primal_group, self.problems )

        # collect results
        for r in results:
            self._primals[ r[0] ] = r[1]

        # report back how long things took (if verbose)
        if self.verbose:
            self._primal_timing = time.time()-tic
            
        return self._primals
    
    def _sharing_update(self):
        """ Combine primal updates in sharing step"""
        if self.verbose:
            tic = time.clock()

        # get group estimates and their mean
        for j in xrange(self.n_groups):
            self._group_ests[:,j] =  np.dot( self.Xs[j], self._primals[j] ).T
        self._est_avg = np.mean( self._group_ests, axis = 1 )
        
        # sharing is caring
        if self.n_groups != 0:
            a = 1./(self.n_groups*self.rho)
        else: 
            a = 1e-8
        old_shared = self._shared
        self._shared = a*( self.y + self.rho*(self._est_avg + self._dual) )
        self._shared_change = self._shared - old_shared
        
        # report back how long things took (if verbose)
        if self.verbose:
            self._shared_timing = time.clock()-tic

        return self._shared

    def _dual_update(self):
        """ Calculate dual variables """
        if self.verbose:
            tic = time.clock()
            
        self._dual = self._dual + self._est_avg - self._shared

        # report back how long things took (if verbose)
        if self.verbose:
            self._dual_timing = time.clock()-tic

        return self._dual

    def _get_results(self):
        """ Collect results of fit """        
        self.results = dict()
        self.results['parameters'] = self.parameters
        self.results['coefficients'] = self.coefficients
        self.results['error'] = self.error
        self.results['estimate'] = self.estimate 

    @property
    def coefficients(self):
        """ Reconstruct the coefficient vector from current primal vectors """
        coefs = np.array([])
        self._sub_coefs = []
        for g in xrange(self.n_groups):
            self._sub_coefs.append( self._primals[g] )
            coefs = np.append( coefs, self._primals[ g ] )
        if self.nonnegative:
            coefs = self._nonneg_constraint(coefs)
        return coefs

    @property
    def _convergence_criterion( self ):
        """
        Criterion for convergence: based on primal and dual residuals
        """
        if self.iters <= 1:
            return False
        else:
            N = self.n_groups
            rho = self.rho
            self._dual_prod = np.dot( self.X.T, self._dual )
            self._primal_tol = ( np.sqrt(self.p)*self._abs_tol +
                                 self._rel_tol*N*np.max( (la.norm( self._est_avg ),
                                                          la.norm( self._shared) )) )
            self._dual_tol = ( np.sqrt(self.n)*self._abs_tol +
                               self._rel_tol*la.norm(self._dual_prod)*rho*N ) 

            self._primal_resid = la.norm( N*(self._est_avg - self._shared) )
            self._dual_resid = la.norm( N*rho*np.dot( self.X.T, self._shared_change ) )
            if self._primal_resid <= self._primal_tol and self._dual_resid <= self._dual_tol:
                return True
            else:
                return False
            
    @property
    def error(self):
        """ Mean squared error of fit """
        # get current approximation to self.y
        self.predict()
        # return mean squared error (Note: this is where to sub in another loss fxn)
        return self._get_current_loss()

    def predict( self ):
        """ Get estimate for current data """
        self.estimate = np.dot( self.X, self.coefficients )

    def _get_current_loss(self, loss_type = 'MSE', *args):
        """ Loss function applied to fit residuals at each iteration """
        if loss_type is 'MSE':
            # get mean squared error
            return la.norm( self.y - self.estimate )**2
        else:
            raise NotImplementedError("Loss function has not yet been added")
    
    def _nonneg_constraint(self, coefs):
        """ Impose non-negativity constraint on the data by projection. """
        for c in xrange(self.p):
            if coefs[c] < 0: coefs[c] = 0.
        return coefs

    def _get_lambda_max(self):
        """ 
        Find the value of lambda at which all coefficients are set to zero
        by finding the minimum value such that 0 is in the subdifferential
        and the coefficients are all zero.
        """
        smooth_grads = np.zeros((self.n_groups,))
        for g in xrange(self.n_groups):
            smooth_grads[g] = la.norm( np.dot( self.Xs[g].T, self.y ) )
        return np.max( np.fabs( smooth_grads ) )

    def _clean_up( self ):
        """ Clean up workers, etc. """
        self.pool.close()

    def _print_to_screen(self):
        """ Print per iteration diagnostics to screen """
        trunc = 6
        if self.iters == 0:
            print " IT  ERR   PRIM_R   PRIM_TOL   DUAL_R   DUAL_TOL   PRIM_T   SHARE_T   DUAL_T"
            print "----------------------------------------------------------------------------------"
            print ( str( self.iters ) + '   ' + 
                    str( self.error )[0:trunc] + '   ' + 
                    str( self._primal_resid )[0:trunc] + '   ' + 
                    str( self._primal_tol )[0:trunc] + '   ' + 
                    str( self._dual_resid )[0:trunc] + '   ' +
                    str( self._dual_tol )[0:trunc] + '   ' + 
                    str( self._primal_timing )[0:trunc] + '   ' + 
                    str( self._shared_timing )[0:trunc] + '   ' +  
                    str( self._dual_timing )[0:trunc] ) 
        else:
            print ( str( self.iters ) + '   ' + 
                    str( self.error )[0:trunc] + '   ' + 
                    str( self._primal_resid )[0:trunc] + '   ' + 
                    str( self._primal_tol )[0:trunc] + '   ' + 
                    str( self._dual_resid )[0:trunc] + '   ' +
                    str( self._dual_tol )[0:trunc] + '   ' + 
                    str( self._primal_timing )[0:trunc] + '   ' + 
                    str( self._shared_timing )[0:trunc] + '   ' +  
                    str( self._dual_timing )[0:trunc] ) 

if __name__ == "__main__":
    pass
