import numpy as np

def isscalar(x): return isinstance(x, (int, long, float))
def islogical(x): return isinstance(x, (bool))

def ftvcs_al_opts(opts):
    '''
    Returns: opts

    Set default options.
    Written by: Chengbo Li
    '''

    if opts.has_key('mu'):
        if not isscalar(opts['mu']) or opts['mu'] <0:
            raise ValueError('opts.mu must be positive.');
        elif opts['mu'] > np.power(2,13) or opts['mu'] < np.power(2,4):
            print 'WARING: Users may want to choose opts.mu between 2^4 and 2^13 as a priority.';
    else:
        opts['mu'] = np.power(2,8);

    # mu is mainly decided by noise level. Set mu big when b is noise-free
    # whereas set mu small when b is very noisy.

    if opts.has_key('beta'):
        if not isscalar(opts['beta']) or opts['beta'] <0:
            raise ValueError('opts.beta must be positive.');
        elif opts['beta'] > np.power(2,13) or opts['beta'] < np.power(2,4):
            print 'WARNING: uusers may want to choose opts.beta  between 2^4 and 2^13 as a priority.';

    else:
        opts['beta'] = np.power(2,5);

    # outer loop tolerence
    if opts.has_key('tol'):
        if not isscalar(opts['tol']) or opts['tol'] <= 0:
            raise ValueError('opts.tol should be a positive small number.');
    else:
        opts['tol'] = 1.e-6;


    # inner loop tolerence
    if opts.has_key('tol_inn'):
        if not isscalar(opts['tol_inn']) or opts['tol_inn'] <= 0:
            raise ValueError('opts.tol_inn should be a positive small number.');
    else:
        opts['tol_inn'] = 1.e-3;


    if opts.has_key('maxcnt'):
        if not isscalar(opts['maxcnt']) or opts['maxcnt'] <= 0:
            raise ValueError('opts.maxcnt should be a positive integer.');
    else:
        opts['maxcnt'] = 10;



    if opts.has_key('maxit'):
        if not isscalar(opts['maxit']) or opts['maxit'] <= 0:
            raise ValueError('opts.maxit should be a positive integer.');
    else:
        opts['maxit'] = 1025;



    if opts.has_key('init'):
        if length(opts['init']) != 1:
            print('User has supplied opts.init as initial guess matrix......');
        elif not isinInterval(opts['init'],0,1,True) or opts['init'] != floor(opts['init']):
            raise ValueError('opts.init should be either 0/1 or an initial guess matrix.');
    else:
        opts['init'] = 1;



    if opts.has_key('disp'):
        if not islogical(opts['disp']):
            raise ValueError('opts.disp should be true or false.');
    else:
        opts['disp'] = False;



    if opts.has_key('scale_A'):
        if not islogical(opts['scale_A']):
            raise ValueError('opts.scale_A should be true or false.');
    else:
        opts['scale_A'] = True;



    if opts.has_key('scale_b'):
        if not islogical(opts['scale_b']):
            raise ValueError('opts.scale_b should be true or false.');
    else:
        opts['scale_b'] = True;



    if opts.has_key('consist_mu'):
        if not islogical(opts['consist_mu']):
            raise ValueError('opts.consist_mu should be true or false.');
    else:
        opts['consist_mu'] = False;

    # consist_mu decides if mu should be accordingly scaled while scaling A and
    # b. Strongly recomm  setting as 'false' if one try to recover a signal
    # or image instead of solving an exact minimization problem.

    if opts.has_key('mu0'):
        if not isscalar(opts['mu0']) or opts['mu0'] <= 0:
            raise ValueError('opts.mu0 is should be a positive number which is no bigger than beta.');
    else:
        opts['mu0'] = opts['mu'];  

    # initial mu

    if opts.has_key('beta0'):
        if not isscalar(opts['beta0']) or opts['beta0'] <= 0:
            raise ValueError('opts.beta0 is should be a positive number which is no bigger than beta.');

    else:
        opts['beta0'] = opts['beta']; 

# initial beta


    if opts.has_key('rate_ctn'):
        if not isscalar(opts['rate_ctn']) or opts['rate_ctn'] <= 1:
            raise ValueError('opts.rate_ctn is either not a scalar or no bigger than one.');

    else:
        opts['rate_ctn'] = 2;

# continuation parameter for both mu and beta


    if opts.has_key('c'):
        if not isscalar(opts['c']) or opts['c'] <= 0 or opts['c'] > 1:
            raise ValueError('opts.c should be a scalar between 0 and 1.');

    else:
        opts['c'] = 1.e-5;



    if opts.has_key('gamma'):
        if not isscalar(opts['gamma']) or opts['gamma'] <= 0 or opts['gamma'] > 1:
            raise ValueError('opts.gamma should be a scalar between 0 and 1.');

    else:
        opts['gamma'] = .6;



    if opts.has_key('gam'):
        if not isscalar(opts['gam']) or opts['gam'] <= 0 or opts['gam'] > 1:
            raise ValueError('opts.gam should be a scalar between 0 and 1.');

    else:
        opts['gam'] = .9995;

    # Control the degree of nonmonotonicity. 0 corresponds to monotone line search.
    # The best convergence is obtained by using values closer to 1 when the iterates
    # are far from the optimum, and using values closer to 0 when near an optimum.

    if opts.has_key('rate_gam'):
        if not isscalar(opts['rate_gam']) or opts['rate_gam'] <= 0 or opts['rate_gam'] > 1:
            raise ValueError('opts.rate_gam should be a scalar between 0 and 1.');
    else:
        opts['rate_gam'] = .9;

    # shrinkage rate of gam


    if opts.has_key('TVnorm'):
        if opts['TVnorm'] != 1 and opts['TVnorm'] != 2:
            raise ValueError('opts.TVnorm should be either 1(TV/L1 model) or 2(TV/L2 model).');
    else:
        opts['TVnorm'] = 2;



    if opts.has_key('nonneg'):
        if not islogical(opts['nonneg']):
            raise ValueError('opts.nonneg should be true or false.');
    else:
        opts['nonneg'] = False;



    if opts.has_key('isreal'):
        if not islogical(opts['isreal']):
            raise ValueError('opts.isreal should be true or false.');
    else:
        opts['isreal'] = False;



    if opts.has_key('TVL2'):
        if not islogical(opts['TVL2']):
            raise ValueError('opts.TVL2 should be true or false.');
    else:
        opts['TVL2'] = False;

    # Decide the model: TV or TV/L2. The default is TV model, which is recomm ed.


    if opts.has_key('tau'):
        if not isscalar(opts['tau']) or opts['tau'] <= 0:
            raise ValueError('opts.tau is not positive scalar.');
    else:
        opts['tau'] = 1.8;

    return opts
