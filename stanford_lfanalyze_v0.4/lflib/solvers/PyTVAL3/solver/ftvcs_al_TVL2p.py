import numpy as np
from PyTVAL3.solver.ftvcs_al_opts import ftvcs_al_opts
from scipy.sparse.linalg.interface import aslinearoperator

# ----------------------------------------------------------

def D(U):
    '''
    Returns: [Dux,Duy]
    '''

    # [ux,uy] = D u
    temp = U[:,0] - U[:,-1]
    Dux = np.hstack((np.diff(U,1,1), np.reshape(temp, (temp.shape[0],1),'f')))
    temp = U[0,:] - U[-1,:]
    Duy = np.vstack((np.diff(U,1,0), np.reshape(temp, (1,temp.shape[0]),'f')))
    return (Dux, Duy)

def D3(U):
    '''
    Returns: [Dux,Duy,Duz]
    '''

    # [ux,uy,uz] = D u
    temp = U[:,0,:] - U[:,-1,:]
    Dux = np.hstack((np.diff(U,1,1), np.reshape(temp, (temp.shape[0],1,1),'f')))
    temp = U[0,:,:] - U[-1,:,:]
    Duy = np.vstack((np.diff(U,1,0), np.reshape(temp, (1,temp.shape[0],1),'f')))
    temp = U[:,:,0] - U[:,:,-1]
    Duz = np.vstack((np.diff(U,1,2), np.reshape(temp, (1,1,temp.shape[0]),'f')))
    return (Dux, Duy, Duz)


def Dt(X,Y):
    '''
    Returns: DtXY

    DtXY = D_1.T * X + D_2.T * Y
    '''

    temp = X[:,-1] - X[:, 0]
    DtXY =        np.hstack( (np.reshape(temp, (temp.shape[0],1)), -np.diff(X,1,1)) );
    temp = Y[-1,:] - Y[0, :]
    DtXY = DtXY + np.vstack( (np.reshape(temp, (1,temp.shape[0])), -np.diff(Y,1,0)) );
    DtXY = np.reshape(DtXY, (np.prod(DtXY.shape), 1),'f')
    return DtXY

def Dt3(X,Y,Z):
    '''
    Returns: DtXYZ

    DtXYZ = D_1.T * X + D_2.T * Y + D_3.T * Z
    '''

    temp = X[:,-1,:] - X[:, 0,:]
    DtXYZ =        np.hstack( (np.reshape(temp, (temp.shape[0],1,1)), -np.diff(X,1,1)) );
    temp = Y[-1,:,:] - Y[0, :, :]
    DtXYZ = DtXYZ + np.vstack( (np.reshape(temp, (1,temp.shape[0],1)), -np.diff(Y,1,0)) );
    temp = Z[:,:,-1] - Z[:, :, 0]
    DtXYZ = DtXYZ + np.vstack( (np.reshape(temp, (1,1,temp.shape[0])), -np.diff(Y,1,2)) );
    DtXYZ = np.reshape(DtXYZ, (np.prod(DtXY.shape), 1),'f')
    return DtXYZ

# ----------------------------------------------------------

class ScaledLinearOperator(object):
    def __init__(self, A, scale_factor):
        self.A = A
        self.scale_factor = scale_factor

    def matvec(self, x ):
        return self.A.matvec(x) * self.scale_factor

    def rmatvec(self, x ):
        return self.A.rmatvec(x) * self.scale_factor


def ScaleA(n,mu,A,b,option):
    '''
    Returns: [mu,A,b]

    Scales mu, A and f so that the largest eigenvalue of A.T*A is 1 and the
    new problem
 
    min sum_i (||wi|| + beta/2 ||Diu - wi||^2) + mu/2 ||Au - b||^2
    
    is equivalent to the old one.  
 
    If option is assigned, mu will be scaled accordingly.
    
    Written by: Chengbo Li
    '''
    tol = .05;
    max_iterations = 10;

    # PORT : Handle real complex matrices here!
    #
    # if not isreal(A(rand(n,1),1)):
    #     eopts.isreal = false;

    # Compute largest eigenvalue by power iteration
    x = np.random.rand(A.shape[1],1)
    iteration = 0
    norm_z = tol + 1;
    while (iteration < max_iterations) and (norm_z > tol):
        z = A.rmatvec(A.matvec(x))
        norm_z = np.linalg.norm(z)
        print norm_z
        x = z / norm_z
        iteration += 1

    if np.real(norm_z) > 1 + 1e-10:
        if option:
            mu = mu*norm_z;
        b = b/np.sqrt(norm_z);
        A = ScaledLinearOperator(A, 1.0/np.sqrt(norm_z))

    return (mu, A, b)



def Scaleb(mu,b,option):
    '''
    Returns: [mu,b,scl]
    Scales mu and f so that the finite difference of f is neither too small 
    nor too large.

    If option is assigned, mu will be scaled accordingly.
    
    Written by: Chengbo Li
    '''

    threshold1 = .5;      # threshold is chosen by experience.
    threshold2 = 1.5;
    scl = 1.0;
    print b.shape
    b_dif = np.abs(np.max(b) - np.min(b));

    if b_dif < threshold1:
        scl = threshold1/b_dif;
        b = scl*b;
        if option:
            mu = mu/scl;
    elif b_dif > threshold2:
        scl = threshold2/b_dif;
        b = scl*b;
        if option:
            mu = mu/scl;

    return (mu, b, scl)

# -----------------------------------------------------------------------

def ftvcs_al_TVL2p(A,b,dims,opts):
    '''
    Returns: [U, out]
    Goal: solve    min sum ||D_i u|| + mu/2||Au-b||_2^2    (with or without 
        the constraint u>=0) to recover image/signal u from encoded b,

        which is equivalent to solve     min sum ||w_i|| + mu/2||Au-b||_2^2
            s.t. D_i u = w_i
    ftvcs_al solves Augmented Lagrangian function:

        min_{u,w} sum ||w_i|| - sigma^T(Du-w) 
           + beta/2 ||Du-w||_2^2 + mu/2||Au-b||_2^2 ,
 
    by an alternating algorithm:
    i)  while norm(up-u)/norm(up) > tol_inn
      1) Fix w^k, do Gradient Descent to 
        - sigma^T(Du-w^k) + beta/2||Du-w^k||^2 + mu/2||Au-f||^2;
        u^k+1 is determined in the following way:
          a) compute step length tau > 0 by BB formula
          b) determine u^k+1 by
                   u^k+1 = u^k - alpha*g^k,
             where g^k = -D^T*sigma + beta D.T(Du^k - w^k) + mu A.T(Au^k-f), 
             and alpha is determined by Amijo-like nonmonotone line search;
      2) Given u^k+1, compute w^k+1 by shrinkage
                  w^k+1 = shrink(Du^k+1-sigma/beta, 1/beta);
      end
    ii) update the Lagrangian multiplier by
              sigma^k+1 = sigma^k - beta(Du^k+1 - w^k+1).
    iii)accept current u as the initial guess to run the loop again
 
    Inputs:
        A        : either an matrix representing the measurement or a struct 
                   with 2 function handles:
                            A(x,1) defines @(x) A*x;
                            A(x,2) defines @(x) A^T*x;
        b        :  either real or complex input vector representing the noisy observation of a
                    grayscale image
        dims      :  size of the problem (tuple)
        opts     :  structure to restore parameters
 
 
    variables in this code:
 
    lam1 = sum ||wi||
    lam2 = ||Du-w||^2 (at current w).
    lam3 = ||Au-f||^2
    lam4 = sigma^T(Du-w)

        f  = lam1 + beta/2 lam2 + mu/2 lam3 - lam4

        g  = A^T(Au-f)
        g2 = D^T(Du-w) (coefficients beta and mu are not included)


    Numerical tests illustrate that this solver requirs large beta.


    Written by: Chengbo Li
    Advisor: Prof. Yin Zhang and Wotao Yin
    Computational and Applied Mathematics department, Rice University
    May. 12, 2009
    '''

    A = aslinearoperator(A)
    b = np.reshape(b, (b.shape[0], 1), 'f')

    # PORT: FINISH
    #
    #    global D Dt
    #    [D,Dt] = defDDt;

    # unify implementation of A
    #
    # PORT: Assume function handle for now...
    #
    #if ~isa(A,'function_handle'):
    #    A = @(u,mode) f_handleA(A,u,mode); 

    # get or check opts
    opts = ftvcs_al_opts(opts);

    # Create output struct
    out = {}

    # problem dimension
    n = np.prod(dims);

    # mark important constants
    mu = opts['mu']
    beta = opts['beta']
    tol_inn = opts['tol_inn']
    tol_out = opts['tol']
    gam = opts['gam']

    # check if A*A'=I
    #tmp = np.random.rand(b.shape[0],1);
    #if np.linalg.norm( A.matvec(A.rmatvec(tmp)) - tmp,1) / np.linalg.norm(tmp,1) < 1e-3:
    #    opts['scale_A'] = false;

    # check scaling A
    #
    if opts['scale_A']:
       [mu,A,b] = ScaleA(n,mu,A,b,opts['consist_mu']); 

    # check scaling b
    #
    scl = 1.0;
    if opts['scale_b']:
       (mu,b,scl) = Scaleb(mu,b,opts['consist_mu']);

    # calculate A'*b
    Atb = A.rmatvec(b);

    # initialize U, beta
    muf = mu;
    betaf = beta;     # final beta
    (U, mu, beta) = ftvcs_al_init(dims,Atb,scl,opts);    # U: sized according to dimss
    if mu > muf:
        mu = muf
    if beta > betaf:
        beta = betaf
    muDbeta = mu/beta;               # muDbeta: constant
    rcdU = U;

    # initialize multiplers
    sigmax = np.zeros(dims);        # sigmax, sigmay: p*q 
    sigmay = np.zeros(dims);

    # initialize D^T sigma
    DtsAtd = np.zeros((np.prod(dims),1)); 

    # initialize out.n2re
    if opts.has_key('Ut'):
        Ut = opts['Ut']*scl        #true U, just for computing the error
        nrmUt = np.linalg.norm(Ut,'fro')
    else:
        Ut = None

    if Ut != None:
        out['n2re'] = [np.linalg.norm(U - Ut,'fro')/nrmUt]; 

    # prepare for iterations
    out['mus'] = [mu]; out['betas'] = [beta];
    out['res'] = []; out['itrs'] = []; out['f'] = []; out['obj'] = []; out['reer'] = [];
    out['lam1'] = []; out['lam2'] = []; out['lam3'] = []; out['lam4'] = [];
    out['itr'] = np.Inf;
    out['tau'] = []; out['alpha'] = []; out['C'] = []; gp = None;
    out['cnt'] = [];

    (Ux,Uy) = D(U);                   # Ux, Uy: p*q

    opts['TVnorm'] = 2
    if opts['TVnorm'] == 1:
        Wx = np.maximum(np.abs(Ux) - 1.0/beta, 0.0)*np.sign(Ux);
        Wy = np.maximum(np.abs(Uy) - 1.0/beta, 0.0)*np.sign(Uy);
        lam1 = np.sum(np.sum(np.abs(Wx) + np.abs(Wy)));
    else:
        V = np.sqrt(Ux*np.conj(Ux) + Uy*np.conj(Uy));        # V: p*q
        V[np.nonzero(V==0)] = 1.0;
        S = np.maximum(V - 1.0/beta, 0.0)/V;        # S: p*q
        Wx = S*Ux;                              # Wx, Wy: p*q
        Wy = S*Uy;
        lam1 = np.sum(np.sum(np.sqrt(Wx*np.conj(Wx) + Wy*np.conj(Wy))));
    (lam2,lam3,lam4,f,g2,Au,g) = get_g(U,Ux,Uy,Wx,Wy,lam1,beta,mu,A,b,Atb,sigmax,sigmay);
    #lam, f: constant      g2: pq        Au: m         g: pq

    # compute gradient
    d = g2 + muDbeta*g - DtsAtd;

    count = 1;
    Q = 1; C = f;                     # Q, C: costant
    out['f'].append(f)
    out['C'].append(C)
    out['lam1'].append(lam1)
    out['lam2'].append(lam2)
    out['lam3'].append(lam3)
    out['lam4'].append(lam4)

    for ii in range(opts['maxit']):
        if opts['disp']:
            print('outer iter = %d, total iter = %d, normU = %4.2e;' % (count,ii,np.linalg.norm(U,'fro')))

        # compute tau first
        if gp != None:
            dg = g - gp;                        # dg: pq
            dg2 = g2 - g2p;                     # dg2: pq
            ss = np.dot(uup.T,uup);                      # ss: constant
            sy = np.dot(uup.T,(dg2 + muDbeta*dg));       # sy: constant
            # sy = uup'*((dg2 + g2) + muDbeta*(dg + g));
            # compute BB step length
            eps = np.spacing(1);  # Funny way of getting epsilon values in numpy
            tau = np.abs(ss/np.maximum(sy,eps))[0,0];               # tau: constant
            print 'tau = %f' % (tau)
        
            fst_itr = False;
        else:
            # do Steepest Descent at the 1st ieration
            #d = g2 + muDbeta*g - DtsAtd;         # d: pq
            print d.shape
            [dx,dy] = D(np.reshape(d,dims,'f'));
            dDd = np.square(np.linalg.norm(dx,'fro')) + np.square(np.linalg.norm(dy,'fro'))
            Ad = A.matvec(d)

            # compute Steepest Descent step length
            tau = np.abs( np.dot(d.T, d) / (dDd + muDbeta*np.dot(Ad.T,Ad)) )[0,0];

            # mark the first iteration 
            fst_itr = True;
    
        # keep the previous values
        Up = U; gp = g; g2p = g2; Aup = Au; 
        Uxp = Ux; Uyp = Uy;

        #############################
        # ONE-STEP GRADIENT DESCENT #
        #############################
        taud = tau*d;
        U = np.reshape(U, (np.prod(U.shape), 1),'f') - taud;
        # projected gradient method for nonnegtivity
        if opts['nonneg']:
            U = np.maximum(np.real(U),0);
        U = np.reshape(U,dims,'f');                    # U: p*q (still)
        (Ux,Uy) = D(U);                                # Ux, Uy: p*q

        (lam2,lam3,lam4,f,g2,Au,g) = get_g(U,Ux,Uy,Wx,Wy,lam1,beta,mu,A,b,Atb,sigmax,sigmay);

        # Nonmonotone Line Search
        alpha = 1;
        du = U - Up;                          # du: p*q
        const = opts['c']*beta*np.dot(d.T,taud);

        # Unew = Up + alpha*(U - Up)
        cnt = 0; flag = True;
    
        while f > C - alpha*const:
            if cnt == 5:
                # shrink gam
                gam = opts['rate_gam']*gam;

                # give up and take Steepest Descent step
                if opts['disp']:
                    print('    count of back tracking attains 5 ');

                #d = g2p + muDbeta*gp - DtsAtd;
                (dx,dy) = D(np.reshape(d,dims,'f'))
                dDd = np.square(np.linalg.norm(dx,'fro')) + np.square(np.linalg.norm(dy,'fro'))
                Ad = A.matvec(d)
                tau = np.abs(np.dot(d.T,d)/(dDd + muDbeta*np.dot(Ad.T,Ad)))
                U = np.reshape(Up, (np.prod(Up.shape), 1),'f') - tau*d;
                # projected gradient method for nonnegtivity
                if opts['nonneg']:
                    U = np.maximum(np.real(U),0);
                    U = np.reshape(U,dims,'f');
                (Ux, Uy) = D(U);
                Uxbar = Ux - sigmax/beta;
                Uybar = Uy - sigmay/beta;
                if opts['TVnorm'] == 1:
                    # ONE-DIMENSIONAL SHRINKAGE STEP
                    Wx = np.maximum(np.abs(Uxbar) - 1.0/beta, 0.0)*np.sign(Uxbar);
                    Wy = np.maximum(np.abs(Uybar) - 1.0/beta, 0.0)*np.sign(Uybar);
                    lam1 = np.sum(np.sum(np.abs(Wx) + np.abs(Wy)));
                else:
                    # TWO-DIMENSIONAL SHRINKAGE STEP
                    V = np.sqrt(Uxbar*np.conj(Uxbar) + Uybar*np.conj(Uybar)); # V: p*q
                    V[np.nonzero(V==0)] = 1.0;
                    S = np.maximum(V - 1.0/beta, 0.0)/V;                         # S: p*q
                    Wx = S*Uxbar;
                    Wy = S*Uybar;
                    lam1 = np.sum(np.sum(np.sqrt(Wx*np.conj(Wx) + Wy*np.conj(Wy))));
                (lam2,lam3,lam4,f,g2,Au,g) = get_g(U,Ux,Uy,Wx,Wy,lam1,beta,mu,A,b,Atb,sigmax,sigmay);
                alpha = 0; # remark the failure of back tracking
                break;
            if flag:
                dg = g - gp;
                dg2 = g2 - g2p;
                dAu = Au - Aup;                 # dAu: m
                dUx = Ux - Uxp;
                dUy = Uy - Uyp;
                flag = False;
            alpha = alpha*opts['gamma'];
            (U,lam2,lam3,lam4,f,Ux,Uy,Au,g,g2) = update_g(dims,lam1,alpha,beta,mu,Up,du,gp,dg,
                                                          g2p,dg2,Aup,dAu,Wx,Wy,Uxp,dUx,
                                                          Uyp,dUy,b,sigmax,sigmay);
            cnt = cnt + 1;
    
        # if back tracking is succeceful, then recompute
        if alpha != 0:
            Uxbar = Ux - sigmax/beta;
            Uybar = Uy - sigmay/beta;
            if opts['TVnorm'] == 1:
                # ONE-DIMENSIONAL SHRINKAGE STEP
                Wx = np.maximum(np.abs(Uxbar) - 1.0/beta, 0.0)*np.sign(Uxbar);
                Wy = np.maximum(np.abs(Uybar) - 1.0/beta, 0.0)*np.sign(Uybar);
            else:
                # TWO-DIMENSIONAL SHRINKAGE STEP
                V = np.sqrt(Uxbar*np.conj(Uxbar) + Uybar*np.conj(Uybar));
                V[np.nonzero(V==0)] = 1.0;
                S = np.maximum(V - 1.0/beta, 0.0)/V;
                Wx = S*Uxbar;
                Wy = S*Uybar;

            # update parameters related to Wx, Wy
            (lam1,lam2,lam4,f,g2) = update_W(beta,Wx,Wy,Ux,Uy,sigmax,
                                             sigmay,lam1,lam2,lam4,f,opts['TVnorm']);
    
        # update reference value
        Qp = Q; Q = gam*Qp + 1; C = (gam*Qp*C + f)/Q;
        uup = U - Up; uup = np.reshape(uup, (np.prod(uup.shape),1),'f');           # uup: pq
        nrmuup = np.linalg.norm(uup,'fro');                   # nrmuup: constant
    
        out['res'].append(nrmuup);
        out['f'].append(f); out['C'].append(C); out['cnt'].append(cnt);
        out['lam1'].append(lam1); out['lam2'].append(lam2); 
        out['lam3'].append(lam3); out['lam4'].append(lam4);
        out['tau'].append(tau); out['alpha'].append(alpha);

        if Ut != None:
            out['n2re'].append(norm(U - Ut,'fro')/norm(Ut,'fro')); 

        nrmup = np.linalg.norm(Up,'fro');
        RelChg = nrmuup/nrmup;

        # recompute gradient
        d = g2 + muDbeta*g - DtsAtd;
    
        if RelChg < tol_inn and not fst_itr:
            count = count + 1;
            RelChgOut = np.linalg.norm(U-rcdU,'fro')/nrmup;
            out['reer'].append(RelChgOut)
            rcdU = U;
            out['obj'].append(f + lam4)
            if len(out['itrs']) == 0:
                out['itrs'].append(ii);
            else:
                out['itrs'].append(ii - np.sum(out['itrs']));

            # stop if already reached final multipliers
            if RelChgOut < tol_out or count > opts['maxcnt']:
                if opts['isreal']:
                    U = np.real(U);
                U = U/scl;
                out['itr'] = ii;
                print 'Number of total iterations is %d. \n' % (out['itr'])
                return (U, out)
        
            # update multipliers
            (sigmax,sigmay,lam4,f) = update_mlp(beta,Wx,Wy,Ux,Uy,sigmax,sigmay,lam4,f);
        
            # update penality parameters for continuation scheme
            beta0 = beta;
            beta = beta*opts['rate_ctn'];
            mu = mu*opts['rate_ctn'];
            if beta > betaf: beta = betaf
            if mu > muf: mu = muf
            muDbeta = mu/beta;
            out['mus'].append(mu); out['betas'].append(beta)


            # update function value, gradient, and relavent constant
            f = lam1 + beta/2*lam2 + mu/2*lam3 - lam4;
            DtsAtd = -(beta0/beta)*d;     # DtsAtd should be divded by new beta instead of the old one for consistency!  
            d = g2 + muDbeta*g - DtsAtd;

            #initialize the constants
            gp = None;
            gam = opts['gam']; Q = 1; C = f;

    if opts['isreal']:
        U = np.real(U);
    print 'Attained the maximum of iterations %d. \n' % (opts['maxit'])
    U = U/scl;

    return (U, out)

# ----------------------------------------------------------

def ftvcs_al_init(dims,Atb,scl,opts):
    '''
    Returns: [U,mu,beta]
    '''

    # initialize mu beta
    if opts.has_key('mu0'):
        mu = opts['mu0']
    else:
        raise ValueError('Initial mu is not provided.')

    if opts.has_key('beta0'):
        beta = opts['beta0']
    else:
        raise ValueError('Initial beta is not provided.')

    # initialize U
    if opts['init'] == 0:
        U = zeros(dims);
    elif opts['init'] == 1:
        U = np.reshape(Atb,dims,'f')
    else:
        U = opts['init']*scl;
        if dims != U.size():
            raise ValueError('Input initial guess has incompatible size! Switch to the default initial guess. \n');
    return (U, mu, beta)

def get_g(U,Ux,Uy,Wx,Wy,lam1,beta,mu,A,b,Atb,sigmax,sigmay):
    '''
    Returns: [lam2,lam3,lam4,f,g2,Au,g]
    '''

    # A*u 
    Au = A.matvec(np.reshape(U, (np.prod(U.shape),1), 'f'))

    # g
    g = A.rmatvec(Au) - Atb;

    # lam2
    Vx = Ux - Wx;
    Vy = Uy - Wy;
    lam2 = np.sum(np.sum(Vx*np.conj(Vx) + Vy*np.conj(Vy)));

    # g2 = D'(Du-w)
    g2 = Dt(Vx,Vy);

    # lam3
    Aub = Au-b;
    lam3 = np.square(np.linalg.norm(Aub,'fro'))

    #lam4
    lam4 = np.sum(np.sum(np.conj(sigmax)*Vx + np.conj(sigmay)*Vy));

    # f
    f = lam1 + beta/2*lam2 + mu/2*lam3 - lam4;

    return (lam2, lam3, lam4, f, g2, Au, g)


def update_g(dims,lam1,alpha,beta,mu,Up,du,gp,dg,g2p,dg2,Aup,dAu,Wx,Wy,Uxp,dUx,Uyp,dUy,b,sigmax,sigmay):
    '''
    Returns: [U,lam2,lam3,lam4,f,Ux,Uy,Au,g,g2]
    '''

    g = gp + alpha*dg;
    g2 = g2p + alpha*dg2;
    U = Up + alpha*np.reshape(du,dims,'f');
    Au = Aup + alpha*dAu;
    Ux = Uxp + alpha*dUx;
    Uy = Uyp + alpha*dUy;

    Vx = Ux - Wx;
    Vy = Uy - Wy;
    lam2 = np.sum(np.sum(Vx*np.conj(Vx) + Vy*np.conj(Vy)));
    Aub = Au-b;
    lam3 = np.square(np.linalg.norm(Aub,'fro'));
    lam4 = np.sum(np.sum(np.conj(sigmax)*Vx + np.conj(sigmay)*Vy));
    f = lam1 + beta/2*lam2 + mu/2*lam3 - lam4;
    
    return (U,lam2,lam3,lam4,f,Ux,Uy,Au,g,g2)


def update_W(beta,Wx,Wy,Ux,Uy,sigmax,sigmay,lam1,lam2,lam4,f,option):
    '''
    Returns: [lam1,lam2,lam4,f,g2]
    '''

    # update parameters because Wx, Wy were updated
    tmpf = f -lam1 - beta/2*lam2 + lam4;
    if option == 1:
        lam1 = np.sum(np.sum(np.abs(Wx) + np.abs(Wy)));
    else:
        lam1 = np.sum(np.sum(np.sqrt(np.square(Wx) + np.square(Wy))));
    Vx = Ux - Wx;
    Vy = Uy - Wy;
    g2 = Dt(Vx,Vy);
    lam2 = np.sum(np.sum(Vx*np.conj(Vx) + Vy*np.conj(Vy)));
    lam4 = np.sum(np.sum(np.conj(sigmax)*Vx + np.conj(sigmay)*Vy));
    f = tmpf +lam1 + beta/2*lam2 - lam4;

    return (lam1,lam2,lam4,f,g2)

def update_mlp(beta, Wx,Wy,Ux,Uy,sigmax,sigmay,lam4,f):
    '''
    Returns: [sigmax,sigmay,lam4,f]
    '''

    Vx = Ux - Wx;
    Vy = Uy - Wy;
    sigmax = sigmax - beta*Vx;
    sigmay = sigmay - beta*Vy;

    tmpf = f + lam4;
    lam4 = np.sum(np.sum(np.conj(sigmax)*Vx + np.conj(sigmay)*Vy));
    f = tmpf - lam4;

    return (sigmax,sigmay,lam4,f)

# --------------------------------------------------------
