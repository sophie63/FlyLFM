from PyTVAL3.solver.ftvcs_al_TVL2p import ftvcs_al_TVL2p

def TVAL3(A,b,dim,opts):
    '''
    Returns: [U,out]

    Accordingly choose the code based on the model selected by the user.

    1) TV model:        min sum ||D_i u||. 
        s.t. Au = b
    2) TV/L2 model:     min sum ||D_i u|| + mu/2||Au-b||_2^2 

    Please use the default one if the user does not have a specific model to
    solver.

    Written by: Chengbo Li
    Advisor: Prof. Yin Zhang and Wotao Yin
    Computational and Applied Mathematics department, Rice University
    May. 15, 2009
    '''
    
    if not opts.has_key('TVL2'):
        opts['TVL2'] = False

    if opts['TVL2']:
        return ftvcs_al_TVL2p(A,b,dim,opts)
    else: 
        return ftvcs_alp(A,b,dim,opts)
