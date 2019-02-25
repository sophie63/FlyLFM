# __BEGIN_LICENSE__
#
# Copyright (C) 2010-2012 Stanford University.
# All rights reserved.
#
# __END_LICENSE__

'''
This small collection of projective geometry routines are useful for
computeng the intersection of arbitrary planes in 3-space.  These are
useful for computing rayspreads, among other things.
'''

import numpy as np
from numpy.linalg import *

# ---------------------------------------------------------------------------------

def ProjectivePlane(a, b, c, x0, y0, z0):
    return np.array([a, b, c, -a*x0-b*y0-c*z0])

def null_space(A, eps=1e-3):
  (u,s,vh) = svd(A,full_matrices=1,compute_uv=1)
  if len(s) < vh.shape[1]:  # Pad singular values if matrix is non-square
      ss = np.zeros(vh.shape[1])
      ss[0:len(s)] = s
      s = ss
  null_space = np.compress(s <= eps, vh, axis=0)
  return null_space.T

def plane_intersect2(p1, p2, p3):
    result = null_space(np.array([p1, p2, p3]))
    return result.flatten() / result[3]

def plane_intersect(p1, p2, p3):
    '''
    Compute the intersection of three projective planes in 3-space.
    
    This closed-form solution was derived using Mathematica.
    '''
    a1 = p1[0]; b1 = p1[1]; c1 = p1[2]; d1 = p1[3];
    a2 = p2[0]; b2 = p2[1]; c2 = p2[2]; d2 = p2[3];
    a3 = p3[0]; b3 = p3[1]; c3 = p3[2]; d3 = p3[3];

    result = np.array([ (b3*c2*d1 - b2*c3*d1 - b3*c1*d2 + b1*c3*d2 + b2*c1*d3 - b1*c2*d3) /
                        (-a3*b2*c1 + a2*b3*c1 + a3*b1*c2 - a1*b3*c2 - a2*b1*c3 + a1*b2*c3),

                        (a3*c2*d1 - a2*c3*d1 - a3*c1*d2 + a1*c3*d2 + a2*c1*d3 - a1*c2*d3) /
                        (a3*b2*c1 - a2*b3*c1 - a3*b1*c2 + a1*b3*c2 + a2*b1*c3 - a1*b2*c3),

                        (a3*b2*d1 - a2*b3*d1 - a3*b1*d2 + a1*b3*d2 + a2*b1*d3 - a1*b2*d3) /
                        (-a3*b2*c1 + a2*b3*c1 + a3*b1*c2 - a1*b3*c2 - a2*b1*c3 + a1*b2*c3),

                        1.0 ])
    return result
