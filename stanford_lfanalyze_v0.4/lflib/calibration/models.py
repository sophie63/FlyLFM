import numpy as np
import cv  # OpenCV
import math


# ------------------------------------------------------------------
#                     ISOMETRY TRANSFORM
#
# Allows for translation & rotation only.  Computed by solving the
# orthogonal Procrustes problem:
#
#   http://en.wikipedia.org/wiki/Orthogonal_Procrustes_problem
#
# ------------------------------------------------------------------

class IsometryWarp(object):

    def __init__(self, lenslet_pitch, pixel_size):
        self.forwardCoefficients = np.identity(3, dtype=np.float32)
        self.reverseCoefficients = np.identity(3, dtype=np.float32)
        self.lenslet_pitch = lenslet_pitch
        self.pixel_size = pixel_size

        
    def fit(self, data):
        n = data.shape[0];

        X1 = data[:,0:2] * self.pixel_size         # (Scaled) chief rays
        X2 = data[:,2: ] * self.lenslet_pitch      # (Scaled) lenslet centers

        # Compute column mean & matrices with mean removed.
        X1_mean = np.mean(X1, axis = 0);
        X2_mean = np.mean(X2, axis = 0);
        X1_hat = X1 - np.tile(X1_mean, (n, 1));
        X2_hat = X2 - np.tile(X2_mean, (n, 1));

        # Compute procrustes transform using the SVD
        (U,D,V) = np.linalg.svd(np.dot(X1_hat.T, X2_hat));

        R = np.dot(U, V.T);
        mu = X2_mean - np.dot(X1_mean, R);

        S1 = np.array( [[ self.pixel_size, 0.0       , 0.0],
                        [ 0.0       , self.pixel_size, 0.0],
                        [ 0.0       , 0.0       , 1.0]] )

        S2 = np.array( [[ 1.0/self.lenslet_pitch, 0.0       , 0.0],
                        [ 0.0       , 1.0/self.lenslet_pitch, 0.0],
                        [ 0.0       , 0.0                   , 1.0]] )

        self.forwardCoefficients = np.array([[ R[0,0], R[1,0], mu[0] ],
                                             [ R[0,1], R[1,1], mu[1] ], 
                                             [ 0.0   , 0.0   , 1.0   ]]);
        self.forwardCoefficients = np.dot(self.forwardCoefficients, S1)
        self.forwardCoefficients = np.dot(S2, self.forwardCoefficients)

        self.reverseCoefficients = np.linalg.inv(self.forwardCoefficients)
    

    def get_error(self, data):

        X1 = np.hstack( (data[:,0:2], np.ones((data.shape[0], 1)) ) )       # (Scaled) chief rays
        X2 = data[:,2: ]       # (Scaled) lenslet centers
        
        X1_prime = np.dot(self.forwardCoefficients, X1.T).T


        err_per_point = np.sqrt(np.sum(( X2 - X1_prime[:,0:2] )**2, axis=1)) # sum squared error per row
        return err_per_point

    # OUTPUT AND DIGANOSTIC ROUTINES

    def __str__(self):
        s = "\t\tTranslation: ["
        tmp = self.reverseCoefficients[0:2,2]
        for t in tmp:
            s+= ('%0.2f, ' % t)
        s+= "]\n\t\tRotation: ["
        s+= str(self.reverseCoefficients[0,0]) + ' ' + str(self.reverseCoefficients[0,1]) + '; '
        s+= str(self.reverseCoefficients[1,0]) + ' ' + str(self.reverseCoefficients[1,1]) + ']'
        return s

    def check_fit(self, data):

        X1 = np.hstack( (data[:,0:2], np.ones((data.shape[0], 1)) ) )       # (Scaled) chief rays
        X2 = data[:,2: ]       # (Scaled) lenslet centers
        
        X1_prime = np.dot(self.forwardCoefficients, X1.T).T

        diff = X2 - X1_prime[:,0:2]
        err = np.sqrt(np.power(diff[:,0],2) + np.power(diff[:,1],2))
        print ('\t--> Fit Statistics - Mean: %0.2f  Median: %0.2f  Min: %0.2f  Max: %0.2f (units: normalized lenslets)' %
               (err.mean(), np.median(err, axis=0), err.min(), err.max() ))


    def eval(self, data, direction='r'):

        X1 = np.hstack( (data[:,0:2], np.ones((data.shape[0], 1)) ) )       # (Scaled) chief rays

        if (direction == 'R' or direction == 'r'):
            X1_prime = np.dot(self.reverseCoefficients, X1.T)
        elif (direction == 'F' or direction == 'f'):
            X1_prime = np.dot(self.forwardCoefficients, X1.T)
        else:
            raise ValueError("Unknown direction in IsometryWarp.eval()")
        return X1_prime[0:2,:].T

    def eval_point(self, data, direction='r'):

        X1 = np.array([data[0], data[1], 1.0])
        if (direction == 'R' or direction == 'r'):
            X1_prime = np.dot(self.reverseCoefficients, X1.T)
        elif (direction == 'F' or direction == 'f'):
            X1_prime = np.dot(self.forwardCoefficients, X1.T)
        else:
            raise ValueError("Unknown direction in AffineWarp.eval()")
        return X1_prime[0:2]

    def forward(self, coords):
        result = np.dot(self.forwardCoefficients, np.array([coords[0], coords[1], 1.0]))
        return result[0:2]
    
    def reverse(self, coords):
        result = np.dot(self.reverseCoefficients, np.array([coords[0], coords[1], 1.0]))
        return result[0:2]

    # Warp an image using this isometry transform.
    #
    # This function expects input_image to be a numpy float32 array.
    #
    # the cropToInside option can be used to crop the output image to
    # the inside of the valid lenslets.  Note that cropping to the
    # inside will (currently) shift the image output relative to the
    # calibration matrix, so you should be cautious about employing
    # this option where you hope to use the calibration matrix for
    # further computations beyond this point.
    def warp_image(self, input_image, output_pixels_per_lenslet, direction="R",
                   cropToInside = False,  lenslet_offset = None, output_size = None):
        im = cv.fromarray(input_image)

        ul = self.eval_point([0, 0], 'f').flatten()
        ur = self.eval_point([0, im.cols], 'f').flatten()
        ll = self.eval_point([im.rows, 0], 'f').flatten()
        lr = self.eval_point([im.rows, im.cols], 'f').flatten()

        leftbound = np.ceil(max(ul[1], ll[1]))
        rightbound = np.floor(min(ur[1], lr[1]))
        topbound = np.ceil(max(ul[0], ur[0]))
        bottombound = np.floor(min(ll[0], lr[0]))

        # Don't crop left of lenslets (0, y) or above (x, 0)
        leftbound = max(leftbound, 0)
        topbound = max(topbound,0)

        if output_size != None:
            putative_output_size = output_size
        else:
            nt = int(np.floor(bottombound - topbound))
            ns = int(np.floor(rightbound - leftbound))
            putative_output_size = (nt*output_pixels_per_lenslet, ns*output_pixels_per_lenslet)

        # Create the output image
        output_image = cv.CreateMat(putative_output_size[0], putative_output_size[1], im.type)

        # Apply the transform.
        if (direction == 'f' or direction == 'F'):
            coeff = cv.CreateMat(2,3,cv.CV_32FC1)
            coeff[0,0] = self.forwardCoefficients[0,0] * output_pixels_per_lenslet
            coeff[0,1] = self.forwardCoefficients[0,1] * output_pixels_per_lenslet
            coeff[0,2] = self.forwardCoefficients[0,2]
            coeff[1,0] = self.forwardCoefficients[1,0] * output_pixels_per_lenslet
            coeff[1,1] = self.forwardCoefficients[1,1] * output_pixels_per_lenslet
            coeff[1,2] = self.forwardCoefficients[1,2]
            cv.WarpAffine(im, output_image, coeff,
                          flags=cv.CV_INTER_LINEAR+cv.CV_WARP_FILL_OUTLIERS+cv.CV_WARP_INVERSE_MAP)
        else:
            coeff = cv.CreateMat(2,3,cv.CV_32FC1)
            coeff[0,0] = self.reverseCoefficients[0,0] / output_pixels_per_lenslet
            coeff[0,1] = self.reverseCoefficients[0,1] / output_pixels_per_lenslet
            coeff[0,2] = self.reverseCoefficients[0,2]
            coeff[1,0] = self.reverseCoefficients[1,0] / output_pixels_per_lenslet
            coeff[1,1] = self.reverseCoefficients[1,1] / output_pixels_per_lenslet
            coeff[1,2] = self.reverseCoefficients[1,2]
            cv.WarpAffine(im, output_image, coeff,
                          flags=cv.CV_INTER_LINEAR+cv.CV_WARP_FILL_OUTLIERS+cv.CV_WARP_INVERSE_MAP)

        result = np.asarray(output_image)

        # Don't crop left of lenslets (0, y) or above (x, 0)
        leftbound = max(leftbound, 0)
        topbound = max(topbound,0)

        if cropToInside:
            return result[topbound * output_pixels_per_lenslet:bottombound * output_pixels_per_lenslet,
                          leftbound * output_pixels_per_lenslet:rightbound * output_pixels_per_lenslet]
        
        else:
            return result



# ------------------------------------------------------------------
#                          AFFINE WARP
#
# the model has 6 DOF capable of encoding rotation, anisotropic
# scaling, shearing, and translation.
#
#   x' = a + cx + ey 
#   y' = b + dx + fy 
#
# ------------------------------------------------------------------
class AffineWarp(object):

    def __init__(self):
        self.forwardCoefficients = np.zeros((2,3), dtype='float32')
        self.reverseCoefficients = np.zeros((2,3), dtype='float32')
	self.forwardCoefficients[0,1] = 1.0
	self.forwardCoefficients[1,2] = 1.0
	self.reverseCoefficients[0,1] = 1.0
	self.reverseCoefficients[1,2] = 1.0

    def eval(self, data, direction='r'):

        x = data[:,0]
        y = data[:,1]
        composite = np.vstack( (np.ones(x.shape), x, y) )
        if (direction == 'R' or direction == 'r'):
            return np.dot(composite.transpose(), self.reverseCoefficients.transpose())
        elif (direction == 'F' or direction == 'f'):
            return np.dot(composite.transpose(), self.forwardCoefficients.transpose())
        else:
            raise ValueError("Unknown direction in AffineWarp.eval()")

    def eval_point(self, data, direction='r'):

        x = np.mat(data[0])
        y = np.mat(data[1])
        composite = np.vstack( (np.mat(np.ones(x.shape)), np.mat(x), np.mat(y) ) )
        if (direction == 'R' or direction == 'r'):
            return np.array(np.dot(composite.transpose(), self.reverseCoefficients.transpose()))[0]
        elif (direction == 'F' or direction == 'f'):
            return np.array(np.dot(composite.transpose(), self.forwardCoefficients.transpose()))[0]
        else:
            raise ValueError("Unknown direction in AffineWarp.eval()")

    def fit(self, data):
        datalen = data.shape[0]
        A = np.zeros( (datalen,3) )
        bx = np.zeros( (datalen) )
        by = np.zeros( (datalen) )
        Ar = np.zeros( (datalen,3) )
        bxr = np.zeros( (datalen) )
        byr = np.zeros( (datalen) )
        
        for i in range(datalen):
            px = float(data[i,0])
            py = float(data[i,1])
            lx = float(data[i,2])
            ly = float(data[i,3])

            A[i,0] = 1.0;          Ar[i,0] = 1.0;
            A[i,1] = px;           Ar[i,1] = lx;
            A[i,2] = py;           Ar[i,2] = ly;

            bx[i] = lx;            bxr[i] = px
            by[i] = ly;            byr[i] = py

        self.forwardCoefficients = np.zeros((2,3), dtype='float32')
        self.reverseCoefficients = np.zeros((2,3), dtype='float32')

        import numpy.linalg as linalg
        self.forwardCoefficients[0,:] = linalg.lstsq(A,bx)[0]
        self.forwardCoefficients[1,:] = linalg.lstsq(A,by)[0]
        self.reverseCoefficients[0,:] = linalg.lstsq(Ar,bxr)[0]
        self.reverseCoefficients[1,:] = linalg.lstsq(Ar,byr)[0]

    def get_error(self, data, direction='f'):

        x = data[:,0]
        y = data[:,1]
        composite = np.vstack( (np.ones(x.shape), x, y ) )
        if (direction == 'R' or direction == 'r'):
            B_fit = np.dot(composite.transpose(), self.reverseCoefficients.transpose())
        elif (direction == 'F' or direction == 'f'):
            B_fit = np.dot(composite.transpose(), self.forwardCoefficients.transpose())
        else:
            raise ValueError("Unkown direction in AffineWarp.get_error()")

        B = data[:,2:4]
        err_per_point = np.sqrt(np.sum((B-B_fit)**2,axis=1)) # sum squared error per row
        return err_per_point

    def forwardCoefficients(self):
        return self.forwardCoefficients

    def reverseCoefficients(self):
        return self.reverseCoefficients

    def forward(self, coords):
        return np.dot(self.forwardCoefficients, np.array([1, coords[0], coords[1]]))

    def reverse(self, coords):
        return np.dot(self.reverseCoefficients, np.array([1, coords[0], coords[1]]))

    def save(self, filename):
        a = np.vstack((self.forwardCoefficients, self.reverseCoefficients))
        np.savetxt(filename, a, fmt="%12.6G") 
        
    def load(self, filename):
        a = np.loadtxt(filename)
        self.forwardCoefficients = a[0:2,:].astype(np.float32)
        self.reverseCoefficients = a[2:4,:].astype(np.float32)

    # OUTPUT AND DIGANOSTIC ROUTINES

    def __str__(self):
        s = "\t\tTranslation: ["
        tmp = self.reverseCoefficients[:,0];
        for t in tmp:
            s+= ('%0.2f, ' % t)
        s+= "]\n\t\tScaling and rotation: ["
        s+= str(self.reverseCoefficients[0,1]) + ' ' + str(self.reverseCoefficients[0,2]) + '; '
        s+= str(self.reverseCoefficients[1,1]) + ' ' + str(self.reverseCoefficients[1,2]) + ']'
        return s

    def check_fit(self, data):

        projected_camlenses = self.eval(data[:,0:2], 'f')
        diff = data[:,2:4] - projected_camlenses
        err = np.sqrt(np.power(diff[:,0],2) + np.power(diff[:,1],2))
        print ('\t--> Fit Statistics - Mean: %0.2f  Median: %0.2f  Min: %0.2f  Max: %0.2f (units: normalized lenslets)' %
               (err.mean(), np.median(err, axis=0), err.min(), err.max() ))

    # Warp an image using this affine transform.
    #
    # This function expects input_image to be a numpy float32 array.
    #
    # the cropToInside option can be used to crop the output image to
    # the inside of the valid lenslets.  Note that cropping to the
    # inside will (currently) shift the image output relative to the
    # calibration matrix, so you should be cautious about employing
    # this option where you hope to use the calibration matrix for
    # further computations beyond this point.
    def warp_image(self, input_image, output_pixels_per_lenslet, direction="R", cropToInside = False):
        im = cv.fromarray(input_image)

        ul = self.eval_point([0, 0], 'f').flatten()
        ur = self.eval_point([0, im.cols], 'f').flatten()
        ll = self.eval_point([im.rows, 0], 'f').flatten()
        lr = self.eval_point([im.rows, im.cols], 'f').flatten()

        leftbound = np.ceil(max(ul[1], ll[1]))
        rightbound = np.floor(min(ur[1], lr[1]))
        topbound = np.ceil(max(ul[0], ur[0]))
        bottombound = np.floor(min(ll[0], lr[0]))

        # Don't crop left of lenslets (0, y) or above (x, 0)
        leftbound = max(leftbound, 0)
        topbound = max(topbound,0)

        nt = int(np.floor(bottombound - topbound))
        ns = int(np.floor(rightbound - leftbound))
        outputsize = (nt*output_pixels_per_lenslet, ns*output_pixels_per_lenslet)

        # Create the output image
        output_image = cv.CreateMat(outputsize[0], outputsize[1], im.type)

        # Apply the transform.
        if (direction == 'f' or direction == 'F'):
            coeff = cv.CreateMat(2,3,cv.CV_32FC1)
            coeff[0,0] = self.forwardCoefficients[0,1] * output_pixels_per_lenslet
            coeff[0,1] = self.forwardCoefficients[0,2] * output_pixels_per_lenslet
            coeff[0,2] = self.forwardCoefficients[0,0]
            coeff[1,0] = self.forwardCoefficients[1,1] * output_pixels_per_lenslet
            coeff[1,1] = self.forwardCoefficients[1,2] * output_pixels_per_lenslet
            coeff[1,2] = self.forwardCoefficients[1,0]
            cv.WarpAffine(im, output_image, coeff,
                          flags=cv.CV_INTER_LINEAR+cv.CV_WARP_FILL_OUTLIERS+cv.CV_WARP_INVERSE_MAP)
        else:
            coeff = cv.CreateMat(2,3,cv.CV_32FC1)
            coeff[0,0] = self.reverseCoefficients[0,1] / output_pixels_per_lenslet
            coeff[0,1] = self.reverseCoefficients[0,2] / output_pixels_per_lenslet
            coeff[0,2] = self.reverseCoefficients[0,0]
            coeff[1,0] = self.reverseCoefficients[1,1] / output_pixels_per_lenslet
            coeff[1,1] = self.reverseCoefficients[1,2] / output_pixels_per_lenslet
            coeff[1,2] = self.reverseCoefficients[1,0]
            cv.WarpAffine(im, output_image, coeff,
                          flags=cv.CV_INTER_LINEAR+cv.CV_WARP_FILL_OUTLIERS+cv.CV_WARP_INVERSE_MAP)

        result = np.asarray(output_image)

        # Don't crop left of lenslets (0, y) or above (x, 0)
        leftbound = max(leftbound, 0)
        topbound = max(topbound,0)
    
        if cropToInside:
            return result[topbound * output_pixels_per_lenslet:bottombound * output_pixels_per_lenslet,
                          leftbound * output_pixels_per_lenslet:rightbound * output_pixels_per_lenslet]
        
        else:
            return result
            
# ------------------------------------------------------------------
#                          CUBIC WARP
# ------------------------------------------------------------------
#
# the model has 20 DOF
# x' = a + cx + ey + gx^2 + ixy + ky^2 + mx^3 + ox^2y + qxy^2 + sy^3
# y' = b + dx + fy + hx^2 + jxy + ly^2 + nx^3 + px^2y + rxy^2 + ty^3
class CubicWarp(object):

    def __init__(self):
        self.forwardCoefficients = np.zeros((2,10), dtype='float32')
        self.reverseCoefficients = np.zeros((2,10), dtype='float32')

    def eval(self, data, direction='r'):
        x = data[:,0]
        y = data[:,1]
        composite = np.vstack( (np.ones(x.shape), x, y, x*x, x*y, y*y, x*x*x, x*x*y, x*y*y, y*y*y ) )
        if (direction == 'R' or direction == 'r'):
            return np.dot(composite.transpose(), self.reverseCoefficients.transpose())
        elif (direction == 'F' or direction == 'f'):
            return np.dot(composite.transpose(), self.forwardCoefficients.transpose())
        else:
            raise ValueError("Unkown direction in CubicWarp.eval()")

    def eval_point(self, data, direction='r'):
        x = np.mat(data[0])
        y = np.mat(data[1])
        composite = np.vstack( (np.mat(np.ones(x.shape)), np.mat(x), np.mat(y), np.mat(x*x), np.mat(x*y), np.mat(y*y), np.mat(x*x*x), np.mat(x*x*y), np.mat(x*y*y), np.mat(y*y*y) ) )
        if (direction == 'R' or direction == 'r'):
            return np.array(np.dot(composite.transpose(), self.reverseCoefficients.transpose()))[0]
        elif (direction == 'F' or direction == 'f'):
            return np.array(np.dot(composite.transpose(), self.forwardCoefficients.transpose()))[0]
        else:
            raise ValueError("Unkown direction in CubicWarp.eval()")

    def fit(self, data):
        datalen = data.shape[0]
        A = np.zeros((datalen,10))
        bx = np.zeros(datalen)
        by = np.zeros(datalen)
        Ar = np.zeros((datalen,10) )
        bxr = np.zeros(datalen)
        byr = np.zeros(datalen)
      
        for i in range(datalen):
            px = float(data[i,0])
            py = float(data[i,1])
            lx = float(data[i,2])
            ly = float(data[i,3])

            A[i,0] = 1.0;          Ar[i,0] = 1.0;
            A[i,1] = px;           Ar[i,1] = lx;
            A[i,2] = py;           Ar[i,2] = ly;
            A[i,3] = px*px;        Ar[i,3] = lx*lx;
            A[i,4] = px*py;        Ar[i,4] = lx*ly;
            A[i,5] = py*py;        Ar[i,5] = ly*ly;
            A[i,6] = px*px*px;     Ar[i,6] = lx*lx*lx;
            A[i,7] = px*px*py;     Ar[i,7] = lx*lx*ly;
            A[i,8] = px*py*py;     Ar[i,8] = lx*ly*ly;
            A[i,9] = py*py*py;     Ar[i,9] = ly*ly*ly;

            bx[i] = lx;            bxr[i] = px
            by[i] = ly;            byr[i] = py

        self.forwardCoefficients = np.zeros((2,10), dtype='float32')
        self.reverseCoefficients = np.zeros((2,10), dtype='float32')


        import numpy.linalg as linalg

        self.forwardCoefficients[0,:] = linalg.lstsq(A,bx)[0]
        self.forwardCoefficients[1,:] = linalg.lstsq(A,by)[0]
        self.reverseCoefficients[0,:] = linalg.lstsq(Ar,bxr)[0]
        self.reverseCoefficients[1,:] = linalg.lstsq(Ar,byr)[0]

    def get_error(self, data, direction='f'):

        x = data[:,0]
        y = data[:,1]
        composite = np.vstack( (np.ones(x.shape),
                                x, y, x*x, x*y, y*y, x*x*x,
                                x*x*y, x*y*y, y*y*y ) )
        if (direction == 'R' or direction == 'r'):
            B_fit = np.dot(composite.transpose(), self.reverseCoefficients.transpose())
        elif (direction == 'F' or direction == 'f'):
            B_fit = np.dot(composite.transpose(), self.forwardCoefficients.transpose())
        else:
            raise ValueError("Unkown direction in CubicWarp.get_error()")

        B = data[:,2:4]
        err_per_point = np.sqrt(np.sum((B-B_fit)**2,axis=1)) # sum squared error per row
        return err_per_point

    def forwardCoefficients(self):
        return self.forwardCoefficients

    def reverseCoefficients(self):
        return self.reverseCoefficients

    def forward(self, coords):
        return np.dot(self.forwardCoefficients, np.array([1, coords[0], coords[1], coords[0]*coords[0], coords[0]*coords[1], coords[1]*coords[1], coords[0]*coords[0]*coords[0], coords[0]*coords[0]*coords[1], coords[0]*coords[1]*coords[1], coords[1]*coords[1]*coords[1]]))

    def reverse(self, coords):
        return np.dot(self.reverseCoefficients, np.array([1, coords[0], coords[1], coords[0]*coords[0], coords[0]*coords[1], coords[1]*coords[1], coords[0]*coords[0]*coords[0], coords[0]*coords[0]*coords[1], coords[0]*coords[1]*coords[1], coords[1]*coords[1]*coords[1]]))

    def save(self, filename):
        a = np.vstack((self.forwardCoefficients, self.reverseCoefficients))
        np.savetxt(filename, a, fmt="%12.6G") 
      
    def load(self, filename):
        a = np.loadtxt(filename)
        self.forwardCoefficients = a[0:2,:].astype(np.float32)
        self.reverseCoefficients = a[2:4,:].astype(np.float32)

    # OUTPUT AND DIGANOSTIC ROUTINES

    def __str__(self):
        s = "\t\tTranslation: ["
        tmp = self.reverseCoefficients[:,0];
        for t in tmp:
            s+= ('%0.2f, ' % t)
        s+= "]\n\t\tScaling and rotation: ["
        s+= str(self.reverseCoefficients[0,1]) + ' ' + str(self.reverseCoefficients[0,2]) + '; '
        s+= str(self.reverseCoefficients[1,1]) + ' ' + str(self.reverseCoefficients[1,2]) + ']'
        s+= '\n\t\tHigher Order terms (x): ['
        tmp = self.reverseCoefficients[0,3:];
        for t in tmp:
            s+= ('%0.2g, ' % t)
        s+= ']\n\t\tHigher Order terms (y): [' 
        tmp = self.reverseCoefficients[1,3:];
        for t in tmp:
            s+= ('%0.2g, ' % t)
        s+= ']'
        return s

    def check_fit(self, data):

        projected_camlenses = self.eval(data[:,0:2], 'f')
        diff = data[:,2:4] - projected_camlenses
        err = np.sqrt(np.power(diff[:,0],2) + np.power(diff[:,1],2))
        print ('\t--> Fit Statistics - Mean: %0.2f  Median: %0.2f  Min: %0.2f  Max: %0.2f' %
               (err.mean(), np.median(err, axis=0), err.min(), err.max() ))

    # Warp an image.
    #
    # This function expects im to be a numpy float32 array.  Returns 
    def warp_image(self, input_image, output_pixels_per_lenslet, direction="R",
                   cropToInside = False, lenslet_offset = None, output_size = None):
        im = cv.fromarray(input_image)

        ul = self.eval_point([0, 0], 'f')
        ur = self.eval_point([0, im.cols], 'f')
        ll = self.eval_point([im.rows, 0], 'f')
        lr = self.eval_point([im.rows, im.cols], 'f')

        leftbound = np.ceil(max(ul[1], ll[1]))
        rightbound = np.floor(min(ur[1], lr[1]))
        topbound = np.ceil(max(ul[0], ur[0]))
        bottombound = np.floor(min(ll[0], lr[0]))

        # Don't crop left of lenslets (0, y) or above (x, 0)
        leftbound = max(leftbound, 0)
        topbound = max(topbound,0)

        if output_size != None:
            putative_output_size = output_size
        else:
            nt = int(np.floor(bottombound - topbound))
            ns = int(np.floor(rightbound - leftbound))
            putative_output_size = (nt*output_pixels_per_lenslet, ns*output_pixels_per_lenslet)
        
        # Create the output image
        output_image = cv.CreateMat(putative_output_size[0], putative_output_size[1], im.type)

        # Apply the transform.
        scaled_shift = (0.0, 0.0)
        if (direction == 'f' or direction == 'F'):
            coeff = np.copy(self.forwardCoefficients)
            coeff[:,1:] *= output_pixels_per_lenslet
            coeff[:,3:] *= output_pixels_per_lenslet
            coeff[:,6:] *= output_pixels_per_lenslet
            if lenslet_offset != None:
                scaled_shift = (lenslet_offset[0] / output_pixels_per_lenslet,
                                lenslet_offset[1] / output_pixels_per_lenslet)
        else:
            coeff = np.copy(self.reverseCoefficients)
            coeff[:,1:] /= output_pixels_per_lenslet
            coeff[:,3:] /= output_pixels_per_lenslet
            coeff[:,6:] /= output_pixels_per_lenslet
            if lenslet_offset != None:
                scaled_shift = (lenslet_offset[0] * output_pixels_per_lenslet,
                                lenslet_offset[1] * output_pixels_per_lenslet)

        (x, y) = np.meshgrid(np.arange(putative_output_size[1], dtype='float32'), np.arange(putative_output_size[0], dtype='float32'))

        ix_array = (np.ones((putative_output_size[0], putative_output_size[1])).astype('float32') * coeff[0,0] + x * coeff[0,1] + y * coeff[0,2] + x * x * coeff[0,3] + x * y * coeff[0,4] + y * y * coeff[0,5] + x * x * x * coeff[0,6] + x * x * y * coeff[0,7] + x * y * y * coeff[0,8] + y * y * y * coeff[0,9] + scaled_shift[0])
        iy_array = (np.ones((putative_output_size[0], putative_output_size[1])).astype('float32') * coeff[1,0] + x * coeff[1,1] + y * coeff[1,2] + x * x * coeff[1,3] + x * y * coeff[1,4] + y * y * coeff[1,5] + x * x * x * coeff[1,6] + x * x * y * coeff[1,7] + x * y * y * coeff[1,8] + y * y * y * coeff[1,9] + scaled_shift[1])
        ix = cv.fromarray(ix_array.astype(np.float32))
        iy = cv.fromarray(iy_array.astype(np.float32))
        cv.Remap(im, output_image, ix, iy, flags=cv.CV_INTER_LINEAR+cv.CV_WARP_FILL_OUTLIERS+cv.CV_WARP_INVERSE_MAP)

        result = np.asarray(output_image)

        # Don't crop left of lenslets (0, y) or above (x, 0)
        leftbound = max(leftbound, 0)
        topbound = max(topbound,0)
    
        if cropToInside:
            return result[topbound * output_pixels_per_lenslet:bottombound * output_pixels_per_lenslet,
                          leftbound * output_pixels_per_lenslet:rightbound * output_pixels_per_lenslet]
        
        else:
            return result
