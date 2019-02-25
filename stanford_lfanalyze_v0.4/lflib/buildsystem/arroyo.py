from SCons.Script import *

import os

def setup_arroyo_environment(conf):

    if "win32" == os.sys.platform:
        print 'WARNING: Arroyo is not available for Windows.  Physical optics models will not be available.'

    else: # Linux and MacOS
        conf.env.Append( CPPPATH=[ '/usr/local/include/arroyo' ] )
        conf.env.Append( LIBPATH=[ '/usr/local/lib' ] )  
        conf.env.Append(LIBS=["arroyo", "cfitsio", "fftw3", "fftw3f", "lapack"])

        if not conf.CheckLibWithHeader('arroyo', 'arroyo.h', 'c++'):
            print 'ERROR: Could not locate arroyo library. Physical optics models will not be available.'
            

