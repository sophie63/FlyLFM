from SCons.Script import *

import os

def find_opencv():
    opencvDir = "C:/opencv/build"
    if os.path.exists( opencvDir ):
        return opencvDir
    
    for x in ('', ' (x86)'):
        for v in ('2.1','2.2', '2.3'):
            opencvDir = "C:/Program Files" + x + "/OpenCV" + v
            if os.path.exists( opencvDir ):
                return opencvDir
    return None

def setup_opencv_environment(conf, opencv_libraries):

    if "win32" == os.sys.platform:
        opencv_dir = find_opencv()
        if opencv_dir is None:
            print( "ERROR: can't find opencv" )
            Exit(1)

        print( "[OpenCV] : '" + opencv_dir + "'" )
        conf.env.Append( CPPPATH=[ opencv_dir + '/include' ] )
        conf.env.Append( LIBPATH=[ opencv_dir + '/x64/vc10/lib' ] )

        # Check for headers
        if not conf.CheckCXXHeader( "opencv/cv.h" ):
            print( "Can't find OpenCV headers" )
            Exit(1)

        # Check for libraries
        for lib in opencv_libraries:
            fullname = 'opencv_' + lib + '231'
            if not conf.CheckLib(fullname):
                print 'ERROR: Could not locate ', fullname
                sys.exit(1)

    else: # Linux and MacOS
        for lib in opencv_libraries:
            fullname = 'opencv_' + lib
            if not conf.CheckLibWithHeader(fullname, 'opencv/cv.h', 'c++'):
                print 'ERROR: Could not locate', fullname, 'library'
                sys.exit(1)

