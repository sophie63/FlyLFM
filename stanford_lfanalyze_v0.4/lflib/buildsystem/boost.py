from SCons.Script import *

import os

def find_boost():
    for x in ('', ' (x86)'):	
        boostDir = "C:/Program Files" + x + "/boost/latest"
        if os.path.exists( boostDir ):
            return boostDir
        for bv in reversed( range(33,50) ):
            for extra in ('', '_0', '_1'):
                boostDir = "C:/Program Files" + x + "/boost/boost_1_" + str(bv) + extra
                if os.path.exists( boostDir ):
                    return boostDir
    if os.path.exists( "C:/boost" ):
        return "C:/boost"
    if os.path.exists( "/boost" ):
        return "/boost"
    print( "can't find boost" )
    Exit(1)

def find_boost_version():
    for x in ('', ' (x86)'):	
        boostDir = "C:/Program Files" + x + "/boost/latest"
        if os.path.exists( boostDir ):
            return boostDir
        for bv in reversed( range(33,50) ):
            for extra in ('', '_0', '_1'):
                boostDir = "C:/Program Files" + x + "/boost/boost_1_" + str(bv) + extra
                if os.path.exists( boostDir ):
                    return str(bv) + extra
    print( "can't find boost" )
    Exit(1)

def find_zlib():
    for x in ('', ' (x86)'):	
        zlibDir = "C:/Program Files" + x + "/GnuWin32"
        if os.path.exists( zlibDir ):
            return zlibDir

    print( "can't find zlib" )
    Exit(1)

def setup_boost_environment(conf, boost_libraries):

    if "win32" == os.sys.platform:
        boost_dir = find_boost()
        if boost_dir is None:
            print( "ERROR: can't find boost" )
            Exit(1)

        print( "[Boost] : '" + boost_dir + "'" )
        conf.env.Append( CPPPATH=[ boost_dir ] )
        conf.env.Append( LIBPATH=[ boost_dir + os.path.sep + 'lib64' ] )  # download 64 bit pre-build dll's here: http://boost.teeks99.com/
        conf.env.Append( CPPDEFINES=['BOOST_ALL_DYN_LINK'] )

        # Check for headers
        if not conf.CheckCXXHeader( "boost/shared_ptr.hpp" ):
            print( "Can't find boost headers" )
            Exit(1)

        # Check for libraries
        boost_version = find_boost_version()
        for lib in boost_libraries:
            fullname = 'boost_' + lib + '-vc100-mt-1_' + boost_version
            if not conf.CheckLib(fullname):
                print 'ERROR: Could not locate ', fullname
                sys.exit(1)

    else: # Linux and MacOS
        for lib in boost_libraries:
            fullname = 'boost_' + lib + '-mt'
            if not conf.CheckLibWithHeader(fullname, 'boost/shared_ptr.hpp', 'c++'):
                print 'ERROR: Could not locate', fullname, 'library'
                sys.exit(1)

