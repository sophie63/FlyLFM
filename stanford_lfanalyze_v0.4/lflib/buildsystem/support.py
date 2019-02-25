# ------------------------------------------------------------------------------
# This file contains helper functions to augment the scons build system.
#
# Portions of this file were adapted from the Field3D build system.
# http://field3d.googlecode.com/svn-history/r12/trunk/BuildSupport.py
# ------------------------------------------------------------------------------

import os
import sys

from SCons.Script import *

from os.path import join

# ------------------------------------------------------------------------------
# Custom Library and Program Builders
#
# At the moment, these routines add proper use of the
# 'install_name_tool' on OSX, but other custom build behavior could be
# added here too.
# ------------------------------------------------------------------------------


def Library (env, target=None, source=None, **kw):
    if source is None:
        source = target
        target = None


    if (sys.platform == 'win32'):
        l   = env.SharedLibrary (target, source, **kw)
        env.Append ( BUILT_SHLIBS = [[ l[0].abspath ]] )
        return l
    
    l   = env.SharedLibrary (target, source, **kw)
    l_i = env.Install (env['INSTALL_LIB_DIR'], l, **kw)
    env.Append ( BUILT_SHLIBS = [[ l[0].abspath, l_i[0].abspath ]] )

    if sys.platform == 'darwin':
      env.AddPostAction (l, "install_name_tool -id %s %s" % (l[0].abspath, l[0].path))
      env.AddPostAction (l_i, "install_name_tool -id %s %s" % (l_i[0].abspath, l_i[0].path))

    return l

def Program (env, target=None, source=None, **kw):
    if source is None:
        source = target
        target = None

    p   = env.Program(target,source,**kw)

    if (sys.platform == 'win32'):
        return p

    p_i = env.Install(env['INSTALL_BIN_DIR'], p,**kw)
    if sys.platform == 'darwin':
      try:
        for l in env['BUILT_SHLIBS']:
          env.AddPostAction (p_i, "install_name_tool -change %s %s %s" % ( l[0], l[1], p_i[0].path ))
      except KeyError:
        pass

    return p

# ------------------------------------------------------------------------------
# Strings
# ------------------------------------------------------------------------------

arch32         = "m32"
arch64         = "m64"

# ------------------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------------------

systemIncludePaths = {
    "darwin" : { arch32 : ["/usr/local/include",
                           "/opt/local/include"],
                 arch64 : ["/usr/local/include",
                           "/opt/local/include"]},
    "linux2" : { arch32 : ["/usr/local/include"],
                 arch64 : ["/usr/local64/include"]}
}

systemLibPaths = {
    "darwin" : { arch32 : ["/usr/local/lib",
                           "/opt/local/lib"],
                 arch64 : ["/usr/local/lib",
                           "/opt/local/lib"]},
    "linux2" : { arch32 : ["/usr/local/lib"],
                 arch64 : ["/usr/local64/lib"]}
}

systemLibs = {
    "darwin" : [],
    "linux2" : ["dl"]
    }

# ------------------------------------------------------------------------------

def numCPUs():
    if os.sysconf_names.has_key("SC_NPROCESSORS_ONLN"):
        nCPUs = os.sysconf("SC_NPROCESSORS_ONLN")
        if isinstance(nCPUs, int) and nCPUs > 0:
            return nCPUs
    else: 
        return int(os.popen2("sysctl -n hw.ncpu")[1].read())
    if os.environ.has_key("NUMBER_OF_PROCESSORS"):
         nCPUs = int(os.environ["NUMBER_OF_PROCESSORS"]);
         if nCPUs > 0:
             return nCPUs
    return 1
