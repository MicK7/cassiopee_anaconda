#!/usr/bin/env python
from distutils.core import setup, Extension
import os, sys

#=============================================================================
# Post requires:
# C++ compiler
# Fortran compiler
# Numpy
# KCore
#=============================================================================

# Write setup.cfg
import KCore.Dist as Dist
Dist.writeSetupCfg()

# Test if numpy exists =======================================================
(numpyVersion, numpyIncDir, numpyLibDir) = Dist.checkNumpy()

# Test if kcore exists =======================================================
(kcoreVersion, kcoreIncDir, kcoreLibDir) = Dist.checkKCore()

from KCore.config import *
prod = os.getenv("ELSAPROD")
if prod is None:
  variant_dir = "build"
else:
  variant_dir = "build/"+prod


# Setting libraryDirs and libraries ===========================================
libraryDirs = [variant_dir, kcoreLibDir]
libraries = ["post", "kcore"]
(ok, libs, paths) = Dist.checkFortranLibs([], additionalLibPaths)
libraryDirs += paths; libraries += libs
(ok, libs, paths) = Dist.checkCppLibs([], additionalLibPaths)
libraryDirs += paths; libraries += libs

import srcs

# extensions =================================================================
listExtensions = []
listExtensions.append(
    Extension('Post.post',
              sources=["Post/post.cpp"],
              include_dirs=["Post"]+additionalIncludePaths+[numpyIncDir,kcoreIncDir], 
              library_dirs=additionalLibPaths+libraryDirs,
              libraries=libraries+additionalLibs,
              extra_compile_args=Dist.getCppArgs(),
              extra_link_args=Dist.getLinkArgs()
              ) )

# setup ======================================================================
setup(
    name="Post",
    version="2.4",
    description="Post-processing of CFD solutions.",
    author="Onera",
    package_dir={"":"."},
    packages=['Post'],
    ext_modules=listExtensions
    )

# Check PYTHONPATH ===========================================================
Dist.checkPythonPath(); Dist.checkLdLibraryPath()
