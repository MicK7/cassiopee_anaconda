#!/usr/bin/env python
import os, sys
from distutils.core import setup, Extension

#=============================================================================
# Dist2Walls requires:
# C++ compiler
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

# Setting libraryDirs and libraries ===========================================
from KCore.config import *
prod = os.getenv("ELSAPROD")
if prod is None:
  variant_dir = "build"
else:
  variant_dir = "build/"+prod

libraryDirs = [kcoreLibDir, variant_dir]
libraries = ["dist2walls", "kcore"]
from KCore.config import *
(ok, libs, paths) = Dist.checkCppLibs([], additionalLibPaths)
libraryDirs += paths; libraries += libs

# Extensions =================================================================
extensions = [
    Extension('Dist2Walls.dist2walls',
              sources=["Dist2Walls/dist2walls.cpp"],
              include_dirs=["Dist2Walls"]+additionalIncludePaths+[numpyIncDir, kcoreIncDir],
              library_dirs=additionalLibPaths+libraryDirs,
              libraries=libraries+additionalLibs,
              extra_compile_args=Dist.getCppArgs(),
              extra_link_args=Dist.getLinkArgs()
	)
    ]

# Setup ======================================================================
setup(
    name="Dist2Walls",
    version="2.4",
    description="Computation of distance to walls.",
    author="Onera",
    package_dir={"":"."},
    packages=['Dist2Walls'],
    ext_modules=extensions
    )

# Check PYTHONPATH ===========================================================
Dist.checkPythonPath(); Dist.checkLdLibraryPath()
