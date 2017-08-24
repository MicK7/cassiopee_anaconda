#!/usr/bin/env python
from distutils.core import setup, Extension
import os, sys

#=============================================================================
# Connector requires:
# C++ compiler
# Fortran compiler
# Numpy
# KCore
#=============================================================================

# Write setup.cfg file
import KCore.Dist as Dist
Dist.writeSetupCfg()

# Test if numpy exists =======================================================
(numpyVersion, numpyIncDir, numpyLibDir) = Dist.checkNumpy()

# Test if kcore exists =======================================================
(kcoreVersion, kcoreIncDir, kcoreLibDir) = Dist.checkKCore()

# Compilation des fortrans ===================================================
from KCore.config import *
prod = os.getenv("ELSAPROD")
if prod is None:
  variant_dir = "build"
else:
  variant_dir = "build/"+prod

# Setting libraryDirs and libraries ===========================================
libraryDirs = [variant_dir, kcoreLibDir]
libraries = ["connector", "kcore"]
(ok, libs, paths) = Dist.checkFortranLibs([], additionalLibPaths)
libraryDirs += paths; libraries += libs
(ok, libs, paths) = Dist.checkCppLibs([], additionalLibPaths)
libraryDirs += paths; libraries += libs

# setup =======================================================================
listExtensions = []
listExtensions.append(
    Extension('Connector.connector',
              sources=['Connector/connector.cpp'],
              include_dirs=["Connector"]+additionalIncludePaths+[numpyIncDir, kcoreIncDir],
              library_dirs=additionalLibPaths+libraryDirs,
              libraries=libraries+additionalLibs,
              extra_compile_args=Dist.getCppArgs(),
              extra_link_args=Dist.getLinkArgs()
              ) )

# setup ======================================================================
setup(
    name="Connector",
    version="2.4",
    description="Connector for *Cassiopee* modules.",
    author="Onera",
    package_dir={"":"."},
    packages=['Connector'],
    ext_modules=listExtensions)

# Check PYTHONPATH ===========================================================
Dist.checkPythonPath(); Dist.checkLdLibraryPath()
