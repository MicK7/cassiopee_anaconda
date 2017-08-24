# Test KCore installation
import sys, os

# Rm . from PYTHONPATH
try:
  del sys.path[sys.path.index('')]
except:
  pass
try:
  del sys.path[sys.path.index(os.getcwd())]
except:
  pass

try:
    import KCore
    print "KCore correctly installed."
except Exception, inst:
    print "FAILED:",inst
    print "FAILED: KCore badly installed."
