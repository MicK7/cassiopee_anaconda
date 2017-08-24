# Dist2Walls installation test
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
    import Dist2Walls
    print "Dist2Walls correctly installed."
except Exception, inst:
    print "FAILED:",inst
    print "FAILED: Dist2Walls badly installed."
    
