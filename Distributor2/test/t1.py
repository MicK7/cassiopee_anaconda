# Distributor2 installation test
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
    import Distributor2
    print "Distributor2 correctly installed."
except Exception, inst:
    print "FAILED:",inst
    print "FAILED: Distributor2 badly installed."
    
