# Post installation test
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
    import Post
    print "Post correctly installed."
except Exception, inst:
    print "FAILED:",inst
    print "FAILED: Post badly installed."
