# - createElsaHybrid -
import Converter.PyTree as C
import Generator.PyTree as G
import Converter.Internal as Internal
import KCore.test as test

# Cas NGon
a = G.cartNGon((0,0,0), (1,1,1), (10,10,10))
a = C.fillEmptyBCWith(a, 'farfield', 'BCFarfield')
Internal._createElsaHybrid(a)
test.testT(a, 1)
