# - createElsaHybrid -
import Converter.PyTree as C
import Generator.PyTree as G
import Converter.Internal as Internal

# Cas NGon
a = G.cartNGon((0,0,0), (1,1,1), (10,10,10))
a = C.fillEmptyBCWith(a, 'farfield', 'BCFarfield')
Internal._createElsaHybrid(a)
C.convertPyTree2File(a, 'out.cgns')
