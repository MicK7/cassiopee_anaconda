# - groupBCByType (PyTree) -
import Converter.PyTree as C
import Converter.Internal as Internal
import Generator.PyTree as G

a = G.cart((0.,0.,0.),(0.1,0.1,0.1),(10,10,10))

a = C.fillEmptyBCWith(a,'wall','BCWall')

t = C.newPyTree(['Base',a])
Internal._groupBCByBCType(t,'BCWall','FamWall')

C.convertPyTree2File(t,'out.cgns')
