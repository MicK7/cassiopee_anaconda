# - appendBaseName2ZoneName (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Converter.Internal as Internal

a = G.cart((0,0,0), (1,1,1), (10,10,10)); a[0] = 'a'
b = G.cart((11,0,0), (1,1,1), (10,10,10)); b[0] = 'b'
t = C.newPyTree(['Base',a,b])

t = Internal.appendBaseName2ZoneName(t)
C.convertPyTree2File(t, 'out.cgns')
