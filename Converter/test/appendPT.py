# - append (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Converter.Internal as Internal

t = C.newPyTree(['Base', 'Base2'])
a = G.cart((0,0,0), (1,1,1), (10,10,10))

# Append a to 'Base' node of t
t = Internal.append(t, a, 'Base')

C.convertPyTree2File(t, 'out.cgns')
