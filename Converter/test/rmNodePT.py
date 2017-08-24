# - _rmNode (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Converter.Internal as Internal

t = C.newPyTree(['Base', 'Base2'])
a = G.cart((0,0,0), (1,1,1), (10,10,10))
t[2][1][2] += [a]

# rm the node a from t
Internal._rmNode(t, a)

Internal.printTree(t)
