# - printTree (pyTree) -
import Converter.PyTree as C
import Converter.Internal as Internal
import Generator.PyTree as G

a = G.cart((0,0,0), (1,1,1), (10,10,10))
Internal.printTree(a) # can print a zone (or any node)

t = C.newPyTree(['Base',a])
Internal.printTree(t) # can print a tree
