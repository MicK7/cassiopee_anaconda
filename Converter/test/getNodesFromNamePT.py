# - getNodesFromName (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Converter.Internal as Internal

a = G.cart( (0,0,0), (1,1,1), (10,10,10) )
t = C.newPyTree(['Base',a])

# Return nodes named 'cart'
nodes = Internal.getNodesFromName(t, 'cart'); print nodes

# Return the 3 coordinate nodes
nodes = Internal.getNodesFromName(t, 'Coordinate*'); print nodes
