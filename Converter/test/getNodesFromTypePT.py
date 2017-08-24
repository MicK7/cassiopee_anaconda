# - getNodesFromType (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Converter.Internal as Internal

a = G.cart( (0,0,0), (1,1,1), (10,10,10) )
t = C.newPyTree(['Base',a])

# Return nodes of type 'Zone_t'
zones = Internal.getNodesFromType(t, 'Zone_t'); print zones

# Limit search to first level children (faster)
zones = Internal.getNodesFromType1(t, 'Zone_t'); print zones
