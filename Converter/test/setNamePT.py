# - setName (pyTree) -
import Converter.PyTree as C
import Converter.Internal as Internal

node = ['node1', 1., [], 'DataArray_t']
Internal.setName(node, 'myNode'); print node
