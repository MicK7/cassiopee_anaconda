# - createChild (pyTree) -
import Converter.PyTree as C
import Converter.Internal as Internal

node = Internal.createNode('myNode', 'DataArray_t', value=1., children=[])
Internal.createChild(node, 'childName', 'DataArray_t', value=2., children=[])
print node 
