# - newParentElementsPosition (pyTree) -
import Converter.Internal as Internal

# Create a parent elements position node
b = Internal.newParentElementsPosition(value=None); print b

# Attach it to elements
d = Internal.newElements(name='Elements', etype='UserDefined', econnectivity=None, eboundary=0)
Internal.newParentElementsPosition(None, parent=d); print d
