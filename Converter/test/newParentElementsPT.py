# - newParentElements (pyTree) -
import Converter.Internal as Internal

# Create a parent elements node
b = Internal.newParentElements(value=None); print b

# Attach it to elements
d = Internal.newElements('Elements', 'UserDefined', None, 0)
Internal.newParentElements(None, parent=d); print d
