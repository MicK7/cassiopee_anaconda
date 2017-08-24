# - newCGNSBase (pyTree) -
import Converter.Internal as Internal

# Create a base node
b = Internal.newCGNSBase('Base'); print b

# Create a base node and attach it to tree
t = Internal.newCGNSTree()
Internal.newCGNSBase('Base', 3, 3, parent=t); print t
