# - newGridConnectivity1to1 (pyTree) -
import Converter.Internal as Internal

# Create a node
n = Internal.newGridConnectivityType(ctype='Abutting1to1'); print n

# Attach it to a parent node
d = Internal.newGridConnectivity(name='Match', donorName='blk1', ctype=None)
Internal.newGridConnectivityType(ctype='Abutting1to1', parent=d) 
print d
