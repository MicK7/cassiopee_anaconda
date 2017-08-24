# - newGridConnectivity1to1 (pyTree) -
import Converter.Internal as Internal

# Create a node
n = Internal.newGridConnectivityProperty(); print n

# Attach it to a parent node
d = Internal.newGridConnectivity(name='Match', donorName='blk1', ctype='Abutting1to1')
Internal.newGridConnectivityProperty(parent=d) 
print d
