# - newGridConnectivity1to1 (pyTree) -
import Converter.Internal as Internal

# Create a node
n = Internal.newOversetHoles(name='OversetHoles', pointRange=None, pointList=None); print n

# Attach it to a parent node
d = Internal.newZoneGridConnectivity(name='ZoneGridConnectivity')
Internal.newOversetHoles(name='OversetHoles', pointList=[22,1036,101,43], parent=d) 
print d
