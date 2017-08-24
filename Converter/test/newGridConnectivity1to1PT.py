# - newGridConnectivity1to1 (pyTree) -
import Converter.Internal as Internal

# Create a node
n = Internal.newGridConnectivity1to1(name='Match', donorName='blk1', pointRange=[1,1,1,33,1,69], pointRangeDonor=[51,51,1,33,1,69], transform=None); print n

# Attach it to a parent node
d = Internal.newZoneGridConnectivity(name='ZoneGridConnectivity')
Internal.newGridConnectivity1to1(name='Match', donorName='blk1', pointRange=[1,1,1,33,1,69], pointRangeDonor=[51,51,1,33,1,69], transform=None, parent=d) 
print d
