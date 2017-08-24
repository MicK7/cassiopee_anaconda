# - newPointRange (pyTree) -
import Converter.Internal as Internal
import numpy

# Create a node
n = Internal.newPointRange(name='PointRange', value=[22,56]); print n

# Attach it to a parent node
d = Internal.newBC('wall', [1,80,30,30,1,2], 'BCWall') 
Internal.newPointRange('PointRange', [1,71,29,29,1,2], parent=d)  
print d
