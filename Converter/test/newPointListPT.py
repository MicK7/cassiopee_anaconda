# - newPointList (pyTree) -
import Converter.Internal as Internal
import numpy

# Create a Descriptor node
n = Internal.newPointList(name='PointList', value=[101,51,22,1036,2]); print n

# Attach it to a parent node
d = Internal.newBC('wall', [1,80,30,30,1,2], 'BCWall') 
Internal.newPointList('PointList', [101,51,22,1036,2], parent=d)  
print d
