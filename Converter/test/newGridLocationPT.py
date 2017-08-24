# - newGridLocation (pyTree) -
import Converter.Internal as Internal
import numpy

# Create a GridLocation node
n = Internal.newGridLocation(value='CellCenter'); print n

# Attach it to a parent node
d = Internal.newBC('wall', [1,80,30,30,1,2], 'BCWall') 
Internal.newGridLocation('Vertex', parent=d) 
print d
