# - newBCProperty (pyTree) -
import Converter.Internal as Internal

# Create a BC property node
n = Internal.newBCProperty(wallFunction='Null', area='Null'); print n

# Attach it to a parent node
d = Internal.newBC(name='BC', pointList=[22,1036,101,43], btype='BCWall')
Internal.newBCProperty(wallFunction='Null', area='Null', parent=d)
print d
