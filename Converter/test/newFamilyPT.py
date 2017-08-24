# - newZone (pyTree) -
import Converter.Internal as Internal

# Create a zone node
z = Internal.newFamily(name='FamWall'); print z

# Create a zone node and attach it to tree
t = Internal.newCGNSTree()
b = Internal.newCGNSBase('Base', 3, 3, parent=t)
z = Internal.newFamily(name='FamWall', parent=b)
print t
