# - newZone (pyTree) -
import Converter.Internal as Internal

# Create a zone node
z = Internal.newGravity(value=[0.,0.,9.81]); print z

# Create a zone node and attach it to tree
t = Internal.newCGNSTree()
b = Internal.newCGNSBase('Base', 3, 3, parent=t)
z = Internal.newGravity(value=[0.,0.,9.81], parent=b)
print t
