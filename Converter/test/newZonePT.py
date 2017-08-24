# - newZone (pyTree) -
import Converter.Internal as Internal

# Create a zone node
z = Internal.newZone('Zone', zsize=[[10],[2],[0]], ztype='Structured'); print z

# Create a zone node and attach it to tree
t = Internal.newCGNSTree()
b = Internal.newCGNSBase('Base', 3, 3, parent=t)
z = Internal.newZone('Zone', [[10],[2],[0]], 'Structured', parent=b)
print t
