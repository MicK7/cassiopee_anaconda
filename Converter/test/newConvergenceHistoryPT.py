# - newZone (pyTree) -
import Converter.Internal as Internal

# Create a zone node
z = Internal.newConvergenceHistory(name='ZoneConvergenceHistory', value=100); print z

# Create a zone node and attach it to tree
t = Internal.newCGNSTree()
b = Internal.newCGNSBase('Base', 3, 3, parent=t)
z = Internal.newConvergenceHistory(name='GlobalConvergenceHistory', value=100, parent=b)
print t
