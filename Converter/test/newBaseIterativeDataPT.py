# - newZone (pyTree) -
import Converter.Internal as Internal

# Create a zone node
z = Internal.newBaseIterativeData(name='BaseIterativeData', nsteps=100, itype='IterationValues'); print z

# Create a zone node and attach it to tree
t = Internal.newCGNSTree()
b = Internal.newCGNSBase('Base', 3, 3, parent=t)
z = Internal.newBaseIterativeData(name='BaseIterativeData', nsteps=100, itype='IterationValues', parent=b)
print t
