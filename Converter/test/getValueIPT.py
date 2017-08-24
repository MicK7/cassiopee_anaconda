# - getValue of a node (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Converter.Internal as Internal

# Structured array
Ni = 40; Nj = 50; Nk = 20
a = G.cart((0,0,0), (1., 0.5,1.), (Ni,Nj,Nk))

# Get value stored in a zone node
print Internal.getValue(a)

# Get type of a zone (from ZoneType node)
node = Internal.getNodeFromName(a, 'ZoneType')
print Internal.getValue(node)
