# - getZoneDim (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Converter.Internal as Internal

a = G.cart((0,0,0), (1,1,1), (10,10,10))
dim = Internal.getZoneDim(a); print dim

a = G.cartTetra((0,0,0), (1,1,1), (10,10,10))
dim = Internal.getZoneDim(a); print dim
