# - getType (pyTree) -
import Converter.Internal as Internal
import Generator.PyTree as G

a = G.cart( (0,0,0), (1,1,1), (10,10,10) )
print Internal.getType(a)
