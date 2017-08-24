# - normL0 (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G

a = G.cart((0,0,0), (1,1,1), (11,11,11))
a = C.initVars(a, "F", 1.)
print 'normL0 = ', C.normL0(a, "F")
