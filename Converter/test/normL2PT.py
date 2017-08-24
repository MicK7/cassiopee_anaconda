# - normL2 (pyTree) -
import Converter.PyTree as C    
import Generator.PyTree as G

ni = 11; nj = 11; nk = 11
a = G.cart((0,0,0), (1,1,1), (ni,nj,nk))

# Add variable F
a = C.addVars(a, "F")
a = C.initVars(a, "F", 1.)
print 'normL2 = ', C.normL2(a, "F")
