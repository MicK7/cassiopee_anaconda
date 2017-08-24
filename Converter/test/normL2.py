# - normL2 (array) -
import Converter as C    
import Generator as G

ni = 11; nj = 11; nk = 11
a = G.cart( (0,0,0), (1,1,1), (ni,nj,nk) )
a = C.initVars(a, "F", 1.)
print 'normL2 = ', C.normL2(a, "F")

# cellN variable IS taken into account
cellnf = C.array('celln', ni, nj, nk)
cellnf = C.initVars(cellnf, "celln", 1.)

cellnf[1][0][1] = 0.
cellnf[1][0][2] = 0.

a = C.addVars([a, cellnf])
print 'normL2 = ', C.normL2(a, "F")
