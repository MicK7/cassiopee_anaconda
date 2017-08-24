# - diffArrays (array) -
import Converter as C   
import Generator as G

ni = 11; nj = 11; nk = 11
a = G.cart( (0,0,0), (1,1,1), (ni,nj,nk) )
a = C.initVars(a, "F", 1.)
a = C.initVars(a, "Q", 1.2)

b = G.cart( (0,0,0), (1,1,1), (ni,nj,nk) )
b = C.initVars(b, "Q", 2.)
b = C.initVars(b, "F", 3.)

ret = C.diffArrays([a], [b]); print ret
