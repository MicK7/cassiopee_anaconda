# - node2Center (array) -
import Converter as C
import Generator as G

def F(x,y): return 2*x+y

ni = 30; nj = 40; nk = 1
a = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,nk))
a = C.initVars(a, 'ro', F, ['x','y'])
ac = C.node2Center(a)
C.convertArrays2File([a,ac], "out.plt")
