# - center2Node (array) -
import Converter as C
import Generator as G

ni = 30; nj = 40; nk = 10
a = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,nk))
a = C.initVars(a, 'ro', 1.)
an = C.center2Node(a)
C.convertArrays2File([an], "out.plt")
