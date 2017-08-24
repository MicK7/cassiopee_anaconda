# - cpVars (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
a = G.cart((0,0,0),(1,1,1),(10,10,10)); a[0] = 'cart1'
b = G.cart((0,0,0),(1,1,1),(10,10,10)); b[0] = 'cart2'
a = C.initVars(a, 'Density',2.)
a = C.cpVars(a, 'Density', a, 'Density2')
c = C.cpVars(a, 'Density', b, 'Density')
t = C.newPyTree(['Base',a,c])
C.convertPyTree2File(t, 'out.cgns')
