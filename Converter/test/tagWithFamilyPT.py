# - tagWithFamily (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G

a = G.cylinder((0,0,0), 1., 1.5, 0., 360., 1., (80,30,2))
b = G.cart((-0.1,0.9,0), (0.01,0.01,1.), (20,20,2))
a = C.tagWithFamily(a, 'CYLINDER')
b = C.tagWithFamily(b, 'CART')

t = C.newPyTree(['Base',a,b])
t[2][1] = C.addFamily2Base(t[2][1], 'CYLINDER')
t[2][1] = C.addFamily2Base(t[2][1], 'CART')

C.convertPyTree2File(t, 'out.cgns')
