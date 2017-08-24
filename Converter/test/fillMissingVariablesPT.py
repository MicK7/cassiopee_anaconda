# - fillMissingVariables (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import KCore.test as test

a = G.cart((0,0,0), (1,1,1), (10,10,11)); a[0] = 'cart1'
b = G.cart((1,0,0), (2,1,1), (10,10,11)); b[0] = 'cart2'
a = C.addVars(a, 'rou'); a = C.addVars(a, 'rov')
a = C.addVars(a, 'centers:cellN')
a = C.addVars(a, ['Density', 'Hx', 'centers:Hy'])

t = C.newPyTree(['Base',a,b])
t = C.fillMissingVariables(t)
C.convertPyTree2File(t, 'out.cgns')
