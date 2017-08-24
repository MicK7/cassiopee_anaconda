# - convertPyTree2File (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G

a = G.cart((0.,0.,0.),(0.1,0.1,0.1),(11,11,11))
t = C.newPyTree(['Base']); t[2][1][2] += [a]
C.convertPyTree2File(t, 'out.cgns')
C.convertPyTree2File(t, 'out.plt')
