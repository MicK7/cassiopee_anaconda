# - convertArray2Hexa (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G

t = C.newPyTree(['Base',3,'Base2',2])

# 2D: quad
a = G.cart((0.,0.,0.), (0.1,0.1,0.2), (10,10,1))
b = C.convertArray2Hexa(a); t[2][2][2].append(b)

# 3D: hexa
a = G.cart((0.,0.,0.), (0.1,0.1,0.2), (10,10,10))
b = C.convertArray2Hexa(a); t[2][1][2].append(b)
C.convertPyTree2File(t, 'out.cgns')
