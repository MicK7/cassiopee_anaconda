# - convertArray2Tetra (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G

t = C.newPyTree(['Base1',2,'Base2',3])

# 2D : triangles 
a = G.cart((0.,0.,0.), (0.1,0.1,0.2), (10,10,1))
b = C.convertArray2Tetra(a); t[2][1][2].append(b)

# 3D : tetrahedras
a = G.cart((0.,0.,0.), (0.1,0.1,0.2), (10,10,10))
b = C.convertArray2Tetra(a); t[2][2][2].append(b)
C.convertPyTree2File(t, 'out.cgns')
