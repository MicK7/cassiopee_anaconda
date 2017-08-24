# - convertArray2Tetra (array) -
import Converter as C
import Generator as G

# 2D : triangles 
a = G.cart((0.,0.,0.), (0.1,0.1,0.2), (10,10,1))
b = C.convertArray2Tetra(a)
C.convertArrays2File([b], 'new1.plt', 'bin_tp')

# 3D : tetrahedras
a = G.cart((0.,0.,0.), (0.1,0.1,0.2), (10,10,10))
b = C.convertArray2Tetra(a)
C.convertArrays2File([b], 'new2.plt', 'bin_tp')
