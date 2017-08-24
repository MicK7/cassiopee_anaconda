# - convertArray2Tetra (array) -
import Converter as C
import Generator as G

# 2D: quads -> triangles 
a = G.cartHexa((0.,0.,0.), (0.1,0.1,0.2), (10,10,1))
b = C.convertArray2Tetra(a)
C.convertArrays2File([b], 'new1.plt', 'bin_tp')

# 3D: hexa -> tetrahedra
a = G.cartHexa((0.,0.,0.), (0.1,0.1,0.2), (10,10,10))
b = C.convertArray2Tetra(a)
C.convertArrays2File([b], 'new2.plt', 'bin_tp')