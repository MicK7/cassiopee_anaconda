# - deform (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Geom.PyTree as D
import Transform.PyTree as T

a = G.cart((0.,0.,0.),(1.,1.,1.),(10,10,10))
vect = ['hx','hy','hz']; a = C.addVars(a, vect) 
a = C.initVars(a, 'hx', 10.)
b = T.deform(a, vect)
C.convertPyTree2File(a, 'out.cgns')

