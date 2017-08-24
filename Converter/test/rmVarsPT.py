# - rmVars (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
a = G.cart((0,0,0),(1,1,1),(10,10,10))
b = C.addVars(a, 'Density')
b = C.rmVars(b, 'Density')
C.convertPyTree2File(b, 'out.cgns')
