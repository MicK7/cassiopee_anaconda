# - randomizeVar (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
a = G.cart((0,0,0),(1,1,1),(11,11,1))
a = C.initVars(a, 'F=10.')
b = C.randomizeVar(a, 'F', 0.1, 1.)
C.convertPyTree2File(b, "out.cgns")
