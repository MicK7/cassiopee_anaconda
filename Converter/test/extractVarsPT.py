# - extractVars (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C

a = G.cart( (0,0,0), (1,1,1), (10,10,10) )
a = C.addVars(a, ['F', 'G', 'centers:H'])
# Keep only F
a = C.extractVars(a, ['F'])
t = C.newPyTree(['Base',a])
C.convertPyTree2File(t, "out.cgns")
