# - addBase2PyTree (pyTree) -
import Converter.PyTree as C

t = C.newPyTree(['Base', 3]) # must contain volume zones
t = C.addBase2PyTree(t, 'Base2', 2) # must contain surface zones
print t
