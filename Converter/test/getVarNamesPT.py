# - getVarNames (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
a = G.cart((0,0,0),(1,1,1),(10,10,10))
a = C.addVars(a, ['Density', 'centers:cellN'])
print C.getVarNames(a, loc='nodes')
print C.getVarNames(a, loc='centers')
print C.getVarNames(a, excludeXYZ=True, loc='both')
