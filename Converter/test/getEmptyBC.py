# - getEmptyBC (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C

a1 = G.cart((0.,0.,0.), (0.1, 0.1, 0.1), (11, 21, 2)); a1[0] = 'cart1'
a1 = C.addBC2Zone(a1, 'wall1', 'BCWall', 'imin')
a2 = G.cart((1., 0.2, 0.), (0.1, 0.1, 0.1), (11, 21, 2)); a2[0] = 'cart2'
a2 = C.addBC2Zone(a2, 'wall1', 'BCWall', 'imax')
t = C.newPyTree(['Base',a1,a2])
# Returns undefined windows (as range)
wins = C.getEmptyBC(t,2); print wins
# Returns undefined windows (as face list)
t= C.convertArray2NGon(t)
faceList = C.getEmptyBC(t,2); print faceList

C.convertPyTree2File(t, 'out.cgns')
