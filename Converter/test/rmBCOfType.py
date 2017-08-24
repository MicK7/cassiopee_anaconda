# - rmBCOfType (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G

a = G.cylinder((0,0,0), 1., 1.5, 0., 360., 1., (80,30,2))
b = G.cart((-0.1,0.9,0), (0.01,0.01,1.), (20,20,2))

a = C.addBC2Zone(a, 'wall1', 'BCWall', 'jmin')
a = C.addBC2Zone(a, 'match1', 'BCMatch', 'imin', a, 'imax', [1,2,3]) 
a = C.addBC2Zone(a, 'match2', 'BCMatch', 'imax', a, 'imin', [1,2,3]) 
a = C.addBC2Zone(a, 'overlap1', 'BCOverlap', 'jmax')
b = C.addBC2Zone(b, 'wall2', 'BCWall', 'imin')
b = C.addBC2Zone(b, 'loin', 'FamilySpecified:LOIN', 'imax')

t = C.newPyTree(['Base']); t[2][1][2] += [a,b]
t[2][1] = C.addFamily2Base(t[2][1], 'LOIN', bndType='BCFarfield')

t = C.rmBCOfType(t, 'BCWall')
t = C.rmBCOfType(t, 'BCMatch')
t = C.rmBCOfType(t, 'BCFarfield')
t = C.rmBCOfName(t, 'FamilySpecified:LOIN')
C.convertPyTree2File(t, 'out.cgns')
