# - extractBCOfType (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G

a = G.cylinder((0,0,0), 1., 1.5, 360., 0., 1., (100,30,10))
a = C.addBC2Zone(a, 'wall1', 'BCWall', 'jmin')
Z = C.extractBCOfType(a, 'BCWall')
t = C.newPyTree(['Base',3,'Skin',2]); t[2][1][2] += [a]; t[2][2][2] += Z
C.convertPyTree2File(t, 'out.cgns')

