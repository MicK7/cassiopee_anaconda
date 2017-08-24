# - getCurvatureRadius (pyTree) -
import Geom.PyTree as D
import Converter.PyTree as C

a = D.circle((0,0,0), 1, 10, 0, 10)
a = D.getCurvatureRadius(a)
t = C.newPyTree(['Base',1]); t[2][1][2].append(a)
C.convertPyTree2File(t, 'out.cgns')
