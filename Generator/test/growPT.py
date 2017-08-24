# - grow (pyTree)-
import Converter.PyTree as C
import Generator.PyTree as G
import Geom.PyTree as D

a = D.sphere((0,0,0), 1., 50)
a = G.getNormalMap(a)
a = C.center2Node(a, C.__FlowSolutionCenters__)
a = C.rmVars(a, C.__FlowSolutionCenters__)
b = G.grow(a, ['sx','sy','sz'])
t = C.newPyTree(['Base1',2,'Base2',3])
t[2][1][2].append(a); t[2][2][2].append(b)
C.convertPyTree2File(t, 'out.cgns')
