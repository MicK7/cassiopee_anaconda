# - volumeFromCrossSection (pyTree) -
import Converter.PyTree as C
import Geom.PyTree as D
import KCore.test as test

contours = []
a = D.polyline([(0.,0.,0.),(1.,0.,0.),(1.,1.,0.),(0., 1., 0.),(0.,0.,0.)])
contours = [D.polyline([(0.,0.,2.),(1.,0.,2.),(1.,1.,2.),(0., 1., 2.),(0.,0.,2.)])]
contours.append(a)

vol = D.volumeFromCrossSections(contours)
t = C.newPyTree(['Base',2,vol])
test.testT(t, 1)
