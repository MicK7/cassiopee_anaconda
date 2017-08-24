# - volumeFromCrossSection (pyTree) -
import Converter.PyTree as C
import Geom.PyTree as D

contours = []
a = D.polyline([(0.,0.,0.),(1.,0.,0.),(1.,1.,0.),(0., 1., 0.),(0.,0.,0.)])
contours = [D.polyline([(0.,0.,2.),(1.,0.,2.),(1.,1.,2.),(0., 1., 2.),(0.,0.,2.)])]
contours.append(a)

vol = D.volumeFromCrossSections(contours)
t = C.newPyTree(['Base',3,vol])
C.convertPyTree2File(t, 'out.cgns')
