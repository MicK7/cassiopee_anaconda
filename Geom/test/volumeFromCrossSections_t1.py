# - volumeFromCrossSections -
import Converter as C
import Geom as D
import KCore.test as test

contours = []
a = D.polyline([(0.,0.,0.),(1.,0.,0.),(1.,1.,0.),(0., 1., 0.),(0.,0.,0.)])
contours = [D.polyline([(0.,0.,2.),(1.,0.,2.),(1.,1.,2.),(0., 1., 2.),(0.,0.,2.)])]
contours.append(a)

vol = D.volumeFromCrossSections(contours)
test.testA([vol],1)
