# - volumeFromCrossSections (array) -
import Converter as C
import Geom as D

contours = []
a = D.polyline([(0.,0.,0.),(1.,0.,0.),(1.,1.,0.),
                (0., 1., 0.),(0.,0.,0.)])
contours = [D.polyline([(0.,0.,2.),(1.,0.,2.),(1.,1.,2.),
                        (0., 1., 2.),(0.,0.,2.)])]
contours.append(a)

vol = D.volumeFromCrossSections(contours)
C.convertArrays2File([vol]+contours, 'out.plt')
