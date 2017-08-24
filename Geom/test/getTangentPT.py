# - getTangent (PyTree) -
import Geom.PyTree as D
import Converter.PyTree as C

c = D.polyline([(0,0,0),(1,1,0),(2,-1,0)])
a = D.spline(c, order=3, density=10.)  

print D.getTangent(a)

