# - getTangent (array) -
import Geom as D

c = D.polyline([(0,0,0),(1,1,0),(2,-1,0)])
a = D.spline(c, order=3, density=10.)  

print D.getTangent(a)

