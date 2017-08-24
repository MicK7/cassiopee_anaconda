# - getMinValue (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G

a = G.cart((0,0,0), (1.,1.,1.), (11,1,1))
minval = C.getMinValue(a, 'CoordinateX'); print minval
minval = C.getMinValue(a, ['CoordinateX', 'CoordinateY']); print minval
minval = C.getMinValue(a, 'GridCoordinates'); print minval
