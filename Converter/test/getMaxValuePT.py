# - getMaxValue (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G

a = G.cart((0,0,0), (1.,1.,1.), (11,2,2))
maxval = C.getMaxValue(a, 'CoordinateX'); print maxval
maxval = C.getMaxValue(a, ['CoordinateX', 'CoordinateY']); print maxval
maxval = C.getMaxValue(a, 'GridCoordinates'); print maxval
