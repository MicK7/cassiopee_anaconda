# - convertArrays2File (array) -
import Generator as G
import Converter as C

# Create a cartesian mesh
a = G.cart( (0,0,0), (0.1, 0.2, 1.), (11, 11, 2))
C.convertArrays2File([a], 'out.plt', 'bin_tp')
