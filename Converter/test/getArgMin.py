# - getArgMin (array) -
import Converter as C
import Generator as G

a = G.cart((0,0,0), (1.,1.,1.), (10,10,10))
a = C.initVars(a, 'F={x}+{y}+{z}')
argmin = C.getArgMin(a, 'F'); print argmin
