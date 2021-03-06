# - connectMatch 2D (array)-
import Generator as G
import Connector as X
import KCore.test as test
import Transform as T
import Converter as C

a = G.cart((0.,0.,0.), (0.1, 0.1, 0.1), (11, 21, 1))
ca = C.array('cellN',10,20,0)
ca = C.initVars(ca, 'cellN', 1.)
a = C.addVars(a, 'cellN')

# partiellement coincident
a2 = G.cart((1., 0.5, 0.), (0.1, 0.1, 0.1), (11, 21, 1))
ca2 = C.array('cellN',10,20,0)
ca2 = C.initVars(ca2, 'cellN', 1.)
a2 = C.addVars(a2, 'cellN')
res = X.connectMatch( a, a2, dim=2)
test.testO(res,1)
# coincident
a3 = G.cart((1.00001, 0., 0.), (0.1, 0.1, 0.1), (11, 21, 1))
res = X.connectMatch( a,a3,sameZone=0,tol=1e-4,dim=2)
test.testO(res,2)
# raccord 1 sur 2
a3 = T.oneovern(a3,(2,2,1))
res = X.connectMatch( a,a3,sameZone=0,tol=1e-4,dim=2)
test.testO(res,3)
