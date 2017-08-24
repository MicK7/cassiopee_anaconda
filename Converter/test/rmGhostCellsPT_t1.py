# - addGhostCells (pyTree) -
import Generator.PyTree as G
import Converter.PyTree as C
import Converter.Internal as Internal
import Transform.PyTree as T
import Connector.PyTree as X
import KCore.test as test

a = G.cart((1,1,1), (1.,1.,1.), (4,2,3)); a[0]='cart1'
b = G.cart((4,1,1), (0.5,1.,1.), (4,2,3)); b[0]='cart2'
c = G.cart((1,1,-3), (1.,1.,0.5), (4,2,9)); c[0]='cart3'

a = T.reorder(a, (-2,1,3))
b = T.reorder(b, (1,2,3))
c = T.reorder(c, (3,1,2))

# Physical BC (here BCWall)
a = C.addBC2Zone(a, 'wall1', 'BCWall', 'jmin')

# Chimere
#-----------------------------
# initialisation cellNatureField
a = C.addVars(a, 'centers:cellN')
a = C.initVars(a, 'centers:cellN', 1.)
b = C.addVars(b, 'centers:cellN')
b = C.initVars(b, 'centers:cellN', 0.)
#-----------------------------
# initialisation density

# Create a function
def F(x1, x2):
    return 3.*x1+2.*x2
a = C.addVars(a, 'nodes:Density')
a = C.initVars(a, 'Density', F, ['CoordinateX','CoordinateY'])
a = C.addVars(a, 'nodes:StagnationPressure')
a = C.initVars(a, 'StagnationPressure', 0.5)

c = C.addVars(c, 'nodes:StagnationPressure')
c = C.initVars(c, 'StagnationPressure', 1.)

t = C.newPyTree(['Base']); t[2][1][2] += [a,b,c]
# Matching BC
t = X.connectMatch(t)
a = t[2][1][2][0]
ag = Internal.addGhostCells(t, a, 2)
t[2][1][2].append(ag)
ag = Internal.rmGhostCells(t, ag, 3)
t[2][1][2].append(ag)
ag = Internal.rmGhostCells(t, ag, 1)
t[2][1][2].append(ag)
ag = Internal.rmGhostCells(t, ag, 1)
t[2][1][2].append(ag)
test.testT(t, 1)
