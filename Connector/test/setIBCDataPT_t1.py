# - setIBCData (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Connector.PyTree as X
import Post.PyTree as P
import Dist2Walls.PyTree as DTW
import Geom.PyTree as D
import Transform.PyTree as T
import KCore.test as test
import numpy as N

a = G.cart((-1,-1,-1),(0.1,0.1,1),(21,21,3))
s = G.cylinder((0,0,-1), 0, 0.4, 360, 0, 4, (31,31,5)) 
s = C.convertArray2Tetra(s); s = T.join(s); s = P.exteriorFaces(s)
t = C.newPyTree(['Base']); t[2][1][2] = [a]
tb = C.newPyTree(['Base']); tb[2][1][2] = [s]

# Blanking
bodies = [[s]]
BM = N.array([[1]],N.int32)
t = X.blankCells(t, bodies, BM, blankingType='center_in')
t = X.setHoleInterpolatedPoints(t, depth=-1)
# Dist2Walls
t = DTW.distance2Walls(t, [s], type='ortho', loc='centers', signed=1)
t = C.center2Node(t, 'centers:TurbulentDistance')
# Gradient de distance localise en centres => normales
t = P.computeGrad(t, 'TurbulentDistance')
t2 = X.setIBCData(t, t, loc='centers', storage='direct')
test.testT(t2,1)
t2 = X.setIBCData(t, t, loc='centers', storage='direct',hi=0.1)
test.testT(t2,2)

# NODES
t = C.newPyTree(['Base']); t[2][1][2] = [a]
# Blanking
bodies = [[s]]
BM = N.array([[1]],N.int32)
t = X.blankCells(t,bodies,BM,blankingType='node_in')
t = X.setHoleInterpolatedPoints(t,depth=-1)
# Dist2Walls
t = DTW.distance2Walls(t, [s], type='ortho', loc='nodes', signed=1)
# Gradient de distance localise en centres => normales
t = P.computeGrad(t, 'TurbulentDistance')
t = C.center2Node(t, 'FlowSolution#Centers')
t2 = X.setIBCData(t, t, loc='nodes', storage='direct')
test.testT(t2,3)
t2 = X.setIBCData(t, t, loc='nodes', storage='direct', hi=0.1)
test.testT(t2,4)

# COMBINATION INSIDE/OUTSIDE
t = C.newPyTree(['Base']); t[2][1][2] = [a]

# Blanking
bodies = [[s]]
BM = N.array([[1]],N.int32)
t = X.blankCells(t,bodies,BM,blankingType='center_in')
t = X.setHoleInterpolatedPoints(t,depth=-1)
t = X.setHoleInterpolatedPoints(t,depth= 2)
# Dist2Walls
t = DTW.distance2Walls(t,[s],type='ortho',loc='centers',signed=1)
t = C.center2Node(t,'centers:TurbulentDistance')
# Gradient de distance localise en centres => normales
t = P.computeGrad(t, 'TurbulentDistance')
t2 = X.setIBCData(t, t, loc='centers', storage='direct',hi=0.,he=0.1*1.4)
test.testT(t2,5)
