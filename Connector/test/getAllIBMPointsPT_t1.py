# - getAllIBMPoints (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Connector.ToolboxIBM as IBM
import Post.PyTree as P
import Converter
import KCore.test as test
N = 21
a = G.cart((0,0,0),(1,1,1),(N,N,N))
C._initVars(a,'centers:TurbulentDistance={centers:CoordinateZ}')
C._initVars(a,'centers:cellN=({centers:CoordinateZ}>0.8)*1.+({centers:CoordinateZ}<0.8)*2.')
a = P.computeGrad(a,'centers:TurbulentDistance')
res = IBM.getAllIBMPoints(a, loc='centers',hi=0.,he=2.,tb=None)
test.testA(res[0]+res[1])
test.testA(res[2],12)
N = 21
a = G.cart((0,0,0),(1,1,1),(N,N,N))
C._initVars(a,'TurbulentDistance={CoordinateZ}')
C._initVars(a,'cellN=({CoordinateZ}>0.8)*1.+({CoordinateZ}<0.8)*2.')
a = P.computeGrad(a,'TurbulentDistance'); a = C.center2Node(a,'FlowSolution#Centers')
res = IBM.getAllIBMPoints(a, loc='nodes',hi=0.,he=2.,tb=None)
test.testA(res[0]+res[1],2)
test.testA(res[2],22)
