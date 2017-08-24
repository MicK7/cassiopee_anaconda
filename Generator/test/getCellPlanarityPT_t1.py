# - getCellPlanarity (pyTree) -
import KCore.test as test
import Generator.PyTree as G
import Geom.PyTree as D
import Converter.PyTree as C

a = D.sphere( (0,0,0), 1., 10)
a = C.initVars(a,'Density',2.); a = C.initVars(a,'centers:cellN',1.)
a = G.getCellPlanarity(a)
test.testT(a,1)

