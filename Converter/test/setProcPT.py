# - setProc (pyTree) -
import Converter.PyTree as C
import Converter.Internal as Internal
import Converter.Mpi as Cmpi
import Distributor2.PyTree as Distributor2
import Generator.PyTree as G
import KCore.test as test

a = G.cart((0,0,0), (1,1,1), (10,10,10))
b = G.cart((12,0,0), (1,1,1), (10,10,10))
t = C.newPyTree(['Base',a,b])
t = Cmpi.setProc(t, 1)

zones = Internal.getNodesFromType(t, 'Zone_t')
for z in zones:
    print z[0]+' -> '+str(Cmpi.getProc(z))
