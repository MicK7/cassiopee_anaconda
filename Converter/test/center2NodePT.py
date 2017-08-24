# - center2Node (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G

def F(x,y): return 2*x+y

# center2Node: create a new zone
ni = 30; nj = 40; nk = 2
a = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,nk))
a = C.initVars(a, 'centers:Density', 1.)
b = C.center2Node(a); b[0] = a[0]+'_nodes'
t = C.newPyTree(['Base1',3,b])
C.convertPyTree2File(t, 'out0.cgns')

# center2Node: modify a variable
ni = 30; nj = 40; nk = 3
a = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,nk))
a = C.initVars(a, 'centers:Density', 1.)
a = C.center2Node(a, 'centers:Density')
t = C.newPyTree(['Base',3,a])
C.convertPyTree2File(t, 'out.cgns')
