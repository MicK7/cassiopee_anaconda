# - node2Center (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G

def F(x,y): return 2*x+y

ni = 30; nj = 40; nk = 3
a = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,nk))
a = C.initVars(a, 'Density=2*{CoordinateX}+{CoordinateY}')

# node2Center : passe une variable en centres (dans la meme zone)
a = C.node2Center(a, 'Density')
C.convertPyTree2File(a, 'out1.cgns')

# node2Center : cree une nouvelle zone contenant les centres
a = G.cart((0,0,0), (10./(ni-1),10./(nj-1),1), (ni,nj,nk))
a = C.initVars(a, 'Density', F, ['CoordinateX','CoordinateY'])
b = C.node2Center(a); b[0] = a[0]+'_centers'
C.convertPyTree2File([a,b], 'out2.cgns')
