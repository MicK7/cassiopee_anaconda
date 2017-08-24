# - convertArray2Node (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G

a = G.cart((0.,0.,0.), (0.1,0.1,0.2), (10,10,1))
a = C.convertArray2Node(a)
t = C.newPyTree(['Base1',2,a])
C.convertPyTree2File(t, 'out.cgns')
