# - addState (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G

a = G.cylinder((0,0,0), 1., 1.5, 0., 360., 1., (80,30,2))
t = C.newPyTree(['Base',a])
t = C.addState(t, adim='adim1', MInf=0.5)
C.convertPyTree2File(t, 'out.cgns')
