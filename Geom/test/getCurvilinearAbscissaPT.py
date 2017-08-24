# - getCurvilinearAbscissa (pyTree)-
import Converter.PyTree as C
import Geom.PyTree as D

a = D.line((0.,0.,0.), (1.,0.,0), 100)
a = D.getCurvilinearAbscissa(a)
t = C.newPyTree(['Base',1]); t[2][1][2].append(a)
C.convertPyTree2File(t, 'out.cgns')
