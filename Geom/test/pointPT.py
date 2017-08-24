# - point (pyTree) -
import Geom.PyTree as D
import Converter.PyTree as C

a = D.point((0,0,0))
b = D.point((1,1,1))
C.convertPyTree2File([a,b], "out.cgns")
