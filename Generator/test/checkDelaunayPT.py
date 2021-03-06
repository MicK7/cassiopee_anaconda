# - checkDelaunay (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Geom.PyTree as D
import Transform.PyTree as T

A = D.text1D('STEPHANIE')
A = C.convertArray2Tetra(A); a = T.join(A)
# Triangulation respecting given contour
tri = G.constrainedDelaunay(a)
res = G.checkDelaunay(a, tri)
C.convertPyTree2File(res, "out.cgns")
