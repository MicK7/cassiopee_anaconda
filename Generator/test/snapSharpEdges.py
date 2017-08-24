# - snapSharpEdges (array) -
import Generator as G
import Converter as C
import Geom as D
import Connector as X
import Transform as T

# polylignes avec angles vifs
s = D.polyline([(0.02,0,0),(1,1,0),(2,1,0),(0.02,0,0)])
# Grille cartesienne (reguliere)
h = 0.1
ni = 30; nj = 20; nk=1
b = G.cart( (-0.5, -0.5, 0), (h, h, 1.), (ni,nj,nk) )

b = G.snapSharpEdges(b, [s], h)

C.convertArrays2File([b,s], 'out.plt')

 
