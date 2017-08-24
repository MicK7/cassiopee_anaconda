# - nearestNodes (array) -
import Converter as C
import Generator as G
import Transform as T
import Post as P

a = G.cart( (0,0,0), (1,1,1), (10,10,10) )
b = T.translate(a,(0.15,0.,0.))
f = P.exteriorFaces(b)

# Enregistre les noeuds de a dans le hook
hook = C.createHook(a, function='nodes')
# Indices des noeuds de a les plus proches des noeuds de f
# et distance correspondante
nodes,dist = C.nearestNodes(hook, f)
print nodes,dist
