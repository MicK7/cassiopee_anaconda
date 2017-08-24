# - getFamilyZoneNames (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G

a = G.cart((0.,0.,0), (0.01,0.01,1.), (20,20,2))
b = G.cart((1.,0.,0), (0.01,0.01,1.), (20,20,2))

a = C.tagWithFamily(a, 'CARTER') 
b = C.tagWithFamily(b, 'CARTER')

t = C.newPyTree(['Base',a,b])
t[2][1] = C.addFamily2Base(t[2][1], 'CARTER')

# Toutes les family zone names de l'arbre
names = C.getFamilyZoneNames(t); print names
