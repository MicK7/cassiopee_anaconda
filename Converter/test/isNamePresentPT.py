# - isNamePresent (PyTree) -
import Converter.PyTree as C
import Generator.PyTree as G

a = G.cart( (0,0,0), (1,1,1), (50,50,50) )
a = C.initVars(a, 'F', 1.)
a = C.initVars(a, 'centers:G', 0.)

b = G.cart( (0,0,0), (1,1,1), (50,50,50) )
b = C.initVars(b, 'F', 2.)
b = C.initVars(b, 'centers:H', 3.)

t = C.newPyTree(['Base',a,b])

print C.getVarNames(a)
print C.getVarNames([a, b])
print C.getVarNames(t[2][1])
print C.getVarNames(t)

print C.isNamePresent(a, 'F')
print C.isNamePresent(a, 'centers:F')
print C.isNamePresent(a, 'centers:G')
print C.isNamePresent([a, b], 'F')
print C.isNamePresent([a, b], 'centers:G')

print C.isNamePresent(t[2][1], 'F')
print C.isNamePresent(t[2][1], 'centers:G')
print C.isNamePresent(t, 'F')
print C.isNamePresent(t, 'centers:G')
