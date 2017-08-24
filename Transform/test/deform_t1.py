# - deform (array) -
import Converter as C
import Transform as T
import KCore.test as test

def F(a):
    b = C.initVars(a, 'dx', 1)
    b = C.initVars(b, 'dy', 0.5)
    b = C.initVars(b, 'dz', 0.5)
    n = C.extractVars(b, ['dx', 'dy', 'dz'])
    a = T.deform(a, n)
    return a

test.stdTestA(F)
