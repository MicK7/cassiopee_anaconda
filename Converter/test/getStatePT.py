# - getState (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G

a = G.cylinder((0,0,0), 1., 1.5, 0., 360., 1., (80,30,2))
t = C.newPyTree(['Base',a])

# Specifie un etat de reference adimensionne par:
# Mach, alpha, Re, MutSMu, TurbRate (adim1)
C._addState(t, adim='adim1', MInf=0.5, alphaZ=0., alphaY=0., 
            ReInf=1.e8, MutSMuInf=0.2, TurbRateInf=1.e-8) 

# Get the ref state
state = C.getState(t)
print state
