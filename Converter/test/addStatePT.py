# - addState (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Converter.Internal as Internal

a = G.cylinder((0,0,0), 1., 1.5, 0., 360., 1., (80,30,2))
t = C.newPyTree(['Base',a])

# Specifie des valeurs
b = Internal.getNodeFromName1(t, 'Base')
C._addState(b, 'EquationDimension', 2)
C._addState(b, 'GoverningEquations', 'Euler')
C._addState(b, 'Mach', 0.6)
C._addState(b, 'Reynolds', 100000)

# Specifie un etat de reference adimensionne par:
# Mach, alpha, Re, MutSMu, TurbRate (adim1)
C._addState(t, adim='adim1', MInf=0.5, alphaZ=0., alphaY=0.,
            ReInf=1.e8, MutSMuInf=0.2, TurbRateInf=1.e-8)

# Specifie un etat de reference adimensionne par:
# Mach, alpha, Re, MutSMu, TurbRate (adim2)
C._addState(t, adim='adim2', MInf=0.5, alphaZ=0., alphaY=0.,
            ReInf=1.e8, MutSMuInf=0.2, TurbRateInf=1.e-8)

# Specifie un etat de reference dimensionne par:
# U, T, P, L, MutSMu, TurbRate (dim1)
C._addState(t, adim='dim1', UInf=35, TInf=294., PInf=101325, LInf=1.,
            alphaZ=0., alphaY=0., MutSMuInf=0.2, TurbRateInf=1.e-8)

# Specifie un etat de reference dimensionne par:
# U, T, Ro, L, MutSMu, TurbRate (dim2)
C._addState(t, adim='dim2', UInf=35, TInf=294., RoInf=1.2, LInf=1.,
            alphaZ=0., alphaY=0., MutSMuInf=0.2, TurbRateInf=1.e-8)

# Specifie un etat de reference dimensionne par:
# U, P, Ro, L, MutSMu, TurbRate (dim2)
C._addState(t, adim='dim3', UInf=35, PInf=101325., RoInf=1.2, LInf=1.,
            alphaZ=0., alphaY=0., MutSMuInf=0.2, TurbRateInf=1.e-8)
