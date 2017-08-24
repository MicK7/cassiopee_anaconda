# - getDistribution (PyTree) -
import Geom.PyTree as D

Foil = D.naca(12., N=49)
print D.getDistribution(Foil)




