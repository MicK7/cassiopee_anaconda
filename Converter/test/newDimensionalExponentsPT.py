# - newDimensionalExponents (pyTree) -
import Converter.Internal as Internal
import numpy

# Create a DimensionalExponents node
# Example for velocity exponents (m/s)
n = Internal.newDimensionalExponents(massExponent=0., lengthExponent=1., timeExponent=-1., temperatureExponent=0., angleExponent=0.); print n

# Attach it to a parent node
d = Internal.newGridCoordinates()
Internal.newDataClass('Dimensional', parent=d)
Internal.newDimensionalExponents(massExponent=0., lengthExponent=1., timeExponent=0., temperatureExponent=0., angleExponent=0., parent=d) 
print d
