# - newDimensionalUnits (pyTree) -
import Converter.Internal as Internal
import numpy

# Create a DimensionalUnits node
n = Internal.newDimensionalUnits(massUnit='Kilogram', lengthUnit='Meter', timeUnit='Second', temperatureUnit='Kelvin', angleUnit='Radian'); print n

# Attach it to a parent node
d = Internal.newGridCoordinates()
Internal.newDataClass('Dimensional', parent=d)
Internal.newDimensionalUnits('Kilogram', 'Meter', 'Second', 'Kelvin', 'Radian', parent=d)
print d
