# - newDataArray (pyTree) -
import Converter.Internal as Internal
import numpy

# Create a DataArray node
n = Internal.newDataArray('CoordinateX', numpy.zeros(10)); print n

# Attach it to a parent node
g = Internal.newGridCoordinates()
Internal.newDataArray('CoordinateX', value=numpy.arange(0,10), parent=g)
Internal.newDataArray('CoordinateY', value=numpy.zeros(10), parent=g)
Internal.newDataArray('CoordinateZ', value=numpy.zeros(10), parent=g)
print g
