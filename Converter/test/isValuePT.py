# - isValue (pyTree) -
import Converter.PyTree as C
import Converter.Internal as Internal
import numpy

# Check a scalar value
node = ['node1', 1., [], 'DataArray_t']
print Internal.isValue(node, 1.)

# Check a numpy array values
node = ['node1', numpy.zeros(10), [], 'DataArray_t']
print Internal.isValue(node, numpy.zeros(10))

# Check a string value
node = ['node1', 'toto', [], 'DataArray_t']
print Internal.isValue(node, 'toto')
