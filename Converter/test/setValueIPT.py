# - setValue (pyTree) -
import Converter.PyTree as C
import Converter.Internal as Internal
import numpy

node = ['node1', 12., [], 'DataArray_t']

# Set a scalar value in node
Internal.setValue(node, 1.); print node

# Set a numpy array in node 
Internal.setValue(node, numpy.zeros(10)); print node

# Set an array as a list
Internal.setValue(node, [1.,12.,13.]); print node
