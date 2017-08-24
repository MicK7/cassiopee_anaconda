# - createNode (pyTree) -
import Converter.PyTree as C
import Converter.Internal as Internal
import numpy

# Create a node named myNode, of type DataArray_t, with one value 1.
node = Internal.createNode('myNode', 'DataArray_t', value=1., children=[])
print node

# Create a node named myNode, of type DataArray_t, with an array valued
# in a list
node = Internal.createNode('myNode', 'DataArray_t', value=[12.,14.,15.], 
                           children=[])
print node

# Create a node named myNode, of type DataArray_t, with a given numpy
a = numpy.zeros( (10) )
a[1] = 1.; a[8] = 2.
node = Internal.createNode('myNode', 'DataArray_t', value=a, children=[])
print node

# Create a node named GoverningEquation, of type 'GoverningEquation_t' and 
# value 'Euler'
node = Internal.createNode('GoverningEquation', 'GoverningEquation_t', 
                           value='Euler')
print node
