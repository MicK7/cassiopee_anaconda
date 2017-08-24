# - newDescriptor (pyTree) -
import Converter.Internal as Internal
import numpy

# Create a Descriptor node
n = Internal.newDescriptor(name='Descriptor', value='Mesh exported from Cassiopee 2.2'); print n

# Attach it to a parent node
b = Internal.newCGNSBase('Base')
Internal.newDescriptor('Descriptor', 'Aircraft with nacelle Mach=0.7 Re=4 million alpha=-2', parent=b)
print b
