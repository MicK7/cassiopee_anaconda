# - newRotatingCoordinates (pyTree) -
import Converter.Internal as Internal

# Create an rotating coordinates node
b = Internal.newRotatingCoordinates(rotationCenter=[0.,0.,0.], rotationRateVector=[0.,0.,0.]); print b

# Attach it to tree
b = Internal.newCGNSBase('Base', 3, 3)
Internal.newRotatingCoordinates([0.,0.,0.], [0.,0.,0.], parent=b); print b
