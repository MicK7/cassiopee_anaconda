# - newAxiSymmetry (pyTree) -
import Converter.Internal as Internal

# Create an axisymmetry node
n = Internal.newAxiSymmetry(referencePoint=[0.,0.,0.], axisVector=[0.,0.,0.])
print n

# Attach it to base
b = Internal.newCGNSBase('Base', 3, 3)
Internal.newAxiSymmetry([0.,0.,0.], [0.,0.,0.], parent=b); print b
