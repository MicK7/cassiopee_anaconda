# - newZone (pyTree) -
import Converter.Internal as Internal

# Create a zone node
z = Internal.newEMMagneticFieldModel(value='Null'); print z

# Create a zone node and attach it to tree
t = Internal.newFlowEquationSet()
z = Internal.newEMMagneticFieldModel(value='Interpolated', parent=t)
print t
