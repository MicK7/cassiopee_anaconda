# - newZone (pyTree) -
import Converter.Internal as Internal

# Create a zone node
z = Internal.newEMElectricFieldModel(value='Null'); print z

# Create a zone node and attach it to tree
t = Internal.newFlowEquationSet()
z = Internal.newEMElectricFieldModel(value='Voltage', parent=t)
print t
