# - newZone (pyTree) -
import Converter.Internal as Internal

# Create a zone node
z = Internal.newThermalRelaxationModel(value='Null'); print z

# Create a zone node and attach it to tree
t = Internal.newFlowEquationSet()
z = Internal.newThermalRelaxationModel(value='ThermalNonequilib', parent=t)
print t
