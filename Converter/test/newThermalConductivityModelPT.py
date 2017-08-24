# - newZone (pyTree) -
import Converter.Internal as Internal

# Create a zone node
z = Internal.newThermalConductivityModel(value='Null'); print z

# Create a zone node and attach it to tree
t = Internal.newFlowEquationSet()
z = Internal.newThermalConductivityModel(value='SutherlandLaw', parent=t)
print t
