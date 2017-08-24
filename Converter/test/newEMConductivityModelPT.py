# - newZone (pyTree) -
import Converter.Internal as Internal

# Create a zone node
z = Internal.newEMConductivityModel(value='Null'); print z

# Create a zone node and attach it to tree
t = Internal.newFlowEquationSet()
z = Internal.newEMConductivityModel(value='Chemistry_LinRessler', parent=t)
print t
