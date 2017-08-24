# - newZone (pyTree) -
import Converter.Internal as Internal

# Create a zone node
z = Internal.newChemicalKineticsModel(value='Null'); print z

# Create a zone node and attach it to tree
t = Internal.newFlowEquationSet()
z = Internal.newChemicalKineticsModel(value='ChemicalNonequilib', parent=t)
print t
