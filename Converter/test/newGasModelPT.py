# - newZone (pyTree) -
import Converter.Internal as Internal

# Create a zone node
z = Internal.newGasModel(value='Ideal'); print z

# Create a zone node and attach it to tree
t = Internal.newFlowEquationSet()
z = Internal.newGasModel(value='Ideal', parent=t)
print t
