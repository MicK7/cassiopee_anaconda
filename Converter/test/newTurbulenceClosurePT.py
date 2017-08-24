# - newZone (pyTree) -
import Converter.Internal as Internal

# Create a zone node
z = Internal.newTurbulenceClosure(value='Null'); print z

# Create a zone node and attach it to tree
t = Internal.newFlowEquationSet()
z = Internal.newTurbulenceClosure(value='ReynoldsStress', parent=t)
print t
