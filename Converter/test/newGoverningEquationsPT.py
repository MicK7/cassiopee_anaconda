# - newZone (pyTree) -
import Converter.Internal as Internal

# Create a zone node
z = Internal.newGoverningEquations(value='Euler'); print z

# Create a zone node and attach it to tree
t = Internal.newFlowEquationSet()
z = Internal.newGoverningEquations(value='NSTurbulent', parent=t)
print t
