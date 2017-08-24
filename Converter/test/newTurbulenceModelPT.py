# - newZone (pyTree) -
import Converter.Internal as Internal

# Create a zone node
z = Internal.newTurbulenceModel(value='TwoEquation_MenterSST'); print z

# Create a zone node and attach it to tree
t = Internal.newFlowEquationSet()
z = Internal.newTurbulenceModel(value='OneEquation_SpalartAllmaras', parent=t)
print t
