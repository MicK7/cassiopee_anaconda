# - newZone (pyTree) -
import Converter.Internal as Internal

# Create a zone node
z = Internal.newUserDefinedData(name='UserDefined', value=None); print z

# Create a zone node and attach it to tree
b = Internal.newFamily(name='FamInjection')
z = Internal.newFamilyBC(value='BCInflow', parent=b)
z = Internal.newUserDefinedData(name='.Solver#BC', value=None, parent=b)
print b
