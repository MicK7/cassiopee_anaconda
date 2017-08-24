# - newZone (pyTree) -
import Converter.Internal as Internal

# Create a zone node
z = Internal.newFamilyBC(value='BCWall'); print z

# Create a zone node and attach it to tree
b = Internal.newFamily(name='FamInjection')
z = Internal.newFamilyBC(value='UserDefined', parent=b)
print b
