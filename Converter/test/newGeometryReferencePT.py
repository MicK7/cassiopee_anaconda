# - newZone (pyTree) -
import Converter.Internal as Internal

# Create a zone node
z = Internal.newGeometryReference(value='ICEM-CFD', file='MyCAD.tin'); print z

# Create a zone node and attach it to tree
b = Internal.newFamily(name='FamWall')
z = Internal.newGeometryReference(value='NASA-IGES', file='MyCAD.iges', parent=b)
print b
