# - newGridCoordinates (pyTree) -
import Converter.Internal as Internal

# Create a GridCoordinates node
n = Internal.newGridCoordinates(); print n

# Create a zone node
z = Internal.newZone('Zone', zsize=[[10],[2],[0]], ztype='Structured')
n = Internal.newGridCoordinates(parent=z)
print z
