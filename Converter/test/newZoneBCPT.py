# - newZoneBC (pyTree) -
import Converter.Internal as Internal

# Create a zoneBC node
n = Internal.newZoneBC(); print n

# Attach it to a parent node
z = Internal.newZone('Zone', zsize=[[10],[2],[0]], ztype='Structured')
Internal.newZoneBC(parent=z) 
print z
