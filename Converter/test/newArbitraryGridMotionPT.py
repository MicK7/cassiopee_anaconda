# - newZoneGridConnectivity (pyTree) -
import Converter.Internal as Internal

# Create a node
n = Internal.newArbitraryGridMotion(name='Motion', value='DeformingGrid'); print n

# Attach it to a parent node
z = Internal.newZone('Zone', zsize=[[10],[2],[0]], ztype='Structured')
Internal.newArbitraryGridMotion(name='Motion', value='NonDeformingGrid', parent=z) 
print z
