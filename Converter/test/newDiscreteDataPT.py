# - newDiscreteData (pyTree) -
import Converter.Internal as Internal

# Create a discrete data node
b = Internal.newDiscreteData(name='DiscreteData'); print b

# Attach it to tree
z = Internal.newZone('Zone', [[10],[2],[0]], 'Structured')
Internal.newDiscreteData('DiscreteData', parent=z); print z
