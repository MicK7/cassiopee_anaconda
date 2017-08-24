# - newElements (pyTree) -
import Converter.Internal as Internal

# Create an elements node
b = Internal.newElements(name='Elements', etype='UserDefined', econnectivity=None, eboundary=0); print b

# Attach it to tree
z = Internal.newZone('Zone', [[10],[2],[0]], 'Unstructured')
Internal.newElements(name='Elements', etype='UserDefined', econnectivity=None, eboundary=0, parent=z); print z
