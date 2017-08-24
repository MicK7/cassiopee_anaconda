# - newOrdinal (pyTree) -
import Converter.Internal as Internal

# Create a ordinal node
b = Internal.newOrdinal(value=1); print b

# Attach it to tree
z = Internal.newZone('Zone', [[10],[2],[0]], 'Structured')
Internal.newOrdinal(value=1, parent=z); print z
