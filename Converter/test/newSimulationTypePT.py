# - newSimulationType (pyTree) -
import Converter.Internal as Internal

# Create a simulation type node
b = Internal.newSimulationType(value='NonTimeAccurate'); print b

# Attach it to parent
b = Internal.newCGNSBase('Base', 3, 3)
Internal.newSimulationType('TimeAccurate', parent=b); print b
