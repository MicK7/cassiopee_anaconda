# - newIntegralData (pyTree) -
import Converter.Internal as Internal

# Create an integral data node
b = Internal.newIntegralData(name='IntegralData'); print b

# Attach it to tree
b = Internal.newCGNSBase('Base', 3, 3)
Internal.newIntegralData('IntegralData', parent=b); print b
