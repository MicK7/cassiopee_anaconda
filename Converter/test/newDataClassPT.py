# - newDataClass (pyTree) -
import Converter.Internal as Internal
import numpy

# Create a DataClass node
n = Internal.newDataClass('Dimensional'); print n

# Attach it to a parent node
d = Internal.newDiscreteData('DiscreteData')
Internal.newDataClass('Dimensional', parent=d)
print d
