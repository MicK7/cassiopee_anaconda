# - newIndexArray (pyTree) -
import Converter.Internal as Internal
import numpy

# Create a node
n = Internal.newIndexArray(name='Index', value=[101,51,22,1024,78]); print n

# Attach it to a parent node
d = Internal.newFlowSolution('FlowSolution', 'Vertex')
Internal.newIndexArray('Index', [101,51,22,1024,78], parent=d)
print d
