# - newFlowSolution (pyTree) -
import Converter.Internal as Internal

# Create a flow solution node
n = Internal.newFlowSolution(name='FlowSolution', gridLocation='Vertex'); print n

# Attach it to a parent node
z = Internal.newZone('Zone', zsize=[[10],[2],[0]], ztype='Structured')
Internal.newFlowSolution(name='FlowSolution', gridLocation='Vertex', parent=z) 
print z
