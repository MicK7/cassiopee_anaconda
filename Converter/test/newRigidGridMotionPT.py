# - newZoneGridConnectivity (pyTree) -
import Converter.Internal as Internal

# Create a node
n = Internal.newRigidGridMotion(name='Motion', origin=[0.,0.,0.], mtype='ConstantRate'); print n

# Attach it to a parent node
z = Internal.newZone('Zone', zsize=[[10],[2],[0]], ztype='Structured')
Internal.newRigidGridMotion(name='Motion', origin=[0.,0.,0.], mtype='ConstantRate', parent=z) 
print z
