# - newGridConnectivity1to1 (pyTree) -
import Converter.Internal as Internal

# Create a node
n = Internal.newPeriodic(rotationCenter=[0.,0.,0.], rotationAngle=[0.,0.,0.], translation=[0.,0.,0.]); print n

# Attach it to a parent node
d = Internal.newGridConnectivityProperty()
Internal.newPeriodic(rotationCenter=[0.,0.,0.], rotationAngle=[0.,0.,0.], translation=[0.,0.,0.], parent=d) 
print d
