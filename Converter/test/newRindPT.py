# - newRind (pyTree) -
import Converter.Internal as Internal

# Create a Rind node (Rind contains the number of ghost cells for a structured block)
n = Internal.newRind([0,0,0,0,1,1]); print n

# Attach it to a parent node
d = Internal.newGridCoordinates()
Internal.newRind([0,0,0,0,1,1], parent=d)
print d
