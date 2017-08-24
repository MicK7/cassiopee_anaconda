# - newBCData (pyTree) -
import Converter.Internal as Internal

# Create a BC data node
n = Internal.newBCData(name='BCData'); print n

# Attach it to a parent node
d = Internal.newBCDataSet(name='BCDataSet', value='UserDefined')
Internal.newBCData('BCNeumann', parent=d)
print d
