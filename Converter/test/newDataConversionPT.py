# - newDataConversion (pyTree) -
import Converter.Internal as Internal
import numpy

#    Data(raw) = Data(nondimensional)*ConversionScale + ConversionOffset

# Create a DataConversion node
n = Internal.newDataConversion(conversionScale=1., conversionOffset=0.) ; print n

# Attach it to a parent node
d = Internal.newDiscreteData('DiscreteData')
Internal.newDataClass('NondimensionalParameter', parent=d)
Internal.newDataConversion(conversionScale=1., conversionOffset=0., parent=d) 
print d
