# Interface pour MPI

try: from Mpi4py import *
except: raise ImportError("Converter:Mpi: requires mpi4py module.") 
import Distributed
import PyTree

#==============================================================================
# IN: t: full/loaded skel/partial
#==============================================================================
def center2Node(t, var='', cellNType=0, graph=None):
    if graph is None: graph = computeGraph(t, type='match')
    tl = addXZones(t, graph)
    tl = convert2PartialTree(tl)
    # print info
    zones = Internal.getZones(tl)
    print 'Rank %d has %d zones.'%(rank, len(zones))
    tl = PyTree.center2Node(tl, var, cellNType)
    tl = rmXZones(tl)
    return tl
