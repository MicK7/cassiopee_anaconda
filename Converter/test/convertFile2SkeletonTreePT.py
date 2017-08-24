# - convertFile2SkeletonTree (pyTree) -
import Converter.PyTree as C
import Generator.PyTree as G
import Converter.Mpi as Cmpi
import Converter

if Cmpi.rank == 0:
    a = G.cart((0.,0.,0.),(0.1,0.1,0.1),(11,11,11))
    t = C.newPyTree(['Base', a])
    C.convertPyTree2File(t, 'in.cgns')
Cmpi.barrier()

t1 = Cmpi.convertFile2SkeletonTree('in.cgns'); print t1
