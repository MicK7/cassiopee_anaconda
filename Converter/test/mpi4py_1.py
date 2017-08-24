# - test de base sur mpi4py -
try: from mpi4py import MPI
except: raise ImportError, 'FAILED mpi4py'

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
print rank, size
