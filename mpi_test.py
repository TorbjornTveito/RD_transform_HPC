from mpi4py import MPI

comm = MPI.COMM_WORLD
mpi_rank = comm.Get_rank()
mpi_size = comm.Get_size()
#mpi_admin = mpi_size - 1



print(mpi_rank)