      program NEKTON

c     Original file
c      include 'mpif.h'
c      integer comm
c      comm = MPI_COMM_WORLD
c
c      call nek_init(comm)
c      call nek_solve()
c      call nek_end()

c     Modified for ADIOS2
      include 'mpif.h'
      integer comm, colour, rank, ierror

      call mpi_init(ierror)
      call mpi_comm_rank(MPI_COMM_WORLD, rank, ierror)
      colour = 0
      call mpi_comm_split(MPI_COMM_WORLD, colour, rank, comm, ierror)

      call nek_init(comm)
      call nek_solve()
      call nek_end()

      end
