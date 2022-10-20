c-------------------------------------------
c       Initialization 
c-------------------------------------------
      subroutine adios2_init()
      include 'SIZE'
      include 'TOTAL'
      include 'RESTART'
      common /nekmpi/ nid_,np_,nekcomm
      nvals = lx1 * ly1 * lz1
      call adios2_setup(nvals,nelv,nelb,nelgv,nelgt,xm1,ym1,zm1,
     & nekcomm)
      end
c-------------------------------------------
c       Write 
c-------------------------------------------
      subroutine adios2_write()
      implicit none
      include 'SIZE'
      include 'INPUT'
      include 'TSTEP'
      include 'SOLN'
      !include 'TOTAL'
      include 'PARALLEL'
      include 'RESTART'
      integer i,nvals, nid_,np_,nekcomm
      common /nekmpi/ nid_,np_,nekcomm
      nvals = lx1 * ly1 * ly1
      end
c-------------------------------------------
c       Finalization 
c-------------------------------------------
      subroutine adios2_end()
      call adios2_finalize()
      end
