c-------------------------------------------
c       Initialization 
c-------------------------------------------
      subroutine async_catalyst_init()
      implicit none
      include 'SIZE'
      include 'INPUT'
      include 'SOLN'
      include 'MASS'
      include 'GEOM'
      include 'TSTEP'
      include 'RESTART'
      include 'PARALLEL'
      integer nid_,np_,nekcomm
      common /nekmpi/ nid_,np_,nekcomm
      call async_catalyst_setup(
     &lx1,ly1,lz1,
     &lx2,ly2,lz2,
     &nelv,nelt,
     &nelb,nelb,
     &nelgv,nelgt,
     &iostep,if3d,
     &time,dt,
     &xm1,ym1,zm1,
     &pr,vx,vy,vz,
     &t,bm1,
     &nekcomm)
      end

c-------------------------------------------
c       Write 
c-------------------------------------------
      subroutine async_catalyst_write()
      implicit none
      include 'SIZE'
      include 'INPUT'
      include 'SOLN'
      include 'MASS'
      include 'GEOM'
      include 'TSTEP'
      if(ifoutfld) call async_catalyst_update(
     &xm1,ym1,zm1,
     &pr,vx,vy,
     &vz,t,bm1)
      end
c-------------------------------------------
c       Finalization 
c-------------------------------------------
      subroutine async_catalyst_end()
      call async_catalyst_finalize()
      end

c-------------------------------------------
c       Copy for OPENMP with SOLN_ADIOS 
c-------------------------------------------
c      subroutine adios2_copy()
c      include 'SIZE'
c      include 'INPUT'
c
c      include 'SOLN'
c      include 'SOLN_ADIOS'
c      vx_adios=vx
c      vy_adios=vy
c      vz_adios=vz
c      t_adios=t
c      pr_adios=pr
c      end

