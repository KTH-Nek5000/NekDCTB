c-------------------------------------------
c       Initialization 
c-------------------------------------------
      subroutine adios2_fides_init()
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
      real pm1 (lx1,ly1,lz1,lelv), xm0(lx1,ly1,lz1,lelt),
     $     ym0(lx1,ly1,lz1,lelt) ! pressure mapping
      common /scrcg/ pm1,xm0,ym0 ! scratch array use
      if (ifsplit) then
            call adios2_fides_setup(
     &lx1,ly1,lz1,
     &nelv,nelb,nelgv,
     &iostep,if3d,
     &time,dt,
     &xm1,ym1,zm1,
     &pr,vx,vy,vz,
     &t,nekcomm)
      else
            call mappr(pm1,pr,xm0,ym0) ! map pressure onto Mesh 1
            call adios2_fides_setup(
     &lx1,ly1,lz1,
     &nelv,nelb,nelgv,
     &iostep,if3d,
     &time,dt,
     &xm1,ym1,zm1,
     &pm1,vx,vy,vz,
     &t,nekcomm)
      endif
      end

c-------------------------------------------
c       Write 
c-------------------------------------------
      subroutine adios2_fides_write()
      implicit none
      include 'SIZE'
      include 'INPUT'
      include 'SOLN'
      include 'MASS'
      include 'GEOM'
      include 'TSTEP'
      real pm1 (lx1,ly1,lz1,lelv), xm0(lx1,ly1,lz1,lelt),
     $     ym0(lx1,ly1,lz1,lelt) ! pressure mapping
      common /scrcg/ pm1,xm0,ym0 ! scratch array use
      if(ifoutfld) then
            if (ifsplit) then
                  call adios2_fides_update(
     &pr,vx,vy,
     &vz,t)
            else
                  call mappr(pm1,pr,xm0,ym0) ! map pressure onto Mesh 1
                  call adios2_fides_update(
     &pm1,vx,vy,
     &vz,t)
            endif
      endif
      end
c-------------------------------------------
c       Finalization 
c-------------------------------------------
      subroutine adios2_fides_end()
      call adios2_fides_finalize()
      end


