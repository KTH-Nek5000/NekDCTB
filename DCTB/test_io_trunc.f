!> @file io_trunc.f
!! @ingroup io_compression
!! @brief Set of I/O related tools for KTH modules
!! @author Adalberto Perez
!! @date Jan 2022
!=======================================================================
!> @brief Control what type of compression will be used. 
!! @ingroup io_compression
      subroutine trunc_test_static()
      implicit none
      include 'SIZE'    ! for size information such as lx1, etc ...
      include 'WZ'      ! for zgm1
      include 'INPUT'   ! for user inputs from par (I think) 
      include 'TSTEP'   ! istep
      include 'MASS'    ! for the mass matrix
      include 'SOLN'

      include 'IOTRUNCD' ! THIS IS MINE, include here the integer a

c     variables that contains the error and norms
      logical ifspec
      real l2err                          !l2 norm of error
      real l2errx                         !l2 norm of error
      real l2erry                         !l2 norm of error
      real l2errz                         !l2 norm of error
      real l2v                            !l2 norm of error
      real l2v_trc                        !l2 norm of error
      real errvec(lx1,ly1,lz1,lelv)       !error vector v_hat-vx
      real errvecx(lx1,ly1,lz1,lelv)      !error vector v_hat-vx
      real errvecy(lx1,ly1,lz1,lelv)      !error vector v_hat-vx
      real errvecz(lx1,ly1,lz1,lelv)      !error vector v_hat-vx
      real errvecpr(lx1,ly1,lz1,lelv)     !error vector v_hat-vx
      real truncvecx(lx1,ly1,lz1,lelv)    !error vector v_hat-vx
      real truncvecy(lx1,ly1,lz1,lelv)    !error vector v_hat-vx
      real truncvecz(lx1,ly1,lz1,lelv)    !error vector v_hat-vx
      real truncvecpr(lx1,ly1,lz1,lelv)   !error vector v_hat-vx
      real truncvec(lx1,ly1,lz1,lelv)     !error vector v_hat-vx
      real normv_specx(nelt), normv_physx(nelt), normv_errx(nelt)
      real normv_specy(nelt), normv_physy(nelt), normv_erry(nelt)
      real normv_specz(nelt), normv_physz(nelt), normv_errz(nelt)
      real normv_specpr(nelt), normv_physpr(nelt), normv_errpr(nelt)

      integer exec_mode, nxyze

      real gl2norm
      real glsum

      exec_mode=1

#ifdef AMR 
      exec_mode=5
#endif

      if (istep.lt.exec_mode) then
        return
      endif 
 
      ! Hard code compression parameters for the test
      nxyze=lx1*ly1*lz1*nelv
      numfile = 1 

      call compress_inputs()

      call decompress_inputs() 

c     calculate the global error norm
      call sub3(errvec,vx_trc,vx_trc_temp,nxyze)
      l2errx = gl2norm(errvec,nxyze)
      call sub3(errvec,vy_trc,vy_trc_temp,nxyze)
      l2erry = gl2norm(errvec,nxyze)
      call sub3(errvec,vz_trc,vz_trc_temp,nxyze)
      l2errz = gl2norm(errvec,nxyze)
      l2v = gl2norm(vz_trc,nxyze)
      l2v_trc = gl2norm(vz_trc_temp,nxyze)
c     write debugging values
      if (nio.eq.0) then
        write(*,*) 'Norm from the reader '
        write(*,*) 'l2 error in x= ', l2errx
        write(*,*) 'l2 error in y= ', l2erry
        write(*,*) 'l2 error in z= ', l2errz
        write(*,*) 'l2 norm vz new= ', l2v
        write(*,*) 'l2 norm vz_old= ', l2v_trc
      endif


!     finish simulation
      LASTEP=1


      return
      end subroutine

!=======================================================================
!> @brief Control what type of compression will be used. 
!! @ingroup io_compression
      subroutine trunc_test_insitu()
      implicit none
      include 'SIZE'    ! for size information such as lx1, etc ...
      include 'WZ'      ! for zgm1
      include 'INPUT'   ! for user inputs from par (I think) 
      include 'TSTEP'   ! istep
      include 'MASS'    ! for the mass matrix
      include 'SOLN'

      include 'IOTRUNCD' ! THIS IS MINE, include here the integer a

c     variables that contains the error and norms
      logical ifspec
      real l2err                          !l2 norm of error
      real l2errx                         !l2 norm of error
      real l2erry                         !l2 norm of error
      real l2errz                         !l2 norm of error
      real l2errpr                         !l2 norm of error
      real l2v                            !l2 norm of error
      real l2v_trc                        !l2 norm of error
      real errvec(lx1,ly1,lz1,lelv)       !error vector v_hat-vx
      real errvecx(lx1,ly1,lz1,lelv)      !error vector v_hat-vx
      real errvecy(lx1,ly1,lz1,lelv)      !error vector v_hat-vx
      real errvecz(lx1,ly1,lz1,lelv)      !error vector v_hat-vx
      real errvecpr(lx1,ly1,lz1,lelv)     !error vector v_hat-vx
      real truncvecx(lx1,ly1,lz1,lelv)    !error vector v_hat-vx
      real truncvecy(lx1,ly1,lz1,lelv)    !error vector v_hat-vx
      real truncvecz(lx1,ly1,lz1,lelv)    !error vector v_hat-vx
      real truncvecpr(lx1,ly1,lz1,lelv)   !error vector v_hat-vx
      real truncvec(lx1,ly1,lz1,lelv)     !error vector v_hat-vx
      real normv_specx(nelt), normv_physx(nelt), normv_errx(nelt)
      real normv_specy(nelt), normv_physy(nelt), normv_erry(nelt)
      real normv_specz(nelt), normv_physz(nelt), normv_errz(nelt)
      real normv_specpr(nelt), normv_physpr(nelt), normv_errpr(nelt)

      integer exec_mode, nxyze

      real gl2norm
      real glsum

      ! Hard code compression parameters for the test
      nxyze=lx1*ly1*lz1*nelv
      iotruncstep = 5 
      ioreadstep  = 6

      call trunc_data()

c     calculate the global error norm
      call sub3(errvec,vx_trc,vx_trc_temp,nxyze)
      l2errx = gl2norm(errvec,nxyze)
      call sub3(errvec,vy_trc,vy_trc_temp,nxyze)
      l2erry = gl2norm(errvec,nxyze)
      call sub3(errvec,vz_trc,vz_trc_temp,nxyze)
      l2errz = gl2norm(errvec,nxyze)
      call sub3(errvec,pr_trc,pr_trc_temp,lx2*lx2*lx2*nelv)
      l2errpr = gl2norm(errvec,nxyze)
      l2v = gl2norm(vz_trc,nxyze)
      l2v_trc = gl2norm(vz_trc_temp,nxyze)

c     write debugging values
      if (nio.eq.0) then
        write(*,*) 'Norm from the reader in z '
        write(*,*) 'l2 error in x= ', l2errx
        write(*,*) 'l2 error in y= ', l2erry
        write(*,*) 'l2 error in z= ', l2errz
        write(*,*) 'l2 error in pr= ', l2errpr
        write(*,*) 'l2 norm vz new= ', l2v
        write(*,*) 'l2 norm vz_old= ', l2v_trc
      endif

!     finish simulation
      if (istep.ge.29) LASTEP=1

      return
      end subroutine
