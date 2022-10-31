### \<casename\>.usr
The ```<casename>.usr``` serves as the driving script. It is here that the toolbox will be initialized, executed and finalized. The following bolck of code must be added to the ```userchk``` subroutine:
```fortran
!     start framework
      if (ISTEP.eq.0) call frame_start

!     monitor simulation
      call frame_monitor

      call trunc_main()

!     finalise framework
      if (ISTEP.eq.NSTEPS.or.LASTEP.eq.1) then
         call frame_end      
      endif
```

Aditionally, the following initialization subroutines must be copied into the ```<casename>.usr``` file:
```fortran
!======================================================================
!> @brief Register user specified modules
      subroutine frame_usr_register
      implicit none

      include 'SIZE'
      include 'FRAMELP'
!-----------------------------------------------------------------------
!     register modules
      call io_register
      call trunc_register

      return
      end subroutine
!======================================================================
!> @brief Initialise user specified modules
      subroutine frame_usr_init
      implicit none

      include 'SIZE'
      include 'FRAMELP'
      include 'SOLN'
!-----------------------------------------------------------------------
!     initialise modules
      call trunc_init

      return
      end subroutine
!======================================================================
!> @brief Finalise user specified modules
      subroutine frame_usr_end
      implicit none

      include 'SIZE'
      include 'FRAMELP'
!-----------------------------------------------------------------------
!     finalise modules
      
      return
      end subroutine
```
