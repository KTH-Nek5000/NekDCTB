!> @file amr_IO.f
!! @ingroup nekamr
!! @brief set of I/O routines for nekamr
!! @author Adam Peplinski
!! @date Mar 7, 2016
!=======================================================================
!> @brief Get free unit number
!! @param[out] iunit     file unit
!! @param[out] ierr      error mark
!! @note Identical routine in Nek_framework/io/io_tools/io_tools.f
      subroutine io_file_freeid(iunit, ierr)
      implicit none

      ! argument list
      integer iunit
      integer ierr

      ! keeep track of max iunit generated
      integer io_iunit_min, io_iunit_max
      common /io_iunit/ io_iunit_min, io_iunit_max

      ! local variables
      logical ifcnnd            ! is unit connected
!-----------------------------------------------------------------------
      ! initialise variables
      ierr=0
      iunit = io_iunit_min

      do
         inquire(unit=iunit,opened=ifcnnd,iostat=ierr)
         if(ifcnnd) then
            iunit = iunit +1
         else
            exit
         endif
      enddo

      if (iunit.gt.io_iunit_max) io_iunit_max = iunit

      return
      end subroutine
!=======================================================================
!> @brief Close opened files
!! @note Identical routine in Nek_framework/io/io_tools/io_tools.f
      subroutine io_file_close()
      implicit none

      ! keeep track of max iunit generated
      integer io_iunit_min, io_iunit_max
      common /io_iunit/ io_iunit_min, io_iunit_max

      ! local variables
      integer iunit, ierr
      logical ifcnnd            ! is unit connected
!-----------------------------------------------------------------------
      do iunit = io_iunit_min, io_iunit_max
         inquire(unit=iunit,opened=ifcnnd,iostat=ierr)
         if(ifcnnd) close(iunit)
      enddo
      io_iunit_max = io_iunit_min

      return
      end subroutine
!=======================================================================
!> @brief Save set of restart files
!! @details Tis is single-step restart routine. Single-step is preferred
!!    for AMR as this way one does not have to take care of refinement time.
!! @remarks This routine uses global scratch space CTMP1
      subroutine amr_restart_save
      implicit none

      include 'SIZE'
      include 'RESTART'
      include 'TSTEP'
      include 'INPUT'
      include 'SOLN'
      include 'AMR'

      ! global scratch arrays
      real tloc(lx1,ly1,lz1,lelt,ldimt)
      common /ctmp1/ tloc

      ! local variables
      character(len=AMR_LSTL_FNM) filename ! file name
      integer ifile, itmp, il

      ! simple timing
      real t1, tmp

      ! file name
      character*6  str
      character*3 prefixl    ! local prefix
      character*(*) kst
      parameter(kst='0123456789abcdefx')

      ! save old parameter values
      real ltime, lp63
      logical lif_full_pres, lifxyo, lifpo, lifvo, lifto,
     $        lifpsco(LDIMT1)

      ! to mark first call
      logical ifcalled
      save ifcalled
      data ifcalled /.true./

      ! functions
      real dnekclock
!-----------------------------------------------------------------------
      ! simple timing
      t1 = dnekclock()

      ! save current tree structure
      call amr_log(AMR_LP_PRD,'Saving restart data.')

      ! get number of snapshots in a set
      if (PARAM(27).lt.0) then
         amr_nsnap = NBDINP
      else
         amr_nsnap = 3
      endif

      if(ifcalled) then
        ! due to restrictions of native nek restart writing routines (rs? files counting)
        !AMR_IORSET = AMR_IOSTART
        AMR_IORSET = 0
        ifcalled=.false.
      endif

      ! prepare file name for ATF (AMR Tree File) file
      filename = 'ATF'//trim(adjustl(SESSION))
      if (len_trim(filename).gt.(AMR_LSTL_FNM-13)) then
         call amr_abort
     $        ('ERROR: amr_end; too long output file name')
      endif
      ifile = mod(AMR_IORSET,AMR_IOSET_MAX) + 1
      write(str,'(i5.5)') ifile
      filename = trim(filename)//'.t'//trim(str)//CHAR(0)
#ifdef P4EST
      call fp4est_tree_save(1,filename)
#endif

      if (ifmvbd) call amr_abort
     $       ('Error: AMR restart does not support moving boundarries')

      ! adjust I/O parameters
      lp63 = param(63)
      param(63)  = 1
      lif_full_pres = IF_FULL_PRES
      IF_FULL_PRES = .true.
      lifxyo = IFXYO
      lifpo= IFPO
      IFPO = .TRUE.
      lifvo= IFVO
      IFVO = .TRUE.
      lifto= IFTO
      IFTO = IFHEAT
      do il=1,NPSCAL
         lifpsco(il)= IFPSCO(il)
         IFPSCO(il) = .TRUE.
      enddo

      ! create prefix and name for DNS
      prefixl(1:2) = 'rs'
      itmp=min(17,amr_ioset_max*amr_nsnap) + 1
      prefixl(3:3)=kst(itmp:itmp)

      ! adjust time
      ltime=time

      if (amr_nsnap.eq.3) then
        ! set time
        time =ltime - dtlag(1) - dtlag(2)
        if (ifheat) then
          ! copy passive scalar data
          itmp = lx1*ly1*lz1*nelt
          do il=1,nfield-1
            call copy(tloc(1,1,1,1,il),tlag(1,1,1,1,lorder-1,il),itmp)
          enddo
          itmp=nfield-1
        else
          itmp=0
        endif
        call outpost2(vxlag(1,1,1,1,2),vylag(1,1,1,1,2),
     $    vzlag(1,1,1,1,2),AMR_PRLAG(1,1,1,1,2),tloc,itmp,prefixl)
      endif

      ! set time
      time =ltime - dtlag(1)
      if (ifheat) then
        ! copy passive scalar data
        itmp = lx1*ly1*lz1*nelt
        do il=1,nfield-1
          call copy(tloc(1,1,1,1,il),tlag(1,1,1,1,1,il),itmp)
        enddo
        itmp=nfield-1
      else
        itmp=0
      endif
      call outpost2(vxlag(1,1,1,1,1),vylag(1,1,1,1,1),vzlag(1,1,1,1,1),
     $    AMR_PRLAG(1,1,1,1,1),tloc,itmp,prefixl)

      ! set time
      time = ltime
      call outpost2(vx,vy,vz,pr,t,itmp,prefixl)

      ! restore I/O parameters
      param(63) = lp63
      IF_FULL_PRES = lif_full_pres
      IFXYO = lifxyo
      IFPO = lifpo
      IFVO = lifvo
      IFTO = lifto
      do il=1,NPSCAL
         IFPSCO(il) = lifpsco(il)
      enddo

      ! adjust file set number
      AMR_IORSET = AMR_IORSET + 1

      ! simple timing
      tmp = dnekclock() - t1
      ! total
      AMR_TC = AMR_TC + tmp
      ! finalisation
      AMR_TCRS = AMR_TCRS + tmp

      return
      end subroutine
!=======================================================================
!> @brief Read set of restart files
!! @details Tis is single-step restart routine. Single-step is preferred
!!    for AMR as this way one does not have to take care of refinement time.
      subroutine amr_restart_read()
      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'TSTEP'
      include 'SOLN'
      include 'AMR'

      ! local variables
      character(len=AMR_LSTL_FNM) filename ! file name
      integer ifile, itmp, il, jl

      ! file name
      character*6  str
      character*3 prefixl    ! local prefix
      character*(*) kst
      parameter(kst='0123456789abcdefx')

      ! simple timing
      real t1, tmp

      ! functions
      real dnekclock
!-----------------------------------------------------------------------
      ! simple timing
      t1 = dnekclock()

      call amr_log(AMR_LP_PRD,'Reading restart data.')

      ! get number of snapshots in a set
      if (PARAM(27).lt.0) then
         amr_nsnap = NBDINP
      else
         amr_nsnap = 3
      endif

      ! create prefix and name for DNS
      prefixl(1:2) = 'rs'
      itmp=min(17,amr_ioset_max*amr_nsnap) + 1
      prefixl(3:3)=kst(itmp:itmp)

      filename = prefixl//trim(adjustl(SESSION))
      if (len_trim(filename).gt.(AMR_LSTL_FNM-13)) then
         call amr_abort
     $        ('ERROR: amr_end; too long output file name')
      endif

      call blank(initc,132*15)
      param(67) = 6.00

      ! read fileds
      do il=1,amr_nsnap-1
         ifile = (AMR_IOSTART-1)*amr_nsnap + il
         write(str,'(i5.5)') ifile
         initc(1) = trim(filename)//'0.f'//trim(str)!//CHAR(0)
         call restart(1)
         ! copy data
         AMR_CHK_TIME(amr_nsnap-il) = time
         itmp=lx1*ly1*lz1*nelv
         call copy(vxlag(1,1,1,1,amr_nsnap-il),vx,itmp)
         call copy(vylag(1,1,1,1,amr_nsnap-il),vy,itmp)
         if (IF3D) call copy(vzlag(1,1,1,1,amr_nsnap-il),vz,itmp)
         itmp=lx2*ly2*lz2*nelv
         call copy(AMR_PRLAG(1,1,1,1,amr_nsnap-il),pr,itmp)
         if (ifheat) then
           itmp = lx1*ly1*lz1*nelt
           do jl=1,nfield-1
             call copy(tlag(1,1,1,1,amr_nsnap-il,jl),t(1,1,1,1,jl),itmp)
           enddo
         endif
      enddo

      ! this flag is not reset, so I do it here
      ! remember is should be reset after the last read as well
      if_full_pres = .false.

      ifile = AMR_IOSTART*amr_nsnap
      write(str,'(i5.5)') ifile
      initc(1) = trim(filename)//'0.f'//trim(str)!//CHAR(0)

      ! simple timing
      tmp = dnekclock() - t1
      ! total
      AMR_TC = AMR_TC + tmp
      ! reading
      AMR_TCRS = AMR_TCRS + tmp

      return
      end  subroutine
!=======================================================================

