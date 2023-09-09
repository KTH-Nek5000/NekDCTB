!> @file io_trunc.f
!! @ingroup io_compression
!! @brief Set of I/O related tools for KTH modules
!! @author Adalberto Perez
!! @date Jan 2022
!=======================================================================
!> @brief Control what type of compression will be used. 
!! @ingroup io_compression
      subroutine stream_main()
      implicit none
      include 'SIZE'    ! for size information such as lx1, etc ...
      include 'WZ'      ! for zgm1
      include 'INPUT'   ! for user inputs from par (I think) 
      include 'TSTEP'   ! istep

      include 'IOSTREAMD' ! THIS IS MINE, include here the integer a

      ! You might need to start in timestep 1 because the intialization
      ! is done after user check is called in step 0. So maybe fix that 
      ! before proceeding to do anything regarding post processing
      if (istep.eq.1) then
              if (ifsstream)   call stream_inputs()
              call sleep(10) 
      endif 

      if (ifinsitustream) call stream_data() 

      return
      end subroutine
!
!=======================================================================
!> @brief Register compression module 
!! @ingroup io_compression
!! @note This subroutine should be called in frame_usr_register
      subroutine stream_register()
      implicit none
      include 'SIZE'    ! for size information such as lx1, etc ...
      include 'FRAMELP'
      include 'WZ'      ! for zgm1
      include 'INPUT'   ! for user inputs from par (I think) 
      include 'TSTEP'   ! istep

      include 'IOSTREAMD' ! THIS IS MINE, include here the integer a

      ! local variables
      integer lpmid, il
      real ltim
      character*2 str

      ! functions
      real dnekclock

!-----------------------------------------------------------------------
      ! timing
      ltim = dnekclock()

      ! check if the current module was already registered
      call mntr_mod_is_name_reg(lpmid,iostream_name)
      if (lpmid.gt.0) then
         call mntr_warn(lpmid,
     $        'module ['//trim(iostream_name)//'] already registered')
         return
      endif
      
      ! find parent module
      call mntr_mod_is_name_reg(lpmid,'FRAME')
      if (lpmid.le.0) then
         lpmid = 1
         call mntr_abort(lpmid,
     $        'parent module ['//'FRAME'//'] not registered')
      endif
      
      ! register module
      call mntr_mod_reg(iostream_id,lpmid,iostream_name,
     $      'IO streaming to python')

      ! register and set active section
      call rprm_sec_reg(iostream_sec_id,iostream_id,
     $     '_'//adjustl(iostream_name),
     $     'Runtime paramere section for io streaming module')
      call rprm_sec_set_act(.true.,iostream_sec_id)

      ! register parameters
      call rprm_rp_reg(filetostream_id,iostream_sec_id,'FILETOSTREAM',
     $     'Filename to stream',rpar_str,10,0.0,.false.,' ')

      call rprm_rp_reg(ifsstream_id,iostream_sec_id,'SSTREAM',
     $     'Stream in post prc mode',rpar_log,10,0.0,.false.,' ')

      call rprm_rp_reg(ifinsitustream_id,iostream_sec_id,
     $     'INSITUSTREAM',
     $     'IN SITU streaming',rpar_log,10,0.0,.false.,' ')

      call rprm_rp_reg(numfile_id,iostream_sec_id,'NUMFILE',
     $     'Number of files to stream',rpar_int,10,0.0,.false.,' ')
 
      call rprm_rp_reg(iostreamstep_id,iostream_sec_id,'STREAMSTEP',
     $     'streaming and output interval',rpar_int,0,0.0,.false.,' ')


      ! set initialisation flag
      iostream_ifinit=.false.

      ! timing
      ltim = dnekclock() - ltim

      return
      end subroutine
!
!=======================================================================
!> @brief Initialize the compression module 
!! @ingroup io_compression
!! @note This subroutine should be called in frame_usr_init
      subroutine stream_init()
      implicit none
      include 'SIZE'    ! for size information such as lx1, etc ...
      include 'FRAMELP'
      include 'WZ'      ! for zgm1
      include 'INPUT'   ! for user inputs from par (I think) 
      include 'TSTEP'   ! istep

      include 'IOSTREAMD' ! THIS IS MINE, include here the integer a

      ! local variables
      integer itmp, il
      real rtmp, ltim
      logical ltmp
      character*20 ctmp

      ! functions
      real dnekclock
!-----------------------------------------------------------------------
      ! check if the module was already initialised
      if (iostream_ifinit) then
         call mntr_warn(iostream_id,
     $        'module ['//trim(iostream_name)//'] already initiaised.')
         return
      endif

      ! timing
      ltim = dnekclock()
      
      ! get runtime parameters
      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,filetostream_id,rpar_str)
      filetostream = ctmp

      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,ifsstream_id,rpar_log)
      ifsstream = ltmp

      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,ifinsitustream_id,rpar_log)
      ifinsitustream = ltmp

      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,numfile_id,rpar_int)
      numfile = itmp

      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,iostreamstep_id,rpar_int)
      iostreamstep = itmp

      write(*,*) 'file to stream is ', filetostream

      !Initialize IO data
      call io_init()

      ! everything is initialised
      iostream_ifinit=.true.

      ! timing
      ltim = dnekclock() - ltim
      
      return
      end subroutine
!
!=======================================================================
!> @brief Compress xxxx.fxxx files already written from nek 
!! @ingroup io_compression
      subroutine stream_inputs()
      implicit none
      include 'SIZE'    ! for size information such as lx1, etc ...
      include 'WZ'      ! for zgm1
      include 'INPUT'   ! for user inputs from par (I think) 
      include 'TSTEP'   ! istep
      include 'SOLN'    ! vx,vy,vz,pr
      include 'PARALLEL'! for nelgv, nelgt
      include 'RESTART' ! nelb (for io operations)

      include 'MASS'      ! for bm1
      include 'IOSTREAMD' ! THIS IS MINE, include here the integer a

c     Declare variables -------------------

c     counter for printing debbuging info
      integer i,e

c     Set constants that determine allocation space. Max poly order: 40
      integer lm,lm2
      parameter (lm=40)
      parameter (lm2=lm*lm)

c     actual polynomial order
      integer nx
      integer nx_pr    !this is the order of the pressure
c     total number of points in the mpirank     
      integer nxyze
c     Logical variables
      logical streamio

c     For the reader
      integer nvals,nid_,np_,nekcomm
      integer ifile
      real pm1(lx1,ly1,lz1,lelv) !Create a pointer that you will reference to the common

c     for the communicator
      common /nekmpi/ nid_,np_,nekcomm 
      common /scrcg/ pm1

c     Declare functions  --------------------

      real gl2norm
      real glsum

c     Reader options
      character*132  fname
      character*132  fname1
      character*6  str
      character*6  str2

c     Set up parameters  --------------------

c     polynomial order, opposed to lm which is used for allocation
      nx=lx1     
      nx_pr=lx2 
c     Total number of points in the mpirank
      nxyze=lx1*ly1*lz1*nelv
      nvals=lx1*lx1*lx1

c     Compress the files recursively 
      do ifile=1,numfile

c       compose the name from the base in the input and the ifile        
        fname1=trim(filetostream)
        write (str,'(i5.5)') ifile
        fname=trim(fname1)//trim(str)

c       read the file 
c       trunc_mfi was designed to read the fields of the trunctated
c       legendre coefficients. Here I am using it to read the actual
c       velocity fields, thus, I must transfer the data from those
c       vectors to the velocity fields in order to start the streaming
c       process. This can be changed in the future.
        call stream_mfi(fname,1) 
c        call copy(vx,vx_hat_trc,nxyze)
c        call copy(vy,vy_hat_trc,nxyze)
c        call copy(vz,vz_hat_trc,nxyze)
c        call copy(pr,pr_hat_trc,nx_pr*nx_pr*nx_pr*nelv)

c       copy the pressure into an array full of zeros (without mapping)
c       in the velocity mesh to prepare for writing
c        call rzero(pm1,lx1*ly1*lz1*lelt)
c        do e=1,nelv
c              call copy  (pm1(1,1,1,e),pr(1,1,1,e),lx2*lx2*lx2)
c        enddo

c       ===============================================
c       Send data to python to be post processed
c       ===============================================
c        if (nid.eq.1) then
c                write (*,*) 'Data for current step'
c                do i=1,nxyze
c                    write(*,*) 'vx= ', vx(i,1,1,1)
c                enddo
c        endif



        !Write the fields into a file
        call adios2_stream(lglel,pm1, vx, 
     &   vy, vz, bm1, t)

        enddo

      return
      end subroutine
!
!=======================================================================
!> @brief In-Situ compression and decompression 
!! @ingroup io_compression
!! @note The read functionality has no real use in a normal run ... it
!!       should be used for debugging purposes mostly.
      subroutine stream_data()
      implicit none
      include 'SIZE'    ! for size information such as lx1, etc ...
      include 'WZ'      ! for zgm1
      include 'INPUT'   ! for user inputs from par (I think) 
      include 'TSTEP'   ! istep
      include 'SOLN'    ! vx,vy,vz,pr
      include 'PARALLEL'! for nelgv, nelgt
      include 'RESTART' ! nelb (for io operations)


      include 'MASS'      ! for bm1
      include 'IOSTREAMD' ! THIS IS MINE, include here the integer a

c     Declare variables -------------------

c     counter for printing debbuging info
      integer i,e

c     Set constants that determine allocation space. Max poly order: 40
      integer lm,lm2
      parameter (lm=40)
      parameter (lm2=lm*lm)

c     actual polynomial order
      integer nx
      integer nx_pr    !this is the order of the pressure
c     total number of points in the mpirank     
      integer nxyze

      logical streamio

c     For the reader
      integer nvals,nid_,np_,nekcomm
      real pm1(lx1,ly1,lz1,lelv) !Create a pointer that you will reference to the common

c     for the communicator
      common /nekmpi/ nid_,np_,nekcomm 
      common /scrcg/ pm1

c     Declare functions  --------------------

      real gl2norm
      real glsum
      !integer eg
      integer vid(2,nelt)

c     Reader options
      character*132  fname
      character*132  fname1
      character*6  str
      character*6  str2
      
c     Set up parameters  --------------------

c     polynomial order, opposed to lm which is used for allocation
      nx=lx1     
      nx_pr=lx2 
c     Total number of points in the mpirank
      nxyze=lx1*ly1*lz1*nelv
      nvals=lx1*lx1*lx1

c     Set up parameters. 
      streamio=.false.
      if(mod(istep,iostreamstep).eq.0.and.istep.ne.0) streamio=.true.
      !if (istep.eq.0) ifile=0

      !Initialize I/O data
      call io_init

c     =================
c     Data Stream
c     =================

c     Apply the operators as needed
      if (streamio) then

c       Copy the pressure into an array full of zeros in the vel mesh    
c       in preparation to write the files
        call rzero(pm1,lx1*ly1*lz1*lelt)
        do e=1,nelv
              call copy  (pm1(1,1,1,e),pr(1,1,1,e),lx2*lx2*lx2)
        enddo

c       =============================================
c       Perform the Streaming with ADIOS2
c       =============================================

        !pm1 is the truncated coefficients mapped into velocity mesh
        call adios2_stream(lglel,pm1, vx, 
     &   vy, vz, bm1, t)

      endif    

      return
      end subroutine
!
!=======================================================================


c     Blank Space. Subroutines after this should be in another file.



!=======================================================================
!> @brief Read xxx.fxxx files into _trunc vectors 
!! @ingroup io_compression
      subroutine stream_mfi(fname_in,ifile)
c
c     (1) Open restart file(s)
c     (2) Check previous spatial discretization 
c     (3) Map (K1,N1) => (K2,N2) if necessary
c
c     nfiles > 1 has several implications:
c
c     i.   For std. run, data is taken from last file in list, unless
c          explicitly specified in argument list of filename
c
c     ii.  For MHD and perturbation cases, 1st file is for U,P,T;
c          subsequent files are for B-field or perturbation fields
c
c
      include 'SIZE'
      include 'TOTAL'
      include 'RESTART'
      include 'IOSTREAMD'

      character*132  hdr
      character*132  fname_in

      character*132  fname
      character*1    fnam1(132)
      equivalence   (fnam1,fname)

      character*1    frontc

      parameter (lwk = 7*lx1*ly1*lz1*lelt)
      common /scrns/ wk(lwk)
      common /scrcg/ pm1(lx1*ly1*lz1,lelv)
      integer e
      logical ifanyc

      integer*8 offs0,offs,nbyte,stride,strideB,nxyzr8


c     explicitly set the I/O flags. In the future implement sioflags()
         IFGETX=.TRUE.
         IF (IF3D) IFGETZ=.TRUE.
         IFANYC=.FALSE.
         DO 400 I=1,NFIELD
            IF (IFADVC(I)) IFANYC=.TRUE.
  400    CONTINUE
         IF (IFFLOW.OR.IFANYC) THEN
            IFGETU=.TRUE.
            IF (IF3D) IFGETW=.TRUE.
         ENDIF
         if (ifflow) ifgetp=.true.
         if (ifheat) ifgett=.true.


      tiostart=dnekclock()

      ! add full path if required
      call blank(fname,132)
      call chcopy(frontc, fname_in, 1)
      if (frontc .ne. '/') then
        lenp = ltrunc(path,132)
        lenf = ltrunc(fname_in,132)
        call chcopy(fnam1(1),path,lenp)
        call chcopy(fnam1(lenp+1),fname_in,lenf)
      else
        lenf = ltrunc(fname_in,132)
        call chcopy(fnam1(1),fname_in,lenf)     
      endif

      call mfi_prepare(fname)       ! determine reader nodes +
                                    ! read hdr + element mapping 

      offs0   = iHeadersize + 4 + isize*nelgr
      nxyzr8  = nxr*nyr*nzr
      strideB = nelBr* nxyzr8*wdsizr
      stride  = nelgr* nxyzr8*wdsizr

      iofldsr = 0

      if (ifgetxr) then      ! if available
         offs = offs0 + ldim*strideB
         call byte_set_view(offs,ifh_mbyte)
         if (ifgetx) then
            if(nid.eq.0) write(6,*) 'Reading mesh'
            call mfi_getv(xm1,ym1,zm1,wk,lwk,.false.)
         else                ! skip the data
            call mfi_getv(xm1,ym1,zm1,wk,lwk,.true.)
         endif
         iofldsr = iofldsr + ldim
      endif

      if (ifgetur) then
         offs = offs0 + iofldsr*stride + ldim*strideB
         call byte_set_view(offs,ifh_mbyte)
         if (ifgetu) then
            if (ifmhd.and.ifile.eq.2) then
c               if(nid.eq.0) write(6,*) 'Reading B field'
               call mfi_getv(bx,by,bz,wk,lwk,.false.)
            else
               if(nid.eq.0) write(6,*) 'Reading velocity field'
c               call mfi_getv(vx_hat_trc,vy_hat_trc,vz_hat_trc,
c     &                       wk,lwk,.false.)
               call mfi_getv(vx,vy,vz,
     &                       wk,lwk,.false.)
            endif
         else
c            call mfi_getv(vx_hat_trc,vy_hat_trc,vz_hat_trc,
c     &                    wk,lwk,.true.)
            call mfi_getv(vx,vy,vz,
     &                    wk,lwk,.true.)
         endif
         iofldsr = iofldsr + ldim
      endif

      !Leave pm1 here, the change is in the interpolation step
      if (ifgetpr) then
         offs = offs0 + iofldsr*stride + strideB
         call byte_set_view(offs,ifh_mbyte)
         if (ifgetp) then
c           if(nid.eq.0) write(6,*) 'Reading pressure field'
            call mfi_gets(pm1,wk,lwk,.false.)
         else
            call mfi_gets(pm1,wk,lwk,.true.)
         endif
         iofldsr = iofldsr + 1
      endif
      if (ifgettr) then
         offs = offs0 + iofldsr*stride + strideB
         call byte_set_view(offs,ifh_mbyte)
         if (ifgett) then
c            if(nid.eq.0) write(6,*) 'Reading temperature field'
            call mfi_gets(t,wk,lwk,.false.)
         else
            call mfi_gets(t,wk,lwk,.true.)
         endif
         iofldsr = iofldsr + 1
      endif

      ierr = 0
      do k=1,ldimt-1
         if (ifgtpsr(k)) then
            offs = offs0 + iofldsr*stride + strideB
            call byte_set_view(offs,ifh_mbyte)
            if (ifgtps(k)) then
c               if(nid.eq.0) write(6,'(A,I2,A)') ' Reading ps',k,' field'
               call mfi_gets(t(1,1,1,1,k+1),wk,lwk,.false.)
            else
               call mfi_gets(t(1,1,1,1,k+1),wk,lwk,.true.)
            endif
            iofldsr = iofldsr + 1
         endif
      enddo
      nbyte = 0
      if(nid.eq.pid0r) nbyte = iofldsr*nelr*wdsizr*nxr*nyr*nzr

      if (ifgtim) time = timer

      if(ifmpiio) then
        if(nid.eq.pid0r) call byte_close_mpi(ifh_mbyte,ierr)
      else
        if(nid.eq.pid0r) call byte_close(ierr)
      endif
      call err_chk(ierr,'Error closing restart file, in mfi.$')
      tio = dnekclock()-tiostart

      dnbyte = nbyte
      nbyte = glsum(dnbyte,1)
      nbyte = nbyte + iHeaderSize + 4 + isize*nelgr

      if (tio.eq.0) tio=1
      if (nio.eq.0) write(6,7) istep,time,
     &             nbyte/tio/1024/1024/10,
     &             nfiler
    7 format(/,i9,1pe12.4,' done :: Read checkpoint data',/,
     &       30X,'avg data-throughput = ',f7.1,'MBps',/,
     &       30X,'io-nodes = ',i5,/)


      if (ifaxis) call axis_interp_ic(pm1)            ! Interpolate to axi mesh
      if (ifgetp) call stream_map_pm1_to_pr(pm1,ifile) ! Interpolate pressure

      return
      end
!=======================================================================
!> @brief map pm1 into pr_hat_trc
!! @ingroup io_compression
!! @note This is only useful if we have written the truncated data into
!a xxx.fxxx file and want to read it. This was used at the intial stages
!of the project to verify compression but now the reads are made with
!ADIOS2 instead.

      subroutine stream_map_pm1_to_pr(pm1,ifile)

      include 'SIZE'
      include 'TOTAL'
      include 'RESTART'
      include 'IOSTREAMD'

      real pm1(lx1*ly1*lz1,lelv)
      integer e

      nxyz2 = lx2*ly2*lz2

      if (ifmhd.and.ifile.eq.2) then
         do e=1,nelv
            if (if_full_pres) then
               call copy  (pm(1,1,1,e),pm1(1,e),nxyz2)
            else
               call map12 (pm(1,1,1,e),pm1(1,e),e)
            endif
         enddo
      elseif (ifsplit) then
         call copy (pr,pm1,lx1*ly1*lz1*nelv)
      else
         do e=1,nelv
            if (if_full_pres) then
               call copy  (pr(1,1,1,e),pm1(1,e),nxyz2)
            else
               call map12 (pr(1,1,1,e),pm1(1,e),e)
            endif
         enddo
      endif
   
      return
      end

!=======================================================================
