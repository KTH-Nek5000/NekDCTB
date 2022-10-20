!> @file io_trunc.f
!! @ingroup io_compression
!! @brief Set of I/O related tools for KTH modules
!! @author Adalberto Perez
!! @date Jan 2022
!=======================================================================
!> @brief Control what type of compression will be used. 
!! @ingroup io_compression
      subroutine trunc_main()
      implicit none
      include 'SIZE'    ! for size information such as lx1, etc ...
      include 'WZ'      ! for zgm1
      include 'INPUT'   ! for user inputs from par (I think) 
      include 'TSTEP'   ! istep
      include 'MASS'    ! for the mass matrix
      include 'SOLN'

      include 'IOTRUNCD' ! THIS IS MINE, include here the integer a

      ! You might need to start in timestep 1 because the intialization
      ! is done after user check is called in step 0. So maybe fix that 
      ! before proceeding to do anything regarding post processing
      if (istep.eq.1) then
              if (ifscompress)   call compress_inputs()
              call sleep(10) 
              if (ifsdecompress) call decompress_inputs() 
              call sleep(10)
	      !Print mass matrix
	      call outpost(BM1,BM1,BM1,pr,t,'bm1') 
      endif 

      if (ifinsitucompress) call trunc_data() 

      return
      end subroutine
!
!=======================================================================
!> @brief Register compression module 
!! @ingroup io_compression
!! @note This subroutine should be called in frame_usr_register
      subroutine trunc_register()
      implicit none
      include 'SIZE'    ! for size information such as lx1, etc ...
      include 'FRAMELP'
      include 'WZ'      ! for zgm1
      include 'INPUT'   ! for user inputs from par (I think) 
      include 'TSTEP'   ! istep

      include 'IOTRUNCD' ! THIS IS MINE, include here the integer a

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
      call mntr_mod_is_name_reg(lpmid,iotrunc_name)
      if (lpmid.gt.0) then
         call mntr_warn(lpmid,
     $        'module ['//trim(iotrunc_name)//'] already registered')
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
      call mntr_mod_reg(iotrunc_id,lpmid,iotrunc_name,
     $      'IO truncation and compression')

      ! register and set active section
      call rprm_sec_reg(iotrunc_sec_id,iotrunc_id,
     $     '_'//adjustl(iotrunc_name),
     $     'Runtime paramere section for io truncation module')
      call rprm_sec_set_act(.true.,iotrunc_sec_id)

      ! register parameters
      call rprm_rp_reg(filetocomp_id,iotrunc_sec_id,'FILETOCOMP',
     $     'Filename to compress',rpar_str,10,0.0,.false.,' ')

      call rprm_rp_reg(ifscompress_id,iotrunc_sec_id,'SCOMP',
     $     'Compress in post prc mode',rpar_log,10,0.0,.false.,' ')

      call rprm_rp_reg(ifsdecompress_id,iotrunc_sec_id,
     $     'SDECOMP',
     $     'DeCompress in post prc mode',rpar_log,10,0.0,.false.,' ')

      call rprm_rp_reg(ifinsitucompress_id,iotrunc_sec_id,
     $     'INSITUCOMP',
     $     'IN SITU compression',rpar_log,10,0.0,.false.,' ')

      call rprm_rp_reg(numfile_id,iotrunc_sec_id,'NUMFILE',
     $     'Number of files to compress',rpar_int,10,0.0,.false.,' ')
 
      call rprm_rp_reg(iotruncstep_id,iotrunc_sec_id,'TRUNCSTEP',
     $     'Truncation and output interval',rpar_int,0,0.0,.false.,' ')

      call rprm_rp_reg(ioreadstep_id,iotrunc_sec_id,'READSTEP',
     $     'Read files interval',rpar_int,0,0.0,.false.,' ')

      call rprm_rp_reg(targeterr_id,iotrunc_sec_id,'TARGETERR',
     $     'Target error for compression',rpar_real,0,0.0,.false.,' ')


      ! set initialisation flag
      iotrunc_ifinit=.false.

      ! timing
      ltim = dnekclock() - ltim

      return
      end subroutine
!
!=======================================================================
!> @brief Initialize the compression module 
!! @ingroup io_compression
!! @note This subroutine should be called in frame_usr_init
      subroutine trunc_init()
      implicit none
      include 'SIZE'    ! for size information such as lx1, etc ...
      include 'FRAMELP'
      include 'WZ'      ! for zgm1
      include 'INPUT'   ! for user inputs from par (I think) 
      include 'TSTEP'   ! istep

      include 'IOTRUNCD' ! THIS IS MINE, include here the integer a

      ! local variables
      integer itmp, il
      real rtmp, ltim
      logical ltmp
      character*20 ctmp

      ! functions
      real dnekclock
!-----------------------------------------------------------------------
      ! check if the module was already initialised
      if (iotrunc_ifinit) then
         call mntr_warn(iotrunc_id,
     $        'module ['//trim(iotrunc_name)//'] already initiaised.')
         return
      endif

      ! timing
      ltim = dnekclock()
      
      ! get runtime parameters
      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,filetocomp_id,rpar_str)
      filetocomp = ctmp

      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,ifscompress_id,rpar_log)
      ifscompress = ltmp

      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,ifsdecompress_id,rpar_log)
      ifsdecompress = ltmp

      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,ifinsitucompress_id,rpar_log)
      ifinsitucompress = ltmp

      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,numfile_id,rpar_int)
      numfile = itmp

      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,iotruncstep_id,rpar_int)
      iotruncstep = itmp

      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,ioreadstep_id,rpar_int)
      ioreadstep = itmp

      call rprm_rp_get(itmp,rtmp,ltmp,ctmp,targeterr_id,rpar_real)
      targeterr = rtmp

      if (nio.eq.0) then
        write(*,*) 'file to comp is ', filetocomp
      endif

      !Initialize IO data
      call io_init()

      ! everything is initialised
      iotrunc_ifinit=.true.

      ! timing
      ltim = dnekclock() - ltim
      
      return
      end subroutine
!
!=======================================================================
!> @brief Compress xxxx.fxxx files already written from nek 
!! @ingroup io_compression
      subroutine compress_inputs()
      implicit none
      include 'SIZE'    ! for size information such as lx1, etc ...
      include 'WZ'      ! for zgm1
      include 'INPUT'   ! for user inputs from par (I think) 
      include 'TSTEP'   ! istep
      include 'SOLN'    ! vx,vy,vz,pr
      include 'PARALLEL'! for nelgv, nelgt
      include 'RESTART' ! nelb (for io operations)
      include 'CTIMER'  ! clock


      include 'IOTRUNCD' ! THIS IS MINE, include here the integer a

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
c     filter transfer function
      real trunc_filter(lm2)                    !filter array
c     filter cutout frequency
      integer kut
c     cut out as a function of the elements      
      integer kutvx(nelt),kutvy(nelt),kutvz(nelt), kutvpr(nelt)
c     compression ratio as a function of elements      
      real crsx(nelt),crsy(nelt),crsz(nelt), crspr(nelt)
c     Logical variables
      logical truncio
      logical readio

c     Variable to determine the type of basis function for transforms
      logical trunc_ifboyd
      logical ifsort

c     average kut
      real avrgkutvx,avrgkutvy, avrgkutvz, avrgkutvpr
      real avrgcrsx,avrgcrsy, avrgcrsz, avrgcrspr
      real numel

c     variables that contains the error and norms
      logical ifspec
      real l2err                          !l2 norm of error
      real l2v                          !l2 norm of error
      real l2v_trc                          !l2 norm of error
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

c     spectral coefficient matrices for velocity and pressure
      real spec_trans_v(lm2)        ! Spectral transform matrix 
      real spec_trans_vt(lm2)       ! Spectral transform matrix V^T 
      real spec_trans_vinv(lm2)     ! Spectral transform matrix V^-1 
      real spec_trans_v_pr(lm2)     ! Spectral transform matrix V 
      real spec_trans_vt_pr(lm2)    ! Spectral transform matrix V^T 
      real spec_trans_vinv_pr(lm2)  ! Spectral transform matrix V^-1 

      real*8 etime_t

c     For the reader
      integer nvals,nid_,np_,nekcomm
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

c     What basis function for transform: .false. is legendre
      trunc_ifboyd = .false.
c     polynomial order, opposed to lm which is used for allocation
      nx=lx1     
      nx_pr=lx2 
c     Total number of points in the mpirank
      nxyze=lx1*ly1*lz1*nelv
      nvals=lx1*lx1*lx1

c     Set up parameters
      ifsort=.true.
      ifile = 1

c     Set up some io parameters      
      call blank(initc,132*15)
      param(67) = 6.00

c     Compress the files recursively 
      do ifile=1,numfile

        etime0 = dnekclock_sync()
c       compose the name from the base in the input and the ifile        
        fname1=trim(filetocomp)
        write (str,'(i5.5)') ifile
        fname=trim(fname1)//trim(str)
        initc(1)=trim(fname)


c       read the file 
c       trunc_mfi was designed to read the fields of the trunctated
c       legendre coefficients. Here I am using it to read the actual
c       velocity fields, thus, I must transfer the data from those
c       vectors to the velocity fields in order to start the truncation
c       process. This can be changed in the future.
c        call trunc_mfi(fname,1) 
c        call copy(vx,vx_hat_trc,nxyze)
c        call copy(vy,vy_hat_trc,nxyze)
c        call copy(vz,vz_hat_trc,nxyze)
c        call copy(pr,pr_hat_trc,nx_pr*nx_pr*nx_pr*nelv)

        call restart(1)

        etime_t = dnekclock_sync() - etime0
        if(nio.eq.0) write(6,'(A,1(1g9.2),A,/)')
     &                   ' done :: read data from disk   ',
     &                    etime_t, ' sec'



c       write the reconstructed fields.        
        call outpost(vx,vy,
     $     vz,pr,t,'rdd')
        if(nio.eq.0) write(6,'(A,1(1g9.2),A,/)')
     &                   ' done :: print read data   ',
     &                    etime_t, ' sec'


c       =======================        
c       Now do the truncation
c       =======================

        etime0 = dnekclock_sync()
c       Build operator each time we write... consider putting this out
        call build_spec_matrix(trunc_ifboyd,nx,spec_trans_v,
     &                         spec_trans_vt,spec_trans_vinv)
        call build_spec_matrix(trunc_ifboyd,nx_pr,spec_trans_v_pr,
     &                         spec_trans_vt_pr,spec_trans_vinv_pr)

        if (nio.eq.0) then
          write(*,*) 'Spectral matrix - V'
          do i=1,lx1*lx1
            write(*,*) 'V(',i,')=', spec_trans_v(i)
          enddo 
        endif

        if (nio.eq.0) then
          write(*,*) 'Spectral matrix - Vinv'
          do i=1,lx1*lx1
            write(*,*) 'V(',i,')=', spec_trans_vinv(i)
          enddo 
        endif
c       transform the field into spectral space u_hat=V^-1*u
cc       For vx
c        call trunc_apply_op(vx_hat,vx,spec_trans_vinv,if3d,nx)
cc       For vy
c        call trunc_apply_op(vy_hat,vy,spec_trans_vinv,if3d,nx)
c       For vz
        call trunc_apply_op(vz_hat,vz,spec_trans_vinv,if3d,nx)
cc       For Pressure
c        call trunc_apply_op(pr_hat,pr,spec_trans_vinv_pr,if3d,nx_pr)

        if (nio.eq.0) then
          write(*,*) 'Field - original'
          do i=1,10
            write(*,*) 'vz = ', vz(i,1,1,60),
     &                 'vzhat =', vz_hat(i,1,1,60)
          enddo 
        endif


c       Truncate the field according to an input error        
c       For vx
c        call trunc_field_l2oferror(vx_trc,vx_hat_trc,normv_physx,
c     &       normv_specx,kutvx,vx_hat,spec_trans_v,spec_trans_vt,
c     &       nx,targeterr,crsx,ifsort)
cc       For vy
c        call trunc_field_l2oferror(vy_trc,vy_hat_trc,normv_physy,
c     &       normv_specy,kutvy,vy_hat,spec_trans_v,spec_trans_vt,
c     &       nx,targeterr,crsy,ifsort)
c       For vz
        call trunc_field_l2oferror(vz_trc,vz_hat_trc,normv_physz,
     &       normv_specz,kutvz,vz_hat,spec_trans_v,spec_trans_vt,
     &       nx,targeterr,crsz,ifsort)
cc       For pr
c        call trunc_field_l2oferror(pr_trc,pr_hat_trc,normv_physpr,
c     &       normv_specpr,kutvpr,pr_hat,spec_trans_v_pr,
c     &       spec_trans_vt_pr, nx_pr,targeterr,crspr,ifsort)

c       write the reconstructed fields.        
        call outpost(vx_trc,vy_trc,
     $     vz_trc,pr_trc,t,'rc1')
        if(nio.eq.0) write(6,'(A,1(1g9.2),A,/)')
     &                   ' done :: print truncated data   ',
     &                    etime_t, ' sec'

        if (nio.eq.0) then
          write(*,*) 'Field - after truncation'
          do i=1,10
            write(*,*) 'vz = ', vz(i,1,1,60),
     &                 'vztrc =', vz_trc(i,1,1,60)
          enddo 
        endif

c       continuity by dss on the filtered array (not in spec coeffs)
#ifdef AMR
        call amr_oph1_proj(vx_trc,vy_trc,vz_trc,lx1,ly1,lz1,nelv)
#else
        call opdssum(vx_trc,vy_trc,vz_trc)       !apply dss to the fields
        call opcolv (vx_trc,vy_trc,vz_trc,VMULT) !divide by the multiplicity
#endif

c       copy the pressure into an array full of zeros (without mapping)
c       in the velocity mesh to prepare for writing
        call rzero(pm1,lx1*ly1*lz1*lelt)
        do e=1,nelv
              call copy  (pm1(1,1,1,e),pr_hat_trc(1,1,1,e),lx2*lx2*lx2)
        enddo

c       calculate the global error norm
        call sub3(errvec,vz_trc,vz,nxyze)
        l2err = gl2norm(errvec,nxyze)
        l2v = gl2norm(vz,nxyze)
        l2v_trc = gl2norm(vz_trc,nxyze)

        etime_t = dnekclock_sync() - etime0
        if(nio.eq.0) write(6,'(A,1(1g9.2),A,/)')
     &                   ' done :: truncation performed   ',
     &                    etime_t, ' sec'
c       ============================================
c       Truncation done, now calculate some metrics
c       ============================================

c       Get average compresion ratio
        avrgkutvx=0
        avrgkutvy=0
        avrgkutvz=0
        avrgkutvpr=0
        avrgcrsx=0
        avrgcrsy=0
        avrgcrsz=0
        avrgcrspr=0
        numel=0
        ! sum in the local process
        do e=1,nelv
            avrgkutvx=avrgkutvx+kutvx(e)
            avrgkutvy=avrgkutvy+kutvy(e)
            avrgkutvz=avrgkutvz+kutvz(e)
            avrgkutvpr=avrgkutvpr+kutvpr(e)
            avrgcrsx=avrgcrsx+crsx(e)
            avrgcrsy=avrgcrsy+crsy(e)
            avrgcrsz=avrgcrsz+crsz(e)
            avrgcrspr=avrgcrspr+crspr(e)
            numel=numel+1
        enddo
        !do a reduction to sum from all processes
        avrgkutvx=glsum(avrgkutvx,1)
        avrgkutvy=glsum(avrgkutvy,1)
        avrgkutvz=glsum(avrgkutvz,1)
        avrgkutvpr=glsum(avrgkutvpr,1)
        avrgcrsx=glsum(avrgcrsx,1)
        avrgcrsy=glsum(avrgcrsy,1)
        avrgcrsz=glsum(avrgcrsz,1)
        avrgcrspr=glsum(avrgcrspr,1)
        numel=glsum(numel,1)
        !get the average
        avrgkutvx=avrgkutvx/numel
        avrgkutvy=avrgkutvy/numel
        avrgkutvz=avrgkutvz/numel
        avrgkutvpr=avrgkutvpr/numel
        avrgcrsx=avrgcrsx/numel
        avrgcrsy=avrgcrsy/numel
        avrgcrsz=avrgcrsz/numel
        avrgcrspr=avrgcrspr/numel


c       write debugging values
        if (nio.eq.0) then 
          write(*,*) 'target error  = ', targeterr
          write(*,*) 'l2 error      = ', l2err
          write(*,*) 'l2 norm vz    = ', l2v
          write(*,*) 'l2 norm vz_trc= ', l2v_trc
          write(*,*) 'AVRG CR in x  = ', avrgcrsx
          write(*,*) 'AVRG CR in y  = ', avrgcrsy
          write(*,*) 'AVRG CR in z  = ', avrgcrsz
          write(*,*) 'AVRG CR in pr = ', avrgcrspr
          write(*,*) 'numel         = ', numel
        endif

        ! Copy the compressed data to compare with what is read
        if (ifile.eq.1) call copy(vx_trc_temp,vx_trc,nxyze)
        if (ifile.eq.1) call copy(vy_trc_temp,vy_trc,nxyze)
        if (ifile.eq.1) call copy(vz_trc_temp,vz_trc,nxyze)
        if (ifile.eq.1) call copy(vx_hat_trc_temp,vx_hat_trc,nxyze)
        if (ifile.eq.1) call copy(vy_hat_trc_temp,vy_hat_trc,nxyze)
        if (ifile.eq.1) call copy(vz_hat_trc_temp,vz_hat_trc,nxyze)


c       ===============================================
c       Perform lossless compression of truncated data
c       ===============================================

c       For debugging
c       write the reconstructed fields.        
        call outpost(vx_trc,vy_trc,
     $     vz_trc,pr_trc,t,'rc2')
        if(nio.eq.0) write(6,'(A,1(1g9.2),A,/)')
     &                   ' done :: print truncated data   ',
     &                    etime_t, ' sec'

        etime0 = dnekclock_sync()

        !Write the fields into a file
        call adios2_update(lglel,pm1, vx_hat_trc, 
     &   vy_hat_trc, vz_hat_trc, t)

        etime_t = dnekclock_sync() - etime0
        if(nio.eq.0) write(6,'(A,1(1g9.2),A,/)')
     &                   ' done :: data written by ADIOS   ',
     &                    etime_t, ' sec'
        enddo

      return
      end subroutine
!
!=======================================================================
!> @brief Decompress lossless compressed files with ADIOS2
!! @ingroup io_compression
      subroutine decompress_inputs()
      implicit none
      include 'SIZE'    ! for size information such as lx1, etc ...
      include 'WZ'      ! for zgm1
      include 'INPUT'   ! for user inputs from par (I think) 
      include 'TSTEP'   ! istep
      include 'SOLN'    ! vx,vy,vz,pr
      include 'PARALLEL'! for nelgv, nelgt
      include 'RESTART' ! nelb (for io operations)
      include 'MASS'    ! for the mass matrix
      include 'CTIMER'  ! clock

      include 'IOTRUNCD' ! THIS IS MINE, include here the integer a

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

c     filter transfer function
      real trunc_filter(lm2)                    !filter array
c     filter cutout frequency
      integer kut
c     cut out as a function of the elements      
      integer kutvx(nelt),kutvy(nelt),kutvz(nelt), kutvpr(nelt)
c     compression ratio as a function of elements      
      real crsx(nelt),crsy(nelt),crsz(nelt), crspr(nelt)
c     Logical variables
      logical truncio
      logical readio

c     Variable to determine the type of basis function for transforms
      logical trunc_ifboyd
      logical ifsort

c     average kut
      real avrgkutvx,avrgkutvy, avrgkutvz, avrgkutvpr
      real avrgcrsx,avrgcrsy, avrgcrsz, avrgcrspr
      real numel

c     variables that contains the error and norms
      logical ifspec
      real l2err                          !l2 norm of error
      real l2errx                          !l2 norm of error
      real l2erry                          !l2 norm of error
      real l2errz                          !l2 norm of error
      real l2v                          !l2 norm of error
      real l2v_trc                          !l2 norm of error
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

c     spectral coefficient matrices for velocity and pressure
      real spec_trans_v(lm2)        ! Spectral transform matrix 
      real spec_trans_vt(lm2)       ! Spectral transform matrix V^T 
      real spec_trans_vinv(lm2)     ! Spectral transform matrix V^-1 
      real spec_trans_v_pr(lm2)     ! Spectral transform matrix V 
      real spec_trans_vt_pr(lm2)    ! Spectral transform matrix V^T 
      real spec_trans_vinv_pr(lm2)  ! Spectral transform matrix V^-1 

      real*8 etime_t

c     For the reader
      integer nvals,nid_,np_,nekcomm
      real pm1(lx1,ly1,lz1,lelv) !Create a pointer that you will reference to the common
      integer vid(2,nelt)

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

c     What basis function for transform: .false. is legendre
      trunc_ifboyd = .false.
c     polynomial order, opposed to lm which is used for allocation
      nx=lx1     
      nx_pr=lx2 
c     Total number of points in the mpirank
      nxyze=lx1*ly1*lz1*nelv
      nvals=lx1*lx1*lx1



c     Build operator matrices for spectral transformations
      call build_spec_matrix(trunc_ifboyd,nx,spec_trans_v,
     &                       spec_trans_vt,spec_trans_vinv)
      call build_spec_matrix(trunc_ifboyd,nx_pr,spec_trans_v_pr,
     &                       spec_trans_vt_pr,spec_trans_vinv_pr)

c     Decompress the files recursively
      do ifile=1,numfile

c       Compose the name of each file. This is currently unused...
c       ... the reader subroutine does the same internally.
        fname1="out.f0000"
        write (str,'(i1.1)') ifile
        str2=".bp"
        fname=trim(fname1)//trim(str)//trim(str2)

c       Some debugging data, uncomment if needed        
        !if (nio.eq.0) then
        !    write(*,*) 'inputs to adios'
        !    write(*,*) 'nvals', nvals
        !    write(*,*) 'elements before rank0', nelb
        !    write(*,*) 'nelgv', nelgv
        !    write(*,*) 'nelgt', nelgt
        !    write(*,*) 'nekcomm', nekcomm
        !    write(*,*) 'adios should get'
        !    write(*,*) 'nelv=', nelv
        !end if

c       ========================================
c       Lossless decompression of out.xxx files
c       ========================================

        etime0 = dnekclock_sync()
c       Fill the pressure array with zeros, just in case       
        call rzero(pm1,lx1*ly1*lz1*lelt)
c       Adios writes the pressure in mesh 1, so read it in pm1 
        call adios2_read(lglelr,pm1,vx_hat_trc,vy_hat_trc,vz_hat_trc,
     &                   nvals,nelv,nelb,
     &                   nelgv,nelgt,nekcomm,trim(fname))

        etime_t = dnekclock_sync() - etime0
        if(nio.eq.0) write(6,'(A,1(1g9.2),A,/)')
     &                   ' done :: ADIOS2 read file   ',
     &                    etime_t, ' sec'

c       Map the pressure to mesh 2 before transforming to phys space
c       copy the corresponding entries from pm1 to prhat
        do e=1,nelv
              call copy  (pr_hat_trc(1,1,1,e),pm1(1,1,1,e),lx2*ly2*lz2)
        enddo

c       =========================================
c       Element redistribution to appropiate nid
c       =========================================

        etime0 = dnekclock_sync()
c       Redistribute the data in spectral space since we need the...
c       ... to be sure that all is correct when transforming back.
        call trunc_vec_transfer(vx_hat_trc, vid, lx1*ly1*lz1,
     &                          nelt,lglelr)
        call trunc_vec_transfer(vy_hat_trc, vid, lx1*ly1*lz1,
     &                          nelt,lglelr)
        call trunc_vec_transfer(vz_hat_trc, vid, lx1*ly1*lz1,
     &                          nelt,lglelr)
        call trunc_vec_transfer(pr_hat_trc, vid, lx2*ly2*lz2,
     &                          nelt,lglelr)
        etime_t = dnekclock_sync() - etime0
        if(nio.eq.0) write(6,'(A,1(1g9.2),A,/)')
     &                   ' done :: data redistributed   ',
     &                    etime_t, ' sec'

c       =============================================
c       Reverse truncation (transform to phys space)
c       =============================================

        etime0 = dnekclock_sync()
c       transform the read field into physical space u=V*u_hat
c       For vx
        call trunc_apply_op(vx_trc,vx_hat_trc,spec_trans_v,if3d,nx)
c       For vy
        call trunc_apply_op(vy_trc,vy_hat_trc,spec_trans_v,if3d,nx)
c       For vz
        call trunc_apply_op(vz_trc,vz_hat_trc,spec_trans_v,if3d,nx)
c       For Pressure
        call trunc_apply_op(pr_trc,pr_hat_trc,spec_trans_v_pr,if3d,
     &                      nx_pr)

c       continuity by dss on the filtered array (not in spec coeffs)
#ifdef AMR
        call amr_oph1_proj(vx_trc,vy_trc,vz_trc,lx1,ly1,lz1,nelv)
#else
        call opdssum(vx_trc,vy_trc,vz_trc)       !apply dss to the fields
        call opcolv (vx_trc,vy_trc,vz_trc,VMULT) !divide by the multiplicity
#endif

c       calculate the global error norm
        call sub3(errvec,vz_trc,vz_trc_temp,nxyze)
        l2err = gl2norm(errvec,nxyze)
        l2v = gl2norm(vz_trc,nxyze)
        l2v_trc = gl2norm(vz_trc_temp,nxyze)
        etime_t = dnekclock_sync() - etime0
        if(nio.eq.0) write(6,'(A,1(1g9.2),A,/)')
     &                   ' done :: data reconstructed   ',
     &                    etime_t, ' sec'

c       write debugging values
        if (nio.eq.0) then
                write(*,*) 'Read file ', fname
        endif

        etime0 = dnekclock_sync()
c       write the reconstructed fields.        
        call outpost(vx_trc,vy_trc,
     $     vz_trc,pr_trc,t,'rct')
        if(nio.eq.0) write(6,'(A,1(1g9.2),A,/)')
     &                   ' done :: data writen to disk   ',
     &                    etime_t, ' sec'

c       =================================
c       Metrics to evaluate read process
c       =================================

        if (ifile.eq.1) then

c       calculate the global error norm
        call sub3(errvec,vx_trc,vx_trc_temp,nxyze)
        l2errx = gl2norm(errvec,nxyze)
        call sub3(errvec,vy_trc,vy_trc_temp,nxyze)
        l2erry = gl2norm(errvec,nxyze)
        call sub3(errvec,vz_trc,vz_trc_temp,nxyze)
        l2errz = gl2norm(errvec,nxyze)
        l2v = gl2norm(vz_trc,nxyze)
        l2v_trc = gl2norm(vz_trc_temp,nxyze)
c         write debugging values
        if (nio.eq.0) then
          write(*,*) 'This only works if you compress and decompress...'
          write(*,*) '... in the same session'
          write(*,*) 'Norm from the reader '
          write(*,*) 'l2 error in x= ', l2errx
          write(*,*) 'l2 error in y= ', l2erry
          write(*,*) 'l2 error in z= ', l2errz
          write(*,*) 'l2 norm vz new= ', l2v
          write(*,*) 'l2 norm vz_old= ', l2v_trc
        endif
        endif


      enddo

      return
      end subroutine
!
!=======================================================================
!> @brief In-Situ compression and decompression 
!! @ingroup io_compression
!! @note The read functionality has no real use in a normal run ... it
!!       should be used for debugging purposes mostly.
      subroutine trunc_data()
      implicit none
      include 'SIZE'    ! for size information such as lx1, etc ...
      include 'WZ'      ! for zgm1
      include 'INPUT'   ! for user inputs from par (I think) 
      include 'TSTEP'   ! istep
      include 'SOLN'    ! vx,vy,vz,pr
      include 'PARALLEL'! for nelgv, nelgt
      include 'RESTART' ! nelb (for io operations)


      include 'IOTRUNCD' ! THIS IS MINE, include here the integer a

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

c     filter transfer function
      real trunc_filter(lm2)                    !filter array
c     filter cutout frequency
      integer kut
c     cut out as a function of the elements      
      integer kutvx(nelt),kutvy(nelt),kutvz(nelt), kutvpr(nelt)
c     compression ratio as a function of elements      
      real crsx(nelt),crsy(nelt),crsz(nelt), crspr(nelt)
c     Logical variables
      logical truncio
      logical readio

c     Variable to determine the type of basis function for transforms
      logical trunc_ifboyd
      logical ifsort

c     average kut
      real avrgkutvx,avrgkutvy, avrgkutvz, avrgkutvpr
      real avrgcrsx,avrgcrsy, avrgcrsz, avrgcrspr
      real numel

c     variables that contains the error and norms
      logical ifspec
      real l2err                          !l2 norm of error
      real l2errx                         !l2 norm of error
      real l2erry                         !l2 norm of error
      real l2errz                         !l2 norm of error
      real l2errpr                         !l2 norm of error
      real l2v                          !l2 norm of error
      real l2v_trc                          !l2 norm of error
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

c     spectral coefficient matrices for velocity and pressure
      real spec_trans_v(lm2)        ! Spectral transform matrix 
      real spec_trans_vt(lm2)       ! Spectral transform matrix V^T 
      real spec_trans_vinv(lm2)     ! Spectral transform matrix V^-1 
      real spec_trans_v_pr(lm2)     ! Spectral transform matrix V 
      real spec_trans_vt_pr(lm2)    ! Spectral transform matrix V^T 
      real spec_trans_vinv_pr(lm2)  ! Spectral transform matrix V^-1 

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

c     What basis function for transform: .false. is legendre
      trunc_ifboyd = .false.
c     polynomial order, opposed to lm which is used for allocation
      nx=lx1     
      nx_pr=lx2 
c     Total number of points in the mpirank
      nxyze=lx1*ly1*lz1*nelv
      nvals=lx1*lx1*lx1

c     Set up parameters. 
      truncio=.false.
      readio=.false.
      ifsort=.true.
      if(mod(istep,iotruncstep).eq.0.and.istep.ne.0) truncio=.true.
      if(mod(istep,ioreadstep).eq.0.and.istep.ne.0) readio=.true.
      if (istep.eq.0) ifile=0

      !Initialize I/O data
      call io_init

c     =================
c     Data compression
c     =================

c     Apply the operators as needed
      if (truncio) then

c       Build operator each time we write... consider putting this out
        call build_spec_matrix(trunc_ifboyd,nx,spec_trans_v,
     &                         spec_trans_vt,spec_trans_vinv)
        call build_spec_matrix(trunc_ifboyd,nx_pr,spec_trans_v_pr,
     &                         spec_trans_vt_pr,spec_trans_vinv_pr)

c       transform the field into spectral space u_hat=V^-1*u
c       For vx
        call trunc_apply_op(vx_hat,vx,spec_trans_vinv,if3d,nx)
c       For vy
        call trunc_apply_op(vy_hat,vy,spec_trans_vinv,if3d,nx)
c       For vz
        call trunc_apply_op(vz_hat,vz,spec_trans_vinv,if3d,nx)
c       For Pressure
        call trunc_apply_op(pr_hat,pr,spec_trans_vinv_pr,if3d,nx_pr)

c       Truncate the field according to an input error        
c       For vx
        call trunc_field_l2oferror(vx_trc,vx_hat_trc,normv_physx,
     &       normv_specx,kutvx,vx_hat,spec_trans_v,spec_trans_vt,
     &       nx,targeterr,crsx,ifsort)
c       For vy
        call trunc_field_l2oferror(vy_trc,vy_hat_trc,normv_physy,
     &       normv_specy,kutvy,vy_hat,spec_trans_v,spec_trans_vt,
     &       nx,targeterr,crsy,ifsort)
c       For vz
        call trunc_field_l2oferror(vz_trc,vz_hat_trc,normv_physz,
     &       normv_specz,kutvz,vz_hat,spec_trans_v,spec_trans_vt,
     &       nx,targeterr,crsz,ifsort)
c       For pr
        call trunc_field_l2oferror(pr_trc,pr_hat_trc,normv_physpr,
     &       normv_specpr,kutvpr,pr_hat,spec_trans_v_pr,
     &       spec_trans_vt_pr, nx_pr,targeterr,crspr,ifsort)

c       continuity by dss on the filtered array (not in spec coeffs)
        call opdssum(vx_trc,vy_trc,vz_trc)       !apply dss to the fields
        call opcolv (vx_trc,vy_trc,vz_trc,VMULT) !divide by the multiplicity

c       Copy the pressure into an array full of zeros in the vel mesh    
c       in preparation to write the files
        call rzero(pm1,lx1*ly1*lz1*lelt)
        do e=1,nelv
              call copy  (pm1(1,1,1,e),pr_hat_trc(1,1,1,e),lx2*lx2*lx2)
        enddo

c       calculate the global error norm
        call sub3(errvec,vz_trc,vz,nxyze)
        l2err = gl2norm(errvec,nxyze)
        l2v = gl2norm(vz,nxyze)
        l2v_trc = gl2norm(vz_trc,nxyze)

c       ============================================
c       Truncation done, now calculate some metrics
c       ============================================

c ---------------------------------------------------------------
c       This piece is for debugging, uncumment if needed        
        if (truncio) then
          !!write(*,*) 'copying values since it is step: ', istep       
          call copy(vx_trc_temp,vx_trc,nxyze)
          call copy(vx_hat_trc_temp,vx_hat_trc,nxyze)
          call copy(vy_trc_temp,vy_trc,nxyze)
          call copy(vy_hat_trc_temp,vy_hat_trc,nxyze)
          call copy(vz_trc_temp,vz_trc,nxyze)
          call copy(vz_hat_trc_temp,vz_hat_trc,nxyze)
          call copy(pr_trc_temp,pr_trc,lx2*lx2*lx2*nelv)
        endif
c ---------------------------------------------------------------        

c       Get the error vectors that can be visualized in visit
c       For the velocity
        do e=1,nelv
          do i=1,nx*nx*nx
            errvecx(i,1,1,e)=normv_physx(e)
            truncvecx(i,1,1,e)=crsx(e)
            errvecy(i,1,1,e)=normv_physy(e)
            truncvecy(i,1,1,e)=crsy(e)
            errvecz(i,1,1,e)=normv_physz(e)
            truncvecz(i,1,1,e)=crsz(e)
          enddo
        enddo
c       For the pressure        
        do e=1,nelv
          do i=1,nx_pr*nx_pr*nx_pr
            errvecpr(i,1,1,e)=normv_physpr(e)
            truncvecpr(i,1,1,e)=crspr(e)
          enddo
        enddo

c       Get average compresion ratio
        avrgkutvx=0
        avrgkutvy=0
        avrgkutvz=0
        avrgkutvpr=0
        avrgcrsx=0
        avrgcrsy=0
        avrgcrsz=0
        avrgcrspr=0
        numel=0
        ! sum in the local process
        do e=1,nelv
            avrgkutvx=avrgkutvx+kutvx(e)
            avrgkutvy=avrgkutvy+kutvy(e)
            avrgkutvz=avrgkutvz+kutvz(e)
            avrgkutvpr=avrgkutvpr+kutvpr(e)
            avrgcrsx=avrgcrsx+crsx(e)
            avrgcrsy=avrgcrsy+crsy(e)
            avrgcrsz=avrgcrsz+crsz(e)
            avrgcrspr=avrgcrspr+crspr(e)
            numel=numel+1
        enddo
        !do a reduction to sum from all processes
        avrgkutvx=glsum(avrgkutvx,1)
        avrgkutvy=glsum(avrgkutvy,1)
        avrgkutvz=glsum(avrgkutvz,1)
        avrgkutvpr=glsum(avrgkutvpr,1)
        avrgcrsx=glsum(avrgcrsx,1)
        avrgcrsy=glsum(avrgcrsy,1)
        avrgcrsz=glsum(avrgcrsz,1)
        avrgcrspr=glsum(avrgcrspr,1)
        numel=glsum(numel,1)
        !get the average
        avrgkutvx=avrgkutvx/numel
        avrgkutvy=avrgkutvy/numel
        avrgkutvz=avrgkutvz/numel
        avrgkutvpr=avrgkutvpr/numel
        avrgcrsx=avrgcrsx/numel
        avrgcrsy=avrgcrsy/numel
        avrgcrsz=avrgcrsz/numel
        avrgcrspr=avrgcrspr/numel


c       write debugging values
        if (nio.eq.0) then 
          write(*,*) 'target error  = ', targeterr
          write(*,*) 'l2 error      = ', l2err
          write(*,*) 'l2 norm vz    = ', l2v
          write(*,*) 'l2 norm vz_trc= ', l2v_trc
          write(*,*) 'AVRG CR in x  = ', avrgcrsx
          write(*,*) 'AVRG CR in y  = ', avrgcrsy
          write(*,*) 'AVRG CR in z  = ', avrgcrsz
          write(*,*) 'AVRG CR in pr = ', avrgcrspr
          write(*,*) 'numel         = ', numel
        endif

c ---------------------------------------------------------
c    --------- To write the different parameters ----------
c ---------------------------------------------------------        
c      param(63) = 1       ! Enforce 64-bit output
c     Write to file
c      call outpost(vx_trc,vy_trc,
c     $     vz_trc,pr_trc,t,'dss')
c
c      call outpost(vx_hat_trc,vy_hat_trc,
c     $     vz_hat_trc,pr_hat_trc,t,'hat')

c      call outpost(errvecx,errvecy,
c     $     errvecz,pr,t,'err')

c      call outpost(truncvecx,truncvecy,
c     $     truncvecz,pr,t,'crt')

c      call outpost(vx,vy,
c     $     vz,pr,t,'fll')
c ----------------------------------------------------------

c       =============================================
c       Perform the lossless compression with ADIOS2
c       =============================================

        !pm1 is the truncated coefficients mapped into velocity mesh
        call adios2_update(lglel,pm1, vx_hat_trc, 
     &   vy_hat_trc, vz_hat_trc, t)

      endif    

c     ===================
c     Data decompression
c     ===================

      if (readio) then

c       Build operator each time we write... consider putting this out
        call build_spec_matrix(trunc_ifboyd,nx,spec_trans_v,
     &                         spec_trans_vt,spec_trans_vinv)
        call build_spec_matrix(trunc_ifboyd,nx_pr,spec_trans_v_pr,
     &                         spec_trans_vt_pr,spec_trans_vinv_pr)


c       Compose the name of the files...
c       currently unused, this is done directly in read routine.
        ifile=ifile+1
        fname1="out.f0000"
        write (str,'(i1.1)') ifile
        str2=".bp"
        fname=trim(fname1)//trim(str)//trim(str2)

c       ========================================
c       Lossless decompression of out.xxx files
c       ========================================

c       fill the pressure array with zeros, just in case       
        call rzero(pm1,lx1*ly1*lz1*lelt)
c       Adios writes the pressure in mesh 1, so read it in pm1 
        call adios2_read(lglelr,pm1,vx_hat_trc,vy_hat_trc,vz_hat_trc,
     &                   nvals,nelv,nelb,
     &                   nelgv,nelgt,nekcomm,trim(fname))

c       Map the pressure to mesh 2 before transforming to phys space
c       copy the corresponding entries from pm1 to prhat
        do e=1,nelv
              call copy  (pr_hat_trc(1,1,1,e),pm1(1,1,1,e),lx2*ly2*lz2)
        enddo

c       =========================================
c       Element redistribution to appropiate nid
c       =========================================

c       Redistribute the read data
        call trunc_vec_transfer(vx_hat_trc, vid, lx1*ly1*lz1,
     &                          nelt,lglelr)
        call trunc_vec_transfer(vy_hat_trc, vid, lx1*ly1*lz1,
     &                          nelt,lglelr)
        call trunc_vec_transfer(vz_hat_trc, vid, lx1*ly1*lz1,
     &                          nelt,lglelr)
        call trunc_vec_transfer(pr_hat_trc, vid, lx2*ly2*lz2,
     &                          nelt,lglelr)

c       ===================
c       Reverse truncation
c       ===================

c       transform the read field into physical space u=V*u_hat
c       For vx
        call trunc_apply_op(vx_trc,vx_hat_trc,spec_trans_v,if3d,nx)
c       For vy
        call trunc_apply_op(vy_trc,vy_hat_trc,spec_trans_v,if3d,nx)
c       For vz
        call trunc_apply_op(vz_trc,vz_hat_trc,spec_trans_v,if3d,nx)
c       For Pressure
        call trunc_apply_op(pr_trc,pr_hat_trc,spec_trans_v_pr,if3d,
     &                      nx_pr)

c       continuity by dss on the filtered array (not in spec coeffs)
        call opdssum(vx_trc,vy_trc,vz_trc)       !apply dss to the fields
        call opcolv (vx_trc,vy_trc,vz_trc,VMULT) !divide by the multiplicity


c       ================================
c       Metrics to evaluate the process
c       ================================

c       calculate the global error norm
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

c       write debugging values
        if (nio.eq.0) then
          write(*,*) 'Norm from the reader in z '
          write(*,*) 'l2 error in x= ', l2errx
          write(*,*) 'l2 error in y= ', l2erry
          write(*,*) 'l2 error in z= ', l2errz
          write(*,*) 'l2 error in pr= ', l2errpr
          write(*,*) 'l2 norm vz new= ', l2v
          write(*,*) 'l2 norm vz_old= ', l2v_trc
        endif

c       ===============================
c       Write out reconstructed fields
c       ===============================

c       write the reconstructed fields for debbuging.        
        call outpost(vx_trc,vy_trc,
     $     vz_trc,pr_trc,t,'rct')

      endif



      return
      end subroutine
!
!=======================================================================
!> @brief Build the transformation matrices from and to lengendre space 
!! @ingroup io_compression
      subroutine build_spec_matrix(trunc_ifboyd,nx,spec_trans_v,
     &           spec_trans_vt,spec_trans_vinv)
      implicit none
      include 'SIZE'    ! for size information such as lx1, etc ...
      include 'WZ'      ! for quadrature weigths

c     for debuggin and loops
      integer i,k,kk

c     Polynomial order in memory (simply for allocation)
      integer lm,lm2
      parameter (lm=40)            !Hard set to 40
      parameter (lm2=lm*lm)
 
c     Transformation matrices
      real spec_trans_v(lm2)       ! Spectral transform matrix V 
      real spec_trans_vt(lm2)      ! Spectral transform matrix V^T 
      real spec_trans_vinv(lm2)    ! Spectral transform matrix V^-1

c     Variable that determines what type of basis function is used
      logical trunc_ifboyd
      logical iforthonorm
c     Polynomial order
      integer nx

c     Variables for matrix inversion
      real rmult(nx)
      integer indr(nx),indc(nx),ipiv(nx)
      integer ierr
 
c     weigth matrix
      real wka(lm2)      !working array
      real spec_w(lm2)

c     initialize the diagonal weigth matrix
      call ident(spec_w,nx)
c     now copy the weights into the main diagonal
      do k=1,nx
         kk = k+nx*(k-1)
         if (nx.eq.lx1) then
           spec_w(kk) = wxm1(k)
         endif
         if (nx.eq.lx2) then
           spec_w(kk) = wxm2(k)
         endif
      enddo      

c     Determine the spectral transformation matrix V
      call spec_coeff_init_varord(spec_trans_v,trunc_ifboyd,nx)

c     scale the V matrix by f=sqrt((2(j+1)-1/)2) to ensure V.T*W*V=I
      iforthonorm = .true.
      if (iforthonorm) then
      call spec_coeff_scale(spec_trans_v,spec_w,nx)
      endif

c     transpose the spectral coefficients      
      call transpose(spec_trans_vt,nx,spec_trans_v,nx)

c     ### Note ### Inverse by gaujordf produce less error than V.T*W
c     Copy the values of the matrix into the array that holds the inv
      call copy(spec_trans_vinv,spec_trans_v,nx*nx)
c     Calculate the inverse V^-1
      call gaujordf(spec_trans_vinv,nx,nx,indr,indc,ipiv,ierr,rmult)

      return
      end subroutine
!
!=======================================================================
!> @brief Scale the transformation matrix to ensure orthogonality wrt. W 
!! @ingroup io_compression
!! @note This subroutine test that V^T*W*V=I
      subroutine spec_coeff_scale(v,w,nx)
      implicit none
      include 'SIZE'    ! for size information such as lx1, etc ...

c     counters
      integer i,j

c     size parameters
      integer nx

c     inputs 
      real w(nx,nx)    !weight matrix                  
      real v(nx,nx)    !coeff matrix                  
      real vinv(nx,nx) !coeff matrix inverse                  
      real vt(nx,nx)   !coeff matrix transposed                 
      real vtwv(nx,nx) !v.t*w*v                  
      real delta(nx)   !scaling factor for legendre polynomials

c     Calculate the nominal scaling factors
      do i=1,nx
        delta(i)=2./(2*(i-1)+1)
      enddo
      
c     In GLL points we need to modify the last entry. Not needed for GL      
      if (nx.eq.lx1) then
        delta(nx)=2./(nx-1)
      endif

c     calculate the inverse to multiply the matrix
      do i=1,nx
        delta(i)=sqrt(1./delta(i))
      enddo
   
c     scale the matrix      
      do i=1,nx
        do j=1,nx
          v(i,j)=v(i,j)*delta(j)
        enddo
      enddo
      
c     Veryfy that the scalinf was okay by checking V.T*W*V=I
      call transpose(vt,nx,v,nx)
c     Calculate the inverse      
      call mxm(vt,nx,w,nx,vinv,nx)
      call mxm(vinv,nx,v,nx,vtwv,nx)
      
      return
      end subroutine
!=======================================================================
!> @brief Apply the operators to the input matrix. 
!! @ingroup io_compression
!! @note This subroutine applies the operator
      subroutine trunc_apply_op(vx_hat,vx,trunc_op,if3d,nx)
      implicit none
      include 'SIZE'    ! for size information such as lx1, etc ...

c     For pressure nx=lx2 and for velocity nx=lx1

c     outputs
      integer nx
      real trunc_op(nx,nx)           !Op modal, filter
      real trunc_opt(nx,nx)          !Op transposed
      real vx_hat(nx,nx,nx,lelv) 
      real vx(nx,nx,nx,lelv) 

c     work variables to apply the filter 
      real wk1  (nx,nx,nx,lelt)              !working array 1
      real wk2  (nx,nx,nx)                   !working awway 2
      integer uumax,vvmax,wwmax,nxyze,i           !monitoring variables

c     type of problem
      logical if3d

c     get the total number of entries in solution array
      nxyze=nx*nx*nx*nelv
c     copy the solution vector into vx_hat to initialize
      call copy(vx_hat,vx,nxyze)

      if (nio.eq.0) then
        write(*,*) 'Field - Copied before transforming'
        do i=1,10
          write(*,*) 'vz = ', vx(i,1,1,60),
     &               'vz_copy =', vx_hat(i,1,1,60)
        enddo 
      endif


c     apply the operator and get solution in v_hat 
      call trunc_filterq(vx_hat,trunc_op,nx,nx,wk1,wk2,
     &                   trunc_opt,if3d,uumax)
c         call filterq(vy,intv,lx1,lz1,wk1,wk2,intt,if3d,vmax)
c         if (if3d)
c     $   call filterq(vz,intv,lx1,lz1,wk1,wk2,intt,if3d,wmax)

      return
      end subroutine
!

!=======================================================================
!> @brief Build a simple filter. Identity matrix with zeros in main diag
!! @ingroup io_compression
      subroutine build_filter_tf(trunc_filter,
     &           kut,nx)
      implicit none
      include 'SIZE'    ! for size information such as lx1, etc ...

c     filter transfer function
c     Allocate arrays that depend on poly order, hard set to ordr 40
      real trunc_filter(nx*nx)                    !filter array

c     filter cutout frequency
      integer kut
c     polynomial order
      integer nx
c     counters
      integer k,kk,k0                       !frecuency variables

c     make the filter an identity matrix (all frequencies are kept)
c     here use nx instead of lm2, the extra allocated space is not used.
      call ident(trunc_filter,nx)

c     Create the actual filter function
      k0=nx-kut
      do k=k0+1,nx
        kk= k+nx*(k-1)
        trunc_filter(kk)=0
      enddo

      return
      end subroutine
!=======================================================================
!> @brief Truncate the fields with l2 norm of error vector
!! @ingroup io_compression
!! @note This subroutine truncates legendre coefficients and evaluate
!the error. It can sort or not sort the coefficients before truncating.
!The truncation is done in a while loop that progressively deletes more
!entries until an error bound is met. THIS USES THE NORM OF THE ERROR
      subroutine trunc_field_l2oferror(u_trc,uhat_trc,normvphys,normv,
     &           kutv,v,spec_v,spec_vt,nx,targeterr,crs,ifsort)
      implicit none
      include 'SIZE'    ! for size information such as lx1, etc ...
      include 'MASS'    ! for the mass matrix
      include 'WZ'
      include 'GEOM'

c     counters
      integer i,j,e,k
      integer nx

c     size parameters
      integer nxyz

c     working arrays
      real sum
      real w1  (nx,nx,nx)
      real w2  (nx,nx,nx)
      real errvec(nx,nx,nx)          !error vector
      real errvecphys(nx,nx,nx)      !error vector
      real ubase(nx,nx,nx,nelt)      !baseline for norm calculation
      real filterx(nx,nx)             !filter array
      real filtery(nx,nx)             !filter array
      real filterz(nx,nx)             !filter array
      real volv(nelt)                   !volume per element
      integer kutx,kuty,kutz          
      real tempnorm
      real uhat_temp(nx,nx,nx)  !transformed spectral coeffs
      real u_temp(nx,nx,nx)     !transformed spectral coeffs

c     inputs 
      real v(nx,nx,nx,lelv)       ! Spectral field to truncate
      real spec_v(nx,nx)             ! Spectral transform matrix V 
      real spec_vt(nx,nx)            ! Spectral transform matrix V^T 
      integer kut                    ! Cutoff frequency
      real targeterr

c     outputs     
      real uhat_trc(nx,nx,nx,nelt)  !transformed spectral coeffs
      real u_trc(nx,nx,nx,nelt)     !transformed spectral coeffs
      real normv(nelt)                 !norm per element
      real normvphys(nelt)                 !norm per element
      integer kutv(nelt)
      real crs(nelt)
      integer cr1,cr2,cr3

c     sorting variables
      real vsort(nx,nx,nx)       !error vector v_hat-vx
      real jacsort(nx,nx,nx)       !error vector v_hat-vx
      real wrksort(nx,nx,nx)       !error vector v_hat-vx
      integer isort(nx,nx,nx)       !error vector v_hat-vx
      integer wrkisort(nx,nx,nx)       !error vector v_hat-vx
      logical ifsort


!!  @note Volume is equal to int(dv) in real space, in computational
!grid it would be int(Jdv) and if we want to numerically integrate it
!it would be W.T*J assuming both are columns. The mass matrix BM1
!already has the multiplication, for the folume we sum the entires

      nxyz=nx*nx*nx
      !Define the apporpiate weights depending on the mesh you are.
      if (nx.eq.lx1) then
      do e=1,nelv
        sum = 0.
        do i=1,nxyz
          sum = sum + BM1(i,1,1,e)
        enddo
        volv(e)= sum 
      enddo
      endif
      if (nx.eq.lx2) then
      do e=1,nelv
        sum = 0.
        do i=1,nxyz
          sum = sum + BM2(i,1,1,e)
        enddo
        volv(e)= sum 
      enddo
      endif


      do e=1,nelv

c       First transform to phys space applying (VxVxV) * u with all
c       the frquenciesi to have a baseline 

        ! copy values for the element into a working array
        call copy(w2,v(1,1,1,e),nxyz)
        ! apply the operator to the x direction, result in w1
        call mxm(spec_v,nx,w2,nx,w1,nx*nx)
        !apply matrix to y direction, result in w2
        do k=1,nx
          call mxm(w1(1,1,k),nx,spec_vt,nx,w2(1,1,k),nx)
        enddo
        !apply matrix to z direction, result in w1
        call mxm (w2,nx*nx,spec_vt,nx,w1,nx)
        !the results are in w1, so copy them to u
        call copy(ubase(1,1,1,e),w1,nxyz)

c       write debugging values
        if (nio.eq.0.and.e.eq.60) then
          write(*,*) 'Field - Entry'
          do i=1,10
            write(*,*) 'vz = ', ubase(i,1,1,e),
     &                 'vz_hat =', v(i,1,1,e)
          end do 
        endif

        if (ifsort) then
c       At this point re order the coefficients! copy v(1,1,1,e) and reorder
          call sortcoeff(vsort,v(1,1,1,e),isort,wrksort,wrkisort,nxyz)
        else
          call copy(vsort,v(1,1,1,e),nxyz)
        endif

c       NOTE:Consider reordering the jacobian in another temporal variable

!! @iterative procedure starts here
        ! initialize some values
        tempnorm=0.
        kut=0

        do while (tempnorm.le.targeterr.and.kut.le.(nx-1)*3)

        !save values from previous iterations since if they enter this, it
        !means that they met tolerance
          kutv(e)=kut
          normv(e)=tempnorm 
          call copy(uhat_trc(1,1,1,e),uhat_temp,nxyz)

        !now advance kut
          kut=kut+1  

!! @note truncate and check the new error

      !the frecuency vector is arrayed as u_alpha_beta_gamma which gamma
      !being the slower, so if I set only the frecuencies to 0 in the z
      !filter I will effectivelly be deleting the highest possible
      !frecuencies

        !set the values acording to kut
          if (kut.eq.0) then
              kutz=0
              kuty=0
              kutx=0
          endif
          if (kut.gt.0.and.kut.le.(nx-1)) then
              kutz=kut
              kuty=0
              kutx=0
          endif
          if (kut.gt.(nx-1).and.kut.le.(nx-1)*2) then
              kutz=nx-1
              kuty=kut-(nx-1)
              kutx=0
          endif
          if (kut.gt.(nx-1)*2.and.kut.le.(nx-1)*3) then
              kutz=nx-1
              kuty=nx-1
              kutx=kut-(nx-1)*2
          endif

      !obtain the filter operator for the current it:
          call build_filter_tf(filterx,kutx,nx)
          call build_filter_tf(filtery,kuty,nx)
          call build_filter_tf(filterz,kutz,nx)

      !Discard Some Frecuencies by applying the filter.
      ! keep in mind that filter=filter^T
      ! copy values for the element into a working array
      ! normaly use v(1,1,1,e) but change it to vsort to have a general
      ! case of sorted and unsorted implementation.
          call copy(w2,vsort(1,1,1),nxyz)
      ! apply the operator to the x direction, result in w1
          call mxm(filterx,nx,w2,nx,w1,nx*nx)
      !apply matrix to y direction, result in w2
          do k=1,nx
             call mxm(w1(1,1,k),nx,filtery,nx,w2(1,1,k),nx)
          enddo
      !apply matrix to z direction, result in w1
          call mxm (w2,nx*nx,filterz,nx,w1,nx)
      !the results are in w1, so copy them to utr
          call copy(uhat_temp,w1,nxyz)
          
          if (ifsort) then
      !reorder the values after truncating them. Ideally sort the
      !jacobian instead of reordering each time.
          call reord(wrksort,isort,nx*nx*nx,uhat_temp)
          call copy(uhat_temp,wrksort,nxyz)
          endif
              
      !get the error between solutions
          call sub3(errvec,v(1,1,1,e),uhat_temp,nxyz)

      !now get the spectral space norm of err vect
          sum = 0.
          do i=1,nxyz
            if (nx.eq.lx1) then
              sum = sum + errvec(i,1,1)*errvec(i,1,1)*JACM1(i,1,1,e)
            endif
            if (nx.eq.lx2) then
              sum = sum + errvec(i,1,1)*errvec(i,1,1)*JACM2(i,1,1,e)
            endif
          enddo
          tempnorm= sqrt(sum)/sqrt(volv(e)) 

        enddo   

        if (kutv(e).eq.0) then
          call copy(uhat_trc(1,1,1,e),v(1,1,1,e),nxyz)
        endif


c       write debugging values
        if (nio.eq.0.and.e.eq.60) then
          write(*,*) 'Field - coefficients'
          do i=1,10
            write(*,*) 'uhat = ', v(i,1,1,e),
     &                 'uhat_trc =', uhat_trc(i,1,1,e)
          end do
        endif

      !transform back to phys space by applying (VxVxV)*u
      ! copy values from the truncated coefficients
        call copy(w2,uhat_trc(1,1,1,e),nxyz)
      ! apply the operator to the x direction, result in w1
        call mxm(spec_v,nx,w2,nx,w1,nx*nx)
      !apply matrix to y direction, result in w2
        do k=1,nx
           call mxm(w1(1,1,k),nx,spec_vt,nx,w2(1,1,k),nx)
        enddo
      !apply matrix to z direction, result in w1
        call mxm (w2,nx*nx,spec_vt,nx,w1,nx)
      !the results are in w1, so copy them to u
        call copy(u_trc(1,1,1,e),w1,nxyz)

c       write debugging values
        if (nio.eq.0.and.e.eq.60) then
          write(*,*) 'Field - after compression'
          do i=1,10
            write(*,*) 'ubase = ', ubase(i,1,1,e),
     &                 'u_trc =', u_trc(i,1,1,e)
          end do 
        endif

      !get the error in phys space
          call sub3(errvecphys,ubase(1,1,1,e),u_trc(1,1,1,e),nxyz)
      ! now get the norm
          sum = 0.
          do i=1,nxyz
            if (nx.eq.lx1) then
            sum = sum + errvecphys(i,1,1)*errvecphys(i,1,1)*BM1(i,1,1,e)
            endif
            if (nx.eq.lx2) then
            sum = sum + errvecphys(i,1,1)*errvecphys(i,1,1)*BM2(i,1,1,e)
            endif
          enddo
          normvphys(e)= sqrt(sum)/sqrt(volv(e)) 

      enddo

c     translate kutv into compression ratios
      !set the values acording to kut
      do e=1,nelv
        kut=kutv(e) 
        cr1=0
        cr2=0
        cr3=0 
        if (kut.le.(nx-1)) then
            cr1=kut*nxyz/nx    
        endif
        if (kut.gt.(nx-1).and.kut.le.(nx-1)*2) then 
            cr1=(nx-1)*nxyz/nx
            cr2=(kut-(nx-1))*nxyz/nx/nx    
        endif
        if (kut.gt.(nx-1)*2.and.kut.le.(nx-1)*3) then
            cr1=(nx-1)*nxyz/nx
            cr2=(nx-1)*nxyz/nx/nx    
            cr3=(kut-(nx-1)*2)*nxyz/nx/nx/nx    
        endif
        kutv(e)=nxyz/(nxyz-(cr1+cr2+cr3))
        crs(e)=real(cr1+cr2+cr3)/nxyz

      enddo
      


      return
      end subroutine
!=======================================================================
!> @brief Create the spectral transformation matrix V 
!! @ingroup io_compression
      subroutine spec_coeff_init_varord(ref_xmap,ifboyd,nx)
c     Initialise spectral coefficients
c     For legendre transform

c      implicit none

      include 'SIZE'
      include 'WZ'

      integer lm, lm2
      parameter (lm=40)
      parameter (lm2=lm*lm)

c     local variables
      integer i, j, k, n, nx, kj
c     Legendre polynomial
      real plegx(lm)
      real z
      real ref_xmap(lm2)
      real pht(lm2)

c     Change of basis
      logical IFBOYD

      kj = 0
      n  = nx-1
      do j=1,nx
        if (nx.eq.lx1) then
          z = ZGM1(j,1)
        endif
        if (nx.eq.lx2) then
          z = ZGM2(j,1)
        endif        
        call legendre_poly(plegx,z,n)
        kj = kj+1
        pht(kj) = plegx(1)
        kj = kj+1
        pht(kj) = plegx(2)

        if (IFBOYD) then        ! change basis to preserve element boundary values
          do k=3,nx
             kj = kj+1
             pht(kj) = plegx(k)-plegx(k-2)
          enddo
        else                    ! legendre basis    
          do k=3,nx
             kj = kj+1
             pht(kj) = plegx(k)
          enddo         
        endif
      enddo

      call transpose (ref_xmap,nx,pht,nx)

      return
      end
!=======================================================================
!> @brief reorder an array b by prestored indices. INVERSE OF SWAP 
!! @ingroup io_compression
!! @note this subroutine reorders array b by prestored indices this is
!the inverse of swap using the same indices
      subroutine reord(b,ind,n,temp)
      real B(1),TEMP(1)
      integer IND(1)

      do I=1,N
         JJ=IND(I)
         B(JJ)=TEMP(I)
      enddo

      return
      end
!=======================================================================
!> @brief Flip vector b and ind 
!! @ingroup io_compression
!! @note This subroutine flips a vector.
      subroutine flipv(b,ind,temp,tempind,n)

      real b(1),temp(1)
      integer ind(1),tempind(1)
      integer n,i,jj

      do i=1,n
         jj=n+1-i
         temp(jj)=b(i)
         tempind(jj)=ind(i)
      enddo
      do i=1,n
         b(i)=temp(i)
         ind(i)=tempind(i)
      enddo

      return
      end
!=======================================================================
!> @brief Sort the spectral coefficient in descending order 
!! @ingroup io_compression
!! @note This subroutine sorts the values of v in descending manner into
!array vsort. The original indices are stored in the isort vector.
      subroutine sortcoeff(vsort,v,isort,wrksort,wrkisort,nxyz)
   
      real vsort(1),v(1),wrksort(1)
      integer isort(1), wrkisort(1)
      integer nxyz      

c     copy absolute values to sort by magnitude
      do i=1,nxyz
         vsort(i)=abs(v(i))
         isort(i)=i
      enddo

c     sort the absolute values of the vectors, here the index is the
c     important part (isort)
      call sort(vsort,isort,nxyz)

c     Flip the indices so they are in a descending order
      call flipv(vsort,isort,wrksort,wrkisort,nxyz)

c     now re order the fields based on the indices
      call copy(vsort,v,nxyz)
      call swap(vsort,isort,nxyz,wrksort)

      return
      end subroutine        
!=======================================================================
!> @brief Truncate fields with the error of the l2 norms 
!! @ingroup io_compression
!! @note This subroutine truncates legendre coefficients and evaluate
!the error. It can sort or not sort the coefficients before truncating.
!The truncation is done in a while loop that progressively deletes more
!entries until an error bound is met. Here the comparison is done after
!separately obtaining the norms. THIS USES THE ERROR OF THE NORM.
      subroutine trunc_field_errorofl2(u_trc,uhat_trc,normvphys,normv,
     &           kutv,v,spec_v,spec_vt,nx,targeterr,crs,ifsort)
      implicit none
      include 'SIZE'    ! for size information such as lx1, etc ...
      include 'MASS'    ! for the mass matrix
      include 'WZ'
      include 'GEOM'

c     counters
      integer i,j,e,k
      integer nx

c     size parameters
      integer nxyz

c     working arrays
      real sum
      real w1  (nx,nx,nx)
      real w2  (nx,nx,nx)
      real errvec(nx,nx,nx)     !error vector
      real errvecphys(nx,nx,nx)      !error vector
      real ubase(nx,nx,nx,nelt)      !baseline for norm calculation
      real filterx(nx,nx)             !filter array
      real filtery(nx,nx)             !filter array
      real filterz(nx,nx)             !filter array
      real volv(nelt)                   !volume per element
      integer kutx,kuty,kutz          
      real tempnorm
      real tempnorm_og
      real uhat_temp(nx,nx,nx)  !transformed spectral coeffs
      real u_temp(nx,nx,nx)     !transformed spectral coeffs

c     inputs 
      real v(nx,nx,nx,lelv)       ! Spectral field to truncate
      real spec_v(nx,nx)             ! Spectral transform matrix V 
      real spec_vt(nx,nx)            ! Spectral transform matrix V^T 
      integer kut                    ! Cutoff frequency
      real targeterr

c     outputs     
      real uhat_trc(nx,nx,nx,nelt)  !transformed spectral coeffs
      real u_trc(nx,nx,nx,nelt)     !transformed spectral coeffs
      real normv(nelt)                 !norm per element
      real normvphys(nelt)                 !norm per element
      integer kutv(nelt)
      real crs(nelt)
      integer cr1,cr2,cr3

c     sorting variables
      real vsort(nx,nx,nx)       !error vector v_hat-vx
      real jacsort(nx,nx,nx)       !error vector v_hat-vx
      real wrksort(nx,nx,nx)       !error vector v_hat-vx
      integer isort(nx,nx,nx)       !error vector v_hat-vx
      integer wrkisort(nx,nx,nx)       !error vector v_hat-vx
      logical ifsort


!!  @note Volume is equal to int(dv) in real space, in computational
!grid it would be int(Jdv) and if we want to numerically integrate it
!it would be W.T*J assuming both are columns. The mass matrix BM1
!already has the multiplication, for the folume we sum the entires

      nxyz=nx*nx*nx
      !Define the apporpiate weights depending on the mesh you are.
      if (nx.eq.lx1) then
      do e=1,nelv
        sum = 0.
        do i=1,nxyz
          sum = sum + BM1(i,1,1,e)
        enddo
        volv(e)= sum 
      enddo
      endif
      if (nx.eq.lx2) then
      do e=1,nelv
        sum = 0.
        do i=1,nxyz
          sum = sum + BM2(i,1,1,e)
        enddo
        volv(e)= sum 
      enddo
      endif


      do e=1,nelv

c       First transform to phys space applying (VxVxV) * u with all
c       the frquenciesi to have a baseline 

        ! copy values for the element into a working array
        call copy(w2,v(1,1,1,e),nxyz)
        ! apply the operator to the x direction, result in w1
        call mxm(spec_v,nx,w2,nx,w1,nx*nx)
        !apply matrix to y direction, result in w2
        do k=1,nx
          call mxm(w1(1,1,k),nx,spec_vt,nx,w2(1,1,k),nx)
        enddo
        !apply matrix to z direction, result in w1
        call mxm (w2,nx*nx,spec_vt,nx,w1,nx)
        !the results are in w1, so copy them to u
        call copy(ubase(1,1,1,e),w1,nxyz)


        if (ifsort) then
c       At this point re order the coefficients! copy v(1,1,1,e) and reorder
          call sortcoeff(vsort,v(1,1,1,e),isort,wrksort,wrkisort,nxyz)
        else
          call copy(vsort,v(1,1,1,e),nxyz)
        endif

c       Consider reordering the jacobian in another temporal variable


!! @iterative procedure starts here
        ! initialize some values
        tempnorm=0.
        kut=0

      !get the norm of the untruncated field
          sum = 0.
          do i=1,nxyz
            if (nx.eq.lx1) then
              sum = sum + v(i,1,1,e)*v(i,1,1,e)*JACM1(i,1,1,e)
            endif
            if (nx.eq.lx2) then
              sum = sum + v(i,1,1,e)*v(i,1,1,e)*JACM2(i,1,1,e)
            endif
          enddo
          tempnorm_og= sqrt(sum)/sqrt(volv(e)) 

        do while (tempnorm.le.targeterr.and.kut.le.(nx-1)*3)

        !save values from previous iterations since if they enter this, it
        !means that they met tolerance
          kutv(e)=kut
          normv(e)=tempnorm 
          call copy(uhat_trc(1,1,1,e),uhat_temp,nxyz)

        !now advance kut
          kut=kut+1  

!! @note truncate and check the new error

      !the frecuency vector is arrayed as u_alpha_beta_gamma which gamma
      !being the slower, so if I set only the frecuencies to 0 in the z
      !filter I will effectivelly be deleting the highest possible
      !frecuencies

        !set the values acording to kut
          if (kut.eq.0) then
              kutz=0
              kuty=0
              kutx=0
          endif
          if (kut.gt.0.and.kut.le.(nx-1)) then
              kutz=kut
              kuty=0
              kutx=0
          endif
          if (kut.gt.(nx-1).and.kut.le.(nx-1)*2) then
              kutz=nx-1
              kuty=kut-(nx-1)
              kutx=0
          endif
          if (kut.gt.(nx-1)*2.and.kut.le.(nx-1)*3) then
              kutz=nx-1
              kuty=nx-1
              kutx=kut-(nx-1)*2
          endif

      !obtain the filter operator for the current it:
          call build_filter_tf(filterx,kutx,nx)
          call build_filter_tf(filtery,kuty,nx)
          call build_filter_tf(filterz,kutz,nx)

      !Discard Some Frecuencies by applying the filter.
      ! keep in mind that filter=filter^T
      ! copy values for the element into a working array
      ! normaly use v(1,1,1,e) but change it to vsort to have a general
      ! case of sorted and unsorted implementation.
          call copy(w2,vsort(1,1,1),nxyz)
      ! apply the operator to the x direction, result in w1
          call mxm(filterx,nx,w2,nx,w1,nx*nx)
      !apply matrix to y direction, result in w2
          do k=1,nx
             call mxm(w1(1,1,k),nx,filtery,nx,w2(1,1,k),nx)
          enddo
      !apply matrix to z direction, result in w1
          call mxm (w2,nx*nx,filterz,nx,w1,nx)
      !the results are in w1, so copy them to utr
          call copy(uhat_temp,w1,nxyz)
          
          if (ifsort) then
      !reorder the values after truncating them. Ideally sort the
      !jacobian instead of reordering each time.
          call reord(wrksort,isort,nx*nx*nx,uhat_temp)
          call copy(uhat_temp,wrksort,nxyz)
          endif


      !now get the spectral space norm 
          sum = 0.
          do i=1,nxyz
            if (nx.eq.lx1) then
            sum = sum + uhat_temp(i,1,1)*uhat_temp(i,1,1)*JACM1(i,1,1,e)
            endif
            if (nx.eq.lx2) then
            sum = sum + uhat_temp(i,1,1)*uhat_temp(i,1,1)*JACM2(i,1,1,e)
            endif
          enddo
          tempnorm= sqrt(sum)/sqrt(volv(e))
          tempnorm=abs(tempnorm-tempnorm_og) 

        enddo   

        if (kutv(e).eq.0) then
          call copy(uhat_trc(1,1,1,e),v(1,1,1,e),nxyz)
        endif

      !transform back to phys space by applying (VxVxV)*u
      ! copy values from the truncated coefficients
        call copy(w2,uhat_trc(1,1,1,e),nxyz)
      ! apply the operator to the x direction, result in w1
        call mxm(spec_v,nx,w2,nx,w1,nx*nx)
      !apply matrix to y direction, result in w2
        do k=1,nx
           call mxm(w1(1,1,k),nx,spec_vt,nx,w2(1,1,k),nx)
        enddo
      !apply matrix to z direction, result in w1
        call mxm (w2,nx*nx,spec_vt,nx,w1,nx)
      !the results are in w1, so copy them to u
        call copy(u_trc(1,1,1,e),w1,nxyz)

      !get the error in phys space
      ! now get the norm
          sum = 0.
          do i=1,nxyz
            if (nx.eq.lx1) then
            sum = sum + ubase(i,1,1,e)*ubase(i,1,1,e)*BM1(i,1,1,e)
            endif
            if (nx.eq.lx2) then
            sum = sum + ubase(i,1,1,e)*ubase(i,1,1,e)*BM2(i,1,1,e)
            endif
          enddo
          tempnorm_og= sqrt(sum)/sqrt(volv(e)) 
c         For truncated
          sum = 0.
          do i=1,nxyz
            if (nx.eq.lx1) then
            sum = sum + u_trc(i,1,1,e)*u_trc(i,1,1,e)*BM1(i,1,1,e)
            endif
            if (nx.eq.lx2) then
            sum = sum + u_trc(i,1,1,e)*u_trc(i,1,1,e)*BM2(i,1,1,e)
            endif
          enddo
          tempnorm= sqrt(sum)/sqrt(volv(e))
          normvphys(e)=abs(tempnorm-tempnorm_og) 

      enddo

c     translate kutv into compression ratios
      !set the values acording to kut
      do e=1,nelv
        kut=kutv(e) 
        cr1=0
        cr2=0
        cr3=0 
        if (kut.le.(nx-1)) then
            cr1=kut*nxyz/nx    
        endif
        if (kut.gt.(nx-1).and.kut.le.(nx-1)*2) then 
            cr1=(nx-1)*nxyz/nx
            cr2=(kut-(nx-1))*nxyz/nx/nx    
        endif
        if (kut.gt.(nx-1)*2.and.kut.le.(nx-1)*3) then
            cr1=(nx-1)*nxyz/nx
            cr2=(nx-1)*nxyz/nx/nx    
            cr3=(kut-(nx-1)*2)*nxyz/nx/nx/nx    
        endif
        kutv(e)=nxyz/(nxyz-(cr1+cr2+cr3))
        crs(e)=real(cr1+cr2+cr3)/nxyz

      enddo
      


      return
      end subroutine
!=======================================================================









c     Blank Space. Subroutines after this should be in another file.











!=======================================================================
!> @brief Read xxx.fxxx files into _trunc vectors 
!! @ingroup io_compression
      subroutine trunc_mfi(fname_in,ifile)
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
      include 'IOTRUNCD'

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
               call mfi_getv(vx_hat_trc,vy_hat_trc,vz_hat_trc,
     &                       wk,lwk,.false.)
c               call mfi_getv(vx,vy,vz,
c     &                       wk,lwk,.false.)
            endif
         else
            call mfi_getv(vx_hat_trc,vy_hat_trc,vz_hat_trc,
     &                    wk,lwk,.true.)
c            call mfi_getv(vx,vy,vz,
c     &                    wk,lwk,.true.)
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
      if (ifgetp) call trunc_map_pm1_to_pr(pm1,ifile) ! Interpolate pressure

      return
      end
!=======================================================================
!> @brief map pm1 into pr_hat_trc
!! @ingroup io_compression
!! @note This is only useful if we have written the truncated data into
!a xxx.fxxx file and want to read it. This was used at the intial stages
!of the project to verify compression but now the reads are made with
!ADIOS2 instead.

      subroutine trunc_map_pm1_to_pr(pm1,ifile)

      include 'SIZE'
      include 'TOTAL'
      include 'RESTART'
      include 'IOTRUNCD'

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
               call copy  (pr_hat_trc(1,1,1,e),pm1(1,e),nxyz2)
            else
               call map12 (pr_hat_trc(1,1,1,e),pm1(1,e),e)
            endif
         enddo
      endif
   
      return
      end

!=======================================================================
!> @brief Redistribute single variable after reading from ADIOS2
!! @param[inout] vr      redistributed vector
!! @param[inout] vi      transfer mark array (work array)
!! @param[in]    rsw     array size
!! @param[inout] lbuff   array size
      subroutine trunc_vec_transfer(vr, vi, rsw, lbuff,lglelr)
      implicit none

      include 'SIZE'
      include 'PARALLEL'
      include 'INPUT'

      ! argument list
      integer rsw, lbuff
      real vr(rsw,lbuff)
      integer vi(2,lbuff)
      integer lglelr(lbuff)

      ! local variables
      ! integer array size
      integer isw
      parameter (isw=2)
      integer*8 vl              ! required by crystal rauter
      integer lnelt,eg
      integer key(1)            ! required by crystal rauter; for sorting

      !integer gllnid
!-----------------------------------------------------------------------
      ! take number of local p4est elements
      lnelt = nelt

      ! single send
      ! pack integer array
      do eg=1,lnelt
         ! global element number
         vi(1,eg) = lglelr(eg)
         ! processor id
         vi(2,eg) = gllnid(lglelr(eg))

c         if (nio.eq.0) then
c             write(*,*) 'this element should be in rank: ', vi(2,eg)
c         endif
      enddo

      ! min aray size
      eg = lbuff
      ! transfer arrays
      call fgslib_crystal_tuple_transfer
     $     (cr_h,lnelt,eg,vi,isw,vl,0,vr,rsw,2)

      ! test local element number
      if (lnelt.ne.nelt) then
         !call amr_abort('Error: v1_transfer; lnelt /= nelt')
         write(*,*) 'Error in transfer, lnelt/=nelt'
      endif

      ! sort elements acording to their global number
      key(1) = 1
      call fgslib_crystal_tuple_sort
     $     (cr_h,lnelt,vi,isw,vl,0,vr,rsw,key,1)

      return
      end subroutine

!=======================================================================
!> @brief Redistribute single variable after reading from ADIOS2
!! @param[inout] vr      redistributed vector
!! @param[inout] vi      transfer mark array (work array)
!! @param[in]    rsw     array size
!! @param[inout] lbuff   array size
      subroutine trunc_filterq(v,f,nx,nz,w1,w2,ft,if3d,dmax)
c
      include 'SIZE'
      include 'TSTEP'

      real v(nx*nx*nz,nelt),w1(1),w2(1)
      logical if3d
c
      real f(nx,nx),ft(nx,nx)
c
      integer e
c
      call transpose(ft,nx,f,nx)
c
      nxyz=nx*nx*nz
      dmax = 0.

      if (nio.eq.0 .and. loglevel.gt.2) write(6,*) 'call filterq',ifield
      nel = nelfld(ifield)

      if (if3d) then
         do e=1,nel
c           Filter
            call copy(w2,v(1,e),nxyz)
            call mxm(f,nx,w2,nx,w1,nx*nx)
            i=1
            j=1
            do k=1,nx
               call mxm(w1(i),nx,ft,nx,w2(j),nx)
               i = i+nx*nx
               j = j+nx*nx
            enddo
            call mxm (w2,nx*nx,ft,nx,w1,nx)
            call sub3(w2,v(1,e),w1,nxyz)
            call copy(v(1,e),w1,nxyz)
            smax = vlamax(w2,nxyz)
            dmax = max(dmax,abs(smax))
         enddo
c
      else
         do e=1,nel
c           Filter
            call copy(w1,v(1,e),nxyz)
            call mxm(f ,nx,w1,nx,w2,nx)
            call mxm(w2,nx,ft,nx,w1,nx)
c
            call sub3(w2,v(1,e),w1,nxyz)
            call copy(v(1,e),w1,nxyz)
            smax = vlamax(w2,nxyz)
            dmax = max(dmax,abs(smax))
         enddo
      endif
c
      return
      end
