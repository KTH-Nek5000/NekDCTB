!> @file getmesh_read.f90
!! @ingroup reatotree
!! @brief Subroutines to read mesh datea from ###.re2 and ###.ma2 files
!! @author Adam Peplinski
!! @date 18 Oct 2019
!===============================================================================
module getmesh_read
  use paramm, only : isp, rdp, rsp, zero, ndim, n_vrts, n_fcs
  use getmesh_tool, only : genconn_log
  implicit none

  private
  public :: open_bin_file, get_dxyz, get_crv_bcs

contains
  !===============================================================================
  !> @brief Open nek5000 .re2 file and read header
  !! @param[in]  session   session name
  !! @param[out] ifbswap   big/little endian test
  !! @param[out] nelgt     global number of temperature elements
  !! @param[out] nelgv     global number of velocity elements
  !! @param[out] wdsizi    real size
  subroutine open_bin_file(session,ifbswap,nelgt,nelgv,wdsizi)
    implicit none

    ! argument list
    character(len=*),intent(in) :: session
    logical, intent(out) :: ifbswap
    integer(isp), intent(out) ::  nelgt, nelgv, wdsizi

    ! local variables
    integer(isp) :: ierr
    integer(isp) :: ndimr
    character(len=80) :: hdr
    character(len=5) :: version
    real(rsp) :: test
    !-------------------------------------------------------------------------------
    ! open the file
    ierr=0
    call byte_open(trim(session)//'.re2'//CHAR(0),ierr)
    if(ierr/=0) call fsc_abort('Error opening file ###.re2'//CHAR(0))

    ! read header
    call byte_read(hdr,20,ierr)
    if(ierr/=0) call fsc_abort('Error reading header from ###.re2'//CHAR(0))

    read (hdr,'(a5,i9,i3,i9)') version,nelgt,ndimr,nelgv
    wdsizi=4
    if(version=='#v002')wdsizi=8
    if(version=='#v003')wdsizi=8

    ! check grid dimension
    if (ndim/=ndimr) call fsc_abort ('Error: dimension missmatch.'//CHAR(0))
    
    ! read test varable to distinguish big/little endian
    call byte_read(test,1,ierr)
    if(ierr/=0) call fsc_abort('Error reading test number from ###.re2'//CHAR(0))
    call if_byte_swap_test(ifbswap, test)
    
    return
  end subroutine open_bin_file
  !===============================================================================
  !> @brief Checking big/little endian
  !! @param[out] ifbswap   big/little endian flag
  !! @param[in]  bytetest  real test pattern
  subroutine if_byte_swap_test(ifbswap, bytetest)
    implicit none

    ! argument list
    logical, intent(out) :: ifbswap
    real(rsp), intent(in) :: bytetest

    ! local variables
    integer(isp) :: ierr
    real(rsp) :: test2, etest
    real(rsp), parameter :: test_pattern=6.54321, eps=0.00020
    !-------------------------------------------------------------------------------
    etest = abs(test_pattern-bytetest)
    ifbswap = .true.
    if (etest<=eps) ifbswap = .false.

    ierr  = 0
    test2 = bytetest
    call byte_reverse(test2,1,ierr)
    if(ierr/=0) call fsc_abort('Error with byte_reverse in if_byte_swap_test '//CHAR(0))

    return
  end subroutine if_byte_swap_test
  !===============================================================================
  !> @brief Read coordinates and set characteristic length
  !! @param[out] dx        characteristic length and coordinates
  !! @param[out] igroup    element group
  !! @param[in]  nelt      global number of elements
  !! @param[in]  wdsizi    real size
  !! @param[in]  ifbswap   big/little endian test
  subroutine get_dxyz(dx,igroup,nelt,wdsizi,ifbswap)
    implicit none

    ! argument list
    integer(isp), intent(in) :: nelt, wdsizi
    real(rdp), dimension(0:,:), intent(out)  :: dx
    integer(isp), dimension(:), intent(out)  :: igroup
    logical, intent(in) :: ifbswap
    
    ! local variables
    ! number of words to read; group + 2x4 for 2d, 3x8 for 3d
    integer(isp), parameter :: nwds = 1 + ndim*(2**ndim)
    integer(isp) :: nwdsl

    integer(isp) :: ierr
    integer(isp) :: el, kl, ll
    integer(isp), dimension(0:54) :: buf
    real(rdp) :: bl

    ! to reorder vertices
    real(rdp), dimension(n_vrts) :: xc,yc,zc

    ! hypercube to strange ordering
    integer(isp), parameter, dimension(8) :: h2s=(/ 1,2,4,3,5,6,8,7 /)
    !-------------------------------------------------------------------------------
    call genconn_log('Reading coordinates')
    ! initialise variables
    nwdsl = nwds*(wdsizi/4)
    bl = huge(zero)
    ll = 1
    ierr = 0
    do el = 1, nelt
       call byte_read(buf,nwdsl,ierr)
       if(ierr/=0) call fsc_abort('Error reading xyz byte data in get_dxyz.'//CHAR(0))
       call buf_to_xyz(xc,yc,zc,igroup(el),buf,wdsizi,ifbswap)
       do kl=1,n_vrts
          dx(0,ll) = bl
          dx(1,ll) = xc(h2s(kl))
          dx(2,ll) = yc(h2s(kl))
          dx(3,ll) = zc(h2s(kl))
          ll = ll + 1
       enddo
    enddo

    ! verify right-handedness of given elements.
    call verify_rh(dx,nelt)

    ! find the characteristic size of each element
    call set_d2(dx,nelt)

    return
  end subroutine get_dxyz
  !===============================================================================
  !> @brief Copy element vertices and group fot read buffer to arrays
  !! @param[out]    xc,yc,zc   vertex coordinates
  !! @param[out]    igroup     element group
  !! @param[inout]  buf        read buffer
  !! @param[in]     wdsizi     real size in words
  !! @param[in]     ifbswap    little/big endian flagg
  subroutine buf_to_xyz(xc,yc,zc,igroup,buf,wdsizi,ifbswap)
    implicit none

    ! argument list
    real(rdp), dimension(:), intent(out)  :: xc,yc,zc
    integer(isp), intent(out) :: igroup
    integer(isp), dimension(0:), intent(inout)  :: buf
    integer(isp), intent(in) :: wdsizi
    logical, intent(in) :: ifbswap
    
    ! local variables
    ! number of words to read; group + 2x4 for 2d, 3x8 for 3d
    integer(isp), parameter :: nwds = 1 + ndim*(2**ndim), nwds2 = nwds*2
    
    integer(isp) :: ierr, il
    real(rsp) :: trrs
    real(rdp) :: trrd
    !-------------------------------------------------------------------------------
!!!!!! ADD SANITY CHECK !!!!!!
    ierr = 0

    if (ifbswap.and.wdsizi==8) then
       call byte_reverse8(buf,nwds2,ierr)
    elseif (ifbswap.and.wdsizi==4) then 
       call byte_reverse (buf,nwds,ierr)
    endif

    if (ierr/=0) call fsc_abort('Error byte_reverse in buf_to_xy '//CHAR(0))

    if(wdsizi==8) then
       igroup = int(transfer(buf(0:1),trrd),isp)
       if (ndim==3) then
          forall(il=1:n_vrts)
             xc(il) = transfer(buf(2*il:2*il+1),trrd)
             yc(il) = transfer(buf(16+2*il:16+2*il+1),trrd)
             zc(il) = transfer(buf(32+2*il:32+2*il+1),trrd)
          end forall
       else
          forall(il=1:n_vrts)
             xc(il) = transfer(buf(2*il:2*il+1),trrd)
             yc(il) = transfer(buf(8+2*il:8+2*il+1),trrd)
             zc(il) = zero
          end forall
       endif
    else 
       igroup = buf(0)
       if (ndim==3) then
          forall(il=1:n_vrts)
             xc(il) = transfer(buf(il),trrs)
             yc(il) = transfer(buf(8+il),trrs)
             zc(il) = transfer(buf(16+il),trrs)
          end forall
       else
          forall(il=1:n_vrts)
             xc(il) = transfer(buf(il),trrs)
             yc(il) = transfer(buf(4+il),trrs)
             zc(il) = zero
          end forall
       endif
    endif

    return
  end subroutine buf_to_xyz
  !===============================================================================
  !> @brief Verify right-handedness of given elements.
  !! @param[in]  dx      vertex coordinates
  !! @param[in]  nelt    element number
  subroutine verify_rh(dx,nelt)
    use getmesh_tool, only : crss2d, volum0
    implicit none

    ! argunent list
    integer(isp), intent(in) :: nelt
    real(rdp), dimension(0:,:), intent(in)  :: dx

    ! local variables
    integer(isp) :: ie, ip
    real(rdp) :: c1

    character(len=100) :: str
    !-------------------------------------------------------------------------------
!!!!!!! ADD SANITY CHECK !!!!!!!
    ! check vector product (in 2D) or element volume (in 3D)
    do ie=1, nelt
       ip = (ie-1)*n_vrts
       c1 = 1.0
#if N_DIM == 2
       ! CRSS2D(A,B,O) = (A-O) X (B-O)
       c1=min(c1,crss2d(dx(1:2,ip+2),dx(1:2,ip+3),dx(1:2,ip+1)))
       c1=min(c1,crss2d(dx(1:2,ip+4),dx(1:2,ip+1),dx(1:2,ip+2)))
       c1=min(c1,crss2d(dx(1:2,ip+1),dx(1:2,ip+4),dx(1:2,ip+3)))
       c1=min(c1,crss2d(dx(1:2,ip+3),dx(1:2,ip+2),dx(1:2,ip+4)))
#else
       ! VOLUM0(A,B,C,O) = (A-O)X(B-O).(C-O)
       c1= min(c1,volum0(dx(1:3,ip+2),dx(1:3,ip+3),dx(1:3,ip+5),dx(1:3,ip+1)))
       c1= min(c1,volum0(dx(1:3,ip+4),dx(1:3,ip+1),dx(1:3,ip+6),dx(1:3,ip+2)))
       c1= min(c1,volum0(dx(1:3,ip+1),dx(1:3,ip+4),dx(1:3,ip+7),dx(1:3,ip+3)))
       c1= min(c1,volum0(dx(1:3,ip+3),dx(1:3,ip+2),dx(1:3,ip+8),dx(1:3,ip+4)))
       c1= min(c1,-volum0(dx(1:3,ip+6),dx(1:3,ip+7),dx(1:3,ip+1),dx(1:3,ip+5)))
       c1= min(c1,-volum0(dx(1:3,ip+8),dx(1:3,ip+5),dx(1:3,ip+2),dx(1:3,ip+6)))
       c1= min(c1,-volum0(dx(1:3,ip+5),dx(1:3,ip+8),dx(1:3,ip+3),dx(1:3,ip+7)))
       c1= min(c1,-volum0(dx(1:3,ip+7),dx(1:3,ip+6),dx(1:3,ip+4),dx(1:3,ip+8)))
#endif
       if (c1<=0.0) then
          write(str,"('ERROR: non-right-handed element',i8)") ie
          call fsc_abort(str//CHAR(0))
       endif
    enddo

    return
  end subroutine verify_rh
  !===============================================================================
  !> @brief Get min edge size for given vertex
  !! @param[inout]  dx      vertex position/min edge length
  !! @param[in]     nelt    element number
  subroutine set_d2(dx,nelt)
    implicit none

    ! argument list
    integer(isp), intent(in) :: nelt
    real(rdp), dimension(0:,:), intent(inout)  :: dx

    ! local variables
    integer(isp) :: el,il,jl,kl,nl,ll

    ! symm. ordering
    integer(isp), parameter, dimension(3,8) :: neigh = reshape(&
         &(/ 2,3,5 , 1,4,6 , 1,4,7 , 2,3,8 , 1,6,7 , 2,5,8 , 3,5,8 , 4,6,7 /),shape(neigh))
    real(rdp) :: d2l
    !-------------------------------------------------------------------------------
    do el = 1,nelt
       ll = (el-1)*n_vrts
       do il = 1,n_vrts
          do kl = 1,ndim
             nl   = neigh(kl,il)
             d2l = zero
             do jl=1,ndim
                d2l = d2l + (dx(jl,nl+ll)-dx(jl,il+ll))**2
             enddo
             dx(0,il+ll) = min (dx(0,il+ll),d2l)
          enddo
       enddo
    enddo
    
    return
  end subroutine set_d2
  !===============================================================================
  !> @brief Get curvature and boundarry data
  !! @param[out]  cbc       character boundarry flag
  !! @param[out]  bc        real boundarry parameters
  !! @param[out]  ccrv      character curvature flag
  !! @param[out]  crv       real curvature parameter
  !! @param[in]   nelt      temperature element number
  !! @param[in]   nelv      velocity element number
  !! @param[in]   wdsizi    real size in words
  !! @param[in]   ifbswap   little/big endian flagg
  subroutine get_crv_bcs(cbc,bc,ccrv,crv,nelt,nelv,wdsizi,ifbswap)
    implicit none

    ! argument list
    integer(isp), intent(in) :: nelt,nelv,wdsizi
    character(len=3), dimension(:,:), intent(out)  :: cbc
    real(rdp), dimension(:,:,:), intent(out)  :: bc
    character(len=1), dimension(:,:), intent(out)  :: ccrv
    real(rdp), dimension(:,:,:), intent(out)  :: crv
    logical, intent(in) ::ifbswap

    ! local variables
    integer(isp) :: ierr, ncurve, nwds, kl, eg, jf
    integer(isp), dimension(55) :: buf
    real(rdp) :: rcurve

    character(len=3), dimension(6) :: cbt
    real(rdp), dimension(5,6) :: bt

    ! return Nekton preprocessor face ordering
    integer(isp), parameter, dimension(6) :: eface = (/ 4 , 2 , 1 , 3 , 5 , 6 /)

    ! return symmetric face ordering
    integer(isp), parameter, dimension(6) :: efaci = (/ 3 , 2 , 4 , 1 , 5 , 6 /)

    !-------------------------------------------------------------------------------
#ifdef DEBUG
    if (size(cbc,1)/=6.or.size(cbc,2)/=nelt) call fsc_abort('Error; get_crv_bcs cbc'//CHAR(0))
    if (size(bc,1)/=5.or.size(bc,2)/=6.or.size(bc,3)/=nelt) call fsc_abort('Error; get_crv_bcs bc'//CHAR(0))
    if (size(ccrv,1)/=12.or.size(ccrv,2)/=nelt) call fsc_abort('Error; get_crv_bcs ccrv'//CHAR(0))
    if (size(crv,1)/=6.or.size(crv,2)/=12.or.size(crv,3)/=nelt) call fsc_abort('Error; get_crv_bcs crv'//CHAR(0))
#endif
    call genconn_log('Reading curvature and BC')

    ! reset curvature data
    crv(:,:,:) = 0.0
    ccrv(:,:) = ' '

    ! reset bc data
    bc(:,:,:) = 0.0
    cbc(:,:) = '   '

    ierr = 0
    ! read curved side data
    if(wdsizi==8) then 
       call byte_read(rcurve,2,ierr)
       if (ifbswap) call byte_reverse8(rcurve,2,ierr)
       ncurve = int(rcurve,isp)
    else
       call byte_read(ncurve,1,ierr)
       if (ifbswap) call byte_reverse(ncurve,1,ierr)
    endif
    if(ierr/=0) call fsc_abort('Error reading ncurve in makemesh '//CHAR(0))

    nwds = (2 + 1 + 5)*(wdsizi/4) !eg+iside+ccurve+curve(6,:,:)
    
    do kl = 1,ncurve
       call byte_read(buf,nwds,ierr)
       if(ierr/=0) call fsc_abort('Error reading curve data in makemesh '//CHAR(0))
       call buf_to_curve(ccrv,crv,buf,nelt,wdsizi,ifbswap)
    enddo
  
    ! read bc for periodicity check; depending on nelv and nelt choose velocity or temperature one
    call rd_bc_bin(cbc,bc,buf,nelt,nelv,wdsizi,ifbswap)

    ! change ordering of faces
    do eg=1,nelt ! SWAP TO PREPROCESSOR NOTATION
       cbt(:) = cbc(:,eg)
       bt(:,:) = bc(:,:,eg) 
       do kl=1,n_fcs
          cbc(kl,eg) = cbt(eface(kl))
          bc(:,kl,eg) = bt(:,eface(kl))
       enddo
       ! renumber periodic faces
       do kl=1,n_fcs
          if (cbc(kl,eg)=='P  ') then
             jf = int(bc(2,kl,eg),isp)
             if (jf>0.and.jf<7) bc(2,kl,eg) = efaci(jf)
          endif
       enddo
    enddo
    
    return
  end subroutine get_crv_bcs
  !===============================================================================
  !> @brief Copy element vertices and group from read buffer to arrays
  !! @param[out]    ccrv       character curvature flag
  !! @param[out]    crv        real curvature parameter
  !! @param[inout]  buf        read buffer
  !! @param[in]     nelt       temperature element number
  !! @param[in]     wdsizi     real size in words
  !! @param[in]     ifbswap    little/big endian flagg
  subroutine buf_to_curve(ccrv,crv,buf,nelt,wdsizi,ifbswap)
    implicit none

    ! argument list
    integer(isp), intent(in) :: nelt, wdsizi
    character(len=1), dimension(:,:), intent(out)  :: ccrv
    real(rdp), dimension(:,:,:), intent(out)  :: crv
    integer(isp), dimension(:), intent(inout)  :: buf
    logical, intent(in) ::ifbswap

    ! local variables
    integer(isp) :: eg, fcs, ierr, nwds, il
    real(rsp) :: trrs
    real(rdp) :: trrd
    character(len=1) :: trch
    !-------------------------------------------------------------------------------
#ifdef DEBUG
    if (size(ccrv,1)/=12.or.size(ccrv,2)/=nelt) call fsc_abort('Error; buf_to_curve ccrv'//CHAR(0))
    if (size(crv,1)/=6.or.size(crv,2)/=12.or.size(crv,3)/=nelt) call fsc_abort('Error; buf_to_curve crv'//CHAR(0))
    if (size(buf)) call fsc_abort('Error; buf_to_curve buf'//CHAR(0))
#endif
    nwds = (2 + 1 + 5)*(wdsizi/4) !eg+iside+ccurve+curve(6,:,:)
    if(wdsizi==8) then
       if(ifbswap) call byte_reverse8(buf,nwds-2,ierr)
       eg = int(transfer(buf(1:2),trrd),isp)
       fcs = int(transfer(buf(3:4),trrd),isp)
       forall(il=1:5) crv(il,fcs,eg) = transfer(buf(3+2*il:4+2*il),trrd)
       trch = transfer(buf(15),trch)
       ccrv(fcs,eg) = trch
    else
       if (ifbswap) call byte_reverse(buf,nwds-1,ierr) ! last is char
       eg  = buf(1)
       fcs = buf(2)
       forall(il=1:5) crv(il,fcs,eg) = transfer(buf(2+il),trrs)
       trch = transfer(buf(8),trch)
       ccrv(fcs,eg) = trch
    endif

    ! sanity check
    if (eg>nelt) call fsc_abort ('Error: too big element number; buf_to_'//CHAR(0))

    return
  end subroutine buf_to_curve
  !===============================================================================
  !> @brief Read boundary data for velocity or temperature mesh
  !! @param[out]    cbc        character boundary flag
  !! @param[out]    bc         real boundary parameter
  !! @param[out]    buf        read buffer
  !! @param[in]     nelt       temperature element number
  !! @param[in]     nelv       velocity element number
  !! @param[in]     wdsizi     real size in words
  !! @param[in]     ifbswap    little/big endian flagg
  subroutine rd_bc_bin(cbc,bc,buf,nelt,nelv,wdsizi,ifbswap)
    implicit none

    ! argument list
    integer(isp), intent(in) :: nelt,nelv,wdsizi
    character(len=3), dimension(:,:), intent(out)  :: cbc
    real(rdp), dimension(:,:,:), intent(out)  :: bc
    integer(isp), dimension(:), intent(out)  :: buf
    logical, intent(in) ::ifbswap

    ! local variables
    integer(isp) :: npass, ierr, kpass, fcs, nwds, nbc_max, n8wds
    real(rdp) :: rbc_max
    !-------------------------------------------------------------------------------
#ifdef DEBUG
    if (size(cbc,1)/=6.or.size(cbc,2)/=nelt) call fsc_abort('Error; rd_bc_bin cbc'//CHAR(0))
    if (size(bc,1)/=5.or.size(bc,2)/=6.or.size(bc,3)/=nelt) call fsc_abort('Error; rd_bc_bin bc'//CHAR(0))
    if (size(buf)) call fsc_abort('Error; rd_bc_bin buf'//CHAR(0))
#endif
    npass = 1
    if (nelt>nelv) npass = 2   ! default to thermal topology (for now)

    ierr = 0
    do kpass = 1,npass

       ! fill up cbc w/ default
       cbc(:,:) = 'E  '

       nwds =(2 + 1 + 5)*(wdsizi/4)   ! eg + iside + cbc + bc(5,:,:)

       if(wdsizi==8) then
          call byte_read(rbc_max,2,ierr)
          if (ifbswap) call byte_reverse8(rbc_max,2,ierr) 
          nbc_max = int(rbc_max,isp)
          do fcs=1,nbc_max
             call byte_read(buf,nwds,ierr)
             n8wds=nwds/2
             if (ifbswap) call byte_reverse8(buf,nwds-2,ierr) ! last is char
             if(ierr/=0) call fsc_abort('Error reading byte bcs '//CHAR(0))
             call buf_to_bc(cbc,bc,buf,nelt,wdsizi)
          enddo
       else
          call byte_read(nbc_max,1,ierr)
          if (ifbswap) call byte_reverse(nbc_max,1,ierr) 
          do fcs=1,nbc_max
             call byte_read(buf,nwds,ierr)
             if (ifbswap) call byte_reverse(buf,nwds-1,ierr) ! last is char
             if(ierr/=0) call fsc_abort('Error reading byte bcs '//CHAR(0))
             call buf_to_bc(cbc,bc,buf,nelt,wdsizi)
          enddo
       endif

    enddo

    return
  end subroutine rd_bc_bin
  !===============================================================================
  !> @brief Copy boundary data from read buffer to arrays
  !! @param[inout]  cbc        character boundary flag
  !! @param[inout]  bc         real boundary parameter
  !! @param[inout]  buf        read buffer
  !! @param[in]     nelt       temperature element number
  !! @param[in]     wdsizi     real size in words
  subroutine buf_to_bc(cbc,bc,buf,nelt,wdsizi)
    implicit none

    ! argument list
    integer(isp), intent(in) :: nelt, wdsizi
    character(len=3), dimension(:,:), intent(inout)  :: cbc
    real(rdp), dimension(:,:,:), intent(inout)  :: bc
    integer(isp), dimension(:), intent(inout)  :: buf

    ! local variables
    integer(isp) :: eg, fcs, il
    real(rsp) :: trrs
    real(rdp) :: trrd
    character(len=3) :: trch
    !-------------------------------------------------------------------------------
#ifdef DEBUG
    if (size(cbc,1)/=6.or.size(cbc,2)/=nelt) call fsc_abort('Error; buf_to_bc cbc'//CHAR(0))
    if (size(bc,1)/=5.or.size(bc,2)/=6.or.size(bc,3)/=nelt) call fsc_abort('Error; buf_to_bc bc'//CHAR(0))
    if (size(buf)) call fsc_abort('Error; buf_to_bc buf'//CHAR(0))
#endif
    if(wdsizi==8) then
       eg = int(transfer(buf(1:2),trrd),isp)
       fcs = int(transfer(buf(3:4),trrd),isp)
       forall(il=1:5) bc(il,fcs,eg) = transfer(buf(3+2*il:4+2*il),trrd)
       trch = transfer(buf(15:16),trch)
       cbc(fcs,eg) = trch
    else
       eg = buf(1)
       fcs = buf(2)
       forall(il=1:5) bc(il,fcs,eg) = transfer(buf(2+il),trrs)
       trch = transfer(buf(8),trch)
       cbc(fcs,eg) = trch
       ! Integer assign of connecting periodic element
       if(nelt>=1000000.and.cbc(fcs,eg)=='P  ') bc(1,fcs,eg) = buf(3)
    endif

    ! sanity check
    if (eg>nelt) call fsc_abort ('Error: too big element number; buf_to_bc'//CHAR(0))
    if (fcs>n_fcs) call fsc_abort ('Error: wrong face number; buf_to_bc'//CHAR(0))
      
    return
  end subroutine buf_to_bc
  !===============================================================================
end module getmesh_read
