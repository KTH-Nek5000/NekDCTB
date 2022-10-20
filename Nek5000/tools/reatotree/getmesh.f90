!> @file getmesh.f90
!! @ingroup reatotree
!! @brief Subroutines to 
!! @author Adam Peplinski
!! @date 18 Oct 2019
!===============================================================================
module getmesh
  use paramm, only : isp, rdp, ndim, n_fcs, n_vrts, n_fvrts, n_edg
  use getmesh_tool, only : genconn_log
  implicit none
  
  private
  public :: makemesh, mesh_dealloc, mesh_type

  type mesh_type
     ! mesh variables
     integer(isp) :: nelv, nelt
     ! vertices
     integer(isp) :: crnk
     integer(isp), dimension(:), allocatable :: cell, crn2cell, crn2vrt, crn2off
     integer(isp), dimension(:), allocatable :: iperc
     integer(isp), dimension(:,:), allocatable :: pcrn
     ! faces
     integer(isp) :: frnk
     integer(isp), dimension(:), allocatable :: face, fcs2cell, fcs2fcs, fcs2off
     integer(isp), dimension(:,:), allocatable :: face_v ! face corners
     ! edges
#if N_DIM == 3
     integer(isp) :: ernk 
     integer(isp), dimension(:), allocatable :: edge, edg2cell, edg2edg, edg2off
     integer(isp), dimension(:,:), allocatable :: edge_v ! edge corners
#endif

     ! vertex coordinates
     real(rdp), dimension(:), allocatable :: vert
  
     ! boundarry conditions
     character(len=3), dimension(:,:), allocatable :: cbc
     real(rdp), dimension(:,:,:), allocatable :: bc

     ! curvature data
     character(len=1), dimension(:,:), allocatable :: ccrv
     real(rdp), dimension(:,:,:), allocatable :: crv

     ! element group
     integer(isp), dimension(:), allocatable :: igroup
  end type mesh_type
    
contains
  !===============================================================================
  !> @brief Read nekton .rea file and make a mesh
  !! @param[inout]  mesh      mesh variable
  !! @param[in]     session   session name
  subroutine makemesh(mesh,session)
    use getmesh_read, only : open_bin_file, get_dxyz, get_crv_bcs
    implicit none

    ! argument list
    type(mesh_type), intent(inout) :: mesh
    character(len=80), intent(in) :: session

    ! local variables
    integer(isp) :: ierr, el, fl
    logical :: ifbswap
    real(rdp) :: qo, qin  ! tolerance for vertex counting
    ! corners
    real(rdp), dimension(:,:), allocatable :: dx

    ! number of words per real
    integer(isp) :: wdsizi
    !-------------------------------------------------------------------------------
    ! open the file
    call genconn_log('Reading ###.re2')
    call open_bin_file(session,ifbswap,mesh%nelt,mesh%nelv,wdsizi)

    ! allocate arrays
    call mesh_alloc(mesh)
    allocate(dx(0:3,(size(mesh%cell))))

    ! read vertex data and get element characteristic length
    call get_dxyz (dx,mesh%igroup,mesh%nelt,wdsizi,ifbswap)

    ! read curved and boundary data
    call get_crv_bcs(mesh%cbc,mesh%bc,mesh%ccrv,mesh%crv,mesh%nelt,mesh%nelv,wdsizi,ifbswap)
    
    call byte_close(ierr)
    if(ierr/=0) call fsc_abort('Error closing ###.re2 in makemesh '//CHAR(0))

    ! get the tolerance
    call genconn_log('Input mesh tolerance (default 0.1):')
    call genconn_log('NOTE: smaller is better, but generous is more forgiving for bad meshes.')
    qo=0.1
    read(*,*) qin
    qin = abs(qin)
    if(qin>0) qo = qin

    ! Compress vertices based on coordinates
    call genconn_log('Unique vertices')
    call unique_vertex(qo,dx,mesh%nelt,mesh%crnk,mesh%cell,mesh%crn2cell,mesh%crn2vrt,mesh%crn2off)

    call genconn_log('Periodic boundary conditioons')
    call periodic_vtx(dx,mesh%nelt,mesh%crnk,mesh%cell,mesh%cbc,mesh%bc,mesh%iperc,mesh%pcrn)

    call genconn_log('Unique faces')
    call unique_face(mesh%cell,mesh%nelt,mesh%frnk,mesh%face,mesh%fcs2cell,mesh%fcs2fcs,mesh%fcs2off,mesh%face_v)

#if N_DIM == 3
    call genconn_log('Unique edges')
    call unique_edge(mesh%cell,mesh%nelt,mesh%ernk,mesh%edge,mesh%edg2cell,mesh%edg2edg,mesh%edg2off,mesh%edge_v)
#endif

    ! move vertex coordinates
    do el=1,mesh%crnk
       do fl= 1,3
          mesh%vert((el-1)*3 + fl) = dx(fl,el)
       enddo
    enddo

    deallocate(dx)

    return
  end subroutine makemesh
  !===============================================================================
  !> @brief Allocate arrays
  !! @param[inout]  mesh   mesh variable
  subroutine mesh_alloc(mesh)
    implicit none

    ! argument list
    type(mesh_type), intent(inout) :: mesh

    ! local variables
    integer(isp) :: vsize
    !-------------------------------------------------------------------------------
    vsize = n_vrts*mesh%nelt
    allocate(mesh%cell(vsize),mesh%crn2cell(vsize))
    allocate(mesh%crn2vrt(vsize),mesh%crn2off(vsize+1))
    allocate(mesh%iperc(vsize),mesh%pcrn(ndim,vsize))


    vsize = n_fcs*mesh%nelt
    allocate(mesh%face(vsize),mesh%fcs2cell(vsize))
    allocate(mesh%fcs2fcs(vsize),mesh%fcs2off(vsize+1))
    allocate(mesh%face_v(n_fvrts,vsize))

#if N_DIM == 3
    vsize = n_edg*mesh%nelt
    allocate(mesh%edge(vsize))
    allocate(mesh%edg2cell(vsize),mesh%edg2edg(vsize),mesh%edg2off(vsize+1))
    allocate(mesh%edge_v(2,vsize))
#endif

    vsize = 3*n_vrts*mesh%nelt
    allocate(mesh%vert(vsize))
    vsize = mesh%nelt
    allocate(mesh%cbc(6,vsize),mesh%bc(5,6,vsize),mesh%ccrv(12,vsize),mesh%crv(6,12,vsize))
    allocate(mesh%igroup(vsize))

    return
  end subroutine mesh_alloc
  !===============================================================================
  !> @brief Dellocate arrays
  !! @param[inout]  mesh   mesh variable
  subroutine mesh_dealloc(mesh)
    implicit none

    ! argument list
    type(mesh_type), intent(inout) :: mesh
    !-------------------------------------------------------------------------------

    deallocate(mesh%cell)
    deallocate(mesh%crn2cell,mesh%crn2vrt,mesh%crn2off)
    deallocate(mesh%iperc,mesh%pcrn)

    deallocate(mesh%face)
    deallocate(mesh%fcs2cell,mesh%fcs2fcs,mesh%fcs2off)
    deallocate(mesh%face_v)

#if N_DIM == 3
    deallocate(mesh%edge)
    deallocate(mesh%edg2cell,mesh%edg2edg,mesh%edg2off)
    deallocate(mesh%edge_v)
#endif
    
    deallocate(mesh%vert)
    deallocate(mesh%cbc,mesh%bc,mesh%ccrv,mesh%crv)
    deallocate(mesh%igroup)
    
    return
  end subroutine mesh_dealloc
  !===============================================================================
  !> @brief Find unique vertices
  subroutine unique_vertex(q,dx,nel,nglb,cell,crn2cell,crn2vrt, crn2off)
    use getmesh_tool, only : sort_tuple, swap_tuple, swap_ipi
    implicit none

    ! argument list
    integer(isp), intent(in) :: nel
    real(rdp), dimension(0:,:), intent(inout)  :: dx
    real(rdp), intent(in) :: q
    integer(isp), intent(out) :: nglb
    integer(isp), dimension(:), intent(out)  :: cell, crn2cell, crn2vrt, crn2off

    ! local variables
    integer(isp) :: n, i, j, lda, nseg, ipass, iseg, ig, ic, icm
    integer(isp), dimension(1) :: key
    integer(isp),dimension(:), allocatable :: ninseg, ind
    real(rdp) ::  qq
    real(rdp), dimension(4) :: dxt,t1,t2
    real(rdp),dimension(:,:), allocatable :: wk
    logical,dimension(:), allocatable :: ifseg
    character(len=100) :: str
    !-------------------------------------------------------------------------------
    ! allocate arrays
    ig = size(cell)
    allocate(ninseg(4*ig), ind(ig))
    allocate(wk(0:3,ig))
    allocate(ifseg(4*ig))
    
    
    n = n_vrts*nel

    write(str,*) 'start locglob_lexico vertex:',nel,n,q
    call genconn_log(str)

    qq = q*q  ! Square of relative tolerance

    do i=1,n
       cell(i)   = i
       ifseg (i) = .false.
    enddo

    ! Sort by directions
    lda         = 4
    nseg        = 1
    ifseg(1)    = .true.
    ninseg(1)   = n

    do ipass=1,2   ! Multiple passes eliminates false positives
       do j=1,ndim       ! Sort within each segment
          i =1
          key(1)=j+1
          do iseg=1,nseg
             call sort_tuple(dx(:,i:i+ninseg(iseg)-1),key,ind(1:ninseg(iseg)),dxt)
             ! Swap position
             call swap_ipi(cell(i:i+ninseg(iseg)-1),ind(1:ninseg(iseg)))
             i  =   i + ninseg(iseg)
          enddo

          do i=2,n
             if ((dx(j,i)-dx(j,i-1))**2>qq*min(dx(0,i),dx(0,i-1))) ifseg(i)=.true.
          enddo

          nseg = 0              !  Count up number of different segments
          do i=1,n
             if (ifseg(i)) then
                nseg = nseg+1
                ninseg(nseg) = 1
             else
                ninseg(nseg) = ninseg(nseg) + 1
             endif
          enddo
       enddo
    enddo

    ! create the corner to cell and vertex map
    crn2cell(:) = 0
    crn2vrt(:) = 0
    crn2off(:) = 0
    ! Assign global node numbers (now sorted lexigraphically!)
    ig  = 0
    ic  = 0
    icm = 0
    do i=1,n
       if (ifseg(i)) then
          ig=ig+1
          icm = max(ic,icm)
          crn2off(ig) = ic
          ic  = 0
       endif
       ic = ic+1     ! count number of instances at present ig
       ind(cell(i)) = ig
       crn2cell(i) = (cell(i)-1)/n_vrts+1
       crn2vrt(i) = mod(cell(i)-1,n_vrts)+1
    enddo
    crn2off(ig+1) = ic
    nglb = ig

    crn2off(1) = 1
    do i=2,nglb+1
       crn2off(i) = crn2off(i-1) + crn2off(i)
    enddo

    ! Unshuffle geometry:
    call swap_tuple(dx,cell,t1,t2)

    ! Reassign cell to hold global index numbering
    cell(1:n) = ind(1:n)
    call self_chk(cell,nel,0)       ! check for not self-ptg.

    ! Reassign geometry to match global index numbering
    ! Retain the geometry that is associated with the smallest bounding radius

    wk(:,1:n) = dx(:,1:n)
    ind(:) = 0 ! flag to see if value already assigned

    do i=1,n
       ig = cell(i)
       if (ind(ig)==0) then
          dx(:,ig) = wk(:,i)
          ind(ig) = 1
       elseif (wk(0,i)<dx(0,ig)) then
          dx(:,ig) = wk(:,i)
       endif
    enddo
    write(str,"(' done locglob_lexico vertex:',4i12)") nseg,nglb,n,icm
    call genconn_log(str)

#ifdef DEBUG
    ! testing
    open(unit=23,file='corners.txt')
    write(23,*) 'ndim, n_vrt, n_unique_pts'
    write(23,*) ndim,n,nglb
    write(23,*) 'dx'
    do i=1,nglb
       write(23,*) (dx(j,i),j=1,ndim)
    enddo
    write(23,*) 'crn2off'
    ic=crn2off(1)
    do i=1,nglb+1
       write(23,*)  crn2off(i)-ic,crn2off(i)
       ic = crn2off(i)
    enddo
    write(23,*) 'crn2cell,crn2vrt'
    do i=1,ic-1
       write(23,*) crn2cell(i),crn2vrt(i)
    enddo
    write(23,*) 'cell'
    do i=1,n
       write(23,*) cell(i)
    enddo
    close(23)
    ! testing end
#endif

    ! deallocate arrays
    deallocate(ninseg, ind)
    deallocate(wk)
    deallocate(ifseg)
    
    return
  end subroutine unique_vertex
  !===============================================================================
  !> @brief check for not self-ptg.
  subroutine self_chk(cell,nel,flag)
    implicit none

    ! argument list
    integer(isp), intent(in) :: nel, flag
    integer(isp), dimension(:), intent(in)  :: cell

    ! local variables
    integer(isp) :: el, il, jl
    character(len=100) :: str
!-------------------------------------------------------------------------------
    do el=1,nel
       do il=1,n_vrts
          do jl=il+1,n_vrts
             if (cell(il+(el-1)*n_vrts)==cell(jl+(el-1)*n_vrts)) then
                call genconn_log('Tighten mesh tolerance!')
                write(str,*) 'ABORT: SELF-CHK ',il,jl,el,flag
                call fsc_abort(str//CHAR(0))
             endif
          enddo
       enddo
    enddo

    return
  end subroutine self_chk
  !===============================================================================
  !> @brief Find corresponding verticies for periodic bounday condition
  !! @details Find corresponding verticies for periodic bounday condition
  !!     we use symmetric face and vertex numbering, cbc and bc were
  !!     corrected during the file reading phase
  subroutine periodic_vtx(dx,nel,irnk,cell,cbc,bc,iper,jmin)
    implicit none

    ! argument list
    integer(isp), intent(in) :: nel, irnk
    integer(isp), dimension(:), intent(in)  :: cell
    real(rdp), dimension(0:,:), intent(in)  :: dx
    character(len=3), dimension(:,:), intent(in)  :: cbc
    real(rdp), dimension(:,:,:), intent(inout)  :: bc
    integer(isp), dimension(:), intent(out)  :: iper
    integer(isp), dimension(:,:), intent(out)  :: jmin
    
    ! local variables
    character(len=3) ::  cb,cj
    integer(isp) :: e, f, je, jf, ke, kf
    character(len=100) :: str
    !-------------------------------------------------------------------------------
!!!!!!!!!! ADD SANITY CHECK !!!!!!!!!!!!!!!!!
    
    write(str,*) 'start periodic vtx:',nel,irnk
    call genconn_log(str)

    iper(:) = 0
    do e=1,irnk                  ! Initial permutation = identity
       do f=1,ndim
          jmin(f,e) = e
       enddo
    enddo

    ! This looks for periodic faces and corners related to them.
    ! Notice it does not give the whole set of periodic corners
    do e=1,nel
       do f=1,n_fcs
          cb = cbc(f,e)
          if (cb=='P  ') then
             je = int(abs(bc(1,f,e)),isp)
             jf = int(bc(2,f,e),isp)

             cj = cbc(jf,je)
             ke = int(abs(bc(1,jf,je)),isp)
             kf = int(bc(2,jf,je),isp)

             if (bc(1,f,e)>0 .and. bc(1,jf,je)>0) then
                if (ke/=e .or. kf/=f .or. cj/='P  ') then
                   call genconn_log('PERIODIC MISMATCH 1:')
                   write(str,6)   e,f,cb,' ie '
                   call genconn_log(str)
                   write(str,6) je,jf,cb,' je '
                   call genconn_log(str)
                   write(str,6) ke,kf,cj,' ke '
                   call genconn_log(str)
6                  format(i12,i3,1x,a3,1x,a4)
                   call fsc_abort('ABORT: PERIODIC MISMATCH'//CHAR(0))
                endif

                call find_connctd_pairs(iper,jmin,e,f,je,jf,cell,dx)

             elseif (bc(1,f,e)*bc(1,jf,je)<=0) then
                call genconn_log('PERIODIC MISMATCH 2:')
                write(str,6)   e,f,cb,' ie '
                call genconn_log(str)
                write(str,6) je,jf,cj,' je '
                call genconn_log(str)
               call fsc_abort('ABORT: PERIODIC MISMATCH'//CHAR(0))
            endif
         endif
      enddo
   enddo

   ! Reset bc array
   do e=1,nel
      do f=1,n_fcs
         if (cbc(f,e)=='P  ') bc(1,f,e) = abs(bc(1,f,e))
      enddo
   enddo

   write(str,*) 'done periodic vtx'
   call genconn_log(str)

   return
 end subroutine periodic_vtx
  !===============================================================================
  !> @brief Find periodic point pairs
  !! @details It has to be done per dimension, as points at the domain edges and
  !!     corners can correspond to multiple points.
  subroutine find_connctd_pairs(iper,jmin,e,f,je,jf,cell,dx)
    implicit none

    ! argument list
    integer(isp), dimension(:), intent(inout)  :: iper
    integer(isp), dimension(:,:), intent(inout)  :: jmin
    integer(isp), intent(in) :: e, f, je, jf
    integer(isp), dimension(:), intent(in)  :: cell
    real(rdp), dimension(0:,:), intent(in)  :: dx

    ! local variables
    real(rdp) :: x0(0:3,4),x1(0:3,4)

    integer(isp) :: shift, smin, i, j, k, j0, j1, ic, iv, jc, jv
    real(rdp) :: x0m, x0a, x1a, d2min, d2, eps, tol
    character(len=100) :: str

    ! circulant vertices on symm. faces, 3D
#if N_DIM == 3
    ! order ctr-clkws when looking at face
    integer(isp), parameter, dimension(4,6) :: vface = reshape(&
         &(/ 1,5,7,3 , 2,4,8,6 , 1,2,6,5 , 3,7,8,4 , 1,3,4,2 , 5,6,8,7 /),shape(vface))
#else
    ! circulant vertices on symm. faces, 2D
    integer(isp), parameter, dimension(4,6) ::  vface = reshape(&
         &(/ 3,1,0,0 , 2,4,0,0 , 1,2,0,0 , 4,3,0,0 , 0,0,0,0 , 0,0,0,0 /),shape(vface))
#endif
    !-------------------------------------------------------------------------------
!!!!!!!!!! ADD SANITY CHECK !!!!!!!!!!!!!!!!!

    do i=1,n_fvrts     ! Grab geometry for P-P face pair
       j0 = vface (i ,f)
       x0(:,i) = dx(:,cell(j0+(e-1)*n_vrts))
       j1 = vface (i ,jf)
       x1(:,i) = dx(:,cell(j1+(je-1)*n_vrts))
    enddo

    x0m = 0.
    do k=1,ndim           ! Subtract off mean of faces
       x0a = 0.
       x1a = 0.
       do i=1,n_fvrts
          x0a = x0a + x0(k,i)
          x1a = x1a + x1(k,i)
       enddo
       do i=1,n_fvrts
          x0(k,i) = x0(k,i) - x0a/n_fvrts
          x1(k,i) = x1(k,i) - x1a/n_fvrts
          ! this is modification with respect to original genmap, were x0m
          ! is the max of element position in the grid; I guess it is
          ! better that x0m is the max of element size
          x0m = max(x0m,abs(x0(k,i)))
          x0m = max(x0m,abs(x1(k,i)))
       enddo
    enddo

    smin = 0
    d2min = 1.e22
    do shift=0,n_fvrts-1
       d2 = 0.
       do i=1,n_fvrts
          j=i+shift
          if (j>n_fvrts) j=j-n_fvrts
          j=n_fvrts+1-j              ! go backward for j !
          do k=1,ndim
             d2 = d2 + (x0(k,i)-x1(k,j))**2
          enddo
       enddo

       if (d2<d2min) then
          smin  = shift
          d2min = d2
       endif
    enddo
    shift = smin
    if (d2min>0) d2min = sqrt(d2min)
    ! check correctness
    ! eps = 1.e-4
    eps = 1.e-3   !ADAM why eps is so big????????
    tol = eps*x0m
    if (d2min>tol) then
       write(str,6) e , f,shift,eps,x0m
       call genconn_log(str)
       write(str,6) je,jf,    i,tol,d2min
       call genconn_log(str)
       call fsc_abort('ABORT: FACE MISMATCH'//CHAR(0))
6      format(i12,i2,i3,2e16.8,' abort: FACE MATCH FAIL')
    endif

    ! dimension dependent part
    k = (f-1)/2+1
    do i=1,n_fvrts

       j=i+shift
       if (j>n_fvrts) j=j-n_fvrts
       j=n_fvrts+1-j              ! go backward for j to make faces match

       iv = vface(i, f)
       jv = vface(j,jf)

       ic = cell(iv+(e-1)*n_vrts)
       jc = cell(jv+(je-1)*n_vrts)

       ! check boundary condition consistency
       if ((.not.btest(iper(ic),k-1)).and.(.not.btest(iper(jc),k-1))) then
          jmin(k,ic) = jc
          jmin(k,jc) = ic
          ! mark vertices as periodic
          iper(ic) = ibset(iper(ic),k-1)
          iper(jc) = ibset(iper(jc),k-1)
       else
          if(jmin(k,ic)/=jc.or.jmin(k,jc)/=ic) then
             write(str,"(6i8,'ABORT: VRT MISMATCH')") jmin(k,ic),jc,jmin(k,jc),ic,e,f
             call fsc_abort(str//CHAR(0))
          endif
       endif

       ! to remove multiple points
       ! ijmin = jmin(ic)
       ! jjmin = jmin(jc)
       ! jmin(ic) = min(ijmin,jjmin)
       ! jmin(jc) = min(ijmin,jjmin)

    enddo

    return
  end subroutine find_connctd_pairs
 !===============================================================================
 !> @brief Unique face
 subroutine unique_face(cell,nel,nglb,face,fcs2cell,fcs2fcs, fcs2off,face_v)
   use getmesh_tool, only : sorti, swap_ipi, sort_tuple, swap_tuple
   implicit none

   ! argument list
   integer(isp), intent(in) :: nel
   integer(isp), intent(out) :: nglb
   integer(isp), dimension(:), intent(in)  :: cell
   integer(isp), dimension(:), intent(out)  :: face, fcs2cell, fcs2fcs, fcs2off
   integer(isp), dimension(:,:), intent(out)  :: face_v

   ! local variables
   integer(isp) :: nfcs, n   ! face number per element and total
   integer(isp) :: nfcvrt,nvtx   ! number of face and element corners
   integer(isp) :: i, j, j1 , k ,ipass     ! loop index
   integer(isp) :: nseg, iseg, ig, ic, icm
   integer(isp), dimension(1) :: key
   integer(isp), dimension(:,:), allocatable :: wk, indv ! store vertex ordering
   integer(isp), dimension(:), allocatable :: ninseg, ind, fcst, t1, t2
   logical(isp), dimension(:), allocatable :: ifseg ! is vertex sorted
   character(len=100) :: str

   ! face vertices
   integer(isp), parameter, dimension(4,6) ::  vface = reshape(&
        &(/ 1,3,5,7 , 2,4,6,8 , 1,2,5,6 , 3,4,7,8 , 1,2,3,4 , 5,6,7,8 /),shape(vface))
   !-------------------------------------------------------------------------------
   ! allocate arrays
   ig = size(face_v,1)
   ic = size(face_v,2)
   icm = size(face)
   allocate(wk(ig,ic),indv(ig,ic),ninseg(ig*ic),ind(icm),fcst(ig),t1(ig),t2(ig))
   allocate(ifseg(ig*ic))
   

   nfcs = 2*ndim
   n = nfcs*nel
   nvtx = 2**ndim
   nfcvrt = nvtx/2

   write(str,*) 'start locglob_lexico face:',nfcs,nel,n
   call genconn_log(str)

   ! indentify face by its unique corner (no periodicity at this point)
   j1 = 0
   do i=1,nel
      do j=1,nfcs
         j1= j1 +1
         do k=1, nfcvrt
            face_v(k,j1)=cell((i-1)*nvtx+vface(k,j))
         enddo
         ! sort face corners
         call sorti(face_v(:,j1),indv(:,j1))
      enddo
   enddo

   ! intial face numbering
   do i=1,n
      face(i)   = i
      ifseg (i) = .false.
   enddo

   ! Sort faces by corners
   nseg        = 1
   ifseg(1)    = .true.
   ninseg(1)   = n

   do ipass=1, 2  ! Multiple passes eliminates false positives
      do j=1,  nfcvrt    ! Sort within each segment
         i =1
         key(1) = j
         do iseg=1,nseg
            call sort_tuple(face_v(:,i:i+ninseg(iseg)-1),key,ind(1:ninseg(iseg)),fcst)
            ! Swap position
            call swap_ipi(face(i:i+ninseg(iseg)-1),ind(1:ninseg(iseg)))
            i  =   i + ninseg(iseg)
         enddo

         do i=2,n
            if (face_v(j,i)/=face_v(j,i-1)) ifseg(i)=.true.
         enddo

         nseg = 0              !  Count up number of different segments
         do i=1,n
            if (ifseg(i)) then
               nseg = nseg+1
               ninseg(nseg) = 1
            else
               ninseg(nseg) = ninseg(nseg) + 1
            endif
         enddo
      enddo
   enddo

!     testing
!      open(unit=23,file='face_v.txt')
!      write (23,*) n, nfcvrt
!      do i=1,n
!        write(23,*) (face_v(j,i),j=1,nfcvrt), ifseg(i), face(i)
!      enddo
!      close(23)
!     testing

   ! create the corner to cell and vertex map
   fcs2cell(:) = 0
   fcs2fcs(:) = 0
   fcs2off(:) = 0
   ! Assign global node numbers (now sorted lexigraphically!)
   ig  = 0
   ic  = 0
   icm = 0
   do i=1,n
      if (ifseg(i)) then
         ig=ig+1
         icm = max(ic,icm)
         fcs2off(ig) = ic
         ic  = 0
      endif
      ic = ic+1     ! count number of instances at present ig
      ind(face(i)) = ig
      fcs2cell(i) = (face(i)-1)/nfcs+1
      fcs2fcs(i) = mod(face(i)-1,nfcs)+1
   enddo
   fcs2off(ig+1) = ic
   nglb = ig

   fcs2off(1) = 1
   do i=2,nglb+1
      fcs2off(i) = fcs2off(i-1) + fcs2off(i)
   enddo

   ! check max number of connected faces
   if (icm>2) then
      write(str,"('ERROR: more than 2 faces connected',i10)") icm
      call fsc_abort(str//CHAR(0))
   endif

   ! Unshuffle faces:
   call swap_tuple(face_v,face,t1,t2)

   ! put back corners ordering
   do i=1,n
      call swap_ipi(face_v(:,i),indv(:,i))
   enddo

   ! Reassign cell to hold global index numbering
   face(1:n)=ind(1:n)

   ! Reassign geometry to match global index numbering
   ! Retain the geometry that is associated with the smallest bounding radius
   wk(:,1:n)=face_v(:,1:n)
   ind(:) = 0 ! flag to see if value already assigned

   do i=1,n
      ig = face(i)
      if (ind(ig)==0) then
         face_v(:,ig)=wk(:,i)
         ind(ig) = 1
      endif
   enddo
   write(str,"(' done locglob_lexico face:',4i10)") nseg,nglb,n,icm
   call genconn_log(str)

#ifdef DEBUG
   ! testing
   open(unit=23,file='face.txt')
   write(23,*) 'ndim, n_fcs, n_unique_fcs'
   write(23,*) ndim,n,nglb
   write(23,*) 'face_v'
   do i=1,nglb
      write(23,*) (face_v(j,i),j=1,nfcvrt)
   enddo
   write(23,*) 'fcs2off'
   ic=fcs2off(1)
   do i=1,nglb+1
      write(23,*)  fcs2off(i)-ic,fcs2off(i)
      ic = fcs2off(i)
   enddo
   write(23,*) 'fcs2cell,fcs2fcs'
   do i=1,ic-1
      write(23,*) fcs2cell(i),fcs2fcs(i)
   enddo
   write(23,*) 'face'
   do i=1,n
      write(23,*) face(i)
   enddo
   close(23)
#endif

   ! deallocate arrays
   deallocate(wk,indv,ninseg,ind,fcst,t1,t2)
   deallocate(ifseg)
   
   return
 end subroutine unique_face
 !===============================================================================
#if N_DIM == 3
 !> @brief Unique edge
 subroutine unique_edge(cell,nel,nglb,edge,edg2cell,edg2edg, edg2off, edge_v)
   use getmesh_tool, only : sorti, swap_ipi, sort_tuple, swap_tuple
   implicit none

   ! argument list
   integer(isp), intent(in) :: nel
   integer(isp), intent(out) :: nglb
   integer(isp), dimension(:), intent(in)  :: cell
   integer(isp), dimension(:), intent(out)  :: edge, edg2cell, edg2edg, edg2off
   integer(isp), dimension(:,:), intent(out)  :: edge_v

   ! local variables
   integer(isp) :: nfcs, n   ! face number per element and total
   integer(isp) :: nfcvrt,nvtx   ! number of face and element corners
   integer(isp) :: i, j, j1 , k ,ipass     ! loop index
   integer(isp) :: nseg, iseg, ig, ic, icm
   integer(isp), dimension(1) :: key
   integer(isp), dimension(2) :: fcst, t1, t2
   integer(isp), dimension(:), allocatable :: ninseg, ind
   integer(isp), dimension(:,:), allocatable :: wk, indv ! store vertex ordering
   logical(isp), dimension(:), allocatable :: ifseg ! is vertex sorted
   character(len=100) :: str

   ! face vertices
   integer(isp), parameter, dimension(2,n_edg) :: vface = reshape(&
        &(/1,2 ,3,4 ,5,6 ,7,8 ,1,3 ,2,4 ,5,7 ,6,8 ,1,5 ,2,6 ,3,7 ,4,8 /),shape(vface))
   !-------------------------------------------------------------------------------
   ! allocate arrays
   ig = size(edge)
   allocate(wk(2,ig),indv(2,ig),ninseg(2*ig),ind(ig))
   allocate(ifseg(2*ig))
   
   nfcs = n_edg
   n = nfcs*nel
   nvtx = 2**ndim
   nfcvrt = 2

   write(str,*) 'start locglob_lexico edge:',nfcs,nel,n
   call genconn_log(str)

   ! indentify edge by its unique corner (no periodicity at this point)
   j1 = 0
   do i=1,nel
      do j=1,nfcs
         j1= j1 +1
         do k=1, nfcvrt
            edge_v(k,j1)=cell((i-1)*nvtx+vface(k,j))
         enddo
         ! sort edge corners
         call sorti(edge_v(:,j1),indv(:,j1))
      enddo
   enddo

   ! intial edge numbering
   do i=1,n
      edge(i)   = i
      ifseg (i) = .false.
   enddo

   ! Sort faces by corners
   nseg        = 1
   ifseg(1)    = .true.
   ninseg(1)   = n

   do ipass=1, nfcvrt  ! Multiple passes eliminates false positives
      do j=1,  nfcvrt    ! Sort within each segment
         i =1
         key(1) = j
         do iseg=1,nseg
            call sort_tuple(edge_v(:,i:i+ninseg(iseg)-1),key,ind(1:ninseg(iseg)),fcst)
            ! Swap position
            call swap_ipi(edge(i:i+ninseg(iseg)-1),ind(1:ninseg(iseg)))
            i  =   i + ninseg(iseg)
         enddo

         do i=2,n
            if (edge_v(j,i)/=edge_v(j,i-1)) ifseg(i)=.true.
         enddo

         nseg = 0              !  Count up number of different segments
         do i=1,n
            if (ifseg(i)) then
               nseg = nseg+1
               ninseg(nseg) = 1
            else
               ninseg(nseg) = ninseg(nseg) + 1
            endif
         enddo
      enddo
   enddo

!     testing
!      open(unit=23,file='edge_v.txt')
!      write (23,*) n, nfcvrt
!      do i=1,n
!        write(23,*) (edge_v(j,i),j=1,nfcvrt), ifseg(i), edge(i)
!      enddo
!      close(23)
!     testing

   ! create the corner to cell and vertex map
   edg2cell(:) = 0
   edg2edg(:) = 0
   edg2off(:) = 0
   ! Assign global node numbers (now sorted lexigraphically!)
   ig  = 0
   ic  = 0
   icm = 0
   do i=1,n
      if (ifseg(i)) then
         ig=ig+1
         icm = max(ic,icm)
         edg2off(ig) = ic
         ic  = 0
      endif
      ic = ic+1     ! count number of instances at present ig
      ind(edge(i)) = ig
      edg2cell(i) = (edge(i)-1)/nfcs+1
      edg2edg(i) = mod(edge(i)-1,nfcs)+1
   enddo
   edg2off(ig+1) = ic
   nglb = ig

   edg2off(1) = 1
   do i=2,nglb+1
      edg2off(i) = edg2off(i-1) + edg2off(i)
   enddo

   ! Unshuffle faces:
   call swap_tuple(edge_v,edge,t1,t2)

   ! put back corners ordering
   do i=1,n
      call swap_ipi(edge_v(:,i),indv(:,i))
   enddo

   ! Reassign cell to hold global index numbering
   edge(1:n) =ind(1:n)

   ! Reassign geometry to match global index numbering
   ! Retain the geometry that is associated with the smallest bounding radius
   wk(:,1:n)=edge_v(:,1:n)
   ind(:) = 0      ! Flag to see if value already assigned

   do i=1,n
      ig = edge(i)
      if (ind(ig)==0) then
         edge_v(:,ig)=wk(:,i)
         ind(ig) = 1
      endif
   enddo
   write(str,"(' done locglob_lexico edge:',4i10)") nseg,nglb,n,icm
   call genconn_log(str)

#ifdef DEBUG
   ! testing
   open(unit=23,file='edge.txt')
   write(23,*) 'ndim, n_edgs, n_unique_edgs'
   write(23,*) ndim,n,nglb
   write(23,*) 'edge_v'
   do i=1,nglb
      write(23,*) (edge_v(j,i),j=1,nfcvrt)
   enddo
   write(23,*) 'edg2off'
   ic=edg2off(1)
   do i=1,nglb+1
      write(23,*)  edg2off(i)-ic,edg2off(i)
      ic = edg2off(i)
   enddo
   write(23,*) 'edg2cell,edg2edg'
   do i=1,ic-1
      write(23,*) edg2cell(i),edg2edg(i)
   enddo
   write(23,*) 'edge'
   do i=1,n
      write(23,*) edge(i)
   enddo
   close(23)
#endif

   ! deallocate arrays
   deallocate(wk,indv,ninseg,ind)
   deallocate(ifseg)
   
   return
 end subroutine unique_edge
#endif
 !===============================================================================
end module getmesh
  
