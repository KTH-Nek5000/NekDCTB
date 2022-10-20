!> @file getmesh_tool.f90
!! @ingroup reatotree
!! @brief Tools required by getmesh
!! @author Adam Peplinski
!! @date 18 Oct 2019
!===============================================================================
module getmesh_tool
  use paramm, only : isp, rdp, rsp
  implicit none

  private
  public :: genconn_log, sort_tuple, swap_tuple, sorti, swap_ipi, crss2d, volum0

  ! Generic interface
  interface sort_tuple
     module procedure sort_tupler, sort_tuplei
  end interface sort_tuple

  interface swap_tuple
     module procedure swap_ip_tupler, swap_ip_tuplei
  end interface swap_tuple
  
contains
  !===============================================================================
  !> @brief Logging routine
  !! @param[in]  logs   log message
  subroutine genconn_log(logs)
    use paramm, only : isp
    implicit none

    ! argument list
    character(len=*), intent(in) :: logs

    ! local variables
    integer(isp), save :: icalld= 0, pkg_id = 0
    character(len=*),parameter :: code_name='reatotree'//CHAR(0)
    !-------------------------------------------------------------------------------
    if (icalld==0) then
       icalld=1

       ! register package
       call fsc_pkg_reg(pkg_id,-1,code_name)
       ! call fsc_pkg_is_reg(3, ierr)
    endif

    call fsc_log (pkg_id, 1, 7, trim(logs)//CHAR(0))

    return
  end subroutine genconn_log
  !===============================================================================
  !> @brief Heap Sort (p 231 Num. Rec., 1st Ed.)
  subroutine sort_tupler(a,key,ind,aa)
    implicit none

    ! argument list
    real(rdp), dimension(:,:), intent(inout)  :: a
    real(rdp), dimension(:), intent(inout)  :: aa
    integer(isp), dimension(:), intent(out)  ::  ind
    integer(isp), dimension(:), intent(in)  :: key

    ! local variables
    integer(isp) :: il, jl, kl, ll, ir, ii
    integer(isp) :: lda, nn, nkey
    !-------------------------------------------------------------------------------
    ! sanity check
#ifdef DEBUG
    if (size(a,1)/=size(aa).or.size(a,2)/=size(ind)) &
         &call fsc_abort('Error; sort_tupler inconsistent sizes'//CHAR(0))
#endif

    ! array sizes
    lda = size(aa)
    nn = size(ind)
    nkey = size(key)
    
    ! sort
    forall(il=1:nn) ind(il)=il
    if (nn<=1) return
    ll=nn/2+1
    ir=nn
    do
       if (ll>1) then
          ll=ll-1
          aa(:) = a(:,ll)
          ii  = ind(ll)
       else
          aa(:) = a(:,ir)
          ii = ind(ir)
          a(:,ir) = a(:,1)
          ind(ir) = ind(1)
          ir=ir-1
          if (ir==1) then
             a(:,1) = aa(:)
             ind(1) = ii
             return
          endif
       endif
       il=ll
       jl=ll+ll
       if (jl<=ir) then
          do
             if (jl<ir) then
                if (iftuple_altbr(a(:,jl),a(:,jl+1),key)) jl=jl+1
             endif
             if (iftuple_altbr(aa,a(:,jl),key)) then
                forall(kl=1:lda) a(kl,il) = a(kl,jl)
                ind(il) = ind(jl)
                il=jl
                jl=jl+jl
             else
                jl=ir+1
             endif
             if(jl>ir) exit
          enddo
       endif
       a(:,il) = aa(:)
       ind(il) = ii
    enddo
  end subroutine sort_tupler
!===============================================================================
  function iftuple_altbr(a,b,key)   result(iftuple)
    implicit none

    real(rdp), dimension(:), intent(in)  :: a, b
    integer(isp), dimension(:), intent(in)  :: key

    logical iftuple
    integer(isp) :: il, kl
    integer(isp) :: nkey, nn
    !-------------------------------------------------------------------------------
    ! sanity check
#ifdef DEBUG
    if (size(a)/=size(b)) call fsc_abort('Error; iftuple_altbr inconsistent sizes'//CHAR(0))
#endif

    ! array sizes
    nn = size(a)
    nkey = size(key)
    
    do il=1,nkey
       kl=key(il)
       if (a(kl)<b(kl)) then
          iftuple = .true.
          return
       elseif (a(kl)>b(kl)) then
          iftuple = .false.
          return
       endif
    enddo
    iftuple = .false.
    return
  end function iftuple_altbr
  !===============================================================================
  !> @brief Heap Sort (p 231 Num. Rec., 1st Ed.)
  subroutine sort_tuplei(a,key,ind,aa)
    implicit none

    ! argument list
    integer(isp), dimension(:,:), intent(inout)  :: a
    integer(isp), dimension(:), intent(inout)  :: aa
    integer(isp), dimension(:), intent(out)  ::  ind
    integer(isp), dimension(:), intent(in)  :: key

    ! local variables
    integer(isp) :: il, jl, kl, ll, ir, ii
    integer(isp) :: lda, nn, nkey
    !-------------------------------------------------------------------------------
    ! sanity check
#ifdef DEBUG
    if (size(a,1)/=size(aa).or.size(a,2)/=size(ind)) &
         &call fsc_abort('Error; sort_tuplei inconsistent sizes'//CHAR(0))
#endif

    ! array sizes
    lda = size(aa)
    nn = size(ind)
    nkey = size(key)
    
    ! sort
    forall(il=1:nn) ind(il)=il
    if (nn<=1) return
    ll=nn/2+1
    ir=nn
    do
       if (ll>1) then
          ll=ll-1
          aa(:) = a(:,ll)
          ii  = ind(ll)
       else
          aa(:) = a(:,ir)
          ii = ind(ir)
          a(:,ir) = a(:,1)
          ind(ir) = ind(1)
          ir=ir-1
          if (ir==1) then
             a(:,1) = aa(:)
             ind(1) = ii
             return
          endif
       endif
       il=ll
       jl=ll+ll
       if (jl<=ir) then
          do
             if (jl<ir) then
                if (iftuple_altbi(a(:,jl),a(:,jl+1),key)) jl=jl+1
             endif
             if (iftuple_altbi(aa,a(:,jl),key)) then
                forall(kl=1:lda) a(kl,il) = a(kl,jl)
                ind(il) = ind(jl)
                il=jl
                jl=jl+jl
             else
                jl=ir+1
             endif
             if(jl>ir) exit
          enddo
       endif
       a(:,il) = aa(:)
       ind(il) = ii
    enddo
  end subroutine sort_tuplei
!===============================================================================
  function iftuple_altbi(a,b,key)   result(iftuple)
    implicit none

    integer(isp), dimension(:), intent(in)  :: a, b
    integer(isp), dimension(:), intent(in)  :: key

    logical iftuple
    integer(isp) :: il, kl
    integer(isp) :: nkey, nn
    !-------------------------------------------------------------------------------
    ! sanity check
#ifdef DEBUG
    if (size(a)/=size(b)) call fsc_abort('Error; iftuple_altbi inconsistent sizes'//CHAR(0))
#endif

    ! array sizes
    nn = size(a)
    nkey = size(key)
    
    do il=1,nkey
       kl=key(il)
       if (a(kl)<b(kl)) then
          iftuple = .true.
          return
       elseif (a(kl)>b(kl)) then
          iftuple = .false.
          return
       endif
    enddo
    iftuple = .false.
    return
  end function iftuple_altbi
  !===============================================================================
  !> @brief Heap Sort (p 231 Num. Rec., 1st Ed.)
  subroutine sorti(a,ind)
    implicit none

    ! argument list
    integer(isp), dimension(:), intent(inout)  :: a, ind

    ! local variables
    integer(isp) :: aa, il, jl, ll, ir, ii
    integer(isp) :: nn
    !-------------------------------------------------------------------------------
    ! sanity check
#ifdef DEBUG
    if (size(a)/=size(ind)) call fsc_abort('Error; sorti inconsistent sizes'//CHAR(0))
#endif

    ! array sizes
    nn = size(a)
    
    forall(il=1:nn) ind(il)=il

    if (nn<=1) return
    ll=nn/2+1
    ir=nn
    do
       if (ll>1) then
          ll=ll-1
          aa  = a  (ll)
          ii  = ind(ll)
       else
          aa =   a(ir)
          ii = ind(ir)
          a(ir) = a(1)
          ind(ir) = ind(1)
          ir=ir-1
          if (ir==1) then
             a(1) = aa
             ind(1) = ii
             return
          endif
       endif
       il=ll
       jl=ll+ll
       if (jl<=ir) then
          do
             if (jl<ir) then
                if ( a(jl)<a(jl+1) ) jl=jl+1
             endif
             if (aa<a(jl)) then
                a(il) = a(jl)
                ind(il) = ind(jl)
                il=jl
                jl=jl+jl
             else
                jl=ir+1
             endif
             if(jl>ir) exit
          enddo
       endif
       a(il) = aa
       ind(il) = ii
    enddo
  end subroutine sorti
  !===============================================================================
  !> @brief In-place permutation: x' = x(p)
  subroutine swap_ipi(x,p)
    implicit none
    ! argument list
    integer(isp), dimension(:), intent(inout)  :: x, p
    !local variables
    integer(isp) :: j, k, loop_start, last, next, xstart, nn
    character(len=100) :: str
    !-------------------------------------------------------------------------------
    ! sanity check
#ifdef DEBUG
    if (size(x)/=size(p)) call fsc_abort('Error; swap_ipi inconsistent sizes'//CHAR(0))
#endif

    ! array sizes
    nn = size(x)
    
    do k=1,nn
       if (p(k)>0) then   ! not swapped
          xstart     = x(k)
          loop_start = k
          last       = k
          do j=k,nn
             next    = p(last)
             if (next<0) then
                write (str,*) 'Hey! iswap_ip problem.',j,k,nn,next
                call fsc_abort(str//CHAR(0))
             elseif (next==loop_start) then
                x(last) = xstart
                p(last) = -p(last)
                exit
             else
                x(last) = x(next)
                p(last) = -p(last)
                last    = next
             endif
          enddo
       endif
    enddo
    p = -p

    return
  end subroutine swap_ipi
  !===============================================================================
  !> @brief In-place permutation: x'(p) = x
  subroutine swap_ip_tupler(x,p,t1,t2)
    implicit none

    real(rdp), dimension(:,:), intent(inout)  :: x
    real(rdp), dimension(:), intent(inout)  :: t1, t2
    integer(isp), dimension(:), intent(inout)  :: p

    integer(isp) :: j, k, loop_start, next, nextp
    character(len=100) :: str
    integer(isp) :: m, n
    !-------------------------------------------------------------------------------
    ! sanity check
#ifdef DEBUG
    if (size(x,1)/=size(t1).or.size(x,1)/=size(t2).or.size(x,2)/=size(p)) &
         &call fsc_abort('Error; swap_ip_tupler inconsistent sizes'//CHAR(0))
#endif

    ! array sizes
    m = size(t1)
    n = size(p)
    
    do k=1,n
       if (p(k)>0) then   ! not swapped
          loop_start = k
          next = p(loop_start)
          t1(:) = x(:,loop_start)
          do j=1,n
             if (next<0) then
                write(str,*) 'Hey! iswapt_ip problem.',j,k,n,next
                call fsc_abort(str//CHAR(0))
             elseif (next==loop_start) then
                x(:,next) = t1(:)
                p(next) = -p(next)
                exit
             else
                t2(:) = x(:,next)
                x(:,next) = t1(:)
                t1 = t2
                nextp   =  p(next)
                p(next) = -p(next)
                next    =  nextp
             endif
          enddo
       endif
    enddo
    p = -p

    return
  end subroutine swap_ip_tupler
  !===============================================================================
  !> @brief In-place permutation: x'(p) = x
  subroutine swap_ip_tuplei(x,p,t1,t2)
    implicit none

    integer(isp), dimension(:,:), intent(inout)  :: x
    integer(isp), dimension(:), intent(inout)  :: t1, t2
    integer(isp), dimension(:), intent(inout)  :: p

    integer(isp) :: j, k, loop_start, next, nextp
    character(len=100) :: str
    integer(isp) :: m, n
    !-------------------------------------------------------------------------------
    ! sanity check
#ifdef DEBUG
    if (size(x,1)/=size(t1).or.size(x,1)/=size(t2).or.size(x,2)/=size(p)) &
         &call fsc_abort('Error; swap_ip_tuplei inconsistent sizes'//CHAR(0))
#endif

    ! array sizes
    m = size(t1)
    n = size(p)
    
    do k=1,n
       if (p(k)>0) then   ! not swapped
          loop_start = k
          next = p(loop_start)
          t1(:) = x(:,loop_start)
          do j=1,n
             if (next<0) then
                write(str,*) 'Hey! iswapt_ip problem.',j,k,n,next
                call fsc_abort(str//CHAR(0))
             elseif (next==loop_start) then
                x(:,next) = t1(:)
                p(next) = -p(next)
                exit
             else
                t2(:) = x(:,next)
                x(:,next) = t1(:)
                t1 = t2
                nextp   =  p(next)
                p(next) = -p(next)
                next    =  nextp
             endif
          enddo
       endif
    enddo
    p = -p

    return
  end subroutine swap_ip_tuplei
  !===============================================================================
  !> @brief Calculte vector product in 2d
  !! @param[in]   xy1,xy2,xy0  vertex coordinates
  function crss2d(xy1,xy2,xy0) result(crss)
    implicit none

    ! argument list
    real(rdp),dimension(:)  :: xy1, xy2, xy0

    ! local variables
    real(rdp) :: crss
    real(rdp) :: v1x, v2x, v1y, v2y
    !-------------------------------------------------------------------------------
#ifdef DEBUG
    if (size(xy1)/=2.or.size(xy2)/=2.or.size(xy0)/=2) call fsc_abort('Error; cross2d'//CHAR(0))
#endif
    v1x=xy1(1)-xy0(1)
    v2x=xy2(1)-xy0(1)
    v1y=xy1(2)-xy0(2)
    v2y=xy2(2)-xy0(2)
    crss = v1x*v2y - v1y*v2x

    return
  end function crss2d
  !===============================================================================
  !> @brief Calculate element volume
  !! @param[in]  p1,p2,p3,p0   vertex coordinates
  function volum0(p1,p2,p3,p0) result(volum)
    implicit none

    ! argument list
    real(rdp), dimension(:)  :: p1, p2, p3, p0

    !     local variables
    real(rdp) :: volum
    real(rdp) :: u1, u2, u3, v1, v2, v3, w1, w2, w3, cross1, cross2, cross3
    !-------------------------------------------------------------------------------
#ifdef DEBUG
    if (size(p1)/=3.or.size(p2)/=3.or.size(p3)/=3.or.size(p0)/=3) call fsc_abort('Error; volum0'//CHAR(0))
#endif
    u1=p1(1)-p0(1)
    u2=p1(2)-p0(2)
    u3=p1(3)-p0(3)

    v1=p2(1)-p0(1)
    v2=p2(2)-p0(2)
    v3=p2(3)-p0(3)

    w1=p3(1)-p0(1)
    w2=p3(2)-p0(2)
    w3=p3(3)-p0(3)

    cross1 = u2*v3-u3*v2
    cross2 = u3*v1-u1*v3
    cross3 = u1*v2-u2*v1

    volum  = w1*cross1 + w2*cross2 + w3*cross3

    return
  end function volum0
  !===============================================================================
end module getmesh_tool
