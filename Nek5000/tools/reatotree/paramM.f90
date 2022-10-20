!> @file paramM.f90
!! @ingroup reatotree
!! @brief Global parameter module
!! @author Adam Peplinski
!! @date 17 Oct 2019
!===============================================================================
module paramM
  implicit none

  public
  ! variable types
  ! integer types
  integer, parameter :: i1b = 1, i2b = 2, i4b = 4, i8b = 8
  integer, parameter :: idp = i8b    ! double precission integer kind
  integer, parameter :: isp = i4b    ! single precission integer kind
  ! real types
  integer, parameter :: rsp = 4, rdp = 8

  ! simple mathematical constants
  real(rdp), parameter :: pi=3.141592653589793238462643383279502884197_rdp
  real(rdp), parameter :: pio2=1.57079632679489661923132169163975144209858_rdp
  real(rdp), parameter :: twopi=6.283185307179586476925286766559005768394_rdp
  real(rdp), parameter :: sqrt2=1.41421356237309504880168872420969807856967_rdp
  real(rdp), parameter :: eulerg=0.5772156649015328606065120900824024310422_rdp
  real(rdp), parameter :: lnbs=2.7182818284590452353602874713526624977572_rdp
  real(rdp), parameter :: zero=0.0_rdp, one=1.0_rdp, two=2.0_rdp, three=3.0_rdp, four=4.0_rdp
  real(rdp), parameter :: half=0.5_rdp

  ! array sizes for mesh description
  ! grid dimension
  integer(isp), parameter :: ndim = N_DIM
  ! number of faces and vertices
  integer(isp), parameter :: n_fcs = 2*ndim, n_vrts = 2**ndim
  integer(isp), parameter :: n_fvrts=n_vrts/2
  ! number of edges
  integer(isp), parameter :: n_edg=12
  
end module paramM
