!> @file amr_shm.f
!! @ingroup nekamr
!! @brief Set of routines to use splitting of mpi ranks between nodes
!! @author Adam Peplinski
!! @date Jun 1, 2018
!=======================================================================
!> @brief Global nonblocking global vector commutative operation
!! @param[inout]   vs   send buffer
!! @param[inout]   vr   receive buffer
!! @param[in]      op   operator
!! @param[in]      nl   array length
!! @return       amr_glb_igopi  communication request
      integer function amr_glb_igopi(vs,vr,nl,op)
      implicit none

      include 'mpif.h'
      include 'SIZE'
      include 'PARALLEL'
      include 'AMR_SHM'

      ! argument list
      integer nl
      integer vs(nl), vr(nl)
      character*3 op

      ! local variables
      integer imsg, ierr
!-----------------------------------------------------------------------
      if     (op.eq.'+  ') then
        call mpi_iallreduce (vs,vr,nl,mpi_integer,mpi_sum ,
     $   AMR_GLBCOMM,imsg,ierr)
      elseif (op.EQ.'M  ') then
        call mpi_iallreduce (vs,vr,nl,mpi_integer,mpi_max ,
     $   AMR_GLBCOMM,imsg,ierr)
      elseif (op.EQ.'m  ') then
        call mpi_iallreduce (vs,vr,nl,mpi_integer,mpi_min ,
     $   AMR_GLBCOMM,imsg,ierr)
      elseif (op.EQ.'*  ') then
        call mpi_iallreduce (vs,vr,nl,mpi_integer,mpi_prod,
     $   AMR_GLBCOMM,imsg,ierr)
      else
        call amr_abort('Error: amr_glb_igopi; unknown op')
      endif

      amr_glb_igopi = imsg

      return
      end function
!=======================================================================
!> @brief Global nonblocking global vector inclusive scan
!! @param[inout]   vs   send buffer
!! @param[inout]   vr   receive buffer
!! @param[in]      op   operator
!! @param[in]      nl   array length
!! @return       amr_glb_iscani  communication request
      integer function amr_glb_iscani(vs,vr,nl,op)
      implicit none

      include 'mpif.h'
      include 'SIZE'
      include 'PARALLEL'
      include 'AMR_SHM'

      ! argument list
      integer nl
      integer vs(nl), vr(nl)
      character*3 op

      ! local variables
      integer imsg, ierr
!-----------------------------------------------------------------------
      if     (op.eq.'+  ') then
        call mpi_iscan (vs,vr,nl,mpi_integer,mpi_sum ,
     $   AMR_GLBCOMM,imsg,ierr)
      elseif (op.EQ.'M  ') then
        call mpi_iscan (vs,vr,nl,mpi_integer,mpi_max ,
     $   AMR_GLBCOMM,imsg,ierr)
      elseif (op.EQ.'m  ') then
        call mpi_iscan (vs,vr,nl,mpi_integer,mpi_min ,
     $   AMR_GLBCOMM,imsg,ierr)
      elseif (op.EQ.'*  ') then
        call mpi_iscan (vs,vr,nl,mpi_integer,mpi_prod,
     $   AMR_GLBCOMM,imsg,ierr)
      else
        call amr_abort('Error: amr_glb_iscani; unknown op')
      endif

      amr_glb_iscani = imsg

      return
      end function
!=======================================================================
!> @brief Global nonblocking broadcast for integer array
!! @param[in]   vx       communication buffer
!! @param[in]   nl       array length
!! @param[in]   isrc     root rank
!! @return       amr_glb_ibcsi  communication request
      integer function amr_glb_ibcsi(vx,nl,isrc)
      implicit none

      include 'mpif.h'
      include 'SIZE'
      include 'PARALLEL'
      include 'AMR_SHM'

      ! argument list
      integer nl, isrc
      integer vx(nl)

      ! local variables
      integer nvxb, imsg, ierr
!-----------------------------------------------------------------------
      nvxb = nl*ISIZE

      call mpi_ibcast (vx,nvxb,mpi_byte,isrc,AMR_GLBCOMM,imsg,ierr)
      amr_glb_ibcsi = imsg

      return
      end function
!=======================================================================
!> @brief Global nonblocking broadcast for real array
!! @param[in]   vx       communication buffer
!! @param[in]   nl       array length
!! @param[in]   isrc     root rank
!! @return       amr_glb_ibcsr  communication request
      integer function amr_glb_ibcsr(vx,nl,isrc)
      implicit none

      include 'mpif.h'
      include 'SIZE'
      include 'PARALLEL'
      include 'AMR_SHM'

      ! argument list
      integer nl, isrc
      real vx(nl)

      ! local variables
      integer nvxb, imsg, ierr
!-----------------------------------------------------------------------
      nvxb = nl*WDSIZE

      call mpi_ibcast (vx,nvxb,mpi_byte,isrc,AMR_GLBCOMM,imsg,ierr)
      amr_glb_ibcsr = imsg

      return
      end function
!=======================================================================
!> @brief Global nonblocking allgather for integer array
!! @param[in]   vs       send buffer
!! @param[in]   ns       send buffer length
!! @param[in]   vr       receive buffer
!! @param[in]   nr       receive buffer length
!! @return       amr_glb_iagti  communication request
      integer function amr_glb_iagti(vs,ns,vr,nr)
      implicit none

      include 'mpif.h'
      include 'SIZE'
      include 'PARALLEL'
      include 'AMR_SHM'

      ! argument list
      integer ns, nr
      integer vs(ns), vr(nr)

      ! local variables
      integer imsg, ierr
!-----------------------------------------------------------------------
      ! possible check
      !if (nr.lt.ns*NP)

      call mpi_iallgather (vs,ns,mpi_integer,vr,ns,mpi_integer,
     $     AMR_GLBCOMM,imsg,ierr)
      amr_glb_iagti = imsg

      return
      end function
!=======================================================================
!> @brief Internode nonblocking global vector commutative operation
!! @param[inout]   vs   send buffer
!! @param[inout]   vr   receive buffer
!! @param[in]      op   operator
!! @param[in]      nl   array length
!! @return       amr_nds_igopi  communication request
      integer function amr_nds_igopi(vs,vr,nl,op)
      implicit none

      include 'mpif.h'
      include 'SIZE'
      include 'PARALLEL'
      include 'AMR_SHM'

      ! argument list
      integer nl
      integer vs(nl), vr(nl)
      character*3 op

      ! local variables
      integer imsg, ierr
!-----------------------------------------------------------------------
      if     (op.eq.'+  ') then
        call mpi_iallreduce (vs,vr,nl,mpi_integer,mpi_sum ,
     $   AMR_NDSCOMM,imsg,ierr)
      elseif (op.EQ.'M  ') then
        call mpi_iallreduce (vs,vr,nl,mpi_integer,mpi_max ,
     $   AMR_NDSCOMM,imsg,ierr)
      elseif (op.EQ.'m  ') then
        call mpi_iallreduce (vs,vr,nl,mpi_integer,mpi_min ,
     $   AMR_NDSCOMM,imsg,ierr)
      elseif (op.EQ.'*  ') then
        call mpi_iallreduce (vs,vr,nl,mpi_integer,mpi_prod,
     $   AMR_NDSCOMM,imsg,ierr)
      else
        call amr_abort('Error: amr_nds_igopi; unknown op')
      endif

      amr_nds_igopi = imsg

      return
      end function
!=======================================================================
!> @brief Internode nonblocking receive for integer array
!! @param[in]   msgtag   message tag
!! @param[in]   vx       communication buffer
!! @param[in]   nl       array length
!! @param[in]   isrc     source rank
!! @return       amr_nds_ircvi  communication request
      integer function amr_nds_ircvi(msgtag,vx,nl,isrc)
      implicit none

      include 'mpif.h'
      include 'SIZE'
      include 'PARALLEL'
      include 'AMR_SHM'

      ! argument list
      integer msgtag, nl,isrc
      integer vx(nl)

      ! local variables
      integer nvxb, imsg, ierr
!-----------------------------------------------------------------------
      nvxb = nl*ISIZE

      call mpi_irecv (vx,nvxb,mpi_byte,isrc,msgtag
     $       ,AMR_NDSCOMM,imsg,ierr)
      amr_nds_ircvi = imsg

      return
      end function
!=======================================================================
!> @brief Internode nonblocking receive for real array
!! @param[in]   msgtag   message tag
!! @param[in]   vx       communication buffer
!! @param[in]   nl       array length
!! @param[in]   isrc     source rank
!! @return       amr_nds_ircvr  communication request
      integer function amr_nds_ircvr(msgtag,vx,nl,isrc)
      implicit none

      include 'mpif.h'
      include 'SIZE'
      include 'PARALLEL'
      include 'AMR_SHM'

      ! argument list
      integer msgtag, nl,isrc
      real vx(nl)

      ! local variables
      integer nvxb, imsg, ierr
!-----------------------------------------------------------------------
      nvxb = nl*WDSIZE

      call mpi_irecv (vx,nvxb,mpi_byte,isrc,msgtag
     $       ,AMR_NDSCOMM,imsg,ierr)
      amr_nds_ircvr = imsg

      return
      end function
!=======================================================================
!> @brief Internode nonblocking send for integer array
!! @param[in]   msgtag   message tag
!! @param[in]   vx       communication buffer
!! @param[in]   nl       array length
!! @param[in]   idst     destination rank
!! @return       amr_nds_isndi  communication request
      integer function amr_nds_isndi(msgtag,vx,nl,idst)
      implicit none

      include 'mpif.h'
      include 'SIZE'
      include 'PARALLEL'
      include 'AMR_SHM'

      ! argument list
      integer msgtag, nl, idst
      integer vx(nl)

      ! local variables
      integer nvxb, imsg, ierr
!-----------------------------------------------------------------------
      nvxb = nl*ISIZE

      call mpi_isend (vx,nvxb,mpi_byte,idst,msgtag
     $       ,AMR_NDSCOMM,imsg,ierr)
      amr_nds_isndi = imsg

      return
      end function
!=======================================================================
!> @brief Internode nonblocking send for real array
!! @param[in]   msgtag   message tag
!! @param[in]   vx       communication buffer
!! @param[in]   nl       array length
!! @param[in]   idst     destination rank
!! @return       amr_nds_isndr  communication request
      integer function amr_nds_isndr(msgtag,vx,nl,idst)
      implicit none

      include 'mpif.h'
      include 'SIZE'
      include 'PARALLEL'
      include 'AMR_SHM'

      ! argument list
      integer msgtag, nl, idst
      real vx(nl)

      ! local variables
      integer nvxb, imsg, ierr
!-----------------------------------------------------------------------
      nvxb = nl*WDSIZE

      call mpi_isend (vx,nvxb,mpi_byte,idst,msgtag
     $       ,AMR_NDSCOMM,imsg,ierr)
      amr_nds_isndr = imsg

      return
      end function
!=======================================================================
!> @brief Internode nonblocking broadcast for integer array
!! @param[in]   vx       communication buffer
!! @param[in]   nl       array length
!! @param[in]   isrc     root rank
!! @return       amr_nds_ibcsi  communication request
      integer function amr_nds_ibcsi(vx,nl,isrc)
      implicit none

      include 'mpif.h'
      include 'SIZE'
      include 'PARALLEL'
      include 'AMR_SHM'

      ! argument list
      integer nl, isrc
      integer vx(nl)

      ! local variables
      integer nvxb, imsg, ierr
!-----------------------------------------------------------------------
      nvxb = nl*ISIZE

      call mpi_ibcast (vx,nvxb,mpi_byte,isrc,AMR_NDSCOMM,imsg,ierr)
      amr_nds_ibcsi = imsg

      return
      end function
!=======================================================================
!> @brief Internode nonblocking broadcast for real array
!! @param[in]   vx       communication buffer
!! @param[in]   nl       array length
!! @param[in]   isrc     root rank
!! @return       amr_nds_ibcsr  communication request
      integer function amr_nds_ibcsr(vx,nl,isrc)
      implicit none

      include 'mpif.h'
      include 'SIZE'
      include 'PARALLEL'
      include 'AMR_SHM'

      ! argument list
      integer nl, isrc
      real vx(nl)

      ! local variables
      integer nvxb, imsg, ierr
!-----------------------------------------------------------------------
      nvxb = nl*WDSIZE

      call mpi_ibcast (vx,nvxb,mpi_byte,isrc,AMR_NDSCOMM,imsg,ierr)
      amr_nds_ibcsr = imsg

      return
      end function
!=======================================================================
!> @brief Internode MPI barrier
      subroutine amr_nds_gsync
      implicit none

      include 'SIZE'
      include 'AMR_SHM'

      ! local variables
      integer ierr
!-----------------------------------------------------------------------
      call mpi_barrier(AMR_NDSCOMM,ierr)

      return
      end subroutine
!=======================================================================
!> @brief Intranode nonblocking global vector commutative operation
!! @param[inout]   vs   send buffer
!! @param[inout]   vr   receive buffer
!! @param[in]      op   operator
!! @param[in]      nl   array length
!! @return       amr_shm_igopi  communication request
      integer function amr_shm_igopi(vs,vr,nl,op)
      implicit none

      include 'mpif.h'
      include 'SIZE'
      include 'PARALLEL'
      include 'AMR_SHM'

      ! argument list
      integer nl
      integer vs(nl), vr(nl)
      character*3 op

      ! local variables
      integer imsg, ierr
!-----------------------------------------------------------------------
      if     (op.eq.'+  ') then
        call mpi_iallreduce (vs,vr,nl,mpi_integer,mpi_sum ,
     $   AMR_SHMCOMM,imsg,ierr)
      elseif (op.EQ.'M  ') then
        call mpi_iallreduce (vs,vr,nl,mpi_integer,mpi_max ,
     $   AMR_SHMCOMM,imsg,ierr)
      elseif (op.EQ.'m  ') then
        call mpi_iallreduce (vs,vr,nl,mpi_integer,mpi_min ,
     $   AMR_SHMCOMM,imsg,ierr)
      elseif (op.EQ.'*  ') then
        call mpi_iallreduce (vs,vr,nl,mpi_integer,mpi_prod,
     $   AMR_SHMCOMM,imsg,ierr)
      else
        call amr_abort('Error: amr_nds_igopi; unknown op')
      endif

      amr_shm_igopi = imsg

      return
      end function
!=======================================================================
!> @brief Intranode nonblocking receive for integer array
!! @param[in]   msgtag   message tag
!! @param[in]   vx       communication buffer
!! @param[in]   nl       array length
!! @param[in]   isrc     source rank
!! @return       amr_shm_ircvi  communication request
      integer function amr_shm_ircvi(msgtag,vx,nl,isrc)
      implicit none

      include 'mpif.h'
      include 'SIZE'
      include 'PARALLEL'
      include 'AMR_SHM'

      ! argument list
      integer msgtag, nl,isrc
      integer vx(nl)

      ! local variables
      integer nvxb, imsg, ierr
!-----------------------------------------------------------------------
      nvxb = nl*ISIZE

      call mpi_irecv (vx,nvxb,mpi_byte,isrc,msgtag
     $       ,AMR_SHMCOMM,imsg,ierr)
      amr_shm_ircvi = imsg

      return
      end function
!=======================================================================
!> @brief Intranode nonblocking receive for real array
!! @param[in]   msgtag   message tag
!! @param[in]   vx       communication buffer
!! @param[in]   nl       array length
!! @param[in]   isrc     source rank
!! @return       amr_shm_ircvr  communication request
      integer function amr_shm_ircvr(msgtag,vx,nl,isrc)
      implicit none

      include 'mpif.h'
      include 'SIZE'
      include 'PARALLEL'
      include 'AMR_SHM'

      ! argument list
      integer msgtag, nl,isrc
      real vx(nl)

      ! local variables
      integer nvxb, imsg, ierr
!-----------------------------------------------------------------------
      nvxb = nl*WDSIZE

      call mpi_irecv (vx,nvxb,mpi_byte,isrc,msgtag
     $       ,AMR_SHMCOMM,imsg,ierr)
      amr_shm_ircvr = imsg

      return
      end function
!=======================================================================
!> @brief Intranode nonblocking send for integer array
!! @param[in]   msgtag   message tag
!! @param[in]   vx       communication buffer
!! @param[in]   nl       array length
!! @param[in]   idst     destination rank
!! @return       amr_shm_isndi  communication request
      integer function amr_shm_isndi(msgtag,vx,nl,idst)
      implicit none

      include 'mpif.h'
      include 'SIZE'
      include 'PARALLEL'
      include 'AMR_SHM'

      ! argument list
      integer msgtag, nl, idst
      integer vx(nl)

      ! local variables
      integer nvxb, imsg, ierr
!-----------------------------------------------------------------------
      nvxb = nl*ISIZE

      call mpi_isend (vx,nvxb,mpi_byte,idst,msgtag
     $       ,AMR_SHMCOMM,imsg,ierr)
      amr_shm_isndi = imsg

      return
      end function
!=======================================================================
!> @brief Intranode nonblocking send for real array
!! @param[in]   msgtag   message tag
!! @param[in]   vx       communication buffer
!! @param[in]   nl       array length
!! @param[in]   idst     destination rank
!! @return       amr_shm_isndr  communication request
      integer function amr_shm_isndr(msgtag,vx,nl,idst)
      implicit none

      include 'mpif.h'
      include 'SIZE'
      include 'PARALLEL'
      include 'AMR_SHM'

      ! argument list
      integer msgtag, nl, idst
      real vx(nl)

      ! local variables
      integer nvxb, imsg, ierr
!-----------------------------------------------------------------------
      nvxb = nl*WDSIZE

      call mpi_isend (vx,nvxb,mpi_byte,idst,msgtag
     $       ,AMR_SHMCOMM,imsg,ierr)
      amr_shm_isndr = imsg

      return
      end function
!=======================================================================
!> @brief Intranode nonblocking broadcast for integer array
!! @param[in]   vx       communication buffer
!! @param[in]   nl       array length
!! @param[in]   isrc     root rank
!! @return       amr_shm_ibcsi  communication request
      integer function amr_shm_ibcsi(vx,nl,isrc)
      implicit none

      include 'mpif.h'
      include 'SIZE'
      include 'PARALLEL'
      include 'AMR_SHM'

      ! argument list
      integer nl, isrc
      integer vx(nl)

      ! local variables
      integer nvxb, imsg, ierr
!-----------------------------------------------------------------------
      nvxb = nl*ISIZE

      call mpi_ibcast (vx,nvxb,mpi_byte,isrc,AMR_SHMCOMM,imsg,ierr)
      amr_shm_ibcsi = imsg

      return
      end function
!=======================================================================
!> @brief Intranode nonblocking broadcast for real array
!! @param[in]   vx       communication buffer
!! @param[in]   nl       array length
!! @param[in]   isrc     root rank
!! @return       amr_shm_ibcsr  communication request
      integer function amr_shm_ibcsr(vx,nl,isrc)
      implicit none

      include 'mpif.h'
      include 'SIZE'
      include 'PARALLEL'
      include 'AMR_SHM'

      ! argument list
      integer nl, isrc
      real vx(nl)

      ! local variables
      integer nvxb, imsg, ierr
!-----------------------------------------------------------------------
      nvxb = nl*WDSIZE

      call mpi_ibcast (vx,nvxb,mpi_byte,isrc,AMR_SHMCOMM,imsg,ierr)
      amr_shm_ibcsr = imsg

      return
      end function
!=======================================================================
!> @brief Intranode nonblocking gather for integer array
!! @param[in]   vs       send buffer
!! @param[in]   ns       send buffer length
!! @param[in]   vr       receive buffer
!! @param[in]   nr       receive buffer length
!! @param[in]   isrc     root rank
!! @return       amr_shm_igti  communication request
      integer function amr_shm_igti(vs,ns,vr,nr,isrc)
      implicit none

      include 'mpif.h'
      include 'SIZE'
      include 'PARALLEL'
      include 'AMR_SHM'

      ! argument list
      integer ns, nr, isrc
      integer vs(ns), vr(nr)

      ! local variables
      integer imsg, ierr
!-----------------------------------------------------------------------
      ! possible check
      !if (nr.lt.ns*AMR_SHMNP)

      call mpi_igather (vs,ns,mpi_integer,vr,ns,mpi_integer,isrc,
     $     AMR_SHMCOMM,imsg,ierr)
      amr_shm_igti = imsg

      return
      end function
!=======================================================================
!> @brief Intranode nonblocking allgather for integer array
!! @param[in]   vs       send buffer
!! @param[in]   ns       send buffer length
!! @param[in]   vr       receive buffer
!! @param[in]   nr       receive buffer length
!! @return      amr_shm_iagti  communication request
      integer function amr_shm_iagti(vs,ns,vr,nr)
      implicit none

      include 'mpif.h'
      include 'SIZE'
      include 'PARALLEL'
      include 'AMR_SHM'

      ! argument list
      integer ns, nr
      integer vs(ns), vr(nr)

      ! local variables
      integer imsg, ierr
!-----------------------------------------------------------------------
      ! possible check
      !if (nr.lt.ns*NP)

      call mpi_iallgather (vs,ns,mpi_integer,vr,ns,mpi_integer,
     $     AMR_SHMCOMM,imsg,ierr)
      amr_shm_iagti = imsg

      return
      end function
!=======================================================================
!> @brief Intranode nonblocking gatherv for integer array.
!! @param[in]   vs       send buffer
!! @param[in]   ns       send buffer length
!! @param[in]   vr       receive buffer
!! @param[in]   nr       receive buffer length
!! @param[in]   cr       receive count array
!! @param[in]   dr       receive displacement array
!! @param[in]   isrc     root rank
!! @return       amr_shm_igtvi  communication request
      integer function amr_shm_igtvi(vs,ns,vr,nr,cr,dr,isrc)
      implicit none

      include 'mpif.h'
      include 'SIZE'
      include 'PARALLEL'
      include 'AMR_SHM'

      ! argument list
      integer ns, nr, isrc
      integer vs(ns), vr(nr), cr(AMR_SHM_MAX), dr(AMR_SHM_MAX)

      ! local variables
      integer imsg, ierr
!-----------------------------------------------------------------------
      ! possible check
      !if (nr.lt.ns*AMR_SHMNP)

      call mpi_igatherv (vs,ns,mpi_integer,vr,cr,dr,mpi_integer,isrc,
     $     AMR_SHMCOMM,imsg,ierr)
      amr_shm_igtvi = imsg

      return
      end function
!=======================================================================
!> @brief Intranode nonblocking allgatherv for integer array.
!! @param[in]   vs       send buffer
!! @param[in]   ns       send buffer length
!! @param[in]   vr       receive buffer
!! @param[in]   nr       receive buffer length
!! @param[in]   cr       receive count array
!! @param[in]   dr       receive displacement array
!! @return       amr_shm_iagtvi  communication request
      integer function amr_shm_iagtvi(vs,ns,vr,nr,cr,dr)
      implicit none

      include 'mpif.h'
      include 'SIZE'
      include 'PARALLEL'
      include 'AMR_SHM'

      ! argument list
      integer ns, nr
      integer vs(ns), vr(nr), cr(AMR_SHM_MAX), dr(AMR_SHM_MAX)

      ! local variables
      integer imsg, ierr
!-----------------------------------------------------------------------
      ! possible check
      !if (nr.lt.ns*AMR_SHMNP)

      call mpi_iallgatherv (vs,ns,mpi_integer,vr,cr,dr,mpi_integer,
     $     AMR_SHMCOMM,imsg,ierr)
      amr_shm_iagtvi = imsg

      return
      end function
!=======================================================================
!> @brief Intranode MPI barrier
      subroutine amr_shm_gsync
      implicit none

      include 'SIZE'
      include 'AMR_SHM'

      ! local variables
      integer ierr
!-----------------------------------------------------------------------
      call mpi_barrier(AMR_SHMCOMM,ierr)

      return
      end subroutine
!=======================================================================
!> @brief Wait for finalisation of single nonblocking operation
!! @param[in]    irqs   message request
      subroutine amr_wait(irqs)
      implicit none

      include 'mpif.h'

      ! argument list
      integer irqs

      ! local variables
      integer ierr
!-----------------------------------------------------------------------
      call mpi_wait (irqs,MPI_STATUS_IGNORE,ierr)

      return
      end subroutine
!=======================================================================
!> @brief Wait for finalisation of multiple nonblocking operation
!! @param[in]    irqs   message request
!! @param[in]    nl     array length
      subroutine amr_waitall(irqs,nl)
      implicit none

      include 'mpif.h'

      ! argument list
      integer nl
      integer irqs(nl)

      ! local variables
      integer ierr
!-----------------------------------------------------------------------
      call mpi_waitall (nl,irqs,MPI_STATUSES_IGNORE,ierr)

      return
      end subroutine
!=======================================================================
!> @brief Split mpi communicator taking into account shared memory nodes
      subroutine amr_node_split()
      implicit none

      include 'mpif.h'
      include 'SIZE'
      include 'AMR'
      include 'AMR_SHM'

      include 'PARALLEL'
      include 'SOLN'

      ! local variables
      integer ierr      ! error flag
      integer imsg(2)   ! message request
      integer itmp, itmpv(2)
      integer igri(AMR_NDS_MAX)

      ! functions
      integer amr_glb_ibcsi, amr_glb_iscani,
     $        amr_glb_igopi
      real dnekclock
!-----------------------------------------------------------------------
      ! duplicate communicator
      call create_comm(AMR_GLBCOMM)

!!#define SNGNDE
#undef SNGNDE
#ifdef SNGNDE
      ! for testing on single node only
      itmp = nid/32
      call mpi_comm_split(AMR_GLBCOMM, itmp, 0, AMR_SHMCOMM, ierr)
#else
      ! split communicator between nodes
      call mpi_comm_split_type(AMR_GLBCOMM, MPI_COMM_TYPE_SHARED,
     $         0, MPI_INFO_NULL, AMR_SHMCOMM, ierr)
#endif
#undef SNGNDE

      ! get local id and process number for intranode communicator
      call mpi_comm_size (AMR_SHMCOMM, AMR_SHMNP, ierr)
      call mpi_comm_rank (AMR_SHMCOMM, AMR_SHMID, ierr)

      ! check array size
      ierr = 0
      if (AMR_SHMNP.gt.AMR_SHM_MAX) then
         ierr = 1
      endif
      call amr_chk_abort(ierr,
     $     'ERROR; amr_shm_split: too small AMR_SHM_MAX')

      ! mark node masters
      if (AMR_SHMID.eq.AMR_SHMST) then
         AMR_IFNODEM = .true.
      else
         AMR_IFNODEM = .false.
      endif

      ! get communicator for node masters
      ! get number of node masters
      if (AMR_IFNODEM) then
         itmpv(1) = 1
      else
         itmpv(1) = 0
      endif
      imsg(1) = amr_glb_iscani(itmpv(1),itmpv(2),1,'+  ')
      call amr_wait(imsg(1))
      itmpv(1) = itmpv(2)
      call mpi_comm_size (AMR_GLBCOMM, itmp, ierr)
      itmp = itmp -1
      imsg(1) = amr_glb_ibcsi(itmpv(1),1,itmp)
      call amr_wait(imsg(1))

      ! check array size
      ierr = 0
      if (itmpv(1).gt.AMR_NDS_MAX) then
         ierr = 1
      endif
      call amr_chk_abort(ierr,
     $     'ERROR; amr_shm_split: too small AMR_NDS_MAX')

      ! collect set of node masters
      call izero(igri,AMR_NDS_MAX)
      if (AMR_IFNODEM) igri(itmpv(2)) = nid

      imsg(1) = amr_glb_igopi(igri,AMR_NODE_GOFF,itmpv(1),'+  ')
      ! create group and communicators
      call mpi_comm_group(AMR_GLBCOMM,AMR_GLBGROUP,ierr)

      call amr_wait(imsg(1))

      call mpi_group_incl(AMR_GLBGROUP,itmpv(1),AMR_NODE_GOFF,
     $     AMR_NDSGROUP,ierr)
      call mpi_comm_create(AMR_GLBCOMM,AMR_NDSGROUP,AMR_NDSCOMM,ierr)

      ierr = 0
      ! get size and rank id of internode communicator
      if (AMR_IFNODEM) then
         call mpi_comm_size (AMR_NDSCOMM, AMR_NDSNP, ierr)
         call mpi_comm_rank (AMR_NDSCOMM, AMR_NDSID, ierr)

         ! check consistency
         if (AMR_NDSNP.ne.itmpv(1)) then
            ierr = 1
            call amr_log(AMR_LP_ERR,
     $       'ERROR; amr_shm_split: inconsistent node number')
         endif

         if ((AMR_NDSID+1).ne.itmpv(2)) then
            ierr = 1
            call amr_log(AMR_LP_ERR,
     $       'ERROR; amr_shm_split: inconsistent node id')
         endif
      endif

      call amr_chk_abort(ierr,
     $     'ERROR; amr_shm_split: error in node split')

      ! get node number
      AMR_NNODE = itmpv(1)

      ! get array position
      AMR_NODE_NR = itmpv(2)

      ! fill last entrance in offset array
      AMR_NODE_GOFF(AMR_NNODE+1) = NP

      ! check if numberring of all shm mpi ranks is continuous
      ierr = 0
      if (NID.lt.AMR_NODE_GOFF(AMR_NODE_NR).or.
     $    NID.gt.AMR_NODE_GOFF(AMR_NODE_NR+1)) ierr = 1
      call amr_chk_abort(ierr,
     $     'ERROR; amr_shm_split: shm numberring not continuous.')

      ! chek node offset and intranode communicator size
      ierr = 0
      if ((AMR_NODE_GOFF(AMR_NODE_NR+1)-AMR_NODE_GOFF(AMR_NODE_NR))
     $    .ne.AMR_SHMNP) ierr = 1
      call amr_chk_abort(ierr,
     $       'ERROR; amr_shm_split: inconsistent node offset')

      if (AMR_NNODE.gt.1) then
         AMR_IF2LEV = .true.
      else
         AMR_IF2LEV = .false.
      endif

      ! for testing only
      !AMR_IF2LEV = .false.

      return
      end subroutine
!=======================================================================
!> @brief Free new mpi communicators
      subroutine amr_node_free()
      implicit none

      include 'mpif.h'
      include 'SIZE'
      include 'AMR'
      include 'AMR_SHM'

      ! local variables
      integer ierr      ! error flag
!-----------------------------------------------------------------------
      ! free groups
      call mpi_group_free(AMR_NDSGROUP,ierr)
      call mpi_group_free(AMR_GLBGROUP,ierr)
      ! free communicators
      call mpi_comm_free(AMR_SHMCOMM,ierr)
!      call mpi_comm_free(AMR_NDSCOMM,ierr)
!      call mpi_comm_free(AMR_GLBCOMM,ierr)
      
      return
      end subroutine
!=======================================================================

