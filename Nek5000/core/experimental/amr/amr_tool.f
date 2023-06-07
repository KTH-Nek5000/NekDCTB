!> @file amr_tool.f
!! @ingroup nekamr
!! @brief Set of tools for nekamr
!! @author Adam Peplinski
!! @date Feb 26, 2016
! This file requires number of constant defined in include file
#include "nekamr.h"
!=======================================================================
!> @brief Find threshold for given element number and error
!! @param[in] el_ref          num. of elements to refine/coarsen
!! @param[in] el_dlt          accuracy of element count
!! @param[in] erri            error value per element
!! @param[in] errl            logical flag to exclude element element
!! @param[in] lnel            local element number (size of erri and errl)
!! @param[in] work1,work2     work arrays
!! @param[in] nint            interval number (size of work arrays)
!! @param[in] it_max          max number of iterations
!! @param[in] ifcrs           count for coarsening/refinement threshold
!! @return    amr_thrsh_get
      real function amr_thrsh_get(el_ref,el_dlt,erri,errl,lnel,
     $     work1,work2,nint,it_max,ifcrs)
      implicit none
      include 'SIZE'
      include 'AMR'

      ! argument list
      integer el_ref, el_dlt, lnel, nint, it_max
      integer work1(nint),work2(nint)
      real erri(lnel)
      logical errl(lnel)
      logical ifcrs

      ! local variables
      integer il, jl
      integer ierr
      integer iint, icall, istart, iend, istp
      real lmax, lmin, ldlt
      character(len=AMR_LSTL_LOG) str

      ! functions
      real glmax, glmin
!-----------------------------------------------------------------------
      ! sanity check
      ierr = 0
      if(lnel.lt.1.or.nint.lt.1) ierr = 1
      call amr_chk_abort(ierr,'Error: ref_thrsh_get; wrong array size')

      ! get global error max/min and interval
      lmax = glmax(erri,lnel)
      lmin = glmin(erri,lnel)
      ldlt = (lmax-lmin)/real(nint)
      ! sanity check
      ierr = 0
      if(ldlt.eq.0.0) ierr = 1
      call amr_chk_abort(ierr,'Error: ref_thrsh_get; lmax = lmin')

      if (el_ref.le.0) then
         if (ifcrs) then
            amr_thrsh_get = lmin*0.9
         else
            amr_thrsh_get = lmax*1.1
         endif
      else
         ! get count per interval
         call izero(work1,nint)
         do il=1,lnel
            if (errl(il)) then
               iint = int((erri(il)-lmin)/ldlt)+1
               iint = min(iint,nint)
               work1(iint) = work1(iint) + 1
            endif
         enddo
         call igop(work1,work2,'+  ',nint)

         if (ifcrs) then
            istart = 1
            iend = nint
            istp = 1
         else
            istart = nint
            iend = 1
            istp = -1
         endif

         ! big loop
         icall = 0
         do jl=1, it_max
            ! get approximate interval
            do il=istart,iend,istp
               icall = icall + work1(il)
               if (icall.ge.el_ref.or.il.eq.iend) goto 30
            enddo
 30         continue

            ! did we converged?
            ! exact solution or whole domain within a range
            if (icall.le.el_ref) then
               amr_thrsh_get = lmin + ldlt*(il+min(0,istp))
               goto 90
            else
               ! is the error smaller than assumed delta
               if (work1(il).le.el_dlt) then
                  il = il-istp
                  amr_thrsh_get = lmin + ldlt*(il+min(0,istp))
                  goto 90
               else
                  ! correct current element number
                  icall = icall - work1(il)
                  ! update min/max values
                  lmax = lmin + ldlt*il
                  lmin = lmax - ldlt
                  ldlt = (lmax-lmin)/real(nint)

                  ! update count per interval
                  call izero(work1,nint)
                  do il=1,lnel
                     if (errl(il)) then
                        iint = int((erri(il)-lmin)/ldlt)+1
                        if (iint.gt.0.and.iint.le.nint) then
                           work1(iint) = work1(iint) + 1
                        endif
                     endif
                  enddo
                  call igop(work1,work2,'+  ',nint)
               endif
            endif
         enddo

 90      continue
      endif

      amr_thrsh_get = max(lmin,min(lmax,amr_thrsh_get))

      write(str,120) el_ref,el_dlt,ifcrs
      call amr_log(AMR_LP_PRD,trim(str))
      write(str,121) jl,icall,amr_thrsh_get
      call amr_log(AMR_LP_PRD,trim(str))
 120  format('Threshold search: el_ref=',I6,',el_dlt=',I5,',ifcrs=',L1)
 121  format('   iter=',I4,',elem. count=',I6,',threshold=',1p3E13.4)

      return
      end function
!=======================================================================
!> @brief Write log messages
!! @param[in] priority  log priority
!! @param[in] logs      log body
      subroutine amr_log(priority,logs)
      implicit none

      include 'SIZE'
      include 'AMR'

      ! argument list
      integer priority
      character(len=*) logs

      character*(*) nek_name
      parameter(nek_name='nekamr'//CHAR(0))

      ! local variables
      character(len=AMR_LSTL_LOG) llogs
      integer slen
      logical ifcalled
      save ifcalled
      data ifcalled /.FALSE./
!-----------------------------------------------------------------------
      if (.not.ifcalled) then
        ifcalled=.TRUE.

        ! register nekamr package
#ifdef P4EST
        call fsc_pkg_reg(AMR_PKG_ID,priority,nek_name)
#endif
      endif

      ! check log length
      slen = len_trim(logs)
      if (slen.ge.AMR_LSTL_LOG) then
         llogs = 'Warning; too long log string, shortening'//CHAR(0)
#ifdef P4EST
         call fsc_log (AMR_PKG_ID, 1, priority,trim(llogs))
#endif
         llogs = logs(1:AMR_LSTL_LOG-1)//CHAR(0)
      else
         llogs = trim(logs)//CHAR(0)
      endif

#ifdef P4EST
      call fsc_log (AMR_PKG_ID, 1, priority, trim(llogs))
#endif
      return
      end
!=======================================================================
!> @brief SC based abort function
!! @param[in] logs      log body
      subroutine amr_abort(logs)
      implicit none

      include 'SIZE'
      include 'AMR'

      ! argument list
      character(len=*) logs

      ! local variables
      character(len=AMR_LSTL_LOG) llogs
      integer slen
!-----------------------------------------------------------------------
      ! check log length
      slen = len_trim(logs)
      if (slen.ge.AMR_LSTL_LOG) then
         llogs = 'Warning; too long log string, shortening'//CHAR(0)
#ifdef P4EST
         call fsc_log (AMR_PKG_ID, 1, AMR_LP_ERR,trim(llogs))
#endif
         llogs = logs(1:AMR_LSTL_LOG-1)//CHAR(0)
      else
         llogs = trim(logs)//CHAR(0)
      endif
#ifdef P4EST
      call fsc_abort(trim(llogs))
#endif
      return
      end
!=======================================================================
!> @brief SC based check abort function
!! @param[in] ierr      error indicator
!! @param[in] logs      log body
      subroutine amr_chk_abort(ierr,logs)
      implicit none

      include 'SIZE'
      include 'AMR'

      ! argument list
      integer ierr
      character(len=*) logs

      ! local variables
      character(len=AMR_LSTL_LOG) llogs
      integer slen, ierrl

      ! functions
      integer iglsum
!-----------------------------------------------------------------------
      ! check log length
      slen = len_trim(logs)
      if (slen.ge.AMR_LSTL_LOG) then
         llogs = 'Warning; too long log string, shortening'//CHAR(0)
#ifdef P4EST
         call fsc_log (AMR_PKG_ID, 1, AMR_LP_ERR,trim(llogs))
#endif
         llogs = logs(1:AMR_LSTL_LOG-1)//CHAR(0)
      else
         llogs = trim(logs)//CHAR(0)
      endif

      ierrl = iglsum(ierr,1)
      if (ierrl.gt.0) then
#ifdef P4EST
         call fsc_abort(trim(llogs))
#endif
      else
         return
      endif

      return
      end
!=======================================================================
!> @brief Perform refinement/coarsening consistency check
!! @details This routine should remove all coarsening marks that would
!!  lead to unbalanced mesh. Right now it is quite conservative reseting
!!  all coarsening marks at the level jumps (any neighbour at higher
!!  refinement level) irrespectively of coarsening mark on neighbour
!!  elements.
      subroutine amr_mark_check()
      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'TOPOL'
      include 'PARALLEL'   ! debugging only
      include 'AMR'
      include 'AMR_TOPOL'

      integer mid,mp,nekcomm,nekgroup,nekreal
      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal

      ! local variables
      ! gather-scatter handle
      ! vertices
      integer gs_handle_vrt
      ! faces
      integer gs_handle_fcs
      ! edges
      integer gs_handle_edg

      ! array length
      integer lnoden

      ! loop index
      integer iel, ivr, idim

      ! temporary variables
      integer itmp
      real rtmp

      ! inverse of number of children faces
      real rincf
      parameter (rincf=1.0/(2.0**(LDIM-1)))

      ! communication arrays
      integer ivar_vrt(AMR_NVRT,LELT), ivar_fcs(AMR_NFCS,LELT,2)
      real rvar_fcs(AMR_NFCS,LELT)

      ! faces attached to given vertex and its oposite faces
      integer ivface(3,8), ivfaceo(3,8)
      save ivface, ivfaceo
      data ivface  /1,3,5,  2,3,5,  1,4,5,  2,4,5,
     $              1,3,6,  2,3,6,  1,4,6,  2,4,6/
      data ivfaceo /2,4,6,  1,4,6,  2,3,6,  1,3,6,
     $              2,4,5,  1,4,5,  2,3,5,  1,3,5/

!#define DEBUG
#ifdef DEBUG
      character*3 str1, str2
      integer iunit, ierr
      ! call number
      integer icalld
      save icalld
      data icalld /0/
#endif
!-----------------------------------------------------------------------
      ! check if we refine/derefine elements on AMR_LMAX/0
      do iel=1,NELT
         if ((AMR_MARK(iel).eq.AMR_RM_H_REF.and.
     $        AMR_LEVEL(iel).eq.AMR_LMAX).or.
     $        (AMR_MARK(iel).eq.AMR_RM_H_CRS.and.
     $        AMR_LEVEL(iel).eq.0)) then
            AMR_MARK(iel) = AMR_RM_NONE
         endif
      enddo

      ! set gather-scatter communicators
      lnoden = AMR_NVRT*NELT
      call fgslib_gs_setup(gs_handle_vrt,AMR_GLBNR_VRT,lnoden,
     $     nekcomm,mp)
      lnoden = AMR_NFCS*NELT
      call fgslib_gs_setup(gs_handle_fcs,AMR_GLBNR_FCS,lnoden,
     $     nekcomm,mp)
!      if (IF3D) then
!         lnoden = AMR_NEDG*NELT
!         call fgslib_gs_setup(gs_handle_edg,AMR_GLBNR_EDG,lnoden,
!     $        nekcomm,mp)
!      endif

#if 0
!     There is something wrong with this check, and this is just a hack to fill in holes
!     in levels
      ! distribute element refinement level
      do iel=1,NELT
         do ivr = 1, AMR_NFCS
            ivar_fcs(ivr,iel,1) = AMR_LEVEL(iel)
         enddo
      enddo

      ! find max element levels across faces
      call fgslib_gs_op(gs_handle_fcs,ivar_fcs,2,4,0)

      ! chck neighbours along z-axis
      do iel=1,NELT
         if (ivar_fcs(5,iel,1).gt.AMR_LEVEL(iel).or.
     $       ivar_fcs(6,iel,1).gt.AMR_LEVEL(iel))
     $        AMR_MARK(iel) = AMR_RM_H_REF
      enddo
#endif

      ! distribute element refinement level
      do iel=1,NELT
         do ivr = 1, AMR_NVRT
            ivar_vrt(ivr,iel) = AMR_LEVEL(iel)
         enddo
      enddo

      ! find max element levels across vertices (sufficient for mesh connectivity)
      call fgslib_gs_op(gs_handle_vrt,ivar_vrt,2,4,0)

      ! do not coarsen if neighbour level is higher
      do iel=1,NELT
         itmp = AMR_LEVEL(iel)
         do ivr = 1, AMR_NVRT
            itmp = max(itmp,ivar_vrt(ivr,iel))
         enddo
         itmp = itmp - AMR_LEVEL(iel)

         if (itmp.gt.1) then
            ! consistency check; mesh should be balanced
            call amr_abort
     $             ('ERROR: amr_mark_check; unbalanced mesh')
         elseif (itmp.eq.1) then
            ! jump in refinemnet; do not coarsen
            if (AMR_MARK(iel).lt.AMR_RM_NONE) then
                AMR_MARK(iel)=AMR_RM_NONE
            endif
         endif
      enddo

      ! distribute h-type refinemen mark
      lnoden = AMR_NVRT*NELT
      call izero(ivar_vrt,lnoden)
      do iel = 1, NELT
         if (AMR_MARK(iel).eq.AMR_RM_H_REF) then
            do ivr = 1, AMR_NVRT
               ivar_vrt(ivr,iel) = AMR_LEVEL(iel)
            enddo
         endif
      enddo

      ! find max level mark of refining elements across vertices
      ! (sufficient for mesh connectivity)
      call fgslib_gs_op(gs_handle_vrt,ivar_vrt,2,4,0)

      ! do not coarsen if neighbour is refining
      do iel=1,NELT
         itmp = 0
         do ivr = 1, AMR_NVRT
            itmp = max(itmp,ivar_vrt(ivr,iel))
         enddo
         if (itmp.eq.1) then
            ! jump in refinemnet do not coarsen
            if (AMR_MARK(iel).lt.AMR_RM_NONE) then
                AMR_MARK(iel)=AMR_RM_NONE
            endif
         endif
      enddo

      ! Correct coarsening falg in the family. The family should coarsen
      ! only if all family members are marked for coarsening.
      lnoden = AMR_NVRT*NELT
      call izero(ivar_vrt,lnoden)
      do iel=1,NELT
         if (AMR_FML_MARK(1,iel).gt.0) then
            ivar_vrt(AMR_FML_MARK(2,iel),iel) = AMR_MARK(iel)
         endif
      enddo

      ! find max refinement mark across family vertices
      call fgslib_gs_op(gs_handle_vrt,ivar_vrt,2,4,0)

      do iel=1,NELT
         if (AMR_FML_MARK(1,iel).gt.0.and.
     $       AMR_MARK(iel).eq.AMR_RM_H_CRS.and.
     $       ivar_vrt(AMR_FML_MARK(2,iel),iel).ne.AMR_RM_H_CRS)
     $       AMR_MARK(iel) = AMR_RM_NONE
      enddo

      ! exclude single element surrounded by refined elements (this seems
      ! to cause the rise in pressure iterations)
      ! Notice, due to coarsening I need more communication as I have to
      ! check external faces of the whole family not just single element.


      ! repeat whole operations a few times
      do idim=1,ndim
         ! get predicted neighbours' levels after refinement.
         ! This operation is far from exact as mesh balancing is not taken
         ! into account.
         do iel = 1, NELT
            rtmp = real(AMR_LEVEL(iel))
            if (AMR_MARK(iel).eq.AMR_RM_H_REF) then
               rtmp = rtmp + 1.0
            elseif (AMR_MARK(iel).eq.AMR_RM_H_CRS) then
               rtmp = rtmp - 1.0
            endif
            do ivr = 1, AMR_NFCS
               rvar_fcs(ivr,iel) = rtmp
               ! take into account nonconforming faces; notice change of face numbering
               if (AMR_HNG_FCS(EFACE(ivr),iel).gt.-1)
     $             rvar_fcs(ivr,iel) = rvar_fcs(ivr,iel)*rincf
            enddo
         enddo

         ! sum level mark  across faces
         call fgslib_gs_op(gs_handle_fcs,rvar_fcs,1,1,0)

         ! extract neighbours levels
         do iel=1,NELT
            do ivr = 1, AMR_NFCS
               ivar_fcs(ivr,iel,1) = nint(rvar_fcs(ivr,iel)) -
     $                 AMR_LEVEL(iel)
               if (AMR_MARK(iel).eq.AMR_RM_H_REF) then
                  ivar_fcs(ivr,iel,1)=ivar_fcs(ivr,iel,1) - 1
               elseif (AMR_MARK(iel).eq.AMR_RM_H_CRS) then
                  ivar_fcs(ivr,iel,1)=ivar_fcs(ivr,iel,1) + 1
               endif
            enddo
         enddo

         ! transfer mark within families
         lnoden = AMR_NFCS*NELT
         call izero(ivar_fcs(1,1,2),lnoden)
         do iel=1,NELT
            if (AMR_FML_MARK(1,iel).gt.0) then
               do ivr=1,NDIM
                  ivar_fcs(ivface(ivr,AMR_FML_MARK(2,iel)),iel,2) =
     $            ivar_fcs(ivfaceo(ivr,AMR_FML_MARK(2,iel)),iel,1)
               enddo
            endif
         enddo

         ! sum level mark of refining elements across family faces
         call fgslib_gs_op(gs_handle_fcs,ivar_fcs(1,1,2),2,1,0)

         ! put exchanged neighbours levels back
         do iel=1,NELT
            if (AMR_FML_MARK(1,iel).gt.0) then
               do ivr=1,NDIM
                  ivar_fcs(ivface(ivr,AMR_FML_MARK(2,iel)),iel,1) =
     $            ivar_fcs(ivface(ivr,AMR_FML_MARK(2,iel)),iel,2) -
     $            ivar_fcs(ivfaceo(ivr,AMR_FML_MARK(2,iel)),iel,1)
               enddo
            endif
         enddo

         ! perform checks of oposite face neighbours
         do iel=1,NELT
            do ivr = 1, AMR_NFCS, 2
               itmp = AMR_LEVEL(iel)
               if (AMR_MARK(iel).eq.AMR_RM_H_REF) then
                  itmp = itmp + 1
               elseif (AMR_MARK(iel).eq.AMR_RM_H_CRS) then
                  itmp = itmp - 1
               endif
               if ((ivar_fcs(ivr,iel,1).gt.itmp).and.
     $               (ivar_fcs(ivr+1,iel,1).gt.itmp)) then
                  if (AMR_MARK(iel).eq.AMR_RM_NONE) then
                        AMR_MARK(iel)=AMR_RM_H_REF
                  elseif (AMR_MARK(iel).eq.AMR_RM_H_CRS) then
                     AMR_MARK(iel)=AMR_RM_NONE
                  endif
               endif
            enddo
         enddo
      enddo

      ! free communicators
      call fgslib_gs_free (gs_handle_vrt)
      call fgslib_gs_free (gs_handle_fcs)
!      if (IF3D) call fgslib_gs_free (gs_handle_edg)

#ifdef DEBUG
      ! for testing
      ! to output refinement
      icalld = icalld+1
      call io_file_freeid(iunit, ierr)
      write(str1,'(i3.3)') NID
      write(str2,'(i3.3)') icalld
      open(unit=iunit,file='mchk.txt'//str1//'i'//str2)

      do iel=1,NELT
         write(iunit,*) iel, LGLEL(iel), AMR_MARK(iel)
         write(iunit,*) (ivar_fcs(ivr+1,iel,1),ivr=1,AMR_NFCS)
         write(iunit,*) (ivar_fcs(ivr+1,iel,2),ivr=1,AMR_NFCS)
      enddo

      close(iunit)
#endif

      return
      end subroutine
#undef DEBUG
!=======================================================================
!> @brief Fix internal element faces and edges by averaging coordinates
!! @todo Treat properly periodic bc
      subroutine amr_fix_geom()
      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'TSTEP'
      include 'GEOM'
      include 'PARALLEL'
      include 'WZ'
      include 'AMR'
      include 'AMR_TOPOL'

      include 'SOLN'

      ! tmporary storage
      integer lt
      parameter (lt = lx1*ly1*lz1)
      real xb(lt,lelt),yb(lt,lelt),zb(lt,lelt),w1(lt),w2(lt)
      real tmsk(lt,lelt)!, tmskj(lt,lelt)
      common /scrns/ xb, yb, zb, w1, w2, tmsk!, tmaskj

      ! local variables
      integer ntot, kpass, iel, ifc, il, jl, nel
      real xm,ym,zm,xx,yx,zx
      real utest(lx1,lz1)
      real epsl
      parameter (epsl=0.0)
      character(len=AMR_LSTL_LOG) str
      character*3 cb

      ! just temporary hack to return for periodic bc
      integer iperiodicl

      ! functions
      integer iglmax
      real glamax
!-----------------------------------------------------------------------
      ! array sizes
      ntot = lt*nelt
      ifield = 1                ! velocity field
      nel = nelv
      if (nelgv.ne.nelgt .or. .not.ifflow) then    ! temperature field
         ifield = 2
         nel = nelt
      endif

      ! mark periodic bc
      ! just temporary hack to return for periodic bc
      iperiodicl=0
      call rone(tmsk,ntot)
      do iel=1,nel                 ! fill mask where bc is periodic
        do ifc=1,AMR_NFCS          ! so we don't translate periodic bcs (z only)
          cb =cbc(ifc,iel,ifield)
          if (cb.eq.'P  ') then
            call facev (tmsk,iel,ifc,0.0,lx1,ly1,lz1)
            iperiodicl=1
          endif
        enddo
      enddo

      iperiodicl = iglmax(iperiodicl,1)
      if (iperiodicl.gt.0) return

#if 0
      ! take care of not face related periodic edges
      call dsop(tmsk,'*  ',lx1,ly1,lz1)

      ! take care of nonconforming faces touching external boundarries
      call copy(tmskj,tmsk,ntot)
      call amr_apply_j(tmskj,lx1,ly1,lz1,nel)

      do iel=1,nel
        ! hanging element
        if (AMR_HNG_ELM(iel).eq.1) then
          do ifc=1,AMR_NFCS
            ! hanging face
            if (AMR_HNG_FCS(ifc,iel).ne.-1) then
              call ftovecl(utest,tmskj(1,iel),ifc,lx1,ly1,lz1)
              do jl=1,lz1
                do il=1,lx1
                   if (utest(il,jl).gt.epsl) then
                      utest(il,jl) = 1.0
                   else
                      utest(il,jl) = 0.0
                   endif
                enddo
              enddo
              call vectofl(tmskj(1,iel),utest,ifc,lx1,ly1,lz1)
            endif
          enddo
        endif
      enddo
      ! remember to comment common block to use outpost
      !call outpost(tmsk,tmskj,vz,pr,t,'   ')
#endif

      ! loop over dimensions
      do kpass = 1,ldim+1
         ! copy coordinates
         call copy(xb,xm1,ntot)
         call copy(yb,ym1,ntot)
         if (IF3D) call copy(zb,zm1,ntot)

         !call outpost(xb,yb,zb,pr,t,'   ')

         ! average faces
         call h1_proj(xb,lx1,ly1,lz1)
         call h1_proj(yb,lx1,ly1,lz1)
         if (IF3D) call h1_proj(zb,lx1,ly1,lz1)

         ! average faces
         !call amr_h1_proj_msk(xb,tmsk,lx1,ly1,lz1,nel)
         !call amr_h1_proj_msk(yb,tmsk,lx1,ly1,lz1,nel)
         !if (IF3D) call amr_h1_proj_msk(zb,tmsk,lx1,ly1,lz1,nel)

         !call outpost(xb,yb,zb,pr,t,'   ')

         ! get shift
         call sub2(xb,xm1,ntot)
         call sub2(yb,ym1,ntot)
         if (IF3D) call sub2(zb,zm1,ntot)

         !call outpost(xb,yb,zb,pr,t,'   ')

         ! remove periodic bc
         !call col2(xb,tmskj,ntot)
         !call col2(yb,tmskj,ntot)
         !if (IF3D) call col2(zb,tmskj,ntot)

         !call outpost(xb,yb,zb,pr,t,'   ')
         !call exitt0

         xm = glamax(xb,ntot)
         ym = glamax(yb,ntot)
         if (IF3D) zm = glamax(zb,ntot)

         ! perform Gordon-Hall
         do iel = 1,nelt
            if (kpass.le.ldim) then
               call gh_face_extend(xb(1,iel),zgm1,lx1,kpass,w1,w2)
               call gh_face_extend(yb(1,iel),zgm1,lx1,kpass,w1,w2)
               if (IF3D) call gh_face_extend(zb(1,iel),zgm1,lx1,
     $          kpass,w1,w2)
            endif
         enddo

         ! apply corrections
         if (kpass.le.ldim) then
            call add2(xm1,xb,ntot)
            call add2(ym1,yb,ntot)
            if (IF3D) call add2(zm1,zb,ntot)
         endif

         xx = glamax(xb,ntot)
         yx = glamax(yb,ntot)
         if (IF3D) then
           zx = glamax(zb,ntot)
         else
           zx = 0.0
         endif

         write(str,1) xm,ym,zm,xx,yx,zx,kpass
    1    format(1p6e12.4,' xyz repair',i2)
         call amr_log(AMR_LP_PRD,trim(str))
      enddo

      ! mark deformed mesh
      param(59) = 1.       ! ifdef = .true.
      call geom_reset(1)   ! reset metrics, etc.

      return
      end subroutine
!=======================================================================
!> @brief Save missing pressure history
      subroutine amr_lagpr
      implicit none

      include 'SIZE'
      include 'SOLN'
      include 'TSTEP'
      include 'AMR'

      ! local variables
      integer ntot
!-----------------------------------------------------------------------
      ! array sizes
      ntot=nx2*ny2*nz2*nelv
      if (NBD.EQ.3) then
         call copy(AMR_PRLAG(1,1,1,1,2),AMR_PRLAG(1,1,1,1,1),ntot)
      endif
      call copy(AMR_PRLAG(1,1,1,1,1),PR,ntot)

      return
      end
!=======================================================================
!> @brief Divide right hand sides by BM1 before refining/coarsening
      subroutine amr_rhs_rmv_bm
      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'SOLN'
      include 'MASS'

      ! local variables
      integer ntotv,ntott,il
!-----------------------------------------------------------------------
      ntotv = NX1*NY1*NZ1*NELV
      ntott = NX1*NY1*NZ1*NELT
      ! arrays for time integration
      if (IFTRAN) then
        if (IFFLOW) then
            call invcol2(ABX1,BM1,ntotv)
            call invcol2(ABY1,BM1,ntotv)
            if(IF3D) call invcol2(ABZ1,BM1,ntotv)

            call invcol2(ABX2,BM1,ntotv)
            call invcol2(ABY2,BM1,ntotv)
            if(IF3D) call invcol2(ABZ2,BM1,ntotv)
        endif


        ! temperature and passive scalars
        if (IFHEAT) then
            ! varialbe T; temperature nad passive scalars
            do il=2,NFIELD
                ! arrays time integration VGRADT[12]
                call invcol2(VGRADT1(1,1,1,1,il-1),BM1,ntott)
                call invcol2(VGRADT2(1,1,1,1,il-1),BM1,ntott)
            enddo                 ! passive scalar loop
        endif                     ! IFHEAT
      endif                       ! IFTRAN

      return
      end subroutine
!=======================================================================
!> @brief Multiply right hand sides by BM1 after refining/coarsening
      subroutine amr_rhs_mlt_bm
      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'SOLN'
      include 'MASS'

      ! local variables
      integer ntotv,ntott,il
!-----------------------------------------------------------------------
      ntotv = NX1*NY1*NZ1*NELV
      ntott = NX1*NY1*NZ1*NELT
      ! arrays for time integration
      if (IFTRAN) then
        if (IFFLOW) then
            call col2(ABX1,BM1,ntotv)
            call col2(ABY1,BM1,ntotv)
            if(IF3D) call col2(ABZ1,BM1,ntotv)

            call col2(ABX2,BM1,ntotv)
            call col2(ABY2,BM1,ntotv)
            if(IF3D) call col2(ABZ2,BM1,ntotv)
        endif

        ! temperature and passive scalars
        if (IFHEAT) then
            ! varialbe T; temperature nad passive scalars
            do il=2,NFIELD
                ! arrays time integration VGRADT[12]
                call col2(VGRADT1(1,1,1,1,il-1),BM1,ntott)
                call col2(VGRADT2(1,1,1,1,il-1),BM1,ntott)
            enddo                 ! passive scalar loop
        endif                     ! IFHEAT
      endif                       ! IFTRAN

      return
      end subroutine
!> @brief Make all the interpolated fields continuous
      subroutine amr_h1avg
      implicit none

      include 'SIZE'
      include 'INPUT'           ! IF3D
      include 'SOLN'
      include 'TSTEP'           ! NBDINP
      include 'GEOM'            ! IFGEOM
      include 'AMR'        ! AMR_PRLAG

      ! local variables
      integer nxyz1, ntott, ntotv
      integer il, jl, k, fs, fe
!-----------------------------------------------------------------------
      ! array sizes
      nxyz1=nx1*ny1*nz1
      ntott=nelt*nxyz1
      ntotv=nelv*nxyz1

      ! fix coordinate arrays
      call fix_geom

      ! density (VTRANS) and viscosity (VDIFF) are set every step by
      ! setprop in nek_advance so not needed here

      ifield = 1
      ! velocity
      call amr_oph1_proj(vx,vy,vz,lx1,ly1,lz1,nelv)

      if(IFFLOW.or.IFMHD) then
        ! lag arrays velocity
        do il =1, NBDINP-1
            call amr_oph1_proj(VXLAG (1,1,1,1,il),VYLAG (1,1,1,1,il),
     $           VZLAG (1,1,1,1,il),lx1,ly1,lz1,nelv)
        enddo

         ! pressure; not really sure; check for pn-pn implementation
         if (ifsplit) call h1_proj(PR,nx1,ny1,nz1)  ! continuous pressure

        ! arrays for time integration; not really sure; check for pn-pn implementation
        if (IFTRAN) then
            if (ifsplit) then  ! continuous pressure
              do il =1, NBDINP-1
                call h1_proj(AMR_PRLAG(1,1,1,1,il),nx1,ny1,nz1)
              enddo
            endif

        endif
      endif

      ! temperature and passive scalars
      if (IFHEAT) then
        do il=2,NFIELD
            ifield = il
            call h1_proj(T(1,1,1,1,il-1) ,nx1,ny1,nz1)
        enddo

        ! lag arrays TLAG !!!NOT FINISHED RIGTH NOW!!!
        do il=2,NFIELD
            ifield = il
            do jl =1, NBDINP-1
                call h1_proj(TLAG(1,1,1,1,jl,il-1),nx1,ny1,nz1)
            enddo
        enddo
      endif                     ! IFHEAT

      ! user defined divergence USRDIV; updated in turb_outflow, so not
      ! needed here

      ! mhd
      ! NOT TESTED AND NOT SURE IT IS OK!!!!!!!!!!!!!!!!
      if (IFMHD) then
        ifield = ifldmhd
        ! mhd blocks reside in velocity mesh????
        ! magnetic field
        call amr_oph1_proj(BX,BY,BZ,lx1,ly1,lz1,nelv)

        ! lag arrays
        do il =1, NBDINP-1
            call amr_oph1_proj(BXLAG(1,il),BYLAG(1,il),BZLAG(1,il),
     $       lx1,ly1,lz1,nelv)
        enddo
      endif

      ! it could be place for perturbation, but interpolation would not
      ! give correct base flow structure, so I skip it for now

      return
      end subroutine
!=======================================================================
!> @brief integer8 version of isort
!! @details Use Heap Sort (p 231 Num. Rec., 1st Ed.)
!! @param[inout] as   sorted array
!! @param[out]   ind  permutation
!! @param[in]    nl   array size
      subroutine i8sort(as,ind,nl)
      implicit none

      ! argument list
      integer nl
      integer*8 as(nl)
      integer ind(nl)

      ! local variables
      integer im, ie, ii, il, jl
      integer*8 aa
!-----------------------------------------------------------------------
      ! initial permuatation
      do im=1,nl
         ind(im)=im
      enddo

      if (nl.eq.1) return
      ! initial middle and end position
      im=nl/2+1
      ie=nl
      ! infinite loop
      do
         if (im.gt.1) then
            im = im-1
            aa = as(im)
            ii = ind(im)
         else
            aa =  as(ie)
            ii = ind(ie)
            as(ie) = as(1)
            ind(ie)= ind(1)
            ie=ie-1
            if (ie.eq.1) then
                as(1) = aa
               ind(1) = ii
               return
            endif
         endif
         il=im
         jl=im+im
         do while (jl.le.ie)
            if (jl.lt.ie) then
               if (as(jl).lt.as(jl+1) ) jl=jl+1
            endif
            if (aa.lt.as(jl)) then
                as(il) = as(jl)
               ind(il) = ind(jl)
               il=jl
               jl=jl+jl
            else
               jl=ie+1
            endif
         enddo
         as(il) = aa
         ind(il) = ii
      enddo

      end
!=======================================================================
!> @brief modiffied version of iswap_ip
!! @details In-place reverse permutation
!! @param[inout] as   permuted array
!! @param[in]    pi   permutation
!! @param[in]    nl   array size
      subroutine i8swap_rip(as,pi,nl)
      implicit none

      ! argument list
      integer nl
      integer*8 as(nl)
      integer pi(nl)

      ! local variables
      integer il, jl
      integer*8 as1, as2
      integer j, k, loop_start, last, next
      character*100 str
!-----------------------------------------------------------------------
      do il=1,nl
         if (pi(il).gt.0) then   ! not swapped
            as1 = as(il)
            loop_start = il
            next = pi(il)
            loop : do jl=il,nl
               if (next.lt.0) then
                  write (str,*) 'Hey! i8swap_ip problem.',jl,il,nl,next
#ifdef P4EST
                  call fsc_abort(str)
#endif
               elseif (next.eq.loop_start) then
                  as(next) = as1
                  pi(next) = -pi(next)
                  exit loop
               else
                  as2 = as(next)
                  as(next) = as1
                  as1 = as2
                  last = next
                  next = pi(last)
                  pi(last) = -pi(last)
               endif
            enddo loop
         endif
      enddo

      ! reset permutation mark
      do il=1,nl
         pi(il) = -pi(il)
      enddo
      return
      end
!=======================================================================
!     All edge related routines used IXCN and ESKIP, but this caused some
!     problems as these arrays have to be continuosly updated for different
!     levels in pressure solver, so I add this routine.
      subroutine edgind(istart,istop,iskip,iedg,nx,ny,nz)
      implicit none
      include 'SIZE'
      include 'INPUT'
      include 'TOPOL'

!     argument list
      integer istart,istop,iskip,iedg,nx,ny,nz

!     local variables
      integer ivrt, icx, icy, icz
!-----------------------------------------------------------------------
!     find vertex position
!     start
      ivrt = icedg(1,iedg)
      icx = mod(ivrt +1,2)
      icy = mod((ivrt-1)/2,2)
      icz = (ivrt-1)/4
      istart =  1 + (nx-1)*icx + nx*(ny-1)*icy + nx*ny*(nz-1)*icz

!     stop
      ivrt = icedg(2,iedg)
      icx = mod(ivrt +1,2)
      icy = mod((ivrt-1)/2,2)
      icz = (ivrt-1)/4
      istop =  1 + (nx-1)*icx + nx*(ny-1)*icy + nx*ny*(nz-1)*icz

!     find stride
      if (iedg.le.4) then
         iskip = 1
      elseif (iedg.le.8) then
         iskip = nx
      else
         iskip =nx*nx
      endif

      return
      end
!=======================================================================
!     this subroutine fills edge edg (including vertices) with real cons
      subroutine ctoe(vfld,edg,cons,nx,ny,nz)
      implicit none

      include 'SIZE'
      include 'TOPOL'

!     argument list
      real vfld(nx*ny*nz)
      integer edg,nx,ny,nz
      real cons

!     local variables
      integer is, ie, isk, il
!-----------------------------------------------------------------------
!      is  = IXCN(icedg(1,edg))
!      ie  = IXCN(icedg(2,edg))
!      isk = ESKIP(edg,3)
      call edgind(is,ie,isk,edg,nx,ny,nz)
      do il=is,ie,isk
        vfld(il) = cons
      enddo

      return
      end
!=======================================================================
!     this subroutine fills edge edg (including vertices) with real field
      subroutine vectoe(vfld,edg,vec,nx,ny,nz)
      implicit none

      include 'SIZE'
      include 'TOPOL'

!     argument list
      real vfld(nx*ny*nz)
      integer edg,nx,ny,nz
      real vec(nx)

!     local variables
      integer is, ie, isk, il, jl
!-----------------------------------------------------------------------
!      is  = IXCN(icedg(1,edg))
!      ie  = IXCN(icedg(2,edg))
!      isk = ESKIP(edg,3)
      jl = 1
      call edgind(is,ie,isk,edg,nx,ny,nz)
      do il=is,ie,isk
        vfld(il) = vec(jl)
        jl = jl + 1
      enddo

      return
      end
!=======================================================================
!     this subroutine copies edge edg (including vertices) to real field
      subroutine etovec(vec,edg,vfld,nx,ny,nz)
      implicit none

      include 'SIZE'
      include 'TOPOL'

!     argument list
      real vfld(nx*ny*nz)
      integer edg,nx,ny,nz
      real vec(nx)

!     local variables
      integer is, ie, isk, il, jl
!-----------------------------------------------------------------------
!      is  = IXCN(icedg(1,edg))
!      ie  = IXCN(icedg(2,edg))
!      isk = ESKIP(edg,3)
      call edgind(is,ie,isk,edg,nx,ny,nz)
      jl = 1
      do il=is,ie,isk
        vec(jl) = vfld(il)
        jl = jl + 1
      enddo

      return
      end
!=======================================================================
      subroutine vectof_addl(b,a,iface,nx,ny,nz)
C
C     Copy vector A to the face (IFACE) of B
C     IFACE is the input in the pre-processor ordering scheme.
C
      real A(NX*NZ)
      real B(NX,NY,NZ)

      CALL FACIND (KX1,KX2,KY1,KY2,KZ1,KZ2,NX,NY,NZ,IFACE)
      k = 0
      DO 100 IZ=KZ1,KZ2
      DO 100 IY=KY1,KY2
      DO 100 IX=KX1,KX2
        k = k + 1
        B(IX,IY,IZ) = B(IX,IY,IZ) + A(k)
  100 CONTINUE
      return
      END
c-----------------------------------------------------------------------
      subroutine ftovec_0l(a,b,iface,nx,ny,nz)
C
C     Copy the face (IFACE) of B to vector A.
C     IFACE is the input in the pre-processor ordering scheme.
C
      real A(NX*NZ)
      real B(NX,NY,NZ)

      CALL FACIND (KX1,KX2,KY1,KY2,KZ1,KZ2,NX,NY,NZ,IFACE)
      k = 0
      DO 100 IZ=KZ1,KZ2
      DO 100 IY=KY1,KY2
      DO 100 IX=KX1,KX2
        k = k + 1
        A(k)=B(IX,IY,IZ)
        B(IX,IY,IZ)=0.0
  100 CONTINUE
      return
      END
c-----------------------------------------------------------------------
      subroutine ftovecl(a,b,iface,nx,ny,nz)
C
C     Copy the face (IFACE) of B to vector A.
C     IFACE is the input in the pre-processor ordering scheme.
C
      real A(NX*NZ)
      real B(NX,NY,NZ)
      CALL FACIND (KX1,KX2,KY1,KY2,KZ1,KZ2,NX,NY,NZ,IFACE)
      k = 0
      DO 100 IZ=KZ1,KZ2
      DO 100 IY=KY1,KY2
      DO 100 IX=KX1,KX2
        k = k + 1
        A(k)=B(IX,IY,IZ)
  100 CONTINUE
      return
      END
c-----------------------------------------------------------------------
      subroutine vectofl(b,a,iface,nx,ny,nz)
C
C     Copy vector A to the face (IFACE) of B
C     IFACE is the input in the pre-processor ordering scheme.
C
      real A(NX*NZ)
      real B(NX,NY,NZ)
      CALL FACIND (KX1,KX2,KY1,KY2,KZ1,KZ2,NX,NY,NZ,IFACE)
      k = 0
      DO 100 IZ=KZ1,KZ2
      DO 100 IY=KY1,KY2
      DO 100 IX=KX1,KX2
        k = k + 1
        B(IX,IY,IZ) = A(k)
  100 CONTINUE
      return
      END

!=======================================================================
