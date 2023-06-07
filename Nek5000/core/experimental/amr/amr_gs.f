!> @file amr_gs.f
!! @ingroup nekamr
!! @brief Gater-scatter related routines
!! @author Adam Peplinski
!! @date Jun 28, 2016
!=======================================================================
!> @brief Get topology information.
!! @details Get global vertex, face and edge numberring together with
!!  hanging node and orinetation information for faces and edges.
      subroutine amr_topol_get()
      implicit none
!-----------------------------------------------------------------------
      ! global node numberring; hangign node mark
      call amr_node_get()

      ! face, edge orientation
      call amr_algn_get()

      ! mark element families sharing the same parent
      call amr_fml_get()

      ! get vertex to crs base mapping
      call amr_gen_crs_map()

      return
      end subroutine
!=======================================================================
!> @brief Get node information.
!! @details Get global vertex, face and edge numberring together with
!!  hanging node mark.
!! @todo Check strange transpose of vertex numberring caused by data
!!  redistribution. It does not cause problems, but is unexpected.
      subroutine amr_node_get()
      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'
      include 'AMR'
      include 'AMR_TOPOL'

      ! local variables
      integer lnelt ! local number of elements; p4est count
      integer lnoden ! local number of indep. nodes; p4est count
      integer gnoden ! local number of independent nodes not shared with other proc.
      integer vnode ! number of nodes; counted vertices, faces, edges
      parameter (vnode = 3**LDIM-1)
      integer*8 lnodes(vnode*LELT) ! global numberring of nodes

      integer el, il, jl ! loop index
      integer itmp
      integer ntot

      ! node splitting and sorting
      integer*8 itmp8, lnvrt8, lnfcs8, lnedg8
      integer*8 glnodes ! global number of independent nodes
      integer prm(vnode*LELT) ! permutation array
      integer*8 unodes(vnode*LELT) ! unique local nodes
      integer ioff(vnode*LELT+1) ! offset position

      ! communicator
      integer gs_handle ! gather-scatter handle

      integer mid,mp,nekcomm,nekgroup,nekreal
      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal

      ! functions
      integer*8 i8glsum, i8gl_running_sum, mod1

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
      call amr_log(AMR_LP_PRD,'Get node numberring.')

#ifdef P4EST
      ! create gll points numbering for topology
      ! I set degree to -LDIM to be able to distinguish between vertices,
      ! edges and faces.
      ! Element interior is discarded.
      itmp = -LDIM
      call fp4est_lnodes_new(itmp)

      ! get mesh topology information to nek
#if N_DIM == 3
      call fp4est_msh_get_tplg(lnelt, lnoden, gnoden, lnodes,
     $      AMR_HNG_ELM, AMR_HNG_FCS, AMR_HNG_EDG)
#else
      call fp4est_msh_get_tplg(lnelt, lnoden, gnoden, lnodes,
     $      AMR_HNG_ELM, AMR_HNG_FCS)
#endif

      ! destroy p4est lnode
      call fp4est_lnodes_del()

      ! test mesh size with respect to p4est
      if (lnelt.ne.AMR_NELT) call amr_abort
     $     ('Error: amr_glbnr_get; lnelt /= AMR_nelt face')

      ! create vertex information
      call fp4est_nodes_new()

      ! get hangign vertex information
      call fp4est_msh_get_node(lnelt, AMR_HNG_VRT)

      ! destroy p4est node
      call fp4est_nodes_del()

      ! test mesh size with respect to p4est
      if (lnelt.ne.AMR_NELT) call amr_abort
     $     ('Error: amr_glbnr_get; lnelt /= AMR_nelt vertex')
#endif
      ! data transfer
      ! this data transfer causes strange transpose in vertex numberring
      ! for now I can't find the reason
      call amr_node_transfer(lnodes,vnode,lnelt)

      ! p4est provides continuous node numberring without distiguisng between
      ! different node type. Split and count nodes in groups: vertices, faces
      ! and edges.

      ! get global number of unique nodes
      itmp8 = gnoden
      glnodes = i8glsum(itmp8,1)

      ! sort nodes
      ntot = vnode*lnelt
      call i8sort(lnodes,prm,ntot)
      ! get unique nodes and offset
      itmp = 1 ! unique node counter
      unodes(1) = lnodes(1)
      ioff(1) = 1
      do il=2,ntot
         if (unodes(itmp).ne.lnodes(il)) then
            itmp = itmp +1
            unodes(itmp) = lnodes(il)
            ioff(itmp) = il
         endif
      enddo
      ioff(itmp+1) = ntot + 1

      ! update local node number
      lnoden = itmp

      ! set gather-scatter communicator
      call fgslib_gs_setup(gs_handle,unodes,lnoden,nekcomm,mp)

      ! mark nodes with process id
      itmp8 = NID
      do il=1,lnoden
         unodes(il) = itmp8
      enddo

      ! find min for given node
      call fgslib_gs_op(gs_handle,unodes,3,3,0)

      ! count local and mark non-local nodes
      AMR_GLBNR_LFCS = 0
      AMR_GLBNR_LEDG = 0
      AMR_GLBNR_LVRT = 0
      do il=1,lnoden
         if (unodes(il).eq.itmp8) then
            itmp = mod1(prm(ioff(il)),vnode)
            if (itmp.le.AMR_NFCS) then ! face
               AMR_GLBNR_LFCS = AMR_GLBNR_LFCS  + 1
               unodes(il) = AMR_GLBNR_LFCS
            else if (itmp.le.(AMR_NFCS + AMR_NEDG)) then ! edge
               AMR_GLBNR_LEDG = AMR_GLBNR_LEDG  + 1
               unodes(il) = AMR_GLBNR_LEDG
            else if (itmp.le.(AMR_NFCS + AMR_NEDG + AMR_NVRT)) then ! vertex
               AMR_GLBNR_LVRT = AMR_GLBNR_LVRT  + 1
               unodes(il) = AMR_GLBNR_LVRT
            else
               call amr_abort
     $             ('ERROR: amr_topol_get; wrong prm')
            endif
         else
            unodes(il) = 0
         endif
      enddo

      ! get global counts
      AMR_GLBNR_GFCS = i8glsum(AMR_GLBNR_LFCS,1)
      AMR_GLBNR_GEDG = i8glsum(AMR_GLBNR_LEDG,1)
      AMR_GLBNR_GVRT = i8glsum(AMR_GLBNR_LVRT,1)

      ! correctness test
      if (glnodes.ne.(AMR_GLBNR_GFCS+AMR_GLBNR_GEDG+AMR_GLBNR_GVRT))
     $ call amr_abort
     $             ('ERROR: amr_topol_get; inconsistent glnodes')

      ! get number of nodes on processes with lower id
      lnfcs8 = i8gl_running_sum(AMR_GLBNR_LFCS) - AMR_GLBNR_LFCS
      lnedg8 = i8gl_running_sum(AMR_GLBNR_LEDG) - AMR_GLBNR_LEDG
      lnvrt8 = i8gl_running_sum(AMR_GLBNR_LVRT) - AMR_GLBNR_LVRT

      ! renumber nodes
      do il=1,lnoden
         if (unodes(il).gt.0) then
            itmp = mod1(prm(ioff(il)),vnode)
            if (itmp.le.AMR_NFCS) then ! face
               unodes(il) = unodes(il) + lnfcs8
            else if (itmp.le.(AMR_NFCS + AMR_NEDG)) then ! edge
               unodes(il) = unodes(il) + lnedg8
            else if (itmp.le.(AMR_NFCS + AMR_NEDG + AMR_NVRT)) then ! vertex
               unodes(il) = unodes(il) + lnvrt8
            else

            endif
         else
            unodes(il) = 0
         endif
      enddo

      ! redistribute node numbers
      call fgslib_gs_op(gs_handle,unodes,3,1,0)

      ! put node numberring back
      do il=1,lnoden
         do jl = ioff(il),ioff(il + 1) - 1
            lnodes(jl) = unodes(il)
         enddo
      enddo

      ! reverse permute lnoden
      call i8swap_rip(lnodes,prm,ntot)

      ! extract numberring
      itmp = 0
      do el=1,lnelt
         ! face
         do il = 1,AMR_NFCS
            itmp = itmp + 1
            AMR_GLBNR_FCS(il,el) = lnodes(itmp)
         enddo
         ! edge
#if N_DIM == 3
         do il = 1,AMR_NEDG
            itmp = itmp + 1
            AMR_GLBNR_EDG(il,el) = lnodes(itmp)
         enddo
#endif
         ! vertex
         do il = 1,AMR_NVRT
            itmp = itmp + 1
            AMR_GLBNR_VRT(il,el) = lnodes(itmp)
         enddo
      enddo

      ! free communicator
      call fgslib_gs_free (gs_handle)

#ifdef DEBUG
      ! for testing
      ! to output refinement
      icalld = icalld+1
      call io_file_freeid(iunit, ierr)
      write(str1,'(i3.3)') NID
      write(str2,'(i3.3)') icalld
      open(unit=iunit,file='topol.txt'//str1//'i'//str2)

      write(iunit,*) lnelt, glnodes
      write(iunit,*) AMR_GLBNR_GVRT, AMR_GLBNR_GFCS, AMR_GLBNR_GEDG
      write(iunit,*) AMR_GLBNR_LVRT, AMR_GLBNR_LFCS, AMR_GLBNR_LEDG
      write(iunit,*) 'FACE'
      do el=1,lnelt
         write(iunit,*) el,LGLEL(el),
     $                 (AMR_GLBNR_FCS(il,el),il=1,AMR_NFCS)
      enddo
#if N_DIM == 3
      write(iunit,*) 'EDGE'
      do el=1,lnelt
         write(iunit,*) el,LGLEL(el),
     $                 (AMR_GLBNR_EDG(il,el),il=1,AMR_NEDG)
      enddo
#endif
      write(iunit,*) 'VERTEX'
      do el=1,lnelt
         write(iunit,*) el,LGLEL(el),
     $                 (AMR_GLBNR_VRT(il,el),il=1,AMR_NVRT)
      enddo
      write(iunit,*) 'VERTEX HANGING'
      do el=1,lnelt
         write(iunit,*) el,LGLEL(el),
     $                 (AMR_HNG_VRT(il,el),il=1,AMR_NVRT)
      enddo
      write(iunit,*) 'FACE HANGING'
      do el=1,lnelt
         write(iunit,*) el,LGLEL(el),
     $                 (AMR_HNG_FCS(il,el),il=1,AMR_NFCS)
      enddo
#if N_DIM == 3
      write(iunit,*) 'EDGE HANGING'
      do el=1,lnelt
         write(iunit,*) el,LGLEL(el),
     $                 (AMR_HNG_EDG(il,el),il=1,AMR_NEDG)
      enddo
#endif
      close(iunit)
#endif
      return
      end subroutine
#undef DEBUG
!=======================================================================
!> @brief Get face, edge alignment.
      subroutine amr_algn_get()
      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'TOPOL'
!#define DEBUG
#ifdef DEBUG
      include 'PARALLEL'
#endif
      include 'AMR'
      include 'AMR_TOPOL'

      ! local variables
      integer lnelt, lnelt2 ! local number of elements; p4est count
      integer el, il, jl, kl ! loop index
      integer itmp1, itmp2, iskip ! dummy variables
      integer degree, vnode ! polynomial degree; number of nodes
      parameter (degree = 4,vnode = degree**LDIM)  ! degree
      integer*8 lnodes(vnode*LELT) ! global numberring of nodes
      integer*8 ledge(degree)

#ifdef DEBUG
      character*3 str1, str2
      integer iunit, ierr
      ! call number
      integer icalld
      save icalld
      data icalld /0/
#endif
!-----------------------------------------------------------------------
      call amr_log(AMR_LP_PRD,'Get element alignment.')

      ! reset arrays
      il = LELT*AMR_NFCS
      call izero(AMR_ALGN_FCS,il)
#if N_DIM == 3
      il = LELT*AMR_NEDG
      call izero(AMR_ALGN_EDG,il)
#endif
#ifdef P4EST
      ! get face orientation
      call fp4est_msh_get_algn(AMR_ALGN_FCS,lnelt)
#endif
#if N_DIM == 3
      ! p4est marks as 'edges' the real edges in the mesh for which the
      ! face-to-face communication does not suffice for data redistribution.
      ! This causes some problem in edge orientation mark. That is why
      ! I let p4est to create simple node numbering and I extract edge
      ! orientation out of it.
      ! create gll points numbering
      itmp1 = degree-1
#ifdef P4EST
      call fp4est_lnodes_new(itmp1)

      ! get lnodes numbering to nekton
      call fp4est_msh_get_lnode(lnelt2, itmp2, lnodes)

      ! destroy p4est lnode
      call fp4est_lnodes_del()
#endif
      ! check calls consistency
      if (vnode.ne.itmp2) then
         call amr_abort('Error: amr_algn_get degree inconsistent')
      endif

      if (lnelt.ne.lnelt2) then
         call amr_abort('Error: amr_algn_get lnelt inconsistent')
      endif

      ! extract edge orientation
      do el=1,lnelt
         do il=1,AMR_NEDG
            ! extract edge
            call edgind(itmp1,itmp2,iskip,il,degree,degree,degree)
            kl=1
            do jl=itmp1,itmp2,iskip
               ledge(kl) = lnodes(jl + (el-1)*vnode)
               kl = kl + 1
            enddo
            ! check orientation
            if (ledge(3).gt.ledge(2)) then
               AMR_ALGN_EDG(il,el) = 0
            else
               AMR_ALGN_EDG(il,el) = 1
            endif
         enddo
      enddo
#endif

      ! redistribute data
      call amr_algn_transfer(lnelt)

#ifdef DEBUG
      ! for testing
      ! to output refinement
      icalld = icalld+1
      call io_file_freeid(iunit, ierr)
      write(str1,'(i3.3)') NID
      write(str2,'(i3.3)') icalld
      open(unit=iunit,file='aligment.txt'//str1//'i'//str2)
      write(iunit,*) lnelt
      write(iunit,*) 'FACE'
      do el=1,lnelt
         write(iunit,*) el,LGLEL(el),
     $         (AMR_ALGN_FCS(il,el),il=1,AMR_NFCS)
      enddo
#if N_DIM == 3
      write(iunit,*) 'EDGE'
      do el=1,lnelt
         write(iunit,*) el,LGLEL(el),
     $         (AMR_ALGN_EDG(il,el),il=1,AMR_NEDG)
      enddo
#endif
      close(iunit)

#endif
      return
      end subroutine
!#undef DEBUG
!=======================================================================
!> @brief Mark element families sharing the same parent
!! @details This routine provides information about sets of elements that
!!  share the same parent and could be destroyed togeter during coarsening
!!  step.
      subroutine amr_fml_get()
      implicit none

      include 'SIZE'
      include 'PARALLEL'   ! debugging only
      include 'TOPOL'      ! debugging only
      include 'AMR'
      include 'AMR_TOPOL'

      ! local variables
      integer ntot, itmp
      integer lnelf_off, lnelf, lfamily(LELT)
      integer il, jl

      ! functions
      integer igl_running_sum

!#define DEBUG
#ifdef DEBUG
      character*3 str1, str2
      integer iunit, ierr, iel
      ! call number
      integer icalld
      save icalld
      data icalld /0/
      integer mid,mp,nekcomm,nekgroup,nekreal
      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal
      integer gs_handle_fcs
      integer lfml_vrt(AMR_NFCS,LELT)
      integer ivface(3,8)
      save ivface
      data ivface /1,3,5,  2,3,5,  1,4,5,  2,4,5,
     $             1,3,6,  2,3,6,  1,4,6,  2,4,6/
#endif
!-----------------------------------------------------------------------
      ! reset family flag
      ntot = 2*LELT
      call izero(AMR_FML_MARK,ntot)
#ifdef P4EST
      ! get list of element families from p4est
      call fp4est_family_get(lfamily,lnelf)
#endif
      ! consistency check
      if (lnelf.gt.AMR_NELT.or.mod(lnelf,AMR_NVRT).ne.0) then
         call amr_abort('Error: amr_fml_get; wrong lnelf')
      endif

      ! get global parent offset
      lnelf = lnelf/AMR_NVRT
      lnelf_off = igl_running_sum(lnelf) - lnelf

      ! fill in AMR_FML_MARK array
      do il=1,lnelf
         lnelf_off = lnelf_off + 1
         do jl=1,AMR_NVRT
            itmp = lfamily(jl + (il-1)*AMR_NVRT) - AMR_nelit
#ifdef DEBUG
            if (itmp.le.0.or.itmp.gt.AMR_nelt)
     $     call amr_abort('Error: amr_fml_get; wrong el. nr.')
#endif
            AMR_FML_MARK(1,itmp) = lnelf_off
            AMR_FML_MARK(2,itmp) = AMR_NVRT + 1 - jl
         enddo
      enddo

#ifdef DEBUG
      ! for testing
      ! to output family mark
      icalld = icalld+1
      call io_file_freeid(iunit, ierr)
      write(str1,'(i3.3)') NID
      write(str2,'(i3.3)') icalld
      open(unit=iunit,file='fmlp4.txt'//str1//'i'//str2)

      do iel=1,AMR_NELT
         write(iunit,*) iel, AMR_nelit+iel, AMR_FML_MARK(:,iel)
      enddo

      close(iunit)
#endif


      ! redistribute family mark
      lnelf = AMR_nelt
      call amr_fml_transfer(lnelf)

#ifdef DEBUG
      ! for testing
      ! check if the family vertex is set correctly
      ! set up communicator
      lnelf = AMR_NFCS*NELT
      call fgslib_gs_setup(gs_handle_fcs,AMR_GLBNR_FCS,lnelf,
     $     nekcomm,mp)

      ! fill own parent number
      do il=1,NELT
         do jl=1,AMR_NFCS
            lfml_vrt(jl,il) = AMR_FML_MARK(1,il)
         enddo
      enddo

      ! sum with the neighbours
      call fgslib_gs_op(gs_handle_fcs,lfml_vrt,2,1,0)

      ! extract neighbour parent number
      do il=1,NELT
         do jl=1,AMR_NFCS
            lfml_vrt(jl,il) = lfml_vrt(jl,il) - AMR_FML_MARK(1,il)
         enddo
      enddo

      ! compare face values connected to the vertex
      do il=1, NELT
         if (AMR_FML_MARK(1,il).gt.0) then
            do jl=1, NDIM
               if (AMR_FML_MARK(1,il).ne.
     $          lfml_vrt(ivface(jl,AMR_FML_MARK(2,il)),il)) then
                  write(*,*) 'FAMILY VERTEX TEST ERROR',
     $             nid,il,jl,AMR_FML_MARK(1,il),
     $             AMR_FML_MARK(2,il),ivface(jl,AMR_FML_MARK(2,il)),
     $             lfml_vrt(ivface(jl,AMR_FML_MARK(2,il)),il)
               endif
            enddo
         endif
      enddo

      call fgslib_gs_free (gs_handle_fcs)

      ! to output family mark
      call io_file_freeid(iunit, ierr)
      write(str1,'(i3.3)') NID
      write(str2,'(i3.3)') icalld
      open(unit=iunit,file='fmlnek.txt'//str1//'i'//str2)

      do iel=1,NELT
         write(iunit,*) iel, LGLEL(iel), AMR_FML_MARK(:,iel)
      enddo

      close(iunit)
#endif

      return
      end subroutine
#undef DEBUG
!=======================================================================
!> @brief Generate global GLL points numbering to generate communicator
!! @details I use symmetric vertex, face and edge numberring described in
!!  setedge. It is consistent with p4est ordering as well. In genral
!! the numberring is consistent with memory aligment but takes into account
!! face and edge orientation. This routine replaces set_vert.
!! @param[out]   glo_num    global node numberring
!! @param[out]   ngv        number of unique nodes
!! @param[in]    nx         number of points in element along single dimension
!! @param[in]    nel        element number
!! @param[in]    ifcenter   do we include element interior
!! @todo Test face and edge orientation in 3D
      subroutine amr_set_vert(glo_num,ngv,nx,nel,ifcenter)
      implicit none

      include 'SIZE'
      include 'PARALLEL'
      include 'AMR'
      include 'AMR_TOPOL'

      ! argument list
      integer nx, nel
      logical ifcenter
      integer*8 glo_num(nx**LDIM*nel),ngv

      ! local variables
      integer nn, mm ! counter
      integer el, il, jl, kl ! loop index
      integer jini, jend, jinc  ! loop limits
      integer kini, kend, kinc  ! loop limits
      logical iftranspose
      integer nxl, els, nstrd(AMR_NMAX), nstrdr(AMR_NMAX) ! data stride
      integer nstrt(AMR_NMAX) ! data start

      integer ny, nz

      ! functions
      integer iglsum

!#define DEBUG
#ifdef DEBUG
!     for testing
      character*3 str1, str2
      integer iunit, ierr
      ! call number
      integer icalld
      save icalld
      data icalld /0/
#endif
!-----------------------------------------------------------------------
      ! exclude vertices
      nxl = nx - 2
      ! reset global numberring
      nn = nx**LDIM*nel
      call i8zero(glo_num,nn)

      ! check if we have to work with vertices only
      if (nx.eq.2) then
         ! vertices only
         call i8copy(glo_num,AMR_GLBNR_VRT,nn)
         ngv = AMR_GLBNR_GVRT
      elseif (nx.gt.2) then
         ! vertices
         nn = 0
         ! data stride
         nstrd(1) = 1
         nstrd(2) = nx-1
         nstrd(3) = nx*(nx-2) + 1
         nstrd(4) = nx -1
#if N_DIM == 3
         nstrd(5) = nx*nx*(nx-2) + 1
         nstrd(6) = nx-1
         nstrd(7) = nx*(nx-2) + 1
         nstrd(8) = nx -1
#endif
         do el=1,nel
            do il=1,AMR_NVRT
               ! vertex posiotin in the element
               nn = nn + nstrd(il)
               glo_num(nn) = AMR_GLBNR_VRT(il,el)
            enddo
         enddo
         ngv = AMR_GLBNR_GVRT
         ! faces
#if N_DIM == 2
         els = nx*nx
         ! data stride
         nstrd(1) = nx
         nstrd(2) = nx
         nstrd(3) = 1
         nstrd(4) = 1
         ! data start
         nstrt(1) = 1
         nstrt(2) = nx
         nstrt(3) = 1
         nstrt(4) = nx*(nx-1) + 1
         do el=1,nel
            ! element start
            mm = (el-1) * els
            do il=1,AMR_NFCS
               ! initial face posiotin in the element
               nn = mm +nstrt(il)
               ! check aligment
               if (AMR_ALGN_FCS(il,el).eq.0) then
                  jini = 1
                  jend = nxl
                  jinc = 1
               else if (AMR_ALGN_FCS(il,el).eq.1) then
                  jini = nxl
                  jend = 1
                  jinc = -1
               else
                  call amr_abort
     $             ('Error: amr_set_vert; wrong face alignment.')
               endif
               do jl=jini,jend,jinc
                  nn = nn + nstrd(il)
                  glo_num(nn) = ngv + nxl*(AMR_GLBNR_FCS(il,el)-1) + jl
               enddo

            enddo
         enddo
         ngv = ngv + AMR_GLBNR_GFCS*nxl
         ! element centre
         if (ifcenter) then
            do el=1,nel
               ! element start + initial centre posiotin in the element
               nn = (el-1) * els + nx + 1

                  do jl=1,nxl
                     do il=1,nxl
                        nn = nn + 1
                        glo_num(nn) = ngv + ((LGLEL(el)-1)*nxl +
     $                                             jl - 1)*nxl + il
                     enddo
                     nn = nn + 2
                  enddo
            enddo
            ! get global sum of counted vertices
            ! I do it because I have no information which mesh (V;T) is used
            mm =  iglsum(nel,1)
            ngv = ngv + mm*nxl*nxl
         endif
#else
         els = nx*nx*nx
         ! data stride; within 1D row
         nstrd(1) = nx
         nstrd(2) = nx
         nstrd(3) = 1
         nstrd(4) = 1
         nstrd(5) = 1
         nstrd(6) = 1
         ! data stride; between 1D rows
         nstrdr(1) = 2*nx
         nstrdr(2) = 2*nx
         nstrdr(3) = nx*(nx-1) + 2
         nstrdr(4) = nx*(nx-1) + 2
         nstrdr(5) = 2
         nstrdr(6) = 2
         ! data start
         nstrt(1) = nx*nx + 1
         nstrt(2) = nx*(nx + 1)
         nstrt(3) = nx*nx + 1
         nstrt(4) = nx*(2*nx-1) + 1
         nstrt(5) = nx + 1
         nstrt(6) = nx*nx*(nx-1) + nx + 1
         do el=1,nel
            ! element start
            mm = (el-1) * els
            do il=1,AMR_NFCS
               ! initial face posiotin in the element
               nn = mm +nstrt(il)
               ! check aligment
               if (AMR_ALGN_FCS(il,el).eq.0) then ! identity
                  iftranspose = .FALSE.
                  jini = 1
                  jend = nxl
                  jinc = 1
                  kini = 1
                  kend = nxl
                  kinc = 1
               else if (AMR_ALGN_FCS(il,el).eq.1) then ! transpose (T)
                  iftranspose = .TRUE.
                  jini = 1
                  jend = nxl
                  jinc = 1
                  kini = 1
                  kend = nxl
                  kinc = 1
               else if (AMR_ALGN_FCS(il,el).eq.2) then ! permutation in x (P_x)
                  iftranspose = .FALSE.
                  jini = nxl
                  jend = 1
                  jinc = -1
                  kini = 1
                  kend = nxl
                  kinc = 1
               else if (AMR_ALGN_FCS(il,el).eq.3) then ! P_x T
                  iftranspose = .TRUE.
                  jini = 1
                  jend = nxl
                  jinc = 1
                  kini = nxl
                  kend = 1
                  kinc = -1
               else if (AMR_ALGN_FCS(il,el).eq.4) then ! P_y T
                  iftranspose = .TRUE.
                  jini = nxl
                  jend = 1
                  jinc = -1
                  kini = 1
                  kend = nxl
                  kinc = 1
               else if (AMR_ALGN_FCS(il,el).eq.5) then ! P_y
                  iftranspose = .FALSE.
                  jini = 1
                  jend = nxl
                  jinc = 1
                  kini = nxl
                  kend = 1
                  kinc = -1
               else if (AMR_ALGN_FCS(il,el).eq.6) then ! P_x P_y T
                  iftranspose = .TRUE.
                  jini = nxl
                  jend = 1
                  jinc = -1
                  kini = nxl
                  kend = 1
                  kinc = -1
               else if (AMR_ALGN_FCS(il,el).eq.7) then ! P_x P_y
                  iftranspose = .FALSE.
                  jini = nxl
                  jend = 1
                  jinc = -1
                  kini = nxl
                  kend = 1
                  kinc = -1
               else
                  call amr_abort
     $             ('Error: amr_set_vert; wrong face alignement.')
               endif

               if (iftranspose) then
                  do kl=kini,kend,kinc
                     do jl=jini,jend,jinc
                        nn = nn + nstrd(il)
                        glo_num(nn) = ngv + (nxl*
     $                         (AMR_GLBNR_FCS(il,el)-1) +jl-1)*nxl + kl
                     enddo
                     nn = nn + nstrdr(il)
                  enddo
               else
                  do kl=kini,kend,kinc
                     do jl=jini,jend,jinc
                        nn = nn + nstrd(il)
                        glo_num(nn) = ngv + (nxl*
     $                         (AMR_GLBNR_FCS(il,el)-1) +kl-1)*nxl + jl
                     enddo
                     nn = nn + nstrdr(il)
                  enddo
               endif
            enddo
         enddo
         ngv = ngv + AMR_GLBNR_GFCS*nxl*nxl
         ! edge
         ! data stride
         nstrd(1)  = 1
         nstrd(2)  = 1
         nstrd(3)  = 1
         nstrd(4)  = 1
         nstrd(5)  = nx
         nstrd(6)  = nx
         nstrd(7)  = nx
         nstrd(8)  = nx
         nstrd(9)  = nx*nx
         nstrd(10) = nx*nx
         nstrd(11) = nx*nx
         nstrd(12) = nx*nx
         ! data start
         nstrt(1)  = 1
         nstrt(2)  = nx*(nx-1) + 1
         nstrt(3)  = nx*nx*(nx-1) + 1
         nstrt(4)  = nx*(nx*nx-1) + 1
         nstrt(5)  = 1
         nstrt(6)  = nx
         nstrt(7)  = nx*nx*(nx-1) + 1
         nstrt(8)  = nx*nx*(nx-1) + nx
         nstrt(9)  = 1
         nstrt(10) = nx
         nstrt(11) = nx*(nx-1) + 1
         nstrt(12) = nx*nx
         do el=1,nel
            ! element start
            mm = (el-1) * els
            do il=1,AMR_NEDG
               ! initial edge posiotin in the element
               nn = mm +nstrt(il)
               ! check aligment
               if (AMR_ALGN_EDG(il,el).eq.0) then
                  jini = 1
                  jend = nxl
                  jinc = 1
               else if (AMR_ALGN_EDG(il,el).eq.1) then
                  jini = nxl
                  jend = 1
                  jinc = -1
               else
                  call amr_abort
     $             ('Error: amr_set_vert; wrong edge alignment.')
               endif
               do jl=jini,jend,jinc
                  nn = nn + nstrd(il)
                  glo_num(nn) = ngv + nxl*(AMR_GLBNR_EDG(il,el)-1) + jl
               enddo
            enddo
         enddo
         ngv = ngv + AMR_GLBNR_GEDG*nxl
         ! element centre
         if (ifcenter) then
            do el=1,nel
               ! element start + initial centre posiotin in the element
               nn = (el-1) * els + nx*(nx+1) + 1
               do kl=1,nxl
                  do jl=1,nxl
                     do il=1,nxl
                        nn = nn + 1
                        glo_num(nn) = ngv + (((LGLEL(el)-1)*nxl +
     $                                  kl-1)*nxl + jl - 1)*nxl + il
                     enddo
                     nn = nn + 2
                  enddo
                  nn = nn + 2*(nx+1)
               enddo
            enddo
            ! get global sum of conted vertices
            ! I do it because I have no information which mesh (V;T) is used
            mm =  iglsum(nel,1)
            ngv = ngv + mm*nxl*nxl*nxl
         endif
#endif
      else ! wrong nx; it cannot be lower thant 2
         call amr_abort
     $   ('Error: amr_set_vert; nx must be >= 2')
      endif

      ! to follow nek5000 user interface
      ny = nx
#if N_DIM == 3
      nz = nx
#else
      nz = 1
#endif
!!!!      if(NIO.eq.0) write(*,*) 'call usrsetvert'
!!!!      call usrsetvert(glo_num,nel,nx,nx,nz)
!!!!      if(NIO.eq.0) write(*,'(A,/)') ' done :: usrsetvert'


#ifdef DEBUG
      ! for testing
      icalld = icalld +1
      call io_file_freeid(iunit, ierr)
      write(str1,'(i3.3)') NID
      write(str2,'(i3.3)') icalld
      open(unit=iunit,file='setvert.txt'//str1//'i'//str2)
      write(iunit,*) nel, nx,  ngv, LDIM
      nn = 0
      do el=1,nel
         write(iunit,*) 'ELEMENT ',el,LGLEL(el)
#if N_DIM == 2
         do jl=1,nx
            write(iunit,*) (glo_num(nn + il),il=1,nx)
            nn = nn + nx
         enddo
#else
         do kl=1,nx
            write(iunit,*) 'K ',kl
            do jl=1,nx
               write(iunit,*) (glo_num(nn + il),il=1,nx)
               nn = nn + nx
            enddo
         enddo
#endif
      enddo
      close(iunit)
#endif
!#undef DEBUG
      return
      end subroutine
!=======================================================================
!> @brief Generate gs setup for given point distributions
!! @details This routine replaces setupds and it is copy of this routine.
!!   In general it is not necessary as modification of setupds with
!!   preprocessing flag switching between set_vert and amr_set_vert
!!   shoud be sufficient, but moving meshes in setup_mesh_dssum
!!   do modify global point numbering playing with vertex numbering. That
!!   is why I keep original version of node numbering separate form
!!   nekamr one.
!! @param[out]   gs_handle  handle to gs setup
!! @param[in]    nx         number of points in element along X dimension
!! @param[in]    ny         number of points in element along Y dimension
!! @param[in]    nz         number of points in element along Z dimension
!! @param[in]    nel        local element number
!! @param[in]    melg       global element number
!! @param[out]   glo_num    global node numberring
      subroutine amr_setupds(gs_handle,nx,ny,nz,nel,melg,glo_num)
      implicit none

      include 'SIZE'
      include 'AMR'

      ! argument list
      integer gs_handle
      integer nx, ny, nz, nel, melg
      integer*8 glo_num(1)

      ! global variables
      integer mid,mp,nekcomm,nekgroup,nekreal
      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal

      ! local variables
      integer ntot
      integer*8 ngv  ! number of unique nodes
      real t0, t1
      character*70 str

      ! functions
      real dnekclock
!-----------------------------------------------------------------------
      t0 = dnekclock()

      ! Global-to-local mapping for gs
      call amr_set_vert(glo_num,ngv,nx,nel,.false.)

      ! Initialize gather-scatter code
      ntot      = nx*ny*nz*nel
      call fgslib_gs_setup(gs_handle,glo_num,ntot,nekcomm,mp)

      t1 = dnekclock() - t0
      write(str,1) t1,gs_handle,nx,ngv,melg
      call amr_log(AMR_LP_PRD,trim(str))
 1    format('setupds time',1pe11.4,' seconds ',2i3,2i12)

      return
      end subroutine
!=======================================================================
!> @brief Get global vertex numbering.
!! @details This routine replaces get_vert_map called in get_vert.
!! @param[out]   vertex    global vertex numberring
      subroutine amr_vertex_get(vertex)
      implicit none

      include 'SIZE'
      include 'AMR'
      include 'AMR_TOPOL'

      ! argument list
      integer vertex (AMR_NVRT*LELT)

      ! locall variables
      integer el, il, itmp
!-----------------------------------------------------------------------
      call amr_log(AMR_LP_PRD,'Get vertex numbering')

      itmp = 0
      do el = 1, NELT
         do il = 1, AMR_NVRT
            itmp = itmp + 1
            vertex(itmp) = AMR_GLBNR_VRT(il,el)
         enddo
      enddo

      return
      end subroutine
!=======================================================================
