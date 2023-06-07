!> @file amr_part.f
!! @ingroup nekamr
!! @brief Partitioning scheme for amr including 2-level partitioning
!! @author Adam Peplinski
!! @date Jun 1, 2018
!
! This file requires some of constant defined in include file
#include "nekamr.h"
!=======================================================================
!> @brief Element-processor mapping
      subroutine amr_get_map()
      implicit none

      include 'SIZE'
      include 'AMR'
!-----------------------------------------------------------------------
      ! partition mesh
#ifdef AMR_PRT_GRPH
      call amr_graph_map()
#else
      call amr_vert_map()
#endif
      ! update mapping count
      AMR_IMAP = AMR_IMAP + 1

      return
      end subroutine
!=======================================================================
!> @brief Original Nek5000 element-processor mapping adopted to AMR scheme
!! @remarks This routine uses global scratch space SCRMG
      subroutine amr_vert_map()
      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'
      include 'AMR'
      include 'AMR_SHM'

      ! local variables
      integer itmp, el, il, jl, ierr
      integer lnelt
      integer iwork(lelt)             ! work array for sorting
      integer ibuf(2)                 !
      integer*8 lnodes(AMR_NVRT,LELT) ! global numberring of vertices
      integer*8 eid8(LELT)            ! global element numbers
      integer imsg(1)                 ! message request
      integer itmpv1(1)               ! dummy array

      ! global variables
      integer mid,mp,nekcomm,nekgroup,nekreal
      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal

      integer vertex(AMR_NVRT,LELT)
      common /ivrtx/ vertex

      ! scratch arrays
      ! current global element offset on p4est side; this could be relatively big array
      integer lelm_goff(0:AMR_SHM_MAX*AMR_NDS_MAX)
      integer gnelt(LELT) ,elmrk(LELT) ,eldst(LELT)
      common /scrmg/ lelm_goff,gnelt,elmrk,eldst

      ! approximate element centre coordinates
      !real qcoordc(LDIM,LELT+1)
      ! approximate element vertex coordinates
      real qcoordv(LDIM,AMR_NVRT,LELT+1)

      ! functions
      integer amr_glb_iagti

!#define DEBUG
#ifdef DEBUG
      character*3 str1, str2
      integer iunit
      ! call number
      integer icalld
      save icalld
      data icalld /0/
#endif
!-----------------------------------------------------------------------
      call amr_log(AMR_LP_PRD,'Get vertex numberring.')

      ! get new p4est element offset
      call izero(lelm_goff,NP+1)
      itmpv1(1) = AMR_NELIT
      imsg(1) = amr_glb_iagti(itmpv1(1),1,lelm_goff,NP)

      ! create vertex numbering for partitioning
      ! This version lets partition library to build the graph
      ! vertices only; interior neglected
#ifdef P4EST
      itmp = 1
      call fp4est_lnodes_new(itmp)

      ! get vertex numberring back to nek
      call fp4est_msh_get_lnode(lnelt, itmp, lnodes)

      ! destroy p4est lnode
      call fp4est_lnodes_del()

      ! check calls consistency
      if (itmp.ne.AMR_NVRT) call amr_abort
     $          ('Error: amr_vert_map degree inconsistent')

      if (lnelt.ne.AMR_NELT) call amr_abort
     $          ('Error: amr_vert_map lnelt inconsistent')

      ! get aproximate element centre coordinates from p4est
      !call fp4est_crd_cnt_get(qcoordc)

      ! get aproximate element vertices coordinates from p4est
      call fp4est_crd_vrt_get(qcoordv)
#endif

#ifdef DEBUG
      ! for testing
      icalld = icalld+1
      call io_file_freeid(iunit, ierr)
      write(str1,'(i3.3)') NID
      write(str2,'(i3.3)') icalld
      open(unit=iunit,file='vertex_map.txt'//str1//'i'//str2)
      write(iunit,*) lnelt
      write(iunit,*) 'Vertex global numberring'
      do el=1,lnelt
         write(iunit,*) el,AMR_NELIT+el,(lnodes(il,el),il=1,AMR_NVRT)
      enddo
      close(iunit)
      open(unit=iunit,file='vertex_map_cnt.txt'//str1//'i'//str2)
      do el=1,lnelt
         write(iunit,*) el,(qcoordc(il,el),il=1,LDIM)
      enddo
      close(iunit)
      open(unit=iunit,file='vertex_map_vrt.txt'//str1//'i'//str2)
      do el=1,lnelt
         do il = 1, AMR_NVRT
           write(iunit,*) el,il,(qcoordv(jl,il,el),jl=1,LDIM)
         enddo
      enddo
      close(iunit)
#endif

      ! fill in global ellement number
      do el=1,lnelt
         eid8(el) = AMR_NELIT+el
      enddo

      call fpartMesh(eid8,lnodes,qcoordv,lelt,lnelt,AMR_NVRT,nekcomm,
     $ meshPartitioner,ierr)
      call amr_chk_abort
     $        (ierr,"Error: amr_get_vert_map partMesh fluid failed")

      NELV = lnelt
      NELT = lnelt

      ! check array sizes
      ierr = 0
      if (nelv .gt. lelv) ierr = 1
      call amr_chk_abort
     $        (ierr,"Error: amr_get_vert_map nelv > lelv")

      ! sort elements
      do el = 1,nelv
         LGLEL(el) = eid8(el)
      enddo
      call isort(LGLEL,iwork,nelv)

      do el = 1,nelv
         call icopy84(vertex(1,el),lnodes(1,iwork(el)),AMR_NVRT)
      enddo

      ! for now velocity mesh only
      if (nelgt.ne.nelgv) call amr_abort
     $       ('Error: amr_get_vert_map; solid elements not supported')

      ! fill in element mapping arrays
      do el = 1,nelt
         itmp = LGLEL(el)
         if (itmp.lt.1 .or. itmp.gt.nelgt)
     $      call amr_chk_abort
     $        (ierr,"Error: amr_get_vert_map invalid global el num.")

         ibuf(1) = el
         ibuf(2) = nid
         call dProcmapPut(ibuf,2,0,itmp)
      enddo

      call icopy48(lnodes,vertex,nelt*AMR_NVRT)
      call printPartStat(lnodes,nelt,AMR_NVRT,nekcomm)

      ! required by AMR
      ! finalise communication
      call amr_wait(imsg(1))
      ! fill last entrance in offset array
      lelm_goff(NP) = AMR_NELGT

      ! mapping of nek5000 to p4est element distribution
      call izero(AMR_NP_NID,LELT)
      do el=1,NELT
        do il=0,NP-1
           if (lelm_goff(il+1).ge.LGLEL(el)) then
               AMR_NP_NID(el) = il
               exit
           endif
        enddo
      enddo

      ! transfer partitioning data back to mesh manager
      do el=1,NELT
        gnelt(el) = LGLEL(el)
        elmrk(el) = el
        eldst(el) = AMR_NP_NID(el)
      enddo
      call amr_int_transfer(lnelt,gnelt,elmrk,eldst,LELT)
      ierr = 0
      call izero(AMR_PN_NID,LELT)
      do el=1,AMR_NELT
        ! sanity check
        if (gnelt(el).ne.(AMR_NELIT+el)) ierr = 1
        AMR_PN_NID(el) = eldst(el)
      enddo
      call amr_chk_abort
     $ (ierr,"Error: amr_get_vert_map; global element number mismatch")

#ifdef P4EST
      ! fill in p4est partitioning data
      call fp4est_egmap_put(gnelt,elmrk,eldst)
#endif

#undef DEBUG
      return
      end subroutine
!=======================================================================
!> @brief Element-processor mapping based on p4est graph
!! @details This routine calculate element-processor mapping storred in
!!  gllnid (PARALLEL). To avoid strong scaling limit for ParMETIS I perform
!!  partitioning in 2 steps. First global partitioning between nodes and
!!  next partitioning within the node.
!! @todo Add edges and vertices to the graph.
!! @todo Add graph partitioning for V and T meshes.
      subroutine amr_graph_map()
      implicit none

      include 'SIZE'
      include 'PARALLEL'
      include 'AMR'
      include 'AMR_SHM'

      ! local variables
      integer imsg(1)             ! message request
      integer itmpv1(1)           ! dummy array
      integer itmp1, itmp2        ! dummy variables
      real    rtmp                ! dummy variable
      integer ierr                ! error mark
      integer il, jl, kl, ll      ! loop index

      integer el_countl   ! vertex number in local graph (p4est element number)
      ! size of graph related arrays. I cannot get them exactly, as edge
      ! and vertex connectivity is case dependent. The better way would
      ! be dynamical allocation of array, but for now I stick to f77
      ! static allocation. In this case there is no guarantee the array
      ! is big enought. It can be guaranteed for face connectivity only.
      ! However, this should be enough in most cases
      ! For now only face connectivity is taken into account.
      integer lgraph, lgoff
      parameter (lgraph=LELT*(max(4,lx1))**LDIM)
      parameter (lgoff = LELT+1) ! ??????? AMR_NVRT*LELT+1
      integer graph(lgraph)        ! graph
      integer graph_offset(lgoff)  ! graph vertex offset
      integer edg_wgt(lgraph)      ! edge weights
      integer vrt_wgt(lgoff)       ! vertex weights
      integer elprt(lgoff)         ! element partitioning
      real qcoord(LDIM*lgoff)      ! approximate vertex coordinates
      real vrt_frac(AMR_NDS_MAX*AMR_SHM_MAX)            ! fraction of vertex weight that should be distributed to each sub-domain
      integer work(lgraph)         ! work array

      ! parameter to scale weights of graph nodes end edges
      ! this is ad-hock number without speciffic investigation
      ! some more general way of edge weighting should be introduced
      integer weight_scale
      parameter (weight_scale = 3)

      ! vertex imbalance tolerance
      real tol_imb

      ! simple timing
      real t1, t2

      ! functions
      real dnekclock
      integer amr_glb_iagti

!#define DEBUG
#ifdef DEBUG
      ! for testing
      character*2 str1, str2
      integer iunit
#endif
!-----------------------------------------------------------------------
      write(*,*) 'THIS IS NOT DONE YET amr_graph_map'
      call exitt0
#ifdef NEKPMETIS
      call amr_log(AMR_LP_PRD,'Runing ParMETIS.')

      ! for now we run it for hydro only
      if (NELGV.eq.NELGT) then
         ! get new p4est element offset
         call izero(AMR_ELM_GOFF,NP+1)
         itmpv1(1) = AMR_NELIT
         imsg(1) = amr_glb_iagti(itmpv1(1),1,AMR_ELM_GOFF,NP)

         ! reset element to process mapping
         call izero(elprt,lgoff)

         ! originally il = nx1, but for bigh element count this
         ! requires parmetis to be compiled with 64bit integer
         il = weight_scale

         ! get graph and coordintes from p4est
         call fp4est_msh_get_graph(el_countl, graph, graph_offset,
     $                             vrt_wgt, edg_wgt,il)
#ifdef AMR_PRT_CRD
         call fp4est_crd_get(qcoord)
#endif

         ! scale node work; it is necessary as weights on nodes and edges
         ! should correspond to amout of time spent at different operations,
         ! and right now they are only related to number of grid points
         ! in the element volume, faces ....
         ! this is just a hack not a proper solution
         do il=1,el_countl
          vrt_wgt(il) = weight_scale*vrt_wgt(il)
         enddo

         ! check graph consistency
         ierr = 0
         if (el_countl.ne.AMR_nelt) ierr = 1
         call amr_chk_abort(ierr,
     $       'ERROR; amr_graph_map: inconsistent node number')

         ! finalise communication
         call amr_wait(imsg(1))

         ! fill last entrance in offset array
         AMR_ELM_GOFF(NP) = AMR_nelgt

#ifdef DEBUG
         ! testing; graph with p4est ordering
         call io_file_freeid(iunit, ierr)
         write(str1,'(i2.2)') NID
         write(str2,'(i2.2)') AMR_IMAP
         open(unit=iunit,file='graph_np4.txt'//str2//'i'//str1)
         write(iunit,*) 'OFFSET'
         write(iunit,*) (AMR_ELM_GOFF(il),il=0,NP)
         write(iunit,*) 'GRAPH; C numbering'
         do il=1,AMR_ELM_GOFF(nid+1)-AMR_ELM_GOFF(NID)
            write(iunit,*) il+AMR_ELM_GOFF(NID)-1,
     $      graph_offset(il+1)-graph_offset(il),
     $      (graph(jl),jl=graph_offset(il)+1,graph_offset(il+1))
         enddo
#ifdef AMR_PRT_CRD
         write(iunit,*) 'COORDINATES'
         do il=1,AMR_ELM_GOFF(nid+1)-AMR_ELM_GOFF(NID)
            write(iunit,*) il+AMR_ELM_GOFF(NID)-1,
     &      (qcoord((il-1)*ldim+jl),jl=1,ldim)
         enddo
#endif
         write(iunit,*) 'NODE WEIGHT'
         do il=1,AMR_ELM_GOFF(nid+1)-AMR_ELM_GOFF(NID)
            write(iunit,*) il+AMR_ELM_GOFF(NID)-1,vrt_wgt(il)
         enddo
         write(iunit,*) 'EDGE WEIGHT'
         do il=1,AMR_ELM_GOFF(nid+1)-AMR_ELM_GOFF(NID)
            write(iunit,*) il+AMR_ELM_GOFF(NID)-1,
     &      graph_offset(il+1)-graph_offset(il),
     $      (edg_wgt(jl),jl=graph_offset(il)+1,graph_offset(il+1))
         enddo
         close(iunit)
#endif

         ! get partitioning
         if (mod(AMR_IMAP,AMR_IMAP_MOD).eq.0) then     ! check mapping count
            call amr_log(AMR_LP_PRD,
     $           'Initial partitioning: bisection')
            ! simple timing
            t1 = dnekclock()
            ! I use sub-domain processor coupling related to p4est
            ! partitioning so no initialisation of elprt is necessary
            if (AMR_IF2LEV) then ! 2-level partitioning
               il = AMR_NNODE
               rtmp = 1.0/real(NP)
               do jl = 1,AMR_NNODE
                  vrt_frac(jl) = (AMR_NODE_GOFF(jl+1) -
     $             AMR_NODE_GOFF(jl))*rtmp
               enddo
            else ! single step partitioning on all the nodes
               il = NP
               rtmp = 1.0/real(NP)
               do jl = 1,NP
                  vrt_frac(jl) = rtmp
               enddo
            endif

            ! imbalance tolerance
            tol_imb = 1.05

            call fpmetis_part(AMR_GLBCOMM,AMR_ELM_GOFF,graph_offset,
     $           graph,il,elprt,ndim,qcoord,vrt_wgt,edg_wgt,vrt_frac,
     $           tol_imb)
            t2 = dnekclock()
            AMR_TCPP = AMR_TCPP + t2 - t1
            AMR_IMAP_F = AMR_IMAP_F + 1

#ifdef DEBUG
            if (AMR_IFNODEM) then
               write(*,*) 'Global time: ',nid,t2-t1
            endif
#endif

            ! 2-level partitioning
            if (AMR_IF2LEV) call amr_map_2lev(graph,graph_offset,
     $            edg_wgt,vrt_wgt,qcoord,elprt,work,lgraph,lgoff)

         else ! update of previous partitioning
#if 0
            ! use last optimal mapping and update it
            ! fill old mapping of new elements
            ! reset refinement/coarsening counters
            il=1
            jl=1
            ! loop over local p4est elements
            do kl =1,AMR_NELT
               ! global element number
               itmp1 = AMR_NELIT+kl
               ! check what happend to the element
               ! if no refinement/coarsening
               itmp2 = AMR_GLGL_MAP(kl)
               ! otherwise
               if (itmp2.eq.0) then
                  if (il.le.AMR_RFN_NR.and.
     $               AMR_GLGL_RFN(1,il).eq.itmp1) then
                     ! refinement
                     itmp2 = AMR_GLGL_RFN(2,il)
                     il = il+1
                  elseif (jl.le.AMR_CRS_NR.and.
     $               AMR_GLGL_CRS(1,1,jl).eq.itmp1) then
                     ! coarsening
                     itmp2 = AMR_GLGL_CRS(2,1,jl)
                     jl = jl+1
                  else
                     ! error; some option had to be taken
                     call amr_abort
     $                    ('Error: remap; wrong pointer')
                  endif
               endif
               ! process id
               elprt(kl) = AMR_GLLNID_O(itmp2)

               ! 2-level partitioning
               if (AMR_IF2LEV) then
                  ! find the node
                  itmp1 = -1
                  do ll=1,AMR_NNODE
                     if (elprt(kl).ge.AMR_NODE_GOFF(ll).and.
     $                   elprt(kl).lt.AMR_NODE_GOFF(ll+1)) then
                        itmp1 = ll-1
                        exit
                     endif
                  enddo
                  elprt(kl) = itmp1
               endif
            enddo

            call amr_log(AMR_LP_PRD,
     $           'Update partitioning: adaptive')
            ! simple timing
            t1 = dnekclock()
            if (AMR_IF2LEV) then ! 2-level partitioning
               il = AMR_NNODE
               rtmp = 1.0/real(NP)
               do jl = 1,AMR_NNODE
                  vrt_frac(jl) = (AMR_NODE_GOFF(AMR_NODE_NR+1) -
     $             AMR_NODE_GOFF(AMR_NODE_NR))*rtmp
               enddo
            else ! single step partitioning on all the nodes
               il = NP
               rtmp = 1.0/real(NP)
               do jl = 1,NP
                  vrt_frac(jl) = rtmp
               enddo
            endif
            ! ParMetis itr parameter describing the ratio of inter-processor
            ! communication time compared to data redistribution time.
            t2 = 1.0E+04
            ! imbalance tolerance
            tol_imb = 1.01
            call fpmetis_rpart(AMR_GLBCOMM,AMR_ELM_GOFF,graph_offset,
     $           graph,il,elprt,t2,vrt_wgt,edg_wgt,vrt_frac,tol_imb)
            t2 = dnekclock()
            AMR_TCPPR = AMR_TCPPR + t2 - t1
            AMR_IMAP_R = AMR_IMAP_R + 1

            if (AMR_IFNODEM) then
               write(*,*) 'Global time: ',nid,t2-t1
            endif

            ! 2-level partitioning
            if (AMR_IF2LEV) call amr_map_2lev(graph,graph_offset,
     $            edg_wgt,vrt_wgt,qcoord,elprt,work,lgraph,lgoff)
#endif
         endif ! AMR_IMAP.eq.0

#ifdef DEBUG
         ! testing
         il = AMR_ELM_GOFF(nid+1)-AMR_ELM_GOFF(NID)
         write(str2,'(i2.2)') AMR_IMAP
         !call fp4est_vtk_iscalar(elprt,il,"metis_part"//str2//CHAR(0))
         call io_file_freeid(iunit, ierr)
         write(str1,'(i2.2)') NID
         open(unit=iunit,file='graph_part.txt'//str2//'i'//str1)
         write(iunit,*) (AMR_ELM_GOFF(il),il=0,NP)
#ifdef AMR_PRT_CRD
         do il=1,AMR_ELM_GOFF(nid+1)-AMR_ELM_GOFF(NID)
            write(iunit,*) il+AMR_ELM_GOFF(NID)-1,elprt(il),
     $           (qcoord(LDIM*(il-1)+jl),jl=1,LDIM)
         enddo
#else
         do il=1,AMR_ELM_GOFF(nid+1)-AMR_ELM_GOFF(NID)
            write(iunit,*) il+AMR_ELM_GOFF(NID)-1,elprt(il)
         enddo
#endif
         close(iunit)
#endif

         ! fill in GLLNID
         call izero(GLLNID,NELGT)
         do il=1,AMR_nelt
            GLLNID(AMR_NELIT+il) = elprt(il)
         enddo

         ! global sum to exchange data
         ! GLLEL is used as scratch
         call igop(GLLNID,GLLEL,'+  ',NELGT)

         ! mapping of nek5000 to p4est element distribution
         call izero(AMR_NP_NID,el_countl)
         il = 0
         do kl=1,NELGT
            if (GLLNID(kl).eq.NID) then
               il = il+1
               loop : do jl=0,NP-1
                   if (AMR_ELM_GOFF(jl+1).ge.kl) then
                       AMR_NP_NID(il) = jl
                       exit loop
                   endif
               enddo loop
            endif
         enddo
      else
         call amr_abort('Error: ParMETIS requires NELGT=NELGV')
      endif

#else
      call amr_abort('Error: nekamr requires ParMETIS')
#endif

      ! Count number of elements on this processor
      NELT=0
      NELV=0
      do il=1,NELGT
         if (GLLNID(il).eq.NID) then
            if (il.le.NELGV) NELV=NELV+1
            if (il.le.NELGT) NELT=NELT+1
         endif
      enddo

      ! update mapping count
      AMR_IMAP = AMR_IMAP + 1
#undef DEBUG
      return
      end
!=======================================================================
!> @brief Prepare and call itranode partitioning
!! @details This routine transfer data between cores according to node
!!    partitioning and calls routine for graph regeneration and partitioning.
!! @param[inout] graph           graph
!! @param[inout] graph_offset    graph offset
!! @param[inout] edg_wgt         edge weight
!! @param[inout] vrt_wgt         vertex weight
!! @param[inout] qcoord          node coordinates
!! @param[inout] elprt           element partitioning
!! @param[inout] work            work array
!! @param[in]    lgrph           graph array size
!! @param[in]    lgoff           graph offset array size
      subroutine amr_map_2lev(graph,graph_offset,edg_wgt,vrt_wgt,
     $            qcoord,elprt,work,lgrph,lgoff)
      implicit none

      include 'SIZE'
      include 'PARALLEL'
      include 'AMR'
      include 'AMR_SHM'

      ! argument list
      integer lgrph, lgoff
      integer graph(lgrph), graph_offset(lgoff), edg_wgt(lgrph)
      integer vrt_wgt(lgoff), elprt(lgoff)
      real qcoord(LDIM*lgoff)
      integer work(lgrph)

      ! local variables
      integer imsg(1)             ! message request
      integer itmpv1(1)           ! dummy array
      integer itmpv2(AMR_NDS_MAX) ! dummy array
      integer itmpv3(AMR_NDS_MAX) ! dummy array
      integer itmp1, itmp2        ! dummy variables
      integer ierr                ! error mark
      integer il, jl              ! loop index

      integer lgnel
      parameter (lgnel = LELT) ! ???? AMR_NVRT*LELT
      integer gnel(lgnel)          ! list of global element numbers
      integer tmp_elm_goff(0:AMR_SHM_MAX)! dummy array for intranode element offset

      ! functions
      integer ivlmin, ivlmax
      integer amr_shm_iagti, amr_glb_iscani,
     $        amr_glb_ibcsi

!#define DEBUG
#ifdef DEBUG
      ! for testing
      character*2 str1, str2
      integer iunit
#endif
!-----------------------------------------------------------------------
#if 0

#ifdef DEBUG
      ! testing; element distribution before transfer
      write(str2,'(i2.2)') AMR_IMAP
      write(str1,'(i2.2)') NID
      open(unit=iunit,file='graph_ipart.txt'//str2//'i'//str1)
      write(iunit,*) nid, AMR_SHMID, AMR_NODE_NR - 1
      write(iunit,*) (AMR_ELM_GOFF(il),il=0,NP)
#ifdef AMR_PRT_CRD
      do il=1,AMR_ELM_GOFF(nid+1)-AMR_ELM_GOFF(NID)
         write(iunit,*) il+AMR_ELM_GOFF(NID)-1, elprt(il),
     $        (qcoord(LDIM*(il-1)+jl),jl=1,LDIM)
      enddo
#else
      do il=1,AMR_ELM_GOFF(nid+1)-AMR_ELM_GOFF(NID)
         write(iunit,*) il+AMR_ELM_GOFF(NID)-1, elprt(il)
      enddo
#endif
      close(iunit)
#endif

      ! get approximate intranode partitioning
      ! count number of elements per node
      call izero(itmpv2,AMR_NNODE)
      do il=1, AMR_nelt
         itmpv2(elprt(il)+1) = itmpv2(elprt(il)+1) + 1
      enddo

      imsg(1) = amr_glb_iscani(itmpv2,itmpv3,AMR_NNODE,'+  ')
      call amr_wait(imsg(1))

      ! broadcast element-to-node partitioning
      do il=1,AMR_NNODE
         itmp1 = itmpv2(il)
         itmpv2(il) = itmpv3(il)
         itmpv3(il) = itmpv3(il) - itmp1
      enddo
      imsg(1) = amr_glb_ibcsi(itmpv2,AMR_NNODE,NP-1)
      call amr_wait(imsg(1))

      itmp1 = ivlmin(itmpv2,AMR_NNODE)
      itmp2 = ivlmax(itmpv2,AMR_NNODE)

      if (NIO.eq.0) write(*,*) 'Node element imbalance :',
     $       itmp2-itmp1,itmp1,itmp2

      ! correct element-to-process mapping
      ierr = 0
      do il=1, AMR_NNODE
         itmp1 = itmpv2(il)/(AMR_NODE_GOFF(il+1)-AMR_NODE_GOFF(il))
         if (itmp1.eq.0) then
            ierr = 0
            exit
         endif
         if (mod(itmpv2(il),itmp1).gt.0) itmp1 = itmp1+1
         itmpv2(il) = itmp1
      enddo
      call amr_chk_abort(ierr,
     $              'ERROR; amr_graph_map: too few elements')

      do il=1, AMR_nelt
         itmp1 = elprt(il)+1
         elprt(il) = AMR_NODE_GOFF(itmp1) + itmpv3(itmp1)/itmpv2(itmp1)
         itmpv3(itmp1) = itmpv3(itmp1) + 1
      enddo

      ! transfer data
      itmp1 = AMR_NELT
      call amr_graph_transfer(itmp1,gnel,graph,
     $     graph_offset,edg_wgt,elprt,lgnel,lgrph,lgoff)
      itmp2 = AMR_NELT
      call amr_grnd_transfer(itmp2,vrt_wgt,qcoord,elprt,lgoff)

      ! get new element offset
      call izero(tmp_elm_goff,AMR_SHMNP+1)
      itmpv1(1) = itmp1
      imsg(1) = amr_shm_iagti(itmpv1(1),1,tmp_elm_goff(1),
     $          AMR_SHMNP)

      ! finalise communication
      call amr_wait(imsg(1))

      do il =1,AMR_SHMNP
         tmp_elm_goff(il) = tmp_elm_goff(il-1) + tmp_elm_goff(il)
      enddo

#ifdef DEBUG
      ! testing; graph after transfer
      call io_file_freeid(iunit, ierr)
      write(str1,'(i2.2)') NID
      write(str2,'(i2.2)') AMR_IMAP
      open(unit=iunit,file='graph_trns.txt'//str2//'i'//str1)
      write(iunit,*) 'OFFSET'
      write(iunit,*) (tmp_elm_goff(il),il=0,AMR_SHMNP)
      write(iunit,*) 'GRAPH; C numbering'
      do il=1,tmp_elm_goff(AMR_SHMID+1)-tmp_elm_goff(AMR_SHMID)
         write(iunit,*) gnel(il)-1, graph_offset(il+1)-graph_offset(il),
     $      (graph(jl),jl=graph_offset(il)+1,graph_offset(il+1))
      enddo
#ifdef AMR_PRT_CRD
      write(iunit,*) 'COORDINATES'
      do il=1,tmp_elm_goff(AMR_SHMID+1)-tmp_elm_goff(AMR_SHMID)
         write(iunit,*) gnel(il)-1, (qcoord((il-1)*ldim+jl),jl=1,ldim)
      enddo
#endif
      write(iunit,*) 'NODE WEIGHT'
      do il=1,tmp_elm_goff(AMR_SHMID+1)-tmp_elm_goff(AMR_SHMID)
         write(iunit,*) gnel(il)-1,vrt_wgt(il)
      enddo
      write(iunit,*) 'EDGE WEIGHT'
      do il=1,tmp_elm_goff(AMR_SHMID+1)-tmp_elm_goff(AMR_SHMID)
         write(iunit,*) gnel(il)-1, graph_offset(il+1)-graph_offset(il),
     $        (edg_wgt(jl),jl=graph_offset(il)+1,graph_offset(il+1))
      enddo
      close(iunit)
#endif

      ! adjust graph
      call amr_map_shm(tmp_elm_goff,gnel,graph_offset,graph,
     $     qcoord,vrt_wgt,edg_wgt,elprt,work,lgnel,lgrph,lgoff)

      ! possible place for improving partitioning;
      ! However in this case intranode partitioning should be
      ! performed on copy of graph variables. Right now I do not
      ! copy them, so following lines are commented
!      il = NP
!      rtmp = 1.0/real(NP)
!      do jl = 1,NP
!         vrt_frac(jl) = rtmp
!      enddo
!      tol_imb = 1.001
!      call fpmetis_refine(nekcomm,AMR_ELM_GOFF,graph_offset,
!     $     graph,il,elprt,vrt_wgt,edg_wgt,vrt_frac,tol_imb)

#undef DEBUG

#endif
      return
      end
!=======================================================================
!> @brief Get intranode partitioning
!! @details Thkis routine performs renembering of graph nodes and resets
!!    the graph removing internode edges and renumbering intranode ones.
!!    Finally it performa intranode partitioning.
!! @param[in]     elm_goff     global element offset
!! @param[inout]  gnel         global element numbers
!! @param[inout]  graph_offset graph vertex offset
!! @param[inout]  graph        graph
!! @param[inout]  qcoord       vertex coordinates
!! @param[inout]  vrt_wgt      vertex weights
!! @param[inout]  edg_wgt      edge weights
!! @param[inout]  elprt        element partitioning
!! @param[inout]  work         work array
!! @param[in]     lgnel        gnel array size
!! @param[in]     lgrph        graph array size
!! @param[in]     lgoff        graph offset array size
      subroutine amr_map_shm(elm_goff,gnel,graph_offset,graph,
     $            qcoord,vrt_wgt,edg_wgt,elprt,work,lgnel,lgrph,lgoff)
      implicit none

      include 'SIZE'
      include 'PARALLEL'
      include 'AMR'
      include 'AMR_SHM'

      ! argument list
      integer lgnel, lgrph, lgoff
      integer elm_goff(0:AMR_SHM_MAX), gnel(lgnel), graph_offset(lgoff),
     $     graph(lgrph), vrt_wgt(lgoff), edg_wgt(lgrph), elprt(lgoff)
      real qcoord(LDIM*lgoff)
      integer work(lgrph)

      ! local variables
      integer imsg(1)             ! message request
      integer itmpv1(LELT*AMR_SHM_MAX) ! dummy array
      integer itmpv2(LELT*AMR_SHM_MAX) ! dummy array
      integer itmp1, itmp2, itmp3 ! dummy variables
      real    rtmp                ! dummy variable
      integer ierr                ! error mark
      integer il, jl, kl          ! loop index

      ! vertex imbalance tolerance
      real vrt_frac(AMR_SHM_MAX) ! fraction of vertex weight that should be distributed to each sub-domain
      real tol_imb

      ! functions
      integer amr_shm_iagtvi, amr_shm_igopi
      real dnekclock

      ! simple timing
      real t1, t2

!#define DEBUG
#ifdef DEBUG
      ! for testing
      character*2 str1, str2
      integer iunit
#endif
!-----------------------------------------------------------------------
#if 0
      ! renumber elements on a single node
      ! gather intranode global element numbers
      itmp1 = AMR_NODE_GOFF(AMR_NODE_NR) ! master node id
      itmp2 = elm_goff(AMR_SHMID+1)- elm_goff(AMR_SHMID) ! core element number
      itmp3 = elm_goff(AMR_SHMNP) ! node element number
      ! get count array
      do il= 1, AMR_SHMNP
         itmpv2(il) = elm_goff(il) - elm_goff(il-1)
      enddo
      imsg(1) = amr_shm_iagtvi(gnel,itmp2,itmpv1,itmp3,
     $          itmpv2,elm_goff)
      call amr_wait(imsg(1))
      ! sort intranode element numbers acording to global number
      call isort(itmpv1,itmpv2,itmp3)
      ! turn to C to fortran numberring
      do il=1,graph_offset(itmp2+1)
          graph(il) = graph(il) + 1
      enddo

      ! adjust local graph renumbering nodes and remove internode edges
      ! sort graph acording to edge global number
      call isort(graph,work,graph_offset(itmp2+1))
      ! compare two sorted lists (intranode elements and sorted graph)
      ! renumbering nodes and marking internode edges
      il = 1
      jl = 1
      ierr = 0
      do
         if (il.gt.graph_offset(itmp2+1)) then
            exit
         else if (jl.gt.itmp3) then
            graph(il) = - graph(il)
            il = il + 1
         else if(graph(il).eq.itmpv1(jl)) then
            graph(il) = itmpv2(jl)
            il = il + 1
         else if(graph(il).gt.itmpv1(jl)) then
            jl = jl + 1
         else if (graph(il).lt.itmpv1(jl)) then
            graph(il) = - graph(il)
            il = il + 1
         else
            ierr = 1
         endif
      enddo

      call amr_chk_abort(ierr,
     $           'ERROR; amr_graph_adj: wrong renumber operation')

      ! permute graph back
      call iswapt_ip(graph,work,graph_offset(itmp2+1))

      ! calculate number of edges for vertex
      do il = 1, itmp2
         work(il) = graph_offset(il+1)-graph_offset(il)
      enddo
      ! remove internode edges
      il = 1
      do jl = 1, itmp2 ! loop over vertices
         do kl = graph_offset(jl)+1,graph_offset(jl+1) ! loop over edges
            if (graph(kl).lt.0) then
               work(jl) = work(jl) - 1
            else
               if (il.ne.kl) then
                  graph(il) = graph(kl)
                  edg_wgt(il) = edg_wgt(kl)
               endif
               il = il + 1
            endif
         enddo
      enddo

      ! check if there are disconnected elements
      ierr = 0
      do il = 1, itmp2
         if (work(il).lt.1) ierr=1
      enddo

      call amr_chk_abort(ierr,
     $           'ERROR; amr_graph_adj: disconnected element')

      ! recalculate offset
      graph_offset(1) = 0
      do il = 2, itmp2+1
         graph_offset(il) = graph_offset(il-1) + work(il-1)
      enddo

      ! turn to fortran to C numberring
      do il=1,graph_offset(itmp2+1)
          graph(il) = graph(il) - 1
      enddo

#ifdef DEBUG
      ! testing; graph after adjusting
      call io_file_freeid(iunit, ierr)
      write(str1,'(i2.2)') NID
      write(str2,'(i2.2)') AMR_IMAP
      open(unit=iunit,file='graph_adj.txt'//str2//'i'//str1)
      write(iunit,*) 'OFFSET'
      write(iunit,*) (elm_goff(il),il=0,AMR_SHMNP)
      write(iunit,*) 'ELEMENT renumbering'
      do il=1,elm_goff(AMR_SHMNP)
         write(iunit,*) il-1,itmpv1(il)-1
      enddo
      write(iunit,*) 'GRAPH; C numbering'
      do il=1,elm_goff(AMR_SHMID+1)-elm_goff(AMR_SHMID)
         write(iunit,*) gnel(il)-1,
     $      graph_offset(il+1)-graph_offset(il),
     $      (graph(jl),jl=graph_offset(il)+1,graph_offset(il+1))
      enddo
#ifdef AMR_PRT_CRD
      write(iunit,*) 'COORDINATES'
      do il=1,elm_goff(AMR_SHMID+1)-elm_goff(AMR_SHMID)
         write(iunit,*) gnel(il)-1,
     &        (qcoord((il-1)*ldim+jl),jl=1,ldim)
      enddo
#endif
      write(iunit,*) 'NODE WEIGHT'
      do il=1,elm_goff(AMR_SHMID+1)-elm_goff(AMR_SHMID)
         write(iunit,*) gnel(il)-1,vrt_wgt(il)
      enddo
      write(iunit,*) 'EDGE WEIGHT'
      do il=1,elm_goff(AMR_SHMID+1)-elm_goff(AMR_SHMID)
         write(iunit,*) gnel(il)-1,
     &     graph_offset(il+1)-graph_offset(il),
     $     (edg_wgt(jl),jl=graph_offset(il)+1,graph_offset(il+1))
      enddo
      close(iunit)
#endif

      ! intranode bisection
      il = AMR_SHMNP
      rtmp = 1.0/real(AMR_SHMNP)
      do jl = 1,AMR_SHMNP
         vrt_frac(jl) = rtmp
      enddo
      ! imbalance tolerance
      tol_imb = 1.01
      ! simple timing
      t1 = dnekclock()

      call fpmetis_part(AMR_SHMCOMM,elm_goff,graph_offset,
     $     graph,il,elprt,ndim,qcoord,vrt_wgt,
     $     edg_wgt,vrt_frac,tol_imb)

      ! simple timing
      t2 = dnekclock()
      AMR_TCPP2 = AMR_TCPP2 + t2 - t1
      AMR_IMAP_2 = AMR_IMAP_2 + 1

#ifdef DEBUG
      if (AMR_IFNODEM) then
         write(*,*) 'Local time: ',nid,t2-t1
      endif

      ! count number of intranode elements per core
      call izero(itmpv1,AMR_SHMNP)
      do il=1, itmp2
         itmpv1(elprt(il)+1) = itmpv1(elprt(il)+1) + 1
      enddo

      imsg(1) = amr_shm_igopi(itmpv1,itmpv2,AMR_SHMNP,'+  ')

      call amr_wait(imsg(1))
      if (AMR_IFNODEM) then
         il = itmpv2(1)
         jl = il
         do kl=2,AMR_SHMNP
            il=min(il,itmpv2(kl))
            jl=max(jl,itmpv2(kl))
         enddo
         write(*,*) 'Local inbalance',nid,il,jl,jl-il
      endif
#endif

      ! shift destination process id with master id
      do il=1,itmp2
         elprt(il) = elprt(il) + itmp1
      enddo

#ifdef DEBUG
      ! testing
      call io_file_freeid(iunit, ierr)
      write(str1,'(i2.2)') NID
      write(str2,'(i2.2)') AMR_IMAP
      open(unit=iunit,file='graph_shm_part.txt'//str2//'i'//str1)
      write(iunit,*) 'OFFSET'
      write(iunit,*) (elm_goff(il),il=0,AMR_SHMNP)
#ifdef AMR_PRT_CRD
      write(iunit,*) 'PARTITIONING, COORDINATES'
      do il=1,elm_goff(AMR_SHMID+1)-elm_goff(AMR_SHMID)
         write(iunit,*) gnel(il)-1,elprt(il),
     &        (qcoord((il-1)*ldim+jl),jl=1,ldim)
      enddo
#else
      write(iunit,*) 'PARTITIONING'
      do il=1,elm_goff(AMR_SHMID+1)-elm_goff(AMR_SHMID)
         write(iunit,*) gnel(il)-1, elprt(il)
      enddo
      close(iunit)
#endif
#endif

      ! redistribute element mark back
      il = elm_goff(AMR_SHMID+1)-elm_goff(AMR_SHMID)
      call amr_prt_transfer(il,gnel,elprt,LELT)

#ifdef DEBUG
      ! testing
      call io_file_freeid(iunit, ierr)
      write(str1,'(i2.2)') NID
      write(str2,'(i2.2)') AMR_IMAP
      open(unit=iunit,file='graph_glb_part.txt'//str2//'i'//str1)
      write(iunit,*) 'OFFSET'
      write(iunit,*) (AMR_ELM_GOFF(il),il=0,NP)
      write(iunit,*) 'PARTITIONING'
      do il=1,AMR_ELM_GOFF(nid+1)-AMR_ELM_GOFF(nid)
         write(iunit,*) gnel(il)-1, elprt(il)
      enddo
      close(iunit)
#endif

#undef DEBUG

#endif
      return
      end subroutine
!=======================================================================


