!> @file amr.f
!! @ingroup nekamr
!! @brief Main interface for nekamr
!! @author Adam Peplinski
!! @date Feb 26, 2016
!
! This file requires number of constant defined in include file
#include "nekamr.h"
!=======================================================================
!> @brief Initialisation of sc and p4est libraries
!! @param[in] intracomm mpi communicator
      subroutine amr_init(intracomm)
      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'DOMAIN'
      include 'AMR'

      ! argument list
      integer intracomm

      ! local variables
      ! simple timing
      real t1, tmp

      ! functions
      real dnekclock
!-----------------------------------------------------------------------
      ! simple timing
      t1 = dnekclock()

      ! reset AMR timers
      AMR_TC = 0.0
      AMR_TCI = 0.0
      AMR_TCF = 0.0
      AMR_TCE = 0.0
      AMR_TCER = 0.0
      AMR_TCP4 = 0.0
      AMR_TCP5 = 0.0
      AMR_TCP6 = 0.0
      AMR_TCP7 = 0.0
      AMR_TCP8 = 0.0
      AMR_TCPN = 0.0
      AMR_TCPS = 0.0
      AMR_TCPD = 0.0
      AMR_TCPM = 0.0
      AMR_TCPP = 0.0
      AMR_TCPPR = 0.0
      AMR_TCPP2 = 0.0
      AMR_TCPR = 0.0
      AMR_TCTD = 0.0
      AMR_TCL = 0.0
      AMR_TCC = 0.0
      AMR_TCS = 0.0
      AMR_TCSC = 0.0
      AMR_TCRS = 0.0
      AMR_TCN = 0.0

      ! reset refinement/coarsening count
      AMR_REF_CNT = 0
      AMR_EER_CNT = 0
      AMR_RIN_CNT = 0

      ! reset mapping count
      AMR_IMAP = 0
      AMR_IMAP_F = 0
      AMR_IMAP_R = 0
      AMR_IMAP_2 = 0
      ! frequency of grph bisection (1 means only graph bisections)
      ! for now remapping is turned off as it gives disjoint graps for
      ! 2-level partitioning
      AMR_IMAP_MOD = 1

      ! set initial default log level
      AMR_LP_DEF=AMR_LP_PRD

#ifdef P4EST
      ! sc init
      call fsc_init(intracomm, 1, 0, AMR_LP_DEF)

      ! p4est init
      call fp4est_init(AMR_LP_DEF)
#endif

      call amr_log(AMR_LP_DEF,'Starting nekamr logs.')
 
      ! check consistency with SIZE
      if (N_DIM.ne.LDIM) then
         call amr_abort('Error: nekamr ldim inconsistent')
      endif

      ! set reset flag
      AMR_IFRESET = .FALSE.

      ! take into account inter- and intra-node communication
      ! critical for 2-level partitioning only, but in one place used for
      ! nek native partitioning as well
      call amr_node_split()

      ! for crs interpolation (vertices only)
      ! check consistency with DOMAIN parameters
      if (AMR_LXC.ne.LXC) then
         call amr_abort('ERROR: amr_init; inconsistent lxc')
      endif

      ! simple timing
      tmp = dnekclock() - t1
      ! total
      AMR_TC = AMR_TC + tmp
      ! initialisation
      AMR_TCI = AMR_TCI + tmp
      return
      end subroutine
!=======================================================================
!> @brief Finalisation of sc and p4est in nekton
      subroutine amr_end()
      implicit none

      include 'SIZE'
      include 'PARALLEL'
      include 'INPUT'
      include 'AMR'
      include 'AMR_SHM'

      ! local variables
      ! simple timing
      real t1, tmp

      ! functions
      real dnekclock
!-----------------------------------------------------------------------
      ! save restart files
      call amr_restart_save()

      ! simple timing
      t1 = dnekclock()

#ifdef P4EST
      ! release memory
      ! delete mesh and ghost cells
      call fp4est_mesh_del()
      call fp4est_ghost_del()
      ! delete tree and connectivity
      call fp4est_tree_del()
      call fp4est_cnn_del()

      ! sc end
      call fsc_pkg_print(AMR_LP_DEF)
      call fsc_finalize()
#endif

      ! free inter- and intra-node communicators
      call amr_node_free()

      ! simple timing
      tmp = dnekclock() - t1
      ! total
      AMR_TC = AMR_TC + tmp
      ! finalisation
      AMR_TCF = AMR_TCF + tmp

      if (NID.eq.0) then
         write(6,*) 
         write(6,'(A)') 'Finalized AMR:'
         write(6,'(A,i6,A)')    'Cumulative time over ',AMR_EER_CNT,
     $      ' calls to amr_refinement'
         write(6,'(A,g13.5,A)') 'Total time in AMR: ',AMR_TC  ,' sec'
         write(6,'(A,g13.5,A)') '      initialisation: ',AMR_TCI ,' sec'
         write(6,'(A,g13.5,A)') '        finalisation: ',AMR_TCF ,' sec'
         write(6,'(A,g13.5,A)') ' restart file saving: ',AMR_TCRS,' sec'
         write(6,'(A,g13.5,A)') '           evolution: ',AMR_TCE ,' sec'
         write(6,'(A,g13.5,A,i6,A)') '        error estimator: ',
     $      AMR_TCER,' sec; ',AMR_EER_CNT,' calls'
         write(6,'(A,g13.5,A,i6,A)') '             p4est part: ',
     $      AMR_TCP4,' sec',AMR_REF_CNT,' calls'
         write(6,'(A,g13.5,A)') '                mesh/ghost: ',AMR_TCP5,
     $      ' sec'
         write(6,'(A,g13.5,A)') '           refinement fill: ',AMR_TCP6,
     $      ' sec'
         write(6,'(A,g13.5,A)') '            refine/balance: ',AMR_TCP7,
     $      ' sec'
         write(6,'(A,g13.5,A)') '                 partition: ',AMR_TCP8,
     $      ' sec'
         write(6,*)
         write(6,'(A,g13.5,A,i6,A)') '        regenarate mesh: ',
     $      AMR_TCPN,' sec',AMR_RIN_CNT,' calls'
         write(6,'(A,g13.5,A)') '                 mesh size: ',AMR_TCPS,
     $      ' sec'
         write(6,'(A,g13.5,A)') '        p4est-nek transfer: ',AMR_TCPD,
     $      ' sec'
         write(6,'(A,g13.5,A)') '   p4est-nek topology data: ',AMR_TCTD,
     $      ' sec'
         write(6,'(A,g13.5,A)') '           element mapping: ',AMR_TCPM,
     $      ' sec'
         write(6,'(A,g13.5,A)') '                ParMetis (bisect): ',
     $      AMR_TCPP,' sec'
         write(6,'(A,g13.5,A)') '                         (bisect): ',
     $      AMR_IMAP_F,' calls'
         write(6,'(A,g13.5,A)') '                ParMetis (remap) : ',
     $      AMR_TCPPR,' sec'
         write(6,'(A,g13.5,A)') '                         (remap) : ',
     $      AMR_IMAP_R,' calls'

         if (AMR_IF2LEV) then
         write(6,'(A,g13.5,A)') '                ParMetis (2 lev) : ',
     $      AMR_TCPP2,' sec'
         write(6,'(A,g13.5,A)') '                         (2 lev) : ',
     $      AMR_IMAP_2,' calls'
         endif

         write(6,'(A,g13.5,A)') '            redistribution: ',AMR_TCPR,
     $      ' sec'
         write(6,*)
         write(6,'(A,g13.5,A)') '         ref/crs  local: ',AMR_TCL ,
     $      ' sec'
         write(6,'(A,g13.5,A)') '       transfer/sorting: ',AMR_TCC ,
     $      ' sec'
         write(6,'(A,g13.5,A)') '         solver restart: ',AMR_TCS ,
     $      ' sec'
         write(6,'(A,g13.5,A)') '               crs restart: ',AMR_TCSC,
     $      ' sec'
         write(6,*)
         write(6,'(A)') 'Averaged values:'
         AMR_TCE = AMR_TCE/real(max(1,AMR_EER_CNT))
         write(6,'(A,g13.5,A)') 'Averaged cycle AMR: ',AMR_TCE ,' sec'
         write(6,'(A,g13.5,A)') 'Averaged cycle nek: ',AMR_TCN ,' sec'
         write(6,*)
      endif

      return
      end subroutine
!=======================================================================
!> @brief Load mesh tree information to mesh manager
      subroutine amr_tree_load()
      implicit none

      include 'SIZE'
      include 'PARALLEL'
      include 'INPUT'
      include 'AMR'

      ! common blocks
      integer nidl,npl,nekcomm,nekgroup,nekreal
      common /nekmpi/ nidl,npl,nekcomm,nekgroup,nekreal

      ! local variables
      character(len=AMR_LSTL_FNM) filename ! file name

      ! file name
      character(len=6)  str

      ! simple timing
      real t1, tmp

      ! functions
      integer iglmax, iglsum, amr_glb_iagti
      real dnekclock
!-----------------------------------------------------------------------
      ! simple timinng
      t1 = dnekclock()

      ! read forest data
      call amr_log(AMR_LP_PRD,'Reading mesh tree data.')

      ! prepare file name for ATF (AMR Tree File) file
      filename = 'ATF'//trim(adjustl(SESSION))
      if (len_trim(filename).gt.(AMR_LSTL_FNM-13)) then
         call amr_abort
     $        ('ERROR: amr_tree_load; too long output file name')
      endif
      write(str,'(i5.5)') AMR_IOSTART
      filename = trim(filename)//'.t'//trim(str)//CHAR(0)
#ifdef P4EST
      call fp4est_tree_load(nekcomm, 1,filename)

      ! create ghost zones
      call fp4est_ghost_new()

      ! get new mesh
      call fp4est_mesh_new()

      ! check boundary condition mark
      call fp4est_bc_check()

#ifdef DEBUG
      ! testing
      call fp4est_vtk_write('mesh_test'//char(0))
#endif

      ! get mesh size
      call fp4est_msh_get_size(AMR_NELGT,AMR_NELIT,AMR_NELT,
     $    AMR_NELV,AMR_MLEV)
#endif

      ! get global max element level
      AMR_MLEV = iglmax(AMR_MLEV,1)

      ! check restart consistency
      if (AMR_MLEV.ne.0.and.AMR_IOSTART.eq.0) call amr_abort
     $          ('Error: amr_tree_load; cannot read #.re2 with LEV>0')

      if (AMR_IOSTART.ne.0) then
        ! global element count
        NELGT = AMR_NELGT
        ! get globalnumber of V elements
        NELGV = iglsum(AMR_NELV,1)

        ! check array sizes for p4est element distribution
        call chk_nel
      endif

      ! to be sure chkParam wont exit a code
      ldimr = LDIM

      ! simple timing
      tmp = dnekclock() - t1
      ! total
      AMR_TC = AMR_TC + tmp
      ! initialisation
      AMR_TCI = AMR_TCI + tmp

      return
      end subroutine
!=======================================================================
!> @brief Load tree information from the file
!! @note This routine can be called after amr_init_par_rea or
!!       amr_init_par_par have been executed.
!! @todo Upgrade curvature data; check setdef and setzer.
      subroutine amr_mesh_init()
      implicit none

      include 'SIZE'
      include 'PARALLEL'
      include 'INPUT'
      include 'GEOM'    ! for bc only
      include 'AMR'

      ! local variables
      integer itmp, lcbc, lrbc, ibcl, nfldl, el, il, jl

      ! tmp array for curvature and boundary surface flags
      integer crvl(AMR_NFCS,LELT),bcl(AMR_NFCS,LELT)

      ! simple timing
      real t1, tmp

      ! functions
      integer iglsum
      real dnekclock
!-----------------------------------------------------------------------
      ! simple timinng
      t1 = dnekclock()

      ! read forest data
      call amr_log(AMR_LP_PRD,'Initialising mesh data.')

      ! check if we generate mesh from .rea or .mesh
      if (AMR_IOSTART.eq.0) then
        ! we use nek5000 standard method to generate the mesh
        ! check consistency of p4est structure and .rea file
        ! T mesh
        if (NELGT.ne.AMR_NELGT) then
           call amr_abort
     $          ('Error: amr_mesh_init; NELGT inconsistent')
        endif
        ! V mesh
        ! get globalnumber of V elements
        itmp = iglsum(AMR_NELV,1)
        if (NELGV.ne.itmp) then
           call amr_abort
     $          ('Error: amr_mesh_init; NELGV inconsistent')
        endif

      endif

      ! initialise arrays
      ! some of them are done in initdat but not all. It would be good to
      ! do it consistently but for now I'm leaving it as it is
      call izero(crvl, AMR_NFCS*LELT)
      call izero(bcl, AMR_NFCS*LELT)

#ifdef P4EST
      ! load mesh data to local arrays
      call fp4est_msh_get_dat(NELGV,LELT,IGROUP,AMR_LEVEL,crvl,bcl)
#endif

      ! prenek to nek face reordering
      ! originally set in setup_topo, but I use them, so I do it here
      call initds

      ! redistribute data
      call amr_tree_transfer(IGROUP,AMR_LEVEL,crvl,bcl)

      ! get b.c. position in bc and cbc array
      ibcl = 2
      if (IFFLOW) ibcl = 1
      nfldl = 1
      if (IFHEAT) nfldl = 2+NPSCAL
      if (IFMHD ) nfldl = 2+NPSCAL+1
!     If p32 = 0.1, there will be no bcs read in
      if (PARAM(32).gt.0) nfldl = ibcl + PARAM(32)-1

      ! partially fill in BC and CBC arrays; only string flag for internal and periodic faces
      lcbc=18*LELT*(LDIMT1 + 1)
      lrbc=30*LELT*(LDIMT1 + 1)
      call rzero(BC, lrbc)
      call blank(CBC,lcbc)
      call blank(cbc_bmap,sizeof(cbc_bmap))
      call izero(boundaryID, size(boundaryID))
      call izero(boundaryIDt, size(boundaryIDt))
      do el = 1, NELT
        do il =1, AMR_NFCS
          if (bcl(il,el).eq.0) then
            do jl=ibcl,nfldl
              CBC(il,el,jl) = 'E  '
            enddo
          else if(bcl(il,el).eq.-1) then
            do jl=ibcl,nfldl
              CBC(il,el,jl) = 'P  '
            enddo
          else
            BC(5,il,el,1) = bcl(il,el)
            ! no difference between V and T mesh for now
            boundaryID(il,el) =  bcl(il,el)
            boundaryIDt(il,el) =  bcl(il,el)
          endif
          BC(4,il,el,1) = crvl(il,el)
        enddo
      enddo

      ! get mesh topology
      call amr_topol_get()

      ! simple timing
      tmp = dnekclock() - t1
      ! total
      AMR_TC = AMR_TC + tmp
      ! initialisation
      AMR_TCI = AMR_TCI + tmp

      return
      end subroutine
!=======================================================================
!> @brief Main interface for mesh refinement/coaresering
      subroutine amr_refinement()
      implicit none

      include 'SIZE'
      include 'PARALLEL'
      include 'TSTEP'
      include 'AMR'
      include 'AMR_REFINE'

      ! local variables
      integer il, jl  ! loop index

      logical ifmesh_ref ! did p4est change the mesh

      ! simple timing
      real t1, t2, tmp, t21, t22
      ! call counter
      integer icalled
      save icalled
      data icalled /0/

      ! functions
      integer iglmin, iglmax, iglsum
      real dnekclock
!-----------------------------------------------------------------------
      ! simple timing
      t1 = dnekclock()
      ! get averaged nek cycle timing
      AMR_T2N = t1
      if (icalled.gt.0)
     $ AMR_TCN = (AMR_TCN*(icalled-1)+AMR_T2N-AMR_T1N)/real(icalled)
      icalled = icalled + 1

      ! get refinement/coarsering mark
      call amr_refinement_mark()
      ! count number of calls
      AMR_EER_CNT = AMR_EER_CNT + 1
      ! simple timing
      AMR_TCER = AMR_TCER + dnekclock() - t1

      ! check if there is any refinement/coarsening mark
      il = iglmin(AMR_MARK,AMR_NELT)
      jl = iglmax(AMR_MARK,AMR_NELT)

      if (il.eq.0.and.jl.eq.0) then
        call amr_log(AMR_LP_PRD,'No refine mark.')
        ! simple timing
        t2 = dnekclock()
        AMR_TC = AMR_TC + t2 - t1
        ! get averaged nek cycle timing
        AMR_T1N = dnekclock()
        return
      endif

      ! stamp logs
      call amr_log(AMR_LP_PRD,'Global refinement; BEGIN.')

      ! save restart files if required
      if (AMR_IOSAVE) call amr_restart_save()

      ! p4est refinement count
      AMR_REF_CNT = AMR_REF_CNT + 1

      ! simple timing
      t21 = dnekclock()
      ! perform operations on p4est side
      call amr_p4est_refine(ifmesh_ref)
      ! simple timing
      t22 = dnekclock()
      AMR_TCP4 = AMR_TCP4 + t22 - t21

      ! if the mesh was modified
      if (ifmesh_ref) then
         ! solver restart count
         AMR_RIN_CNT = AMR_RIN_CNT + 1
         ! simple timing
         t21 = dnekclock()

         ! before refining make sure BM1 is removed from rhs
         call amr_rhs_rmv_bm


         ! regenerate mesh info on nek side
         call amr_mesh_regenerate
         ! simple timing
         t22 = dnekclock()
         AMR_TCPN = AMR_TCPN + t22 - t21

         ! refine, transfer and coarsen data
         call amr_refine_coarsen

         ! set reset flag
         AMR_IFRESET = .TRUE.

         ! simple timing
         t21 = dnekclock()
         ! reinitialise solver
         call amr_nek_reinit

         ! finally make sure rhs is multiplied by BM1
         call amr_rhs_mlt_bm

         ! simple timing
         t22 = dnekclock()
         AMR_TCS = AMR_TCS + t22 - t21

         ! unset reset flag
         AMR_IFRESET = .FALSE.
      else
         call amr_log(AMR_LP_PRD,'Mesh not changed at ref./crs. step')
      endif

      ! stamp logs
      call amr_log(AMR_LP_PRD,'Global refinement; END.')

      ! simple timing
      t2 = dnekclock()
      tmp = t2 - t1
      ! total
      AMR_TC = AMR_TC + tmp
      ! evolution
      AMR_TCE = AMR_TCE + tmp

      ! get averaged nek cycle timing
      AMR_T1N = dnekclock()

      return
      end subroutine
!=======================================================================
!> @brief Get refinement mark
      subroutine amr_refinement_mark
      implicit none

      include 'SIZE'
      include 'AMR'
!-----------------------------------------------------------------------
      call amr_log(AMR_LP_PRD,'Get refinement mark.')

      ! reset refinement mark
      call izero(AMR_MARK,NELT)

      ! user supplied refinement marker
      call user_ref_mark(AMR_MARK,AMR_LEVEL)

      ! check if coarsening does not cause balance problem
      call amr_mark_check

      ! redistribute refinement mark back to p4est element distribution
      call amr_mark_transfer

      return
      end subroutine
!=======================================================================
!> @brief Perform refinement/coaresering of the mesh on p4est side
!! @param[inout] ifp4est_ref    did p4est modify the mesh
      subroutine amr_p4est_refine(ifp4est_ref)
      implicit none

      include 'SIZE'
      include 'AMR'

      ! argument list
      logical ifp4est_ref

      ! local variables
      ! mesh modification marker
      integer imod
      ! arays for global operator
      integer itmp(1), iwrk(1)

      ! partition strategy for p4est
      ! 0 - equal element count
      ! 1 - octants families prepared for coarsening
      integer p4part
      parameter (p4part = 1)

      ! copy paremeter for forest comparison in p4est
      ! 0 - no check of quadrant data
      ! 1 - check if quadrant data are identical
      integer p4test
      parameter (p4test = 1)

      ! simple timing
      real t1, t2

      ! functions
      real dnekclock
!-----------------------------------------------------------------------
#ifdef P4EST
      ! simple timing
      t1 = dnekclock()
      ! destroy old mesh and ghost structures
      call fp4est_mesh_del()
      call fp4est_ghost_del()
      t2 = dnekclock()
      AMR_TCP5 = AMR_TCP5 + t2 - t1

      ! simple timing
      t1 = dnekclock()
      ! transfer refinement mark to p4est
      call fp4est_refm_put(AMR_MARK)
      t2 = dnekclock()
      AMR_TCP6 = AMR_TCP6 + t2 - t1

      ! simple timing
      t1 = dnekclock()
      ! p4estperform local refine/coarsen/balance
      call fp4est_tree_copy(p4test)
      call fp4est_refine(AMR_LMAX)
      call fp4est_coarsen()
      call fp4est_balance()
      call fp4est_tree_check(imod,p4test)
      t2 = dnekclock()
      AMR_TCP7 = AMR_TCP7 + t2 - t1

      ! get global mesh modification test
      itmp(1) = imod
      call igop(itmp,iwrk,'*  ',1)

      ! was mesh modified
      if (itmp(1).eq.0) then
         ifp4est_ref = .TRUE.
      else
         ifp4est_ref = .FALSE.
      endif

      ! only if the mesh was modified
      if (ifp4est_ref) then
         ! simple timing
         t1 = dnekclock()
         ! partition  p4est mesh
         call fp4est_part(p4part)
         t2 = dnekclock()
         AMR_TCP8 = AMR_TCP8 + t2 - t1
      endif

      ! simple timing
      t1 = dnekclock()
      ! create ghost zones
      call fp4est_ghost_new()
      ! get new mesh
      call fp4est_mesh_new()
      t2 = dnekclock()
      AMR_TCP5 = AMR_TCP5 + t2 - t1
#endif

      return
      end subroutine
!=======================================================================
!> @brief Regenerate mesh information on nek side
      subroutine amr_mesh_regenerate()
      implicit none

      include 'SIZE'
      include 'PARALLEL'
      include 'INPUT'
      include 'GEOM'    ! for bc only
      include 'AMR'

      ! local variables
      integer el, il, jl, lcbc, lrbc, ibcl, nfldl
      ! simple timing
      real t1, t2

      ! tmp array for curvature and boundary surface flags
      integer crvl(AMR_NFCS,LELT),bcl(AMR_NFCS,LELT)

      ! logs
      character(len=AMR_LSTL_LOG) str

      ! functions
      integer iglmin, iglmax, iglsum
      real dnekclock

#ifdef DEBUG
!     for testing
      character*2 str1, str2
      integer iunit, ierr
!     call number
      integer icalld
      save icalld
      data icalld /0/
#endif
!-----------------------------------------------------------------------
      ! save old partition sizes on nek side
      AMR_NELGT_O = NELGT
      AMR_NELT_O = NELT
      AMR_NELV_O = NELV

#ifdef DEBUG
      ! testing
      icalld = icalld + 1
      write(str2,'(i2.2)') icalld
#ifdef P4EST
      call fp4est_vtk_write('mesh_test'//str2//char(0))
#endif
#endif

      ! check boundary condition mark
#ifdef P4EST
      call fp4est_bc_check()
#endif

      ! simple timing
      t1 = dnekclock()
      ! get new mesh size
#ifdef P4EST
      call fp4est_msh_get_size(AMR_NELGT,AMR_NELIT,AMR_NELT,
     $    AMR_NELV,AMR_MLEV)
#endif
      t2 = dnekclock()
      AMR_TCPS = AMR_TCPS + t2 - t1

      ! get global max quadrant level
      AMR_MLEV = iglmax(AMR_MLEV,1)

      ! global element count
      NELGT = AMR_NELGT
      ! get globalnumber of V elements
      NELGV = iglsum(AMR_NELV,1)

      ! local element count related to p4est element distribution
      NELT = AMR_NELT
      NELV = AMR_NELV

      ! stamp log
      write(str, *) 'After refinement'
      call amr_log(AMR_LP_PRD,trim(str))
      write(str,11) AMR_LMAX, AMR_MLEV
 11   format('Global max/current refinement level: ',2I4)
      call amr_log(AMR_LP_PRD,trim(str))

      ! check array sizes for p4est element distribution
      call chk_nel

      ! simple timing
      t1 = dnekclock()
      call izero(crvl, AMR_NFCS*LELT)
      call izero(bcl, AMR_NFCS*LELT)
#ifdef P4EST
      ! load mesh data to local arrays
      call fp4est_msh_get_dat(NELGV,LELT,IGROUP,AMR_LEVEL,crvl,bcl)

      ! reset refinement history
      call fp4est_msh_get_hst(AMR_MAP_NR,AMR_RFN_NR,AMR_CRS_NR,
     $                        AMR_GLGL_MAP,AMR_GLGL_RFN,AMR_GLGL_CRS)

#endif
      t2 = dnekclock()
      AMR_TCPD = AMR_TCPD + t2 - t1

      ! simple timing
      t1 = dnekclock()
      ! new element-processor mapping
      call mapelpr
      t2 = dnekclock()
      AMR_TCPM = AMR_TCPM + t2 - t1

      ! simple timing
      t1 = dnekclock()
      ! redistribute new BC and curvature face marks
      call amr_tree_transfer(IGROUP,AMR_LEVEL,crvl,bcl)

      ! get b.c. position in bc and cbc array
      ibcl = 2
      if (IFFLOW) ibcl = 1
      nfldl = 1
      if (IFHEAT) nfldl = 2+NPSCAL
      if (IFMHD ) nfldl = 2+NPSCAL+1
      if (PARAM(32).gt.0) nfldl = ibcl + PARAM(32)-1

      ! partially fill in BC and CBC arrays; only string flag for internal and periodic faces
      lcbc=18*LELT*(LDIMT1 + 1)
      lrbc=30*LELT*(LDIMT1 + 1)
      call rzero(BC, lrbc)
      call blank(CBC,lcbc)
      call blank(cbc_bmap,sizeof(cbc_bmap))
      call izero(boundaryID, size(boundaryID))
      call izero(boundaryIDt, size(boundaryIDt))
      do el = 1, NELT
        do il =1, AMR_NFCS
          if (bcl(il,el).eq.0) then
            do jl=ibcl,nfldl
              CBC(il,el,jl) = 'E  '
            enddo
          else if(bcl(il,el).eq.-1) then
            do jl=ibcl,nfldl
              CBC(il,el,jl) = 'P  '
            enddo
          else
            BC(5,il,el,1) = bcl(il,el)
            ! no difference between V and T mesh for now
            boundaryID(il,el) =  bcl(il,el)
            boundaryIDt(il,el) =  bcl(il,el)
          endif
          BC(4,il,el,1) = crvl(il,el)
        enddo
      enddo

      ! redistribute refinement history data
      call amr_hst_map_transfer
      call amr_hst_rfn_transfer
      call amr_hst_crs_transfer

      ! check refinement consistency and buffer size
      if (AMR_NELT_O.ne.(AMR_MAP_NR+AMR_RFN_NR/AMR_NVRT+AMR_CRS_NR_S))
     $   call amr_abort("Error: mesh_regenerate; wrong el. num.")
      if (AMR_NVRT*LELT.le.AMR_RFN_NR_S)
     $   call amr_abort("Error: mesh_regenerate; increase lelt ref.")
      if (AMR_NVRT*LELT.le.(NELT+AMR_CRS_NR*(AMR_NVRT-1)))
     $   call amr_abort("Error: mesh_regenerate; increase lelt crs.")

#ifdef DEBUG
      ! testing
      call io_file_freeid(iunit, ierr)
      write(str1,'(i2.2)') NID
      write(str2,'(i2.2)') icalld
      open(unit=iunit,file='hist.txt'//str1//'i'//str2)
      write(iunit,*) NELV, NELT, NELGV, NELGT,AMR_NELV,AMR_NELT
      write(iunit,*) AMR_MAP_NR,AMR_RFN_NR,AMR_CRS_NR
      write(iunit,*) 'AMR_GLGL_MAP'
      do jl=1,AMR_NELT
         write(iunit,*) jl, jl + AMR_NELIT, AMR_GLGL_MAP(1,jl),
     $                  AMR_GLGL_MAP(3,jl)
      enddo
      write(iunit,*) 'AMR_GLGL_RFN'
      do jl=1,AMR_RFN_NR
         write(iunit,*) jl, (AMR_GLGL_RFN(il,jl),il=1,3)
      enddo
      write(iunit,*) 'AMR_GLGL_CRS'
      do jl=1,AMR_CRS_NR
         write(iunit,*) jl, (AMR_GLGL_CRS(1:2,il,jl),il=1,AMR_NVRT)
      enddo
      close(iunit)
#endif

      t2 = dnekclock()
      AMR_TCPR = AMR_TCPR + t2 - t1

      ! simple timing
      t1 = dnekclock()
      ! get mesh topology
      call amr_topol_get()
      t2 = dnekclock()
      AMR_TCTD = AMR_TCTD + t2 - t1

      return
      end subroutine
!=======================================================================
!> @brief Reinitialise solver
!! @details Mesh dependent logical variables, BC, global sum communicators
      subroutine amr_nek_reinit
      implicit none

      include 'SIZE'
      include 'INPUT'           ! IFMVBD, NPSCAL
      include 'TSTEP'           ! NELFLD, NSTEPS
      include 'RESTART'         ! NELB
      include 'GEOM'
      include 'AMR'
      include 'VPROJ'           ! variables necessary to reset velocity projection for P_n-P_n-2

      ! global variables
      ! for characteristics
      real ct_vx(0:lorder+1)
      common /cchar/ ct_vx ! time for each slice in c_vx()

      ! variables necessary to reset pressure projection for P_n-P_n-2
      integer nprv(2)
      common /orthbi/ nprv

      ! scratch for uzawa_gmres
      integer intype, icg
      real res   (lx2,ly2,lz2,lelv)
      real h1    (lx1,ly1,lz1,lelv)
      real h2    (lx1,ly1,lz1,lelv)
      real h2inv (lx1,ly1,lz1,lelv)
      common /scruz/ res,h1,h2,h2inv

      ! local variables
      integer il,ibcl,nfldl,istart, iend, igeom
      real kwave2
      logical ifemati

      ! simple timing
      real t1, t2

      ! functions
      integer igl_running_sum
      real dnekclock
!-----------------------------------------------------------------------
      ! for AMR I assume bc conditions are set in usrdat and projection
      ! of gll points in usrdat2 (the last one can require communication)
      ! set bc conditions
      call usrdat

      ! get b.c. position in bc and cbc array
      ibcl = 2
      if (IFFLOW) ibcl = 1
      nfldl = 1
      if (IFHEAT) nfldl = 2+NPSCAL
      if (IFMHD ) nfldl = 2+NPSCAL+1
      if (PARAM(32).gt.0) nfldl = ibcl + PARAM(32)-1

      ! logical arrays and BC conditions
      ! setvar
      ! reset field - element number map
      istart = 1
      if (IFMVBD) istart = 0
      iend = nfldl+(LDIMT-1 - NPSCAL)

      do il = istart, iend
        if (IFTMSH(il)) then
            NELFLD(il) = NELT
        else
            NELFLD(il) = NELV
        endif
      enddo

      ! Setup domain topology; main gs communicator
      ! gs_counter not supported yet !!!
      call setup_topo

      ! set interpolation arrays
      call amr_set_noncon()

      ! deform local elements using projection on the surface
      call usrdat2

      ! this is not critical, and in reality the result is not correct anyhow,
      ! as global number of degrees of freedom does change during the simulation,
      ! so work per grid point is not well defined
      call dofcnt

      ! no need to reinitialise I/O

      ! make sure all the interpolated fields are continuous
      call amr_h1avg

      ! recalculate jacobians...
      ! gengeom calls in addition einit, but this doesn't seem to
      ! set anything mesh dependent
      call geom_reset(1)

      ! initalize mesh dependent logical flags
      ! this operation could be simplified to aviod global communication,
      ! but htis means duplication of code
      call setlog(.false.)

      ! zero out masks corresponding to Dirichlet boundary points.
      call bcmask

      ! Need eigenvalues to set tolerances in presolve
      ! to force eigenvalues writing; just temporary solution
      igeom = 2
      call geneig (igeom)

      ! pressure solver initialization (uses "SOLN" space as scratch)
      ! last sentence is not true for AMR
      ! potentially expensive, so time it separately
      t1 = dnekclock()
      if (ifflow.and.iftran) call prinit
      t2 = dnekclock()
      AMR_TCSC = AMR_TCSC + t2 - t1

      ! reset pressure projection for P_n-P_n-2
      if (int(PARAM(95)).gt.0) then
         PARAM(95) = ISTEP+5
         nprv(1) = 0 ! veloctiy field only
      endif

      ! reset velocity projection for P_n-P_n-2
      ! variable from VPROJ include file
      ! hsolve and hmhzsf use it slightly differently
      if (int(PARAM(94)).gt.0.or.ifprojfld(1)) then
         PARAM(94) = ISTEP+5
         ivproj(2,1) = 0
         ivproj(2,2) = 0
         if (IF3D) ivproj(2,3) = 0
      endif

      ! Reset constant mass flow variables
      call vol_flow

      ! restart uzawa_gmres
      if (.not. ifsplit.and.param(42).eq.0)
     $    call uzawa_gmres(res,h1,h2,h2inv,intype,icg)

      ! mark convection dealiasing for recalculation
      call set_dealias_rx

      ! restart catalyst/libsim
#ifdef VISIT
      ! not done yet
#elif CATALYST
      ! check if coprocessor is active; needed because refinement can
      ! be performed before coprocessor is started
      call testcoproc(il)
      if (il.ne.0) then
         ! reset the mesh
         il = 1
         call creategrid(xm1, ym1, zm1, lx1, ly1, lz1, lelt, ldim, il)
      endif
#endif

      ! recalculate convective velocity
!     in the old version for characteristics I simply reset the number of collected fields
!     to zero, but in general I could rebuilt the whole space based on
!     vlag arrays; to consider for future!!!!!!!
!      if (IFCHAR) ct_vx(0) = 0.0
!      call setup_convect(2)
!     for new version nothing done yet

      return
      end subroutine
!=======================================================================
!> @brief Recalculate right hand side arrays
!! @remarks This routine uses global scratch space SCRNS
      subroutine amr_makef
      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'TSTEP'
      include 'SOLN'
      include 'AMR'

      ! local variables
      integer il
      integer ntot1, ntot2, prstep, ierr
      real timel, dtl

      ! work arrays
      real vxl(lx1,ly1,lz1,lelv),vyl(lx1,ly1,lz1,lelv),
     $     vzl(lx1,ly1,lz1,lelv),prl(lx2,ly2,lz2,lelv)
      common /scrns/ vxl, vyl, vzl, prl
!-----------------------------------------------------------------------
      ntot1 = nx1*ny1*nz1*nelv
      ntot2 = nx2*ny2*nz2*nelv

      ! save current fields
      timel = time
      dtl = dt
      call opcopy(vxl,vyl,vzl,vx,vy,vz)
      call copy(prl,pr,ntot2)

      ! initialise logicals
      call setsolv

      ! fill in previous dt values
      if (nbd.eq.3) dtlag(2) = AMR_CHK_TIME(1)-AMR_CHK_TIME(2)
      dtlag(1) = time - AMR_CHK_TIME(1)

      ! set old condition
      time = AMR_CHK_TIME(nbd-1)
      do il=1,nbd-1
         dt = dtlag(nbd-il)
         call opcopy(vx,vy,vz,vxlag(1,1,1,1,nbd-il),
     $      vylag(1,1,1,1,nbd-il),vzlag(1,1,1,1,nbd-il))
         if (il.ge.2) call copy(prlag,pr,ntot2)
         call copy(pr,AMR_PRLAG(1,1,1,1,nbd-il),ntot2)

         call setup_convect(2) ! Save conv vel
         call setprop

         ! compare with routine makef in navier1.f
         call makeuf
         if (filterType.eq.2)                      call make_hpf
         if (ifexplvis.and.ifsplit)                call makevis
         if (ifnav .and..not.ifchar)               call advab
         if (ifmvbd.and..not.ifchar)               call admeshv
         if (iftran)                               call makeabf
         ! the rest is not here as it does not update ab[xyz][12]

         ! palace for passive scalars (VGRADT[12]) and magnetic field (BB[XYZ][12])
         ! perturbation

         ! update time
         time = time +dt
      enddo

      ! put fields back
      time = timel
      dt = dtl
      call opcopy(vx,vy,vz,vxl,vyl,vzl)
      ! copy lag pressure if necessary
      if (nbd.ge.2) call copy(prlag,pr,ntot2)
      call copy(pr,prl,ntot2)

      return
      end subroutine
!=======================================================================
!> @brief Set oefficients for time integration and get right-hand sides
      subroutine amr_set_abbd()
      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'TSTEP'
      include 'AMR'
!-----------------------------------------------------------------------
      if (AMR_IOSTART.gt.0) then
         ! this is just a hack to set correctly variable EETIME0 for timing
         ! it allows to have proper readings of step timing in log
         ! unfortunately EETIME0 is internal variable in comment, so I have to
         ! call it once even though no step is done
         if (NIO.eq.0) then
            write(*,*) 'Fake call to comment to initialilse EETIME0.'
         endif
         istep  = 1
         call comment

         ! set coefficients for time integration
         ! get number of snapshots in a set
         if (PARAM(27).lt.0) then
            istep = NBDINP
         else
            istep = 3
         endif
         call setordbd
         call rzero(bd,10)
         call setbd(bd,dtlag,nbd)
         call rzero(ab,10)
         call setabbd(ab,dtlag,nab,nbd)
         ! possible place for moving boundarries...; not supported yet

         ! set right hand sides
         call amr_makef

         ! correct istep before time loop
         istep = istep - 1

      else
         istep  = 0
      endif

      return
      end subroutine
!=======================================================================
