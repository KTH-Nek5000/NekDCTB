!> @file rea2tree.f90
!! @ingroup reatotree
!! @brief Program to convert Nek5000 mesh format to p4est one based on .rea, .re2 and .map files
!! @author Adam Peplinski
!! @date 18 Oct 2019
program main
  use mpi
  use paramm, only : isp
  use getmesh_tool, only : genconn_log
  use getconn, only : connectivity_create
  use rea2treeM
  implicit none

  integer(isp) :: myid, nprocs, ierr
  character(len=80) :: session
  !-------------------------------------------------------------------------------
  ! Initialise MPI, sc and p4est packages
  call package_init(myid, nprocs)
  ! this version of the program is still serial, so nprocs==1; MPI is used
  ! because p4est requres it
  if (nprocs>1) call fsc_abort('Must be run on single processor'//CHAR(0))

  ! Get session name
  call genconn_log('Give session name; max 80 characters')
  if (myid==0) read(*,*) session
  call mpi_bcast(session,80,mpi_character,0,mpi_comm_world,ierr)

  BCmode = 1
  call genconn_log('Give BC and curved sides treatment mode:')
  call genconn_log('    1 - gmsh style (##.re2 produced by gmsh2nek)')
  call genconn_log('    2 - standard ##.re2 (BC characterised by string; simple geometry)')
  call genconn_log('    3 - manual (small meshes only; BC and CRV can be given separately)')
  if (myid==0) read(*,*) BCmode
  call mpi_bcast(BCmode,1,mpi_integer,0,mpi_comm_world,ierr)

  ! Create connectivity
  call connectivity_create(conn,mesh,session,myid)

  ! create BC mapping list
  if (BCmode == 2) call BCmap_create(mesh,session,myid)

  ! check connectivity correctness; this check doesn't seem to be very detailed
  call fp4est_cnn_valid(ierr)
  call fsc_check_abort(ierr,'Connectivity not valid'//CHAR(0))
  ! saving connectivity; for testing only
  ! call fp4est_cnn_save(trim(session)//'_connectivity'//CHAR(0))
  ! call fp4est_cnn_load(trim(session)//'_connectivity'//CHAR(0))

  ! Create tree forest
  call fp4est_tree_new(MPI_COMM_WORLD,0)
  ! check tree correctness
  call fp4est_tree_valid(ierr)
  call fsc_check_abort(ierr,'Tree forest not valid'//CHAR(0))
  ! save the tree
  call fp4est_tree_save(1,'ATF'//trim(session)//'.t00000'//CHAR(0))
  ! for debugging only
  ! call fp4est_test_print()
  ! call fp4est_tree_load(MPI_COMM_WORLD, 1,'testt_tree'//CHAR(0))

  ! tree structure visualisation
  ! call fp4est_vtk_write(trim(session)//'_plot'//CHAR(0))

  ! delete tree and connectivity
  call fp4est_tree_del()
  call fp4est_cnn_del()

  ! Finalise packages
  call package_finalize
  
  stop
end program main
!===============================================================================
!> @brief Initialise mesh data in p4est
!! @details Subroutine required by p4est to initialise mesh data.
!!     I store there curvature and flows boundary conditions.
!! @param[in]
!! @param[out]
!! @param[out]
!! @param[out]   bcl    externa boundary condition flag
!! @param[out]   crvl   curvature flag
subroutine nekp4est_init_msh_dat(tree,imsh,igrp,bcl,crvl)
  use paramm, only : isp, ndim, n_fcs, n_vrts, n_fvrts
  use rea2treeM, only : mesh, BCmode, nBCmap, cBCmap, iBCmap
  use getconn, only : conn_type
  use getmesh, only : mesh_type
  implicit none
  
  ! argument list
  integer(isp) :: tree, imsh, igrp
  integer(isp), dimension(n_fcs) :: bcl, crvl

  ! local variablesy
  integer(isp) :: il, jl, ll, ml, treel

  ! face renumbering
  ! return Nekton preprocessor face ordering
  integer(isp), parameter, dimension(6) :: eface = (/ 4 , 2 , 1 , 3 , 5 , 6 /)

  ! edge numbering
  ! face related edges
  integer(isp), parameter, dimension(4,6) :: eedge = reshape(&
       &(/ 1,5,9,10, 2,6,10,11, 3,7,11,12, 4,8,9,12, 1,2,3,4, 5,6,7,8/),shape(eedge))

  ! vertex numbering
  ! face related edges
  integer(isp), parameter, dimension(4,6) :: evert = reshape(&
       &(/ 1,3,5,7, 2,4,6,8, 1,2,5,6, 3,4,7,8, 1,2,3,4, 5,6,7,8/),shape(evert))
  !-------------------------------------------------------------------------------
  treel=tree+1

  ! initialise variables
  crvl(:) = 0
  bcl(:) = 0
  ! mark velocity and temperature meshes
  if (tree<mesh%nelv) then
     imsh=0
  else
     imsh=1
  endif

  ! element group
  igrp = mesh%igroup(treel)

  select case(BCmode)
     case (1) ! gmsh style
        do il=1,n_fcs
           if (mesh%cbc(il,treel) == 'P  ') then
               bcl(il) = -1
           else if (mesh%cbc(il,treel) /= 'E  ') then
               bcl(il) = int(mesh%bc(5,il,treel),isp)
               crvl(il) = bcl(il)
           endif
        enddo
     case (2) ! standard ###.re2; surfaces described by strings only
        do il=1,n_fcs
           if (mesh%cbc(il,treel) == 'P  ') then
               bcl(il) = -1
           else if (mesh%cbc(il,treel) /= 'E  ') then
              do jl=1, nBCmap
                 if (cBCmap(jl) == mesh%cbc(il,treel)) then
                    bcl(il) = iBCmap(jl)
                    crvl(il) = bcl(il)
                 endif
              enddo
           endif
        enddo
     case default ! manual
        do il=1,n_fcs
           if (mesh%cbc(il,treel) == 'P  ') then
               bcl(il) = -1
           else if(mesh%cbc(il,treel) /= 'E  ') then
              write(*,*) ' '
              write(*,*) 'Element: ',treel
              write(*,*) 'Face position:'
              write(*,*) '   (p4est) = ',il
              write(*,*) '   (prenek) = ',eface(il)

              write(*,*) 'Vertices (p4est notation):'
              write(*,*) 'x, y, z'
              do jl=1,n_fvrts
                 ml = mesh%cell(tree*n_vrts + evert (jl,il))
                 write(*,*) (mesh%vert((ml-1)*3+ll),ll=1,3)
              enddo

              write(*,*) 'Curvature type (prenek):'
              write(*,*) '   face : ',mesh%ccrv(eface(il),treel)
              if (ndim==3) then
                 write(*,*) '   edges : ',(mesh%ccrv(eedge(jl,il),treel),jl=1,4)
              endif
              write(*,*) 'Boundary type : ',mesh%cbc(il,treel)
              write(*,*) 'Give positive integer for surface id'
              read(*,*) jl
              crvl(il) = abs(jl)
              write(*,*) 'Give positive integer for boundarry condition id'
              read(*,*) jl
              bcl(il) = abs(jl)
              write(*,*) 'Face marked as: ',crvl(il),bcl(il)
           endif
        enddo
  end select

  return
end subroutine nekp4est_init_msh_dat
!===============================================================================
