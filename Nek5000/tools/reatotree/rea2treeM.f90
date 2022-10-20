!> @file rea2treeM.f90
!! @ingroup reatotree
!! @brief Main data structures and routines for rea2tree tool
!! @author Adam Peplinski
!! @date 17 Oct 2019
!===============================================================================
module rea2treeM
  use paramm, only : isp, n_fcs
  use getmesh_tool, only : genconn_log
  use getconn, only : conn_type, conn_dealloc
  use getmesh, only : mesh_type, mesh_dealloc
  implicit none

  private
  public :: mesh, conn, BCmode, nBCmap, cBCmap, iBCmap
  public :: package_init, package_finalize, BCmap_create

  ! main data structures
  type(mesh_type), save :: mesh
  type(conn_type), save :: conn

  ! flag for BC and curved sides reatment mode
  integer, save :: BCmode
  integer(isp), save :: nBCmap
  character(len=3), dimension(:), allocatable, save :: cBCmap
  integer(isp), dimension(:), allocatable, save :: iBCmap

contains
  !===============================================================================
  !> @brief Fill in connectivity information
  !! @param[out]  myid      mpi rank id
  !! @param[out]  numprocs  number of processors
  subroutine package_init(myid, numprocs)
    use mpi
    implicit none

    ! argument list
    integer(isp), intent(out) :: myid, numprocs

    ! local variables
    integer(isp) :: ierr ! error check
    !-------------------------------------------------------------------------------
    ! MPI init
    call mpi_init(ierr)
    call fsc_check (ierr)
    call mpi_comm_rank(mpi_comm_world, myid, ierr)
    call fsc_check (ierr)
    call mpi_comm_size(mpi_comm_world, numprocs, ierr)
    call fsc_check (ierr)

    ! sc init
    call fsc_set_log(7)
    call fsc_init(mpi_comm_world, 1, 1, 7)

    ! p4est init
    call fp4est_init(7)

    return
  end subroutine package_init
  !===============================================================================
  !> @brief Finalisation of packages
  subroutine package_finalize()
    use mpi
    implicit none

    ! local variables
    integer(isp) :: ierr ! error check
    !-------------------------------------------------------------------------------
    ! sc end
    call fsc_pkg_print(5)
    call fsc_finalize()

    ! clean memory
    call conn_dealloc(conn)
    call mesh_dealloc(mesh)

    ! MPI end
    call mpi_finalize(ierr)

    return
  end subroutine package_finalize
  !===============================================================================
  !> @brief Create BC to integer mapping mapping
  !! @param[in]  mesh      mesh variable
  !! @param[in]  session   session name
  !! @param[in]  myid      mpi rank id
  subroutine BCmap_create(mesh,session,myid)
    use mpi
    implicit none

    ! argument list
    type(mesh_type), intent(inout) :: mesh
    character(len=80), intent(in) :: session
    integer(isp), intent(in) :: myid

    ! local variables
    character(len=3), dimension(:), allocatable :: cbctmp

    integer(isp) :: iel, ifc, il
    integer(isp) :: ierr ! error check
    !-------------------------------------------------------------------------------
    ! initialisation
    nBCmap=0
    allocate(cBCmap(10))

    ! gather unique BC
    ! this part is serial
    do iel = 1,mesh%nelt
       do ifc=1,n_fcs
          extern : if (mesh%cbc(ifc,iel) /= 'P  '.and.mesh%cbc(ifc,iel) /= 'E  ') then
             do il = 1, nBCmap
                if(cBCmap(il) == mesh%cbc(ifc,iel)) exit extern
             enddo
             nBCmap = nBCmap + 1
             if (nBCmap > size(cBCmap)) then
                allocate(cbctmp(nBCmap - 1))
                cbctmp=cBCmap
                deallocate(cBCmap)
                allocate(cBCmap(nBCmap + 9))
                cBCmap(1:nBCmap-1)=cbctmp(1:nBCmap-1)
                deallocate(cbctmp)
             endif
             cBCmap(nbcmap) = mesh%cbc(ifc,iel)
          endif extern
        enddo
    enddo

    ! list BC
    call genconn_log('Unique BC list:')
    do il=1, nBCmap
       call genconn_log(cBCmap(il))
    enddo

    ! get mapping
    allocate(iBCmap(nBCmap))
    call genconn_log('Give BC to integer mapping:')
    do il=1, nBCmap
       call genconn_log(cBCmap(il))
       if (myid==0) read(*,*) iBCmap(il)
    enddo
    call mpi_bcast(iBCmap,nBCmap,mpi_integer,0,mpi_comm_world,ierr)

    ! check for duplicates
    ierr = 1
    do il=1, nBCmap-1
       do ifc =il+1,nBCmap
          if (iBCmap(il) == iBCmap(ifc)) then
             ierr = 0
          endif
       enddo
    enddo
    call fsc_check_abort(ierr,'Multiple BC are mapped to single integer'//CHAR(0))


    return
  end subroutine BCmap_create
  !===============================================================================
end module rea2treeM
  
