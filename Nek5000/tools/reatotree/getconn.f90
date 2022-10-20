!> @file getconn.f90
!! @ingroup reatotree
!! @brief Subroutines to generate connectivity for p4est based on data from getmesh.f90
!! @author Adam Peplinski
!! @date 17 Oct 2019
!===============================================================================
module getconn
  use paramm, only : i1b, isp, rdp, ndim, n_fcs, n_vrts, n_fvrts, n_edg
  use getmesh_tool, only : genconn_log
  use getmesh, only : mesh_type
  implicit none

  private
  public :: conn_type, connectivity_create, connectivity_create_simple, conn_dealloc

  type conn_type
     ! connectivity data
     integer(isp) :: num_vertices, num_trees, num_corners
#if N_DIM == 3
     integer(isp) :: num_edges
#endif
     integer(isp), dimension(:), allocatable :: tree_to_vertex
     integer(isp), dimension(:), allocatable :: tree_to_corner
     
     integer(isp), dimension(:), allocatable :: tree_to_tree
     integer(i1b), dimension(:), allocatable :: tree_to_face

     integer(isp), dimension(:), allocatable :: corner_to_tree
     integer(i1b), dimension(:), allocatable :: corner_to_corner
     integer(isp), dimension(:), allocatable :: ctt_offset     
#if N_DIM == 3
     integer(isp), dimension(:), allocatable :: tree_to_edge
     integer(isp), dimension(:), allocatable :: edge_to_tree
     integer(i1b), dimension(:), allocatable :: edge_to_edge
     integer(isp), dimension(:), allocatable :: ett_offset
#endif
     
     ! periodic faces mark
     integer(isp), dimension(:), allocatable :: pfcs
     
     ! vertex physical position
     real(rdp), dimension(:), allocatable :: vertices
  end type conn_type

  ! connectivity parameter arrays
  ! face vertices
  integer(isp), parameter, dimension(4,6) :: vface = reshape(&
       &(/ 1,3,5,7 , 2,4,6,8 , 1,2,5,6 , 3,4,7,8 , 1,2,3,4 , 5,6,7,8 /),shape(vface))

  ! edge vertices
  integer(isp), parameter, dimension(2,n_edg) :: vedge  = reshape(&
       &(/ 1,2 , 3,4 , 5,6 , 7,8 , 1,3 , 2,4 , 5,7 , 6,8 , 1,5 , 2,6 , 3,7 , 4,8 /),shape(vedge))

  ! edge relatd faces
  integer(isp), parameter, dimension(2,n_edg) :: eface  = reshape(&
       &(/ 3,5 , 4,5 , 3,6 , 4,6 , 1,5 , 2,5 , 1,6 , 2,6 , 1,3 , 2,3 , 1,4 , 2,4 /),shape(eface))

  ! corner relatd faces
  integer(isp), parameter, dimension(3,8) :: cface = reshape(&
       &(/ 1,3,5 , 2,3,5 , 1,4,5 , 2,4,5 , 1,3,6 , 2,3,6 , 1,4,6 , 2,4,6 /),shape(cface))

  ! corner relatd edges
  integer(isp), parameter, dimension(3,8) :: cedge = reshape(&
       &(/ 1,5,9 , 1,6,10 , 2,5,11 , 2,6,12 , 3,7,9 , 3,8,10 , 4,7,11 , 4,8,12 /),shape(cedge))

  ! corner to face corner
  integer(isp), parameter, dimension(6,8) :: cfcrn = reshape(&
       &(/ 1,-1, 1,-1, 1,-1 , -1, 1, 2,-1, 2,-1 ,  2,-1,-1, 1, 3,-1 &
       &, -1, 2,-1, 2, 4,-1 ,  3,-1, 3,-1,-1, 1 , -1, 3, 4,-1,-1, 2 &
       &,  4,-1,-1, 3,-1, 3 , -1, 4,-1, 4,-1, 4 /),shape(cfcrn))

  ! corner to face corner
  integer(isp), parameter, dimension(2,12) :: cecrn = reshape(&
       &(/ 1,2 , 3,4 , 5,6 , 7,8 , 1,3 , 2,4 , 5,7 , 6,8 , 1,5 , 2,6 , 3,7 , 4,8 /),shape(cecrn))

  ! face corners
  integer(isp), parameter, dimension(4,6) :: fcsc = reshape(&
       &(/ 1,3,5,7 , 2,4,6,8 , 1,2,5,6 , 3,4,7,8 , 1,2,3,4 , 5,6,7,8 /),shape(fcsc))

  ! to calculate neighbour face corner
  integer(isp), parameter, dimension(6,6) :: rt =reshape( &
       &(/ 1,2,2,1,1,2 , 3,1,1,2,2,1 , 3,1,1,2,2,1 , 1,3,3,1,1,2 &
       &, 1,3,3,1,1,2 , 3,1,1,3,3,1 /),shape(rt))
  integer(isp), parameter, dimension(4,3) :: qt = reshape(&
       &(/ 2,3,6,7 , 1,4,5,8 , 1,5,4,8 /),shape(qt))
  integer(isp), parameter, dimension(4,8) :: pt = reshape(&
       &(/ 1,2,3,4 , 1,3,2,4 , 2,1,4,3 , 2,4,1,3 , 3,1,4,2 , 3,4,1,2 &
       &, 4,2,3,1 , 4,3,2,1 /),shape(pt))
  
contains
  !===============================================================================
  !> @brief Fill in connectivity information
  !! @details This subroutine creates connecitvity based on the full tree information
  !!      including boundary conditions.
  !! @param[inout]  conn      connectivity variable
  !! @param[inout]  mesh      mesh variable
  !! @param[in]  session   session name
  !! @param[in]  myid      mpi rank id
  subroutine connectivity_create(conn,mesh,session,myid)
    implicit none

    ! argument list
    type(conn_type), intent(inout) :: conn
    type(mesh_type), intent(inout) :: mesh
    character(len=80), intent(in) :: session
    integer(isp), intent(in) :: myid

    ! local variables
!#define DEBUG
#ifdef DEBUG
    integer(isp) :: il, jl, kl
#endif
    !-------------------------------------------------------------------------------
    ! get tree structure
    call conn_get(conn,mesh,session)

    ! testing
#ifdef DEBUG
    if (myid==0) then
       open(unit=1,file='tree_test.txt')
       write(1,*) conn%num_vertices, conn%num_trees, conn%num_corners
       write(1,*) 'vertices'
       do il=1,conn%num_vertices
          write(1,*) il-1,conn%vertices(3*il-2),conn%vertices(3*il-1),conn%vertices(3*il)
       enddo
       write(1,*) 'ttv, ttc '
       do il=1,conn%num_trees
          do jl=1,n_vrts
             kl=(il-1)*n_vrts +jl
             write(1,*) il-1,jl-1,conn%tree_to_vertex(kl), conn%tree_to_corner(kl)
          enddo
       enddo
       write(1,*) 'ttt, ttf'
       do il=1,conn%num_trees
          do jl=1,n_fcs
             kl=(il-1)*n_fcs +jl
             write(1,*) il-1,jl-1,conn%tree_to_tree(kl), conn%tree_to_face(kl)
          enddo
       enddo
       write(1,*) 'ctt, ctc'
       do il=1,conn%ctt_offset(conn%num_corners+1)
          write(1,*) conn%corner_to_tree(il), conn%corner_to_corner(il)
       enddo
       write(1,*) 'ctto'
       do il=1,conn%num_corners+1
          write(1,*) conn%ctt_offset(il)
       enddo
#if N_DIM == 3
       write(1,*) 'tte, nedges = ',conn%num_edges
       do il=1,conn%num_trees
          do jl=1,n_edg
             kl=(il-1)*n_edg + jl
             write(1,*) il-1,jl-1,conn%tree_to_edge(kl)
          enddo
       enddo
       write(1,*) 'ett, ete'
       do il=1,conn%ett_offset(conn%num_edges+1)
          write(1,*) conn%edge_to_tree(il), conn%edge_to_edge(il)
       enddo
       write(1,*) 'etto'
       do il=1,conn%num_edges+1
          write(1,*) conn%ett_offset(il)
       enddo
#endif
       close(1)
    endif
#endif
#undef DEBUG
    ! create new connectivity
#if N_DIM == 3
    call fp4est_cnn_new(conn%num_vertices,conn%num_trees,conn%num_edges,conn%num_corners,conn%vertices,&
         &conn%tree_to_vertex,conn%tree_to_tree,conn%tree_to_face,conn%tree_to_edge,&
         &conn%ett_offset,conn%edge_to_tree,conn%edge_to_edge,&
         &conn%tree_to_corner,conn%ctt_offset,conn%corner_to_tree,conn%corner_to_corner)
#else
    call fp4est_cnn_new(conn%num_vertices,conn%num_trees,conn%num_corners,conn%vertices,&
         &conn%tree_to_vertex,conn%tree_to_tree,conn%tree_to_face,&
         &conn%tree_to_corner,conn%ctt_offset,conn%corner_to_tree,conn%corner_to_corner)
#endif

    return
  end subroutine connectivity_create
  !===============================================================================
  !> @brief Fill in connectivity information
  !! @details This subroutine is inteded for testing. It fills connectivity
  !!     arrays using p4est subroutine, but neglects all boundary cond.
  !! @param[inout]  conn      connectivity variable
  !! @param[inout]  mesh      mesh variable
  !! @param[in]  session   session name
  !! @param[in]  myid      mpi rank id
  subroutine connectivity_create_simple(conn,mesh,session,myid)
    implicit none

    ! argument list
    type(conn_type), intent(inout) :: conn
    type(mesh_type), intent(inout) :: mesh
    character(len=80), intent(in) :: session
    integer(isp), intent(in) :: myid

    ! local variables
    integer(isp) :: il
    !-------------------------------------------------------------------------------
    ! get tree structure
    call conn_get(conn,mesh,session)

    ! set initial (not correct) values of corners and trees interconnection
    conn%num_corners = 0
    conn%ctt_offset(:) = 0
    conn%corner_to_tree(:) = 0
    conn%corner_to_corner(:) = 0
    conn%tree_to_corner(:) = -1

    do il=1,n_fcs*conn%num_trees
       ! ttt adn ttf must be initially set to its own numbers
       conn%tree_to_tree(il) = (il-1)/n_fcs
       conn%tree_to_face(il) = int(mod(il-1,n_fcs),i1b)
    enddo

#if N_DIM == 3
    conn%num_edges = 0

    conn%ett_offset(:) = 0
    conn%edge_to_tree(:) = 0
    conn%edge_to_edge(:) = 0
    conn%tree_to_edge(:) = -1
#endif

    ! create new connectivity; this connectivity is not complete,
    ! as e.g. tree_to_tree, tree_to_face adn corners information
    ! is not set properly.
#if N_DIM == 3
    call fp4est_cnn_new(conn%num_vertices,conn%num_trees,conn%num_edges,conn%num_corners,conn%vertices,&
         &conn%tree_to_vertex,conn%tree_to_tree,conn%tree_to_face,conn%tree_to_edge,&
         &conn%ett_offset,conn%edge_to_tree,conn%edge_to_edge,&
         &conn%tree_to_corner,conn%ctt_offset,conn%corner_to_tree,conn%corner_to_corner)
#else
    call fp4est_cnn_new(conn%num_vertices,conn%num_trees,conn%num_corners,conn%vertices,&
         &conn%tree_to_vertex,conn%tree_to_tree,conn%tree_to_face,&
         &conn%tree_to_corner,conn%ctt_offset,conn%corner_to_tree,conn%corner_to_corner)
#endif

    ! complete the information of internal structure of the grid based
    ! on tree_to_vertex. Corners are computed as well. The boundary
    ! conditions are missing
    call fp4est_cnn_complete

    return
  end subroutine connectivity_create_simple
  !===============================================================================
  !> @brief Fill connectivity arrays based on the mesh informations
  !! @param[inout]  conn      connectivity variable
  !! @param[inout]  mesh      mesh variable
  !! @param[in]     session   session name
  subroutine conn_get(conn,mesh,session)
    use getmesh
    implicit none

    ! argument list
    type(conn_type), intent(inout) :: conn
    type(mesh_type), intent(inout) :: mesh
    character(len=80), intent(in) :: session
    !-------------------------------------------------------------------------------
    ! get data from .rea file
    call makemesh(mesh,session)

    conn%num_trees = mesh%nelt
    conn%num_vertices = mesh%crnk

    call conn_alloc(conn)

    call genconn_log('Getting vertices and trees')
    ! fill vertices
    conn%vertices(1:3*conn%num_vertices) = mesh%vert(1:3*conn%num_vertices)

    ! fill cells
    conn%tree_to_vertex(:) = mesh%cell(:)-1

    call genconn_log('Getting face connectivity')
    ! get tree2tree and tree2face
    call get_tree(conn,mesh)

    call genconn_log('Getting edge connectivity')
#if N_DIM == 3
    call get_edge(conn,mesh)
#endif

    call genconn_log('Getting corner connectivity')
    call get_corner(conn,mesh)
    call genconn_log('Connectivity finished; generating tree structure')
    
    return
  end subroutine conn_get
  !===============================================================================
  !> @brief Allocate connectivity arrays
  !! @param[inout]  conn   connectivity variable
  subroutine conn_alloc(conn)
    implicit none

    ! argument list
    type(conn_type), intent(inout) :: conn

    ! local variables
    integer(isp) :: vsize
    !-------------------------------------------------------------------------------
    vsize = n_vrts*conn%num_trees
    allocate(conn%tree_to_vertex(vsize))
    allocate(conn%tree_to_corner(vsize))
    allocate(conn%corner_to_tree(vsize))
    allocate(conn%ctt_offset(vsize+1))
    allocate(conn%corner_to_corner(vsize))

    vsize = n_fcs*conn%num_trees
    allocate(conn%tree_to_tree(vsize))
    allocate(conn%tree_to_face(vsize))

#if N_DIM == 3
    vsize = n_edg*conn%num_trees
    allocate(conn%tree_to_edge(vsize))
    allocate(conn%edge_to_tree(vsize))
    allocate(conn%edge_to_edge(vsize))
    allocate(conn%ett_offset(vsize+1))
#endif

    vsize = 3*conn%num_vertices
    allocate(conn%vertices(vsize))

    vsize = n_fcs*conn%num_trees
    allocate(conn%pfcs(vsize))

    return
  end subroutine conn_alloc
  !===============================================================================
  !> @brief Dellocate connectivity arrays
  !! @param[inout]  conn   connectivity variable
  subroutine conn_dealloc(conn)
    implicit none

    ! argument list
    type(conn_type), intent(inout) :: conn
    !-------------------------------------------------------------------------------
    deallocate(conn%tree_to_vertex)
    deallocate(conn%tree_to_corner)
    deallocate(conn%corner_to_tree)
    deallocate(conn%corner_to_corner)
    deallocate(conn%ctt_offset)

    deallocate(conn%tree_to_tree)
    deallocate(conn%tree_to_face)

#if N_DIM == 3
    deallocate(conn%tree_to_edge)
    deallocate(conn%edge_to_tree)
    deallocate(conn%edge_to_edge)
    deallocate(conn%ett_offset)
#endif
    
    deallocate(conn%vertices)
    
    deallocate(conn%pfcs)
    
    return
  end subroutine conn_dealloc
  !===============================================================================
  !> @brief Get tree2tree and tree2face
  !! @param[inout]  conn      connectivity variable
  !! @param[in]    mesh      mesh variable
  subroutine get_tree(conn,mesh)
    use getmesh_tool, only : sorti, swap_ipi, sort_tuple
    implicit none

    ! argument list
    type(conn_type), intent(inout) :: conn
    type(mesh_type), intent(inout) :: mesh

    ! local variables
    integer(isp) :: nface, nneigh
    integer(isp) :: il, jl, kl, i1, i2, e1, e2, f1, f2, c0, rl
    integer(isp) :: nseg, ipass, iseg
    integer(isp), dimension(n_fvrts) :: indv, fcst
    integer(isp), dimension(1) :: key
    integer(isp), dimension(:), allocatable :: face_l, ninseg, ind
    logical, dimension(:), allocatable :: ifseg ! is vertex sorted
    character(len=100) :: str
    !-------------------------------------------------------------------------------
    ! allocate arrays
    il = size(conn%tree_to_face)
    allocate(face_l(il),ninseg(n_fvrts*il),ind(il))
    allocate(ifseg(n_fvrts*il))

    nface = n_fcs*conn%num_trees

    ! reset arrays
    conn%tree_to_tree(:) = -1
    conn%tree_to_face(:) = -1
    ifseg(:) = .false.

    ! for periodic bc to identify faces sort elements in face_v
    do il=1,mesh%frnk
       conn%pfcs(il) = il
       face_l(il) = il
       ! sort face corners
       call sorti(mesh%face_v(:,il),indv)
    enddo

    ! Sort faces by corners
    nseg        = 1
    ifseg(1)    = .true.
    ninseg(1)   = mesh%frnk

    do ipass=1, 2  ! Multiple passes eliminates false positives
       do jl=1, n_fvrts   ! Sort within each segment
          il =1
          key(1) = jl
          do iseg=1,nseg
             call sort_tuple(mesh%face_v(:,il:il+ninseg(iseg)-1),key,ind(1:ninseg(iseg)),fcst)
             ! Swap position
             call swap_ipi(face_l(il:il+ninseg(iseg)-1),ind(1:ninseg(iseg)))
             il = il + ninseg(iseg)
          enddo

          do il=2,mesh%frnk
             if (mesh%face_v(jl,il)/=mesh%face_v(jl,il-1)) ifseg(il)=.true.
          enddo

          nseg = 0              !  Count up number of different segments
          do il=1,mesh%frnk
             if (ifseg(il)) then
                nseg = nseg+1
                ninseg(nseg) = 1
             else
                ninseg(nseg) = ninseg(nseg) + 1
             endif
          enddo
       enddo
    enddo

    ! go through all unique faces and create neighbour info
    do il=1,mesh%frnk
       ! neighbour number
       i1=mesh%fcs2off(il)
       i2=mesh%fcs2off(il+1)-1
       nneigh=i2-i1

       if (nneigh==0) then
          ! single face
          e1 = mesh%fcs2cell(i1) -1
          f1 = mesh%fcs2fcs(i1)
          ! am I correctly pointing myself
          if (mesh%face(e1*n_fcs + f1)/=il) then
             write(str,4) e1+1,f1
4            format('ERROR: face not pointing itself',2i5)
             call fsc_abort(str//CHAR(0))
          endif
          ! was I already connected
          if (conn%tree_to_tree(e1*n_fcs + f1)/=-1.or.conn%tree_to_face(e1*n_fcs + f1)/=-1) then
             write(str,5) e1+1,f1
5            format('ERROR: face multiple connections',2i5)
             call fsc_abort(str//CHAR(0))
          endif

          ! check for periodic bc
          rl = 1
          kl = (f1-1)/2 +1
          do jl=1,n_fvrts
             ! face corner
             c0 = mesh%cell(e1*n_vrts + vface(jl,f1))
             if (btest(mesh%iperc(c0),kl-1)) then
                fcst(jl) = mesh%pcrn(kl,c0)
             else
                rl = -1
             endif
          enddo

          if (rl==1) then
             ! periodic face; which face should be connected
             ! sort corners and search
             call sorti(fcst,indv)
             rl=1
             do jl=1,n_fvrts
                do kl = rl, mesh%frnk
                   if (fcst(jl)==mesh%face_v(jl,kl)) then
                      rl=kl
                      exit
                   endif
                enddo
             enddo
             rl = face_l(rl)

             ! is f2 face correct one
             i1=mesh%fcs2off(rl)
             i2=mesh%fcs2off(rl+1)-1
             if (i1/=i2) then
                write(str,6) e1+1,f1, rl
6               format('ERROR: periodic face cannot be internal',3i5)
                call fsc_abort(str//CHAR(0))
             endif
             ! check this face
             e2 = mesh%fcs2cell(i1) -1
             f2 = mesh%fcs2fcs(i1)
             ! am I correctly pointing myself
             if (mesh%face(e2*n_fcs + f2)/=rl) then
                write(str,4) e2+1,f2
                call fsc_abort(str//CHAR(0))
             endif
             ! am I connecting to myself; this is complicated issue, so I leave
             ! it for future investigation; for now just exit
             if (e1==e2) then
                write(str,7) e2+1,f1,f2
7               format('ERROR: cell connected to itsef not supported yet ',3i5)
                call fsc_abort(str//CHAR(0))
             endif
             ! save periodic faces
             conn%pfcs(il) = rl

             i1 = e1*n_fcs + f1
             i2 = f2
             conn%tree_to_tree(i1) = e2
             ! find orientation
             ! primary face
             if (f1>f2) then
                rl=f1
                f1=f2
                f2=rl
                rl=e1
                e1=e2
                e2=rl
             endif

             ! find relative corners
             kl = (f1-1)/2 +1
             c0 = mesh%cell(e1*n_vrts + vface(1,f1))
             c0 = mesh%pcrn(kl,c0)
             rl=-1
             do jl=1,n_fvrts
                if (c0==mesh%cell(e2*n_vrts + vface(jl,f2))) rl=jl-1
             enddo
             ! set arrays
             if (rl==-1) then
                write(str,18) c0,e1+1,e2+1,f1,f2
                call fsc_abort(str//CHAR(0))
             endif
             conn%tree_to_face(i1) = int(n_fcs*rl+i2 - 1,i1b)

             ! if the other face was allready done be sure you have the same orientation
             i2 = conn%tree_to_tree(i1)*n_fcs + i2
             i2 = conn%tree_to_face(i2)
             if (i2/=-1.and.i2/n_fcs/=rl) then
                write(str,8) e1+1,e2+1,f1,f2
8               format('ERROR: orientation erron in periodic bc',4i5)
                call fsc_abort(str//CHAR(0))
             endif

          else
             ! no periodic bc; fill arrays
             conn%tree_to_tree(e1*n_fcs + f1) = e1
             conn%tree_to_face(e1*n_fcs + f1) = int(f1 - 1,i1b)
          endif

       elseif (nneigh==1) then
          ! find orientation
          ! primary face
          if (mesh%fcs2fcs(i1)>mesh%fcs2fcs(i2)) then
             e1 = i1
             i1 = i2
             i2 = e1
          endif
          ! elements and faces
          e1 = mesh%fcs2cell(i1) -1
          e2 = mesh%fcs2cell(i2) -1
          f1 = mesh%fcs2fcs(i1)
          f2 = mesh%fcs2fcs(i2)
          ! am I correctly pointing myself
          if (mesh%face(e1*n_fcs + f1)/=il) then
             write(str,4) e1+1,f1
             call fsc_abort(str//CHAR(0))
          endif
          if (mesh%face(e2*n_fcs + f2)/=il) then
             write(str,4) e2+1,f2
             call fsc_abort(str//CHAR(0))
          endif
          ! was I already connected
          if (conn%tree_to_tree(e1*n_fcs + f1)/=-1.or.conn%tree_to_face(e1*n_fcs + f1)/=-1) then
             write(str,4) e1+1,f1
             call fsc_abort(str//CHAR(0))
          endif
          if (conn%tree_to_tree(e2*n_fcs + f2)/=-1.or.conn%tree_to_face(e2*n_fcs + f2)/=-1) then
             write(str,4) e2+1,f2
             call fsc_abort(str//CHAR(0))
          endif
          ! am I connected to myself; this should be done by periodic bc
          if (e1==e2) then
             write(str,16) e2+1,f1,f2
16           format('ERROR: cell connected to itsef not through periodic BC ',3i5)
             call fsc_abort(str//CHAR(0))
          endif
          ! find relative corners
          c0 = mesh%cell(e1*n_vrts + vface(1,f1))
          rl=-1
          do jl=1,n_fvrts
             if (c0==mesh%cell(e2*n_vrts + vface(jl,f2))) rl=jl-1
          enddo
          ! set arrays
          if (rl==-1) then
             write(str,18) c0,e1+1,e2+1,f1,f2
18           format('ERROR: not recognized face crn',5i5)
             call fsc_abort(str//CHAR(0))
          endif
          conn%tree_to_tree(e1*n_fcs + f1) = e2
          conn%tree_to_face(e1*n_fcs + f1) = int(n_fcs*rl+f2 - 1,i1b)
          conn%tree_to_tree(e2*n_fcs + f2) = e1
          conn%tree_to_face(e2*n_fcs + f2) = int(n_fcs*rl+f1 - 1,i1b)

       else
          ! error
          call fsc_abort('ERROR: wrong neighbour number'//CHAR(0))
       endif
    enddo !loop over unique faces

    ! check if all the faces are finally connected
    do il=1,nface
       if(conn%tree_to_tree(il)==-1.or.conn%tree_to_face(il)==-1) then
          write(str,20) (il-1)/n_fcs+1,mod(il-1,n_fcs)+1
20        format('ERROR: face not connected',2i5)
          call fsc_abort(str//CHAR(0))
       endif
    enddo

    ! check if all periodic faces are corretly linked
    do il=1,mesh%frnk
       e2 = conn%pfcs(il)
       if (il/=conn%pfcs(e2)) then
          write(str,22) il, e2
22        format('ERROR: periodic faces not corretly linked',2i5)
          call fsc_abort(str//CHAR(0))
       endif
    enddo

    ! deallocate arrays
    deallocate(face_l,ninseg,ind,ifseg)
    
    return
  end subroutine get_tree
  !===============================================================================
  !> @brief Get edge information
  !! @details p4est marks as 'edges' the real edges in the mesh for which the
  !!     face-to-face communication does not suffice for data redistribution
  !! @param[inout]  conn      connectivity variable
  !! @param[inout]  mesh      mesh variable
  subroutine get_edge(conn,mesh)
    use getmesh
    use getmesh_tool, only : sorti, swap_ipi, sort_tuple
    implicit none

    ! argument list
    type(conn_type), intent(inout) :: conn
    type(mesh_type), intent(inout) :: mesh

#if N_DIM == 3
    ! local variables
    integer(isp) :: nedge, nneigh
    integer(isp) :: il, i1, i2
    integer(isp) :: jl, kl, ll, e1, e2, f1, f2, c0, c1, d1, rl, ledn, fedn
    integer(isp) :: nseg, ipass, iseg
    integer(isp), parameter :: ledd=20
    integer(isp), dimension(ledd) :: led, iled, fed, ifed
    integer(isp), dimension(2) :: indv, fcst
    integer(isp), dimension(1) :: key
    integer(isp), dimension(:), allocatable :: face_l, iedge, icrn, ninseg, ind, ipere
    integer(isp), dimension(:,:), allocatable :: pedge
    logical, dimension(:), allocatable :: ifseg ! is vertex sorted
    character(len=100) :: str
    !-------------------------------------------------------------------------------
    ! allocate arrays
    il = size(conn%tree_to_edge)
    jl = size(conn%tree_to_vertex)
    allocate(face_l(il),iedge(il),icrn(jl),ninseg(2*il),ind(il),ipere(il),pedge(ndim,il))
    allocate(ifseg(2*il))
    
    ! nothing to do if we are not 3D
    nedge = n_edg*conn%num_trees

    conn%num_edges = 0
    ! reset arrays
    conn%tree_to_edge(:)=-2
    conn%ett_offset(:)=0
    conn%edge_to_tree(:)=0
    conn%edge_to_edge(:)=0
    ifseg(:) = .false.

    ! to find periodic edges sort them
    do il=1,mesh%ernk
       do jl=1,ndim
          pedge(jl,il) = il
       enddo
       ipere(il)=0
       iedge(il)=-1
       face_l(il)=il
       ! sort face corners
       call sorti(mesh%edge_v(:,il),indv)
    enddo
    
    ! Sort edges by corners
    nseg        = 1
    ifseg(1)    = .true.
    ninseg(1)   = mesh%ernk

    do ipass=1, 2  ! Multiple passes eliminates false positives
       do jl=1,  2   ! Sort within each segment
          il =1
          key(1) = jl
          do iseg=1,nseg
             call sort_tuple(mesh%edge_v(:,il:il+ninseg(iseg)-1),key,ind(1:ninseg(iseg)),fcst)
             ! Swap position
             call swap_ipi(face_l(il:il+ninseg(iseg)-1),ind(1:ninseg(iseg)))
             il = il + ninseg(iseg)
          enddo

          do il=2,mesh%ernk
             if (mesh%edge_v(jl,il)/=mesh%edge_v(jl,il-1)) ifseg(il)=.true.
          enddo

          nseg = 0              !  Count up number of different segments
          do il=1,mesh%ernk
             if (ifseg(il)) then
                nseg = nseg+1
                ninseg(nseg) = 1
             else
                ninseg(nseg) = ninseg(nseg) + 1
             endif
          enddo
       enddo
    enddo

#ifdef DEBUG
    ! testing
    open(unit=1,file='edge.tmp')
    write(1,*) nseg
    do il=1,mesh%ernk
       write(1,*) mesh%edge_v(1,il),mesh%edge_v(2,il), face_l(il), il
    enddo
    close(1)
    ! testing
#endif

    ! Go through all unique edges, check pointers and periodicity
    ! This part looks for periodic edges related to periodic faces.
    ! Notice it is not the full set of periodic edges
    do il=1, mesh%ernk
       ! neighbour number
       i1=mesh%edg2off(il)
       i2=mesh%edg2off(il+1)-1
       nneigh=i2-i1

       ! do I point correctly myself
       do jl=i1,i2
          e1 = mesh%edg2cell(jl)-1
          f1 = mesh%edg2edg(jl)
          if (mesh%edge(e1*n_edg+f1)/=il) then
             write(str,4) e1+1,f1
4            format('ERROR: edge not pointing itself',2i5)
             call fsc_abort(str//CHAR(0))
          endif
       enddo

       ! is this edge on periodic bc
       e1 = mesh%edg2cell(i1)-1
       f1 = mesh%edg2edg(i1)
       d1 = (f1-1)/4 +1 !  axis

       do ll=1,3    ! directions
          if (ll/=d1) then ! skip direction parallel to the edge
             rl = 1
             do jl=1,2    ! corner
                ! edge corner
                c0 = mesh%cell(e1*n_vrts + vedge(jl,f1))
                if (btest(mesh%iperc(c0),ll-1)) then
                   fcst(jl) = mesh%pcrn(ll,c0)
                else
                   rl = -1
                endif
             enddo   !cell

             if (rl==1) then
                ! periodic edge; which edge should be connected
                ! sort corners and search
                call sorti(fcst,indv)
                rl=1
                do jl=1,2
                   do kl = rl, mesh%ernk
                      if (fcst(jl)==mesh%edge_v(jl,kl)) then
                         rl=kl
                         exit
                      endif
                   enddo
                enddo
                rl = face_l(rl)
                pedge(ll,il)= rl
                ipere(il)=ipere(il)+1

                ! is this edge pointing itself correctly
                i1=mesh%edg2off(rl)
                i2=mesh%edg2off(rl+1)-1
                do jl=i1,i2
                   e2 = mesh%edg2cell(jl)-1
                   f2 = mesh%edg2edg(jl)
                   if (mesh%edge(e2*n_edg+f2)/=rl) then
                      write(str,4) e2+1,f2
                      call fsc_abort(str//CHAR(0))
                   endif

                   ! am I connecting to myself; this is complicated issue, so I leave
                   ! it for future investigation; for now just exit
                   if (e1==e2) then
                      write(str,6) e2+1,f1,f2
6                     format('ERROR: edge; cell connected to itsef not supported yet ',3i5)
                      call fsc_abort(str//CHAR(0))
                   endif

                enddo
             endif   ! periodic edge
          endif   ! skip directions parallel to edge
       enddo   ! directions

    enddo ! unique edges; periodicity check

    ! check that periodic edges correctly point themselves
    do il=1,mesh%ernk
       do ll=1,ndim
          e1 = pedge(ll,il)
          if (il/=pedge(ll,e1)) then
             write(str,8) il, e1
8            format('ERROR: periodic edges not corretly linked',2i5)
             call fsc_abort(str//CHAR(0))
          endif
       enddo
    enddo

#ifdef DEBUG
    ! testing
    open(unit=1,file='eper.tmp')
    do il=1,mesh%ernk
       write(1,*) il,ipere(il),(pedge(jl,il),jl=1,ndim)
    enddo
    close(1)
    ! testing
#endif

    ! for periodic edge orinetation one needs relation between periodic
    ! corners; not all corners are directly pointed by pcrn, so
    ! build the full cornet to corner relation
    ! !!!!!NOTICE, THIS WILL NOT WORK CORRECTLY WITH THE PERIODIC
    ! BOUNDARY CONDITION CONNECTING THE CELL TO ITSELF!!!!!
    ! TO CHANGE IN THE FUTURE
    do il=1,mesh%crnk
       icrn(il) = il
    enddo
    do kl=1,ndim
       do il=1,mesh%crnk
          rl = min(il,icrn(il))
          do ll=1,ndim
             if (btest(mesh%iperc(il),ll-1)) rl = min(rl,mesh%pcrn(ll,il))
          enddo
          icrn(il) = rl
          do ll=1,ndim
             if (btest(mesh%iperc(il),ll-1)) icrn(mesh%pcrn(ll,il)) = rl
          enddo
       enddo
    enddo

#ifdef DEBUG
    ! testing
    open(unit=1,file='icrn.tmp')
    do il=1,mesh%crnk
       write(1,*) il,icrn(il),(mesh%pcrn(jl,il),jl=1,ndim)
    enddo
    close(1)
    ! testing
#endif

    ! go through all unique edges and find number of unique faces the edge
    ! lies on. If the edge ies on more than 3 unique faces the face to face
    ! communication does not allow to update it correctly, so mark it as
    ! "edge" in p4est meaning.
    ! offset
    conn%ett_offset(1)=0

    do il=1,mesh%ernk
       if (iedge(il)==-1) then ! not set yet
          ! is this face periodic; geather informations on the whole set
          if(ipere(il)==0) then ! no periodicity
             ledn = 1
             led(ledn) = il
          elseif(ipere(il)<=2) then ! periodicity in 1 or 2 directions
             ! not all edges are directly pionted by pedge, so follow pointers
             ledn = 1
             led(ledn) = il
             do ll=1,ndim
                if (il/=pedge(ll,il)) then
                   ledn = ledn + 1
                   led(ledn) = pedge(ll,il)
                   if (iedge(led(ledn))/=-1) then
                      write(str,10) led(ledn)
10                    format('ERROR: edge allready assigned',1i5)
                      call fsc_abort(str//CHAR(0))
                   endif
                   do  ! infinite loop to follow pointers
                      if (ipere(led(ledn))==1) then ! dim of the pointed edge
                         exit
                      elseif (ipere(led(ledn))==2) then
                         do kl=1,ndim
                            rl = pedge(kl,led(ledn))
                            if (iedge(rl)/=-1) then
                               write(str,10) led(ledn)
                               call fsc_abort(str//CHAR(0))
                            endif
                            if(led(ledn)/=rl.and.led(ledn-1)/=rl) then
                               if (il==rl) then
                                  goto 12
                               elseif(led(ledn-1)/=rl)then
                                  ledn = ledn + 1
                                  led(ledn) = rl
                               endif
                            endif
                            if (ledn==ledd) goto 12
                         enddo
                      else ! error no more than 2 periodic directions for single edge
                         write(str,20) led(ledn)
                         call fsc_abort(str//CHAR(0))
                      endif ! dim of the pointed edge
                   enddo ! infinite loop to follow pointers
                endif  ! do I point other edge
             enddo ! directions in the starting point

12           continue

             ! write(*,*) 'periodic edge',il,ledn,(led(jl),jl=1,ledn)

             if (ledn==ledd) then ! to many connections
                write(str,14) ledn
14              format('ERROR: to many periodic edges',1i5)
                call fsc_abort(str//CHAR(0))
             endif

             ! usually there are no more then 4 periodic edges, unless you have
             ! strange periodic conditions; warn the user
             if (ledn>4) then ! to many connections
                write(str,16) ledn
16              format('WARNING: many periodic edges',1i5)
                call genconn_log(str)
             endif

             ! edges should not repeat in the list; check it
             do jl=2,ledn
                call sorti(led(1:ledn),iled(1:ledn))
                if (led(jl-1)==led(jl)) then
                   write(str,18) il
18                 format('ERROR: edges repeatin in the periodic list',1i5)
                   call fsc_abort(str//CHAR(0))
                endif
             enddo

          else    ! error no more than 2 periodic directions for single edge
             write(str,20) il
20           format('ERROR: edges to many periodic directions',1i5)
             call fsc_abort(str//CHAR(0))
          endif   ! is periodic

          ! write(*,*) 'edge set',il,ledn,(led(jl),jl=1,ledn)

          ! we have the full set of periodic edges; mark them as used
          do jl=1,ledn
             iedge(led(jl)) = 1
          enddo

          ! build set of unique faces realted to the edges
          fedn = 0
          do jl=1,ledn
             i1=mesh%edg2off(led(jl))
             i2=mesh%edg2off(led(jl)+1)-1
             do kl=i1,i2
                e1 = mesh%edg2cell(kl)-1
                f1 = mesh%edg2edg(kl)
                do ll=1,2
                   fedn = fedn + 1
                   fed(fedn) = mesh%face(e1*n_fcs + eface(ll,f1))
                   ! is this face periodic
                   fed(fedn) = min(fed(fedn),conn%pfcs(fed(fedn)))
                enddo
             enddo
          enddo
          call sorti(fed(1:fedn),ifed(1:fedn))
          ! count unique faces
          rl=1
          do jl=1,fedn-1
             if (fed(jl)/=fed(jl+1)) rl=rl+1
          enddo

          ! fill connectivity arrays
          if (rl<=3) then !  not an 'edge'
             do jl=1,ledn
                i1=mesh%edg2off(led(jl))
                i2=mesh%edg2off(led(jl)+1)-1
                do kl=i1,i2
                   e1 = mesh%edg2cell(kl)-1
                   f1 = mesh%edg2edg(kl)
                   conn%tree_to_edge(e1*n_edg+f1) = -1
                enddo
             enddo
          else ! mark the edge
             conn%num_edges = conn%num_edges + 1
             ! count number of mesh edges related to this one
             rl = 0
             ! to get orientation pickup single corner of the edge
             !     !!!!! NOTICE; THE WAY PERIODICITY OF CORNERS c0, c1, d1, IS
             !       MARKED IS NOT FULLY GENERAL AND DOES NOT ALLOW FOR CELLS
             !       CONNECTED TO ITSELF
             i1 = mesh%edg2off(led(1))
             e1 = mesh%edg2cell(i1)-1
             f1 = mesh%edg2edg(i1)
             c0 = conn%tree_to_vertex(e1*n_vrts + vedge(1,f1)) +1
             c0 = icrn(c0) ! periodic bc
             c1 = conn%tree_to_vertex(e1*n_vrts + vedge(2,f1)) +1
             c1 = icrn(c1) ! periodic bc
             do jl=1,ledn
                i1=mesh%edg2off(led(jl))
                i2=mesh%edg2off(led(jl)+1)-1
                do kl=i1,i2
                   rl= rl+1
                   e1 = mesh%edg2cell(kl)-1
                   f1 = mesh%edg2edg(kl)
                   conn%tree_to_edge(e1*n_edg+f1) = conn%num_edges -1
                   conn%edge_to_tree(conn%ett_offset(conn%num_edges)+rl) = e1
                   ! check orientation
                   d1 = conn%tree_to_vertex(e1*n_vrts + vedge(1,f1)) +1
                   d1 = icrn(d1) ! periodic bc
                   if (c0==d1) then
                      conn%edge_to_edge(conn%ett_offset(conn%num_edges)+rl) = int(f1-1,i1b)
                   elseif (c1==d1) then
                      conn%edge_to_edge(conn%ett_offset(conn%num_edges)+rl) = int(f1-1 + n_edg,i1b)
                   else !  error
                      write(str,22) e1,f1
22                    format('ERROR: edge does not fit corners',2i5)
                      call fsc_abort(str//CHAR(0))
                   endif
                enddo
             enddo
             conn%ett_offset(conn%num_edges+1)=conn%ett_offset(conn%num_edges)+rl
          endif

       endif   ! is already assigned
    enddo ! unique edges

    ! possible orientation flip
    do il = 1,conn%num_edges
       i1 = conn%ett_offset(il) +1
       i2 = conn%ett_offset(il+1)
       c0=0
       c1=0
       do jl= i1,i2
          if (conn%edge_to_edge(jl)>=n_edg) then
             c0 = c0+1
          else
             c1 = c1 +1
          endif
       enddo
       if (c0>c1) then
          do jl= i1,i2
             if (conn%edge_to_edge(jl)>=n_edg) then
                conn%edge_to_edge(jl) = conn%edge_to_edge(jl) - int(n_edg,i1b)
             else
                conn%edge_to_edge(jl) = conn%edge_to_edge(jl) + int(n_edg,i1b)
             endif
          enddo
       endif
    enddo

    ! final checks
    ! are all edges connected
    do il = 1,nedge
       if (conn%tree_to_edge(il)==-2) then
          write(str,24) il/n_edg + 1,mod(il-1,n_edg) + 1
24        format('ERROR: edge not connected',2i5)
          call fsc_abort(str//CHAR(0))
       endif
    enddo
    ! pointers and orientation
    do il = 1,conn%num_edges
       i1 = conn%ett_offset(il) +1
       i2 = conn%ett_offset(il+1)
       do jl= i1,i2
          c0 = conn%edge_to_tree(jl)
          c1 = conn%edge_to_edge(jl)
          d1 = c1/n_edg +1
          c1 = mod(c1,n_edg) +1
          ! are edges correctly pointing themselves
          if (il-1/=conn%tree_to_edge(c0*n_edg+c1)) then
             write(str,26) il,c0+1,c1
26           format('ERROR: edge not pointing itself',3i5)
             call fsc_abort(str//CHAR(0))
          endif
          ! for orientation check
          ind(jl) = conn%tree_to_vertex(c0*n_vrts + vedge(d1,c1))+1
          ind(jl) = icrn(ind(jl)) ! periodic bc
       enddo
       ! orientation check
       if (i2>i1) then
          do jl=i1+1,i2
             if (ind(i1)/=ind(jl)) then
                write(str,28) il,jl,ind(i1),ind(jl)
28              format('ERROR: edge orientation problem',4i5)
                call fsc_abort(str//CHAR(0))
             endif
          enddo
       endif
    enddo

    ! deallocate arrays
    deallocate(face_l,iedge,icrn,ninseg,ind,ipere,pedge)
    deallocate(ifseg)
#endif

    return
  end subroutine get_edge
  !===============================================================================
  !> @brief Get corner information
  !! @details p4est marks as 'corners' the real corners in the mesh for which the
  !!     face-to-face and 'edge' communication does not suffice for data redistribution
  !! @param[inout]  conn      connectivity variable
  !! @param[inout]  mesh      mesh variable
  subroutine get_corner(conn,mesh)
    use getmesh
    use getmesh_tool, only : sorti
    implicit none

    ! argument list
    ! argument list
    type(conn_type), intent(inout) :: conn
    type(mesh_type), intent(inout) :: mesh

    ! local variables
    integer(isp) :: ncorner
    integer(isp) :: il, jl, kl, ll, ml, rl, ncrn, i1, i2, j1, j2, m1, m2, c0, c1
    integer(isp) :: skip
    integer(isp), dimension(:), allocatable :: icrn, icrnp, icrnoff, ictc, icte
    integer(isp), dimension(ndim) :: iface, ntree, nface, orient, fcorner
#if N_DIM == 3
    integer(isp), dimension(ndim) :: iedge, iwhich, iflip
#endif

    character(len=100) :: str
    !-------------------------------------------------------------------------------
    ! allocate arrays
    il = size(conn%tree_to_corner)
    allocate(icrn(il), icrnp(il), icrnoff(il), ictc(il), icte(il))
    
    ncorner = n_vrts*conn%num_trees

    conn%num_corners = 0
    ! reset arrays
    conn%tree_to_corner(:)=-2
    conn%ctt_offset(:)=0
    conn%corner_to_tree(:)=0
    conn%corner_to_corner(:)=0

    ! Information about periodic bc are already included in tree_to_tree,
    ! tree_to_face and edge information, so we don't recreate it here.
    ! However, we have to reduce the number of points removing the multiple ones.
    ! !!!!! THIS STEP IS DONE IN A SIMPLE WAY THAT DOES NOT SUPPORT ELEMETS CONNECTING TO THEMSELVES !!!!!
    ! TO CHANGE IN THE FUTURE
    do il=1,mesh%crnk
       icrn(il) = il
       icrnp(il) = il
    enddo
    do kl=1,ndim
       do il=1,mesh%crnk
          rl = min(il,icrn(il))
          do ll=1,ndim
             if (btest(mesh%iperc(il),ll-1)) rl = min(rl,mesh%pcrn(ll,il))
          enddo
          icrn(il) = rl
          do ll=1,ndim
             if (btest(mesh%iperc(il),ll-1)) icrn(mesh%pcrn(ll,il)) = rl
          enddo
       enddo
    enddo
    ! order points
    ictc(1:mesh%crnk) = icrn(1:mesh%crnk)
    call sorti(ictc(1:mesh%crnk),icrnp(1:mesh%crnk))
    i1=ictc(1)
    icrnoff(1)=1
    ncrn = 1
    do il=2,mesh%crnk
       if (i1/=ictc(il)) then
          i1=ictc(il)
          ncrn = ncrn +1
          icrnoff(ncrn) = il
       endif
    enddo
    icrnoff(ncrn+1) = mesh%crnk+1

    ! testing
#ifdef DEBUG
    open(unit=1,file='icrn2.tmp')
    do il=1,mesh%crnk
       write(1,*) il,icrn(il),ictc(il),icrnp(il),icrnoff(il)
    enddo
    close(1)
#endif
    skip = 0

    do il =1,ncrn  ! loop over all unique points
       i1 = icrnoff(il)
       i2 = icrnoff(il+1) - 1
       ! gather information for all points related to that one
       rl = 0
       do jl=i1,i2
          j1 = mesh%crn2off(icrnp(jl))
          j2 = mesh%crn2off(icrnp(jl)+1) - 1
          do kl=j1,j2
             rl=rl+1
             ictc(rl) = mesh%crn2vrt(kl)
             icte(rl) = mesh%crn2cell(kl)
          enddo
       enddo

       ! write(*,*) 'corner ',il,rl,(ictc(kl),kl=1,rl)
       ! write(*,*) 'element',il,rl,(icte(kl),kl=1,rl)

       ! is it a real 'corner' in p4est meaning
       do jl=1,rl    ! all vertices related to single point
          ! identify all element corners connected through faces
          j1 = ictc(jl)
          j2 = icte(jl)-1
          do ll=1,ndim
             ! corner related faces
             iface(ll) = cface(ll,j1)
             ntree(ll) = conn%tree_to_tree(j2*n_fcs + iface(ll))
             c0 = conn%tree_to_face(j2*n_fcs + iface(ll))
             if (ntree(ll)==j2.and.iface(ll)==(c0+1)) then
                ntree(ll) = -1
                nface(ll) = -1
                orient(ll) = -1
                fcorner(ll) = -1
             else
                nface(ll) = mod(c0,n_fcs)+1
                orient(ll) = c0/n_fcs +1
                fcorner(ll) = cfcrn(iface(ll),j1)
             endif
             ! corner lrelated edges
#if N_DIM == 3
             c0 = cedge(ll,j1)
             ! which edge
             iedge(ll) = conn%tree_to_edge(j2*n_edg + c0)+1
             ! which edge corner
             if (iedge(ll)>0) then
                if (cecrn(1,c0)==j1) then
                   iwhich(ll) = 0
                elseif(cecrn(2,c0)==j1) then
                   iwhich(ll) = 1
                else   ! error
                   write(str,4) il,jl
4                  format('ERROR: edge corner not mathing',2i5)
                   call fsc_abort(str//CHAR(0))
                endif
                ! orientation
                m1 = conn%ett_offset(iedge(ll)) +1
                m2 = conn%ett_offset(iedge(ll)+1)
                do ml=m1,m2
                   if(j2==conn%edge_to_tree(ml)) then
                      iflip(ll) = conn%edge_to_edge(ml)/n_edg
                   endif
                enddo
             else
                iwhich(ll) = -1
                iflip(ll) = -1
             endif
#endif
          enddo

          ! write(*,*) il,jl,(iface(ll),ll=1,ndim),(ntree(ll),ll=1,ndim)
          ! write(*,*) il,jl,(nface(ll),ll=1,ndim),(orient(ll),ll=1,ndim)
          ! write(*,*) il,jl,(fcorner(ll),ll=1,ndim)
          ! write(*,*) il,jl,(iedge(ll),ll=1,ndim),(iwhich(ll),ll=1,ndim)
          ! write(*,*) il,jl,(iflip(ll),ll=1,ndim)

          ! go through all corners and check if all corners are face connected
          do kl=1,rl
             skip = 0
             if (j1==ictc(kl).and.j2==(icte(kl)-1)) then
                skip = 1
             else
                do ll=1,ndim
                   c0 = -1
                   if ((icte(kl)-1)==ntree(ll)) then
                      ! check wether corners match each other
#if N_DIM == 3
                      c0 = rt(nface(ll),iface(ll))
                      c0 = qt(orient(ll),c0)
                      c0 = pt(fcorner(ll),c0)
#else
                      c0 = ieor(orient(ll)-1,fcorner(ll)-1)+1 !!!!! ths for checking !!!!!!
#endif
                      c0 = fcsc(c0,nface(ll))
                   endif
                   if (c0==ictc(kl)) then
                      skip = 1
                      exit
                   endif
                enddo
             endif

             ! if skip==0 face communcation is not sufficient to exchange data
             !     check edges
#if N_DIM == 3
             if (skip==0) then
                do ll=1,ndim
                   if (iedge(ll)>0) then
                      m1 = conn%ett_offset(iedge(ll)) +1
                      m2 = conn%ett_offset(iedge(ll)+1)
                      do ml=m1,m2
                         if((icte(kl)-1)==conn%edge_to_tree(ml)) then
                            c1 = conn%edge_to_edge(ml)
                            c0 = c1/n_edg
                            c1 = mod(c1, n_edg) + 1
                            c0 = ieor(c0,iflip(ll))
                            c0 = ieor(c0,iwhich(ll))
                            c0 = cecrn(c0+1,c1)
                            if(c0==ictc(kl)) then
                               skip = 1
                               exit
                            endif
                         endif
                      enddo
                      if(skip==1) exit
                   endif
                enddo
             endif ! skip==0
#endif

             if (skip==0) exit
          enddo

          ! write(*,*) 'crn skip',il,jl,rl,skip

          ! if skip==0 we've found one corner not connected through face nor edge
          ! exit loop
          if (skip==0) exit
       enddo   ! all vertices related to single point

       ! if skip==0 face and edge communication is not sufficient
       ! mark as 'corner'
       if (skip==0) then
          conn%num_corners = conn%num_corners + 1
          do jl=1,rl
             j1 = ictc(jl)
             j2 = icte(jl)-1
             conn%tree_to_corner(j2*n_vrts+j1) = conn%num_corners -1
             conn%corner_to_tree(conn%ctt_offset(conn%num_corners)+jl) = j2
             conn%corner_to_corner(conn%ctt_offset(conn%num_corners)+jl) = int(j1 -1,i1b)
          enddo
          conn%ctt_offset(conn%num_corners+1) = conn%ctt_offset(conn%num_corners)+rl
       else    ! skip==0
          do jl=1,rl
             j1 = ictc(jl)
             j2 = icte(jl)-1
             conn%tree_to_corner(j2*n_vrts+j1) = -1
          enddo
       endif   ! skip==0
    enddo     ! all unique points

    ! final checks
    ! are all corners connected
    do il = 1,ncorner
       if (conn%tree_to_corner(il)==-2) then
          write(str,10) il/n_vrts + 1,mod(il-1,n_vrts) + 1
10        format('ERROR: corner not connected',2i5)
          call fsc_abort(str//CHAR(0))
       endif
    enddo
    ! pointers
    do il = 1,conn%num_corners
       i1 = conn%ctt_offset(il) +1
       i2 = conn%ctt_offset(il+1)
       do jl= i1,i2
          c0 = conn%corner_to_tree(jl)
          c1 = conn%corner_to_corner(jl)
          ! are edges correctly pointing themselves
          if (il-1/=conn%tree_to_corner(c0*n_vrts+c1+1)) then
             write(str,12) il,c0+1,c1
12           format('ERROR: corner not pointing itself',3i5)
             call fsc_abort(str//CHAR(0))
          endif
       enddo
    enddo

    ! deallocate arrays
    deallocate(icrn, icrnp, icrnoff, ictc, icte)
    
    return
  end subroutine get_corner
  !===============================================================================
end module getconn
