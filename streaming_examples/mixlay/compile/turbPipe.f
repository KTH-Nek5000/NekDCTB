#define RAD  1.0         !pipe radius
#define ZLENPIPE  8.0   !pipe length
!-----------------------------------------------------------------------
!
!     user subroutines required by nek5000
!
!     Parameters used by this set of subroutines:
!
!-----------------------------------------------------------------------
      subroutine uservp (ix,iy,iz,ieg)
      include 'SIZE'
      include 'NEKUSE'          ! UDIFF, UTRANS

      UDIFF =0.0
      UTRANS=0.0

      return
      end
!-----------------------------------------------------------------------
      subroutine userf  (ix,iy,iz,ieg)
      
      include 'SIZE'
      include 'NEKUSE'          ! FF[XYZ]
      include 'PARALLEL'

      integer ix,iy,iz,ieg

      ffx = 0.0
      ffy = 0.0
      ffz = 0.0

      return
      end
!-----------------------------------------------------------------------
      subroutine userq  (ix,iy,iz,ieg)
      include 'SIZE'
      include 'NEKUSE'          ! QVOL

      QVOL   = 0.0

      return
      end
!-----------------------------------------------------------------------
      subroutine userchk
!      implicit none
      include 'SIZE'            !
      include 'TSTEP'           ! ISTEP, lastep, time
      include 'INPUT'           ! IF3D, PARAM
      include 'SOLN'   

!for computing tau_wall on the fly
      real ffx_new,ffy_new,ffz_new  
      common /cforce/ ffx_new,ffy_new,ffz_new
c   for torque calculations
      common /ctorq/ dragx(0:maxobj),dragpx(0:maxobj),dragvx(0:maxobj)
     $             , dragy(0:maxobj),dragpy(0:maxobj),dragvy(0:maxobj)
     $             , dragz(0:maxobj),dragpz(0:maxobj),dragvz(0:maxobj)
c
     $             , torqx(0:maxobj),torqpx(0:maxobj),torqvx(0:maxobj)
     $             , torqy(0:maxobj),torqpy(0:maxobj),torqvy(0:maxobj)
     $             , torqz(0:maxobj),torqpz(0:maxobj),torqvz(0:maxobj)
c
     $             , dpdx_mean,dpdy_mean,dpdz_mean
     $             , dgtq(3,4)

      common /scrcg/ pm1(lx1,ly1,lz1,lelv)  ! Mapped pressure
      real pa(lx1,ly2,lz2), pb(lx1,ly1,lz2) ! Working arrays for mappr

      integer il


!     object definition
      integer iobj_wall, bIDs(1)
      save iobj_wall

!     Point for torque calculation
      real x0(3), wall_area, tauw, tauw2
      data x0 /3*0.0/

      real p_wall, glmax, a_wall

      integer iref_step




      ! regenerate objects
!      call reset_obj()          ! this routine available for AMR only
!      bIDs(1) = 3               ! 'W  '
!      call create_obj(iobj_wall,bIDs,1)
      
!     start framework
      if (ISTEP.eq.0) call frame_start

!     monitor simulation
      call frame_monitor

!     save/load files for full-restart
      call chkpt_main

!     for tripping
      call stat_avg

      call trunc_main()

!     finalise framework
      if (ISTEP.eq.NSTEPS.or.LASTEP.eq.1) then
         call frame_end      
      endif
     
      return
      end
!-----------------------------------------------------------------------
      subroutine userbc (ix,iy,iz,iside,eg)
      include 'SIZE'
      include 'NEKUSE'          ! UX, UY, UZ, TEMP, X, Y

!     velocity
      ux =  0.0
      uy =  0.0
      uz =  0.0

      return
      end
!-----------------------------------------------------------------------
      subroutine useric (ix,iy,iz,ieg)
      implicit none
      
      include 'SIZE'
      include 'INPUT'
      include 'TSTEP'           ! PI
      include 'NEKUSE'          ! UX, UY, UZ, TEMP, Z

!     argument list
      integer ix,iy,iz,ieg

!     local variables
      real C, k, kth, kz,alpha,beta,cons,eps,Ri,rR,xL1,zL1
      real rp,ReTau, th, u_axial, u_th, u_rho, kx, dr
      integer meanProfType

c     Define the initial mean profile
      meanProfType=1;   !1: parabolic velocity profile
                        !2: Reichardt law at a given ReTau

 
c     Parabolic profile ub=1.0
      if (meanProfType.eq.1) then
         r   = sqrt(y*y+x*x)
         Ri   = 1.0      !pipe radius
         rR   = (r*r)/(Ri*Ri)
         cons = 2.0
         uz   = cons*(1-rR)
      else
c        Reichardt function
         ReTau = 360.0
         C      = 5.17
         k      = 0.41

         r   = sqrt(y*y+x*x)
         rp = (1.0-r)*ReTau
         uz  = 1/k*log(1+k*rp) + (C - (1/k)*log(k)) *
     $      (1 - exp(-rp/11.0) - rp/11.0*exp(-rp/3))
         uz  = uz * ReTau*param(2)
      endif

c     perturb
      eps  = 3e-2
      kz   = 23
      kth  = 13


      th = atan2(y,x)


      xL1=(2.0*PI)
      zL1=ZLENPIPE
      alpha = kz * 2*PI/zL1
      beta  = kth * 2*PI/xL1
      dr = 0.2

      u_axial  =      eps*beta  * sin(alpha*z)*cos(beta*th)
      u_th     =     - eps  * alpha     * r * cos(alpha*z)*sin(beta*th)
      if (r>dr) then
         u_rho    =     eps*alpha * sin(alpha*z)*sin(beta*th) * (1/r)
      else
         u_rho   = eps*alpha*sin(alpha*z)*sin(beta*th)*(1/dr)
      endif

      ux = cos(th)*u_rho -r*sin(th)*u_th
      uy = sin(th)*u_rho +r*cos(th)*u_th
      uz = uz + u_axial

      temp=0

      return
      end


!----------------------------------------------------------------------
      subroutine usrdat
      include 'SIZE'

      call setbc(3,1,'W  ')   ! wall

      return
      end
!-----------------------------------------------------------------------
      subroutine usrdat2
      include 'SIZE'


      
      return
      end
!-----------------------------------------------------------------------
      subroutine usrdat3
      implicit none
      include 'SIZE'
      include 'INPUT'

      ! To set up forcing in z-direction
      param(54) = -3 !main flow direction: z-direction
      param(55) = 1  !Ubar
      
      return
      end



c-----------------------------------------------------------------------
!@Subroutine to evaluate the wall pressure. 
! INPUT: BCobj identifies the wall, vect is the pressure vector, n is the number of objects. 
! OUTPUT: pwall is wall pressure and wall area is area of the wall objects
      subroutine pw_eval(BCobj,vect,n,wallpr,wallarea)
      implicit none

      include 'SIZE'
      include 'TOTAL'

      ! Input/Output variables
      real vect(lx1,ly1,lz1,lelv), wallpr, wallarea
      integer BCobj(n), n

      ! Local variables
      integer e, f, i, nn, ntot
      real sint, sarea, pwall, sa

      ! functions
      real vlsum                

!     Array containing the area of each object for checking purposes
!     ( + work array): a and wa are used to check the right area value
      real a, wa

      nn = nx1*nz1
      ntot = nx1*ny1*nz1*nelv 

      if (n.gt.1) then
         write(*,*)'pw_eval sub. works only for 1 OBJECT identifier'
      endif
      

      do e=1,nelv
      do f=1,2*ndim
      do i = 1,n
         if (boundaryID(f,e) .eq. BCobj(n)) then

            ! Alternative way to compute pressure(to check surface_int)
            a = a + vlsum(area(1,1,f,e),nn)   

            call surface_int(sint,sarea,vect,e,f)
            call cadd(pwall,sint,1)
            call cadd(sa,sarea,1)

         endif
      enddo
      enddo
      enddo


      call gop(a, wa, '+  ', maxobj)


      wallarea = sa 
      wallpr   = pwall

      if (nid .eq. 0) then
         write (6,*) 'Surface area of object: ',a, wallarea
         write (6,*) 'Wall pressure: ', wallpr
      endif
 2    format(a2,i1,a3,1es18.10e1)



      return
      end
c-----------------------------------------------------------------------

c-----------------------------------------------------------------------
      subroutine normalize_pres
      implicit none
      include 'SIZE'
      include 'SOLN'
      include 'MASS'

      real pmean, glsc2
      integer nxyz2, ntot2

      nxyz2 = nx2*ny2*nz2
      ntot2 = nxyz2*nelv

      pmean = glsc2(pr,bm2,ntot2)/volvm2
      call cadd(pr,-pmean,ntot2)
c      if (lorder.ge.2) call cadd(prlag,-pmean,ntot2)

      if (nid.eq.0) then
         write(6,*)'Mean wall pressure: ',pmean
      endif
      
      return
      end
c-----------------------------------------------------------------------




!======================================================================
!> @brief Register user specified modules
      subroutine frame_usr_register
      implicit none

      include 'SIZE'
      include 'FRAMELP'
!-----------------------------------------------------------------------
!     register modules
      call io_register
      call chkpt_register
      call stat_register
      call trunc_register

      return
      end subroutine
!======================================================================
!> @brief Initialise user specified modules
      subroutine frame_usr_init
      implicit none

      include 'SIZE'
      include 'FRAMELP'
      include 'SOLN'
!-----------------------------------------------------------------------
!     initialise modules
      call chkpt_init
      call stat_init
      call trunc_init

      return
      end subroutine
!======================================================================
!> @brief Finalise user specified modules
      subroutine frame_usr_end
      implicit none

      include 'SIZE'
      include 'FRAMELP'
!-----------------------------------------------------------------------
!     finalise modules
      
      return
      end subroutine
!======================================================================
!> @brief Provide element coordinates and local numbers (user interface)
!! @param[out]  idir              mapping (uniform) direction
!! @param[out]  ctrs              2D element centres
!! @param[out]  cell              local element numberring
!! @param[in]   lctrs1,lctrs2     array sizes
!! @param[out]  nelsort           number of local 3D elements to sort
!! @param[out]  map_xm1, map_ym1  2D coordinates of mapped elements
!! @param[out]  ierr              error flag
      subroutine user_map2d_get(idir,ctrs,cell,lctrs1,lctrs2,nelsort,
     $     map_xm1,map_ym1,ierr)
      implicit none

      include 'SIZE'
      include 'INPUT'           ! [XYZ]C
      include 'GEOM'            ! [XYZ]M1

!     argument list
      integer idir
      integer lctrs1,lctrs2
      real ctrs(lctrs1,lctrs2)  ! 2D element centres  and diagonals 
      integer cell(lctrs2)      ! local element numberring
      integer nelsort           ! number of local 3D elements to sort
      real map_xm1(lx1,lz1,lelt), map_ym1(lx1,lz1,lelt)
      integer ierr              ! error flag

!     local variables
      integer ntot              ! tmp array size for copying
      integer el ,il ,jl        ! loop indexes
      integer nvert             ! vertex number
      real rnvert               ! 1/nvert
      real xmid,ymid            ! 2D element centre
      real xmin,xmax,ymin,ymax  ! to get approximate element diagonal
      integer ifc               ! face number

!     dummy arrays
      real xcoord(8,LELT), ycoord(8,LELT) ! tmp vertex coordinates

!#define DEBUG
#ifdef DEBUG
!     for testing
      character*3 str1, str2
      integer iunit, ierrl
      ! call number
      integer icalldl
      save icalldl
      data icalldl /0/
#endif

!-----------------------------------------------------------------------
!     initial error flag
      ierr = 0
!     set important parameters
!     uniform direction; should be taken as input parameter
!     x-> 1, y-> 2, z-> 3
      idir = 3
      
!     get element midpoints
!     vertex number
      nvert = 2**NDIM
      rnvert= 1.0/real(nvert)

!     eliminate uniform direction
      ntot = 8*NELV
      if (idir.EQ.1) then  ! uniform X
         call copy(xcoord,YC,ntot) ! copy y
         call copy(ycoord,ZC,ntot) ! copy z
      elseif (idir.EQ.2) then  ! uniform Y
         call copy(xcoord,XC,ntot) ! copy x
         call copy(ycoord,ZC,ntot) ! copy z
      elseif (idir.EQ.3) then  ! uniform Z
         call copy(xcoord,XC,ntot) ! copy x
         call copy(ycoord,YC,ntot) ! copy y
      endif

!     set initial number of elements to sort
      nelsort = 0
      call izero(cell,NELT)

!     for every element
      do el=1,NELV
!     element centre
         xmid = xcoord(1,el)
         ymid = ycoord(1,el)

!     element diagonal
         xmin = xmid
         xmax = xmid
         ymin = ymid
         ymax = ymid
         do il=2,nvert
            xmid=xmid+xcoord(il,el)
            ymid=ymid+ycoord(il,el)
            xmin = min(xmin,xcoord(il,el))
            xmax = max(xmax,xcoord(il,el))
            ymin = min(ymin,ycoord(il,el))
            ymax = max(ymax,ycoord(il,el))
         enddo
         xmid = xmid*rnvert
         ymid = ymid*rnvert

!     count elements to sort
            nelsort = nelsort + 1
!     2D position
!     in general this coud involve some curvilinear transform
            ctrs(1,nelsort)=xmid
            ctrs(2,nelsort)=ymid
!     reference distance
            ctrs(3,nelsort)=sqrt((xmax-xmin)**2 + (ymax-ymin)**2)
            if (ctrs(3,nelsort).eq.0.0) then
               ierr = 1
               return
            endif
!     element index
            cell(nelsort) = el
      enddo

!     provide 2D mesh
!     in general this coud involve some curvilinear transform
      if (idir.EQ.1) then  ! uniform X
         ifc = 4
         do el=1,NELV
            call ftovec(map_xm1(1,1,el),ym1,el,ifc,nx1,ny1,nz1)
            call ftovec(map_ym1(1,1,el),zm1,el,ifc,nx1,ny1,nz1)
         enddo
      elseif (idir.eq.2) then  ! uniform y
         ifc = 1
         do el=1,nelv
            call ftovec(map_xm1(1,1,el),xm1,el,ifc,nx1,ny1,nz1)
            call ftovec(map_ym1(1,1,el),zm1,el,ifc,nx1,ny1,nz1)
         enddo
      elseif (idir.eq.3) then  ! uniform z
         ifc = 5
         do el=1,nelv
            call ftovec(map_xm1(1,1,el),xm1,el,ifc,nx1,ny1,nz1)
            call ftovec(map_ym1(1,1,el),ym1,el,ifc,nx1,ny1,nz1)
         enddo
      endif

#ifdef DEBUG
!     testing
      ! to output refinement
      icalldl = icalldl+1
      call io_file_freeid(iunit, ierrl)
      write(str1,'(i3.3)') NID
      write(str2,'(i3.3)') icalldl
      open(unit=iunit,file='map2d_usr.txt'//str1//'i'//str2)
      
      write(iunit,*) idir, NELV, nelsort
      write(iunit,*) 'Centre coordinates and cells'
      do el=1,nelsort
         write(iunit,*) el, ctrs(:,el), cell(el)
      enddo
      write(iunit,*) 'GLL coordinates'
      do el=1,nelsort
         write(iunit,*) 'Element ', el
         write(iunit,*) 'XM1'
         do il=1,nz1
            write(iunit,*) (map_xm1(jl,il,el),jl=1,nx1)
         enddo
         write(iunit,*) 'YM1'
         do il=1,nz1
            write(iunit,*) (map_ym1(jl,il,el),jl=1,nx1)
         enddo
      enddo
      close(iunit)
#endif

      return
      end subroutine
!=======================================================================
!> @brief Provide velocity, deriv. and vort. in required coordinates and normalise pressure
!! @param[out]   lvel             velocity
!! @param[out]   dudx,dvdx,dwdx   velocity derivatives
!! @param[out]   vort             vorticity
!! @param[inout] pres             pressure
      subroutine user_stat_trnsv(lvel,dudx,dvdx,dwdx,vort,pres)
      implicit none

      include 'SIZE'
      include 'SOLN'
      include 'INPUT'               ! if3d
      include 'GEOM'

      ! argument list
      real lvel(LX1,LY1,LZ1,LELT,3) ! velocity array
      real dudx(LX1,LY1,LZ1,LELT,3) ! velocity derivatives; U
      real dvdx(LX1,LY1,LZ1,LELT,3) ! V
      real dwdx(LX1,LY1,LZ1,LELT,3) ! W
      real vort(LX1,LY1,LZ1,LELT,3) ! vorticity
      real pres(LX1,LY1,LZ1,LELT)   ! pressure

      ! local variables
      integer itmp              ! dummy variable
      integer il, jl            ! loop index
      integer ifll              ! field number for object definition
      real vrtmp(lx1*lz1)       ! work array for face
      real vrtmp2(2)            ! work array
      real pmeanl, obj_srfl
 
      ! functions
      real vlsum
!-----------------------------------------------------------------------
      ! Velocity transformation; simple copy
      itmp = NX1*NY1*NZ1*NELV
      call copy(lvel(1,1,1,1,1),VX,itmp)
      call copy(lvel(1,1,1,1,2),VY,itmp)
      call copy(lvel(1,1,1,1,3),VZ,itmp)

      ! Derivative transformation
      ! No transformation
      call gradm1(dudx(1,1,1,1,1),dudx(1,1,1,1,2),dudx(1,1,1,1,3),
     $      lvel(1,1,1,1,1))
      call gradm1(dvdx(1,1,1,1,1),dvdx(1,1,1,1,2),dvdx(1,1,1,1,3),
     $      lvel(1,1,1,1,2))
      call gradm1(dwdx(1,1,1,1,1),dwdx(1,1,1,1,2),dwdx(1,1,1,1,3),
     $      lvel(1,1,1,1,3))

      ! get vorticity
      if (IF3D) then
         ! curlx
         call sub3(vort(1,1,1,1,1),dwdx(1,1,1,1,2),
     $        dvdx(1,1,1,1,3),itmp)
         ! curly
         call sub3(vort(1,1,1,1,2),dudx(1,1,1,1,3),
     $        dwdx(1,1,1,1,1),itmp)
      endif
      ! curlz
      call sub3(vort(1,1,1,1,3),dvdx(1,1,1,1,1),dudx(1,1,1,1,2),itmp)

      ! normalise pressure
      ! in this example I integrate pressure over all faces marked "W"
      ifll = 1     ! I'm interested in velocity bc
      call rzero(vrtmp2,2)  ! zero work array
      obj_srfl = 0.0 ! initialise object surface
      pmeanl = 0.0  ! initialise pressure mean
      itmp = LX1*LZ1
      do il=1,nelv   ! element loop
         do jl=1,2*ldim
            if (cbc(jl,il,ifll).eq.'W  ') then
               vrtmp2(1) = vrtmp2(1) + vlsum(area(1,1,jl,il),itmp)
               call ftovec(vrtmp,pres,il,jl,lx1,ly1,lz1)
               call col2(vrtmp,area(1,1,jl,il),itmp)
               vrtmp2(2) = vrtmp2(2) + vlsum(vrtmp,itmp)
            endif
         enddo
      enddo
      ! global communication
      call gop(vrtmp2,vrtmp,'+  ',2)
      ! missing error chack vrtmp2(1) == 0
      vrtmp2(2) = -vrtmp2(2)/vrtmp2(1)
!     DEBUG
!      write(*,*)'vrtmp2(1)',vrtmp2(1)      
      ! subtract mean pressure
      itmp = NX1*NY1*NZ1*NELV
      call cadd(pres,vrtmp2(2),itmp)
      
      return
      end subroutine
!======================================================================

c automatically added by makenek
      subroutine usrdat0() 

      return
      end

c automatically added by makenek
      subroutine usrsetvert(glo_num,nel,nx,ny,nz) ! to modify glo_num
      integer*8 glo_num(1)

      return
      end

c automatically added by makenek
      subroutine userqtl

      call userqtl_scig

      return
      end
