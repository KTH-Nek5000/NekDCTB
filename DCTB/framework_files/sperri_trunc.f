!> @file sperri.f
!! @ingroup nekamr
!! @brief Spectral error indicator.
!! @details Set of subroutines calculating error indicator based on
!! variable spectra. Adopted from Catherine Mavriplis code.
!! @author Adam Peplinski
!! @author Nicolas Offermans
!! @date Jun 20, 2016
!=======================================================================
!> @brief Get refinement indicator for single variable; main interface
!! @param[out]   eind      error indicator
!! @param[out]   sig       averaged slope
!! @param[in]    var       tested variable
!! @param[in]    lnelt     local number of elements
!! @param[in]    var_name  variable name
!! @note This routine accepts m1 variables only
      subroutine speri_get(eind,sig,var,lnelt)
      implicit none

      include 'SIZE'

      ! argument list
      real eind(lnelt),sig(lnelt)
      real var(LX1,LY1,LZ1,lnelt)
      integer lnelt

      ! local variables
      ! work arrays
      real xa(LX1,LY1,LZ1), xb(LX1,LY1,LZ1)

      ! initalisation
      logical ifcalled
      save ifcalled
      data ifcalled /.FALSE./
!-----------------------------------------------------------------------
!     ! this is not the best way of initialisation, but for now I leave it here
      if(.not.ifcalled) then
        ifcalled=.TRUE.
        ! initialise error estimator
        call speri_init
      endif

      ! zero arrays
      call rzero(eind,lnelt)
      call rzero(sig,lnelt)

      call speri_var(eind,sig,var,lnelt,xa,xb)

      return
      end subroutine
!=======================================================================
!> @brief Initialise error estimator
      subroutine speri_init
      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'SPERRID'

      ! local variables
      integer il, jl, aa
!-----------------------------------------------------------------------
      ! stamp logs
      if (nio.eq.0) write(*,*)'Spectral error indicator init.'

      ! set cutoff parameters
      ! used for values
      SERI_SMALL = 1.e-14
      ! used for ratios
      SERI_SMALLR = 1.e-10
      ! used for gradients
      SERI_SMALLG = 1.e-5
      ! used for sigma and rtmp in error calculations
      SERI_SMALLS = 0.2

      ! number of points in fitting
      SERI_NP = 4
      ! last modes skipped
      SERI_ELR = 0

      ! correctness check
      if (SERI_NP.gt.SERI_NP_MAX) then
        if (nio.eq.0) write(*,*) 'SETI_NP greater than SERI_NP_MAX' 
      endif
      il = SERI_NP+SERI_ELR
      jl = min(LX1,LY1)
      if (IF3D) jl = min(jl,LZ1)
      if (il.gt.jl) then
        if (nio.eq.0) write(*,*) 'SERI_NP+SERI_ELR greater than L?1'
      endif

      ! initalise coefficient mapping
      call speri_set_map

      return
      end subroutine
!=======================================================================
!> @brief Initialise Legendre transforms for computing spectral coef.
!         and factors for computing error estimators.
      subroutine speri_set_map
      implicit none
      
      include 'SIZE'
      include 'INPUT'   ! if3d
      include 'SPERRID' ! SERI_Lj, SERI_Ljt, SERI_NP, SERI_ELR
      include 'WZ'

      real fack, facj, faci
      real work(lx1)
      integer il,jl,kl
!-----------------------------------------------------------------------
      ! check polynomial order and numer of points for extrapolation
      if (min(NX1,NY1).le.(SERI_NP+SERI_ELR)) then
        if (nio.eq.0) write(*,*) 'Error: increase L[XYZ]1'
      endif

      if (IF3D.and.(NZ1.le.(SERI_NP+SERI_ELR))) then
        if (nio.eq.0) write(*,*) 'Error: increase L[XYZ]1'
      endif

      ! Legendre transforms
      call speri_legend_transform(SERI_Lj,SERI_Ljt,nx1,work)

      ! multiplicity factors
      do kl=1,NZ1
        if (NZ1 .eq. 1) then
           fack = 1.0
        elseif (kl .eq. NZ1) then
           fack = (real(NZ1)-1.0)/2.0
        else
           fack = (2*(real(NZ1)-1.0)+1.0)/2.0
        endif
        do jl=1,NY1
           if (jl .eq. NY1) then
              facj = (real(NY1)-1.0)/2.0 
           else
              facj = (2*(real(NY1)-1.0)+1.0)/2.0
           endif
           do il=1,NX1
              if (il .eq. NX1) then
                 faci = (real(NX1)-1.0)/2.0 
              else
                 faci = (2*(real(NX1)-1.0)+1.0)/2.0
              endif
              SERI_FAC(il,jl,kl) = fack*facj*faci
           enddo
        enddo
      enddo

      return
      end subroutine
!=======================================================================
!> @brief Error and sigma for single variable
!! @details Get error estimater and sigma for single variable on a whole
!!  mesh for all directions.
!! @param[out]  est   estimated error
!! @param[out]  sig   estimated exponent
!! @param[in]   var   tested variable
!! @param[in]   nell  element number
!! @param[in]   xa    work array
!! @param[in]   xb    work array
      subroutine speri_var(est,sig,var,nell,xa,xb)
      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'SPERRID'

      ! argument list
      real est(nell), sig(nell)
      real var(LX1,LY1,LZ1,nell)
      integer nell
      real xa(LX1,LY1,LZ1), xb(LX1,LY1,LZ1)

      ! local variables
      integer il, jl, kl, ll, j_st, j_en
      ! polynomial coefficients
      real coeff(LX1,LY1,LZ1)
      ! Legendre coefficients; first value coeff(1,1,1)
      real coef11
      ! copy of last SERI_NP columns of coefficients
      real coefx(SERI_NP_MAX,LY1,LZ1),coefy(SERI_NP_MAX,LX1,LZ1),
     $     coefz(SERI_NP_MAX,LX1,LY1)
      ! estimated error
      real estx, esty, estz
      ! estimated decay rate
      real sigx, sigy, sigz
      real third
      parameter (third = 1.0/3.0)
!-----------------------------------------------------------------------
      ! loop over elements
      do il = 1,nell
        ! go to Legendre space (done in two operations)
        call tensr3(coeff,nx1,var(1,1,1,il),nx1,SERI_Lj,SERI_Ljt,
     $              SERI_Ljt,xa)
        ! square the coefficients
        call vsq(coeff,nx1*ny1*nz1)

        ! lower left corner
        coef11 = coeff(1,1,1)

        ! small value; nothing to od
        if (coef11.ge.SERI_SMALL) then
           ! extrapolate coefficients
           ! X - direction
           ! copy last SERI_NP collumns (or less if NX1 is smaller)
           ! SERI_ELR allows to exclude last row
            j_st = max(1,NX1-SERI_NP+1-SERI_ELR)
            j_en = max(1,NX1-SERI_ELR)
            do ll = 1,NZ1
                do kl = 1,NY1
                    do jl = j_st,j_en
                        coefx(j_en-jl+1,kl,ll) = coeff(jl,kl,ll)
                    enddo
                enddo
            enddo
            ! get extrapolated values
            call speri_extrap(estx,sigx,coef11,coefx,
     $           j_st,j_en,NY1,NZ1)

            ! Y - direction
            ! copy last SERI_NP collumns (or less if NY1 is smaller)
            ! SERI_ELR allows to exclude last row
            j_st = max(1,NY1-SERI_NP+1-SERI_ELR)
            j_en = max(1,NY1-SERI_ELR)
            do ll = 1,NZ1
                do kl = j_st,j_en
                    do jl = 1,NX1
                        coefy(j_en-kl+1,jl,ll) = coeff(jl,kl,ll)
                    enddo
                enddo
            enddo

            ! get extrapolated values
            call speri_extrap(esty,sigy,coef11,coefy,
     $           j_st,j_en,NX1,NZ1)

            if (IF3D) then
                ! Z - direction
                ! copy last SERI_NP collumns (or less if NZ1 is smaller)
                ! SERI_ELR allows to exclude last row
                j_st = max(1,NZ1-SERI_NP+1-SERI_ELR)
                j_en = max(1,NZ1-SERI_ELR)
                do ll = j_st,j_en
                    do kl = 1,NY1
                        do jl = 1,NX1
                            coefz(j_en-ll+1,jl,kl) = coeff(jl,kl,ll)
                        enddo
                    enddo
                enddo

                ! get extrapolated values
                call speri_extrap(estz,sigz,coef11,coefz,
     $              j_st,j_en,NX1,NY1)

                ! average
                est(il) =  sqrt(estx + esty + estz)
                sig(il) =  third*(sigx + sigy + sigz)
            else
                est(il) =  sqrt(estx + esty)
                sig(il) =  0.5*(sigx + sigy)
            endif

        else
            ! for testing
            estx = 0.0
            esty = 0.0
            estz = 0.0
            sigx = -1.0
            sigy = -1.0
            sigz = -1.0
            ! for testing; end

            est(il) =  0.0
            sig(il) = -1.0
        endif

      enddo

      return
      end subroutine
!=======================================================================
!> @brief Get extrapolated values of sigma and error estimator.
!! @details We assume coef(n) = c*exp(-sigma*n) and estiamte sigma and
!!  eror eest = sqrt(2*(coef(N)**2/(2*N+1+\int_(N+1)^\infty coef(n)**2/(n+1) dn)))
!! @param[out]    estx   estimated error
!! @param[out]    sigx   estimated exponent
!! @param[in]     coef11 Legendre coefficients; first base-mode amplitude (coeff(1,1,1))
!! @param[in]     coef   Legendre coefficients; last SERI_NP columns
!! @param[in]     ix_st  first mode in coef
!! @param[in]     ix_en  last mode in coef
!! @param[in]     nyl    number of modes in y direction (aray size)
!! @param[in]     nzl    number of modes in z direction (aray size)
      subroutine speri_extrap(estx,sigx,coef11,coef,
     $           ix_st,ix_en,nyl,nzl)
      implicit none

      include 'SIZE'
      include 'SPERRID'

      ! argument list
      integer ix_st,ix_en,nyl,nzl
      ! Legendre coefficients; last SERI_NP columns
      real coef(SERI_NP_MAX,nyl,nzl)
      ! Legendre coefficients; first value coeff(1,1,1)
      real coef11
      ! estimated error and decay rate
      real estx, sigx

      ! local variables
      integer il, jl, kl, ll  ! loop index
      integer nsigt, pnr, nzlt
      real sigt, smallr, cmin, cmax, cnm, rtmp, rtmp2, rtmp3
      real sumtmp(4), cffl(SERI_NP_MAX)
      real stmp, estt, clog, ctmp, cave, erlog
      logical cuse(SERI_NP_MAX)
!-----------------------------------------------------------------------
      ! initial values
      estx =  0.0
      sigx = -1.0

      ! relative cutoff
      smallr = coef11*SERI_SMALLR

      ! number of points
      pnr = ix_en - ix_st +1

      ! to few points to interpolate
!      if ((ix_en - ix_st).le.1) return

      ! for averaging, initial values
      sigt = 0.0
      nsigt = 0

      ! loop over all face points
      nzlt = max(1,nzl - SERI_ELR) !  for 2D runs
      do il=1,nzlt
        ! weight
        rtmp3 = 1.0/(2.0*(il-1)+1.0)
        do jl=1,nyl - SERI_ELR

            ! find min and max coef along single row
            cffl(1) = coef(1,jl,il)
            cmin = cffl(1)
            cmax = cmin
            do kl =2,pnr
                cffl(kl) = coef(kl,jl,il)
                cmin = min(cmin,cffl(kl))
                cmax = max(cmax,cffl(kl))
            enddo

            ! are coefficients sufficiently big
            if((cmin.gt.0.0).and.(cmax.gt.smallr)) then
                ! mark array position we use in iterpolation
                do kl =1,pnr
                    cuse(kl) = .TRUE.
                enddo
                ! max n for polynomial order
                cnm = real(ix_en)

                ! check if all the points should be taken into account
                ! in original code by Catherine Mavriplis this part is written
                ! for 4 points, so I place if statement first
                if (pnr.eq.4) then
                    ! should we neglect last values
                    if ((cffl(1).lt.smallr).and.
     &                  (cffl(2).lt.smallr)) then
                        if (cffl(3).lt.smallr) then
                            cuse(1) = .FALSE.
                            cuse(2) = .FALSE.
                            cnm = real(ix_en-2)
                        else
                            cuse(1) = .FALSE.
                            cnm = real(ix_en-1)
                        endif
                    else
                        ! should we take stronger gradient
                        if ((cffl(1)/cffl(2).lt.SERI_SMALLG).and.
     $                      (cffl(3)/cffl(4).lt.SERI_SMALLG)) then
                            cuse(1) = .FALSE.
                            cuse(3) = .FALSE.
                            cnm = real(ix_en-1)
                        elseif ((cffl(2)/cffl(1).lt.SERI_SMALLG).and.
     $                          (cffl(4)/cffl(3).lt.SERI_SMALLG)) then
                            cuse(2) = .FALSE.
                            cuse(4) = .FALSE.
                        endif
                    endif
                endif

                ! get sigma for given face point
                do kl =1,4
                    sumtmp(kl) = 0.0
                enddo
                ! find new min and count number of points
                cmin = cmax
                cmax = 0.0
                do kl =1,pnr
                    if(cuse(kl)) then
                        rtmp  = real(ix_en-kl)
                        rtmp2 = log(cffl(kl))
                        sumtmp(1) = sumtmp(1) +rtmp2
                        sumtmp(2) = sumtmp(2) +rtmp
                        sumtmp(3) = sumtmp(3) +rtmp*rtmp
                        sumtmp(4) = sumtmp(4) +rtmp2*rtmp
                        ! find new min and count used points
                        cmin = min(cffl(kl),cmin)
                        cmax = cmax + 1.0
                    endif
                enddo
                ! decay rate along single row
                stmp = (sumtmp(1)*sumtmp(2) - sumtmp(4)*cmax)/
     $                 (sumtmp(3)*cmax - sumtmp(2)*sumtmp(2))
                ! for averaging
                sigt = sigt + stmp
                nsigt = nsigt + 1

                ! get error estimator depending on calculated decay rate
                estt = 0.0
                if (stmp.lt.SERI_SMALLS) then
                    estt = cmin
                else
                    ! get averaged constant in front of c*exp(-sig*n)
                    clog = (sumtmp(1)+stmp*sumtmp(2))/cmax
                    ctmp = exp(clog)
                    ! average exponent
                    cave = sumtmp(1)/cmax
                    ! check quality of approximation comparing is to the constant cave
                    do kl =1,2
                        sumtmp(kl) = 0.0
                    enddo
                    do kl =1,pnr
                        if(cuse(kl)) then
                            erlog = clog - stmp*real(ix_en-kl)
                            sumtmp(1) = sumtmp(1)+
     $                          (erlog-log(cffl(kl)))**2
                            sumtmp(2) = sumtmp(2)+
     $                          (erlog-cave)**2
                        endif
                    enddo
                    rtmp = 1.0 - sumtmp(1)/sumtmp(2)
                    if (rtmp.lt.SERI_SMALLS) then
                        estt = cmin
                    else
                        ! last coefficient is not included in error estimator
                        estt = ctmp/stmp*exp(-stmp*cnm)
                    endif
                endif
                ! add contribution to error estimator; variable weight
                estx = estx + estt/(2.0*(jl-1)+1.0)*rtmp3
            endif  ! if((cmin.gt.0.0).and.(cmax.gt.smallr))
        enddo
      enddo
      ! constant weight
      ! Multiplication by 4 in 2D / 8 in 3D
      ! Normalization of the error by the volume of the reference element
      ! which is equal to 4 in 2D / 8 in 3D
      ! ==> Both operations cancel each other
      estx = estx/(2.0*(ix_en-1)+1.0)

      ! final everaging
      ! sigt = 2*sigma so we divide by 2
      if (nsigt.gt.0) then
        sigx = 0.5*sigt/nsigt
      endif

      return
      end subroutine
!=======================================================================
!> @brief Build matrices for Legendre transform
!! @details Lj(:,j) = legendre_poly(zgll(j)) * wgll(j)
!! @param[out]    Lj     Legendre transform
!! @param[out]    Ljt    transposed Legendre transform
!! @param[in]     nx     number of ponints
!! @param[in]     pleg   work array
      subroutine speri_legend_transform(Lj,Ljt,nx,pleg)
      implicit none
     
      include 'SIZE'
      include 'WZ'

      real    Lj(lx1,lx1),Ljt(lx1,lx1),pleg(lx1)
      integer i,j,n,nx

      n = nx-1
      do j=1,nx
         call legendre_poly(pleg,zgm1(j,1),n)  ! Return Lk(z), k=0,...,n
         do i=1,nx-1
            Lj(i,j) = pleg(i)*wxm1(j)*(2*(i-1)+1)/2
         enddo
         Lj(nx,j) = pleg(i)*wxm1(j)*(nx-1)/2
      enddo
      call transpose (Ljt,nx,Lj,nx)

      return
      end subroutine
!=======================================================================
