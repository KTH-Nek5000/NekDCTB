!> @file amr_refine.f
!! @ingroup nekamr
!! @brief Set of routines to perform interpolation for refinement/coarsening.
!! @author Adam Peplinski
!! @date Jun 07, 2016
!=======================================================================
!> @brief Get interpolation parameters for octal interpolation
      subroutine amr_genwz
      implicit none

      include 'SIZE'
      include 'INPUT'           ! IF3D, IFSPLIT
      include 'WZ'              ! ZGM(123)
      include 'AMR'
      include 'AMR_REFINE'

      ! local variables
      integer il, jl, kl     !     loop indexes
      integer nt2

      ! tmp  array to calculate interpolation positions
      integer nmax
      parameter (nmax=84) ! like in speclib
      real tmpl(nmax)

      ! for crs interpolation
      real zgmc(AMR_LXC), wgmc(AMR_LXC)

#ifdef DEBUG
      ! for testing
      character*2 str
      integer iunit, ierr
#endif
!-----------------------------------------------------------------------
      !zero arrays
      ! Volume interpolation
      nt2 = LX1*LX1*2
      call rzero(IXAMR1CF,nt2)
      call rzero(IXTAMR1CF,nt2)
      call rzero(IXAMR1FC,nt2)
      call rzero(IXTAMR1FC,nt2)
      nt2 = LX2*LX2*2
      call rzero(IXAMR2CF,nt2)
      call rzero(IXTAMR2CF,nt2)
      call rzero(IXAMR2FC,nt2)
      call rzero(IXTAMR2FC,nt2)
      nt2 = LY1*LY1*2
      call rzero(IYAMR1CF,nt2)
      call rzero(IYTAMR1CF,nt2)
      call rzero(IYAMR1FC,nt2)
      call rzero(IYTAMR1FC,nt2)
      nt2 = LY2*LY2*2
      call rzero(IYAMR2CF,nt2)
      call rzero(IYTAMR2CF,nt2)
      call rzero(IYAMR2FC,nt2)
      call rzero(IYTAMR2FC,nt2)
      if (IF3D) then
         nt2 = LZ1*LZ1*2
         call rzero(IZAMR1CF,nt2)
         call rzero(IZTAMR1CF,nt2)
         call rzero(IZAMR1FC,nt2)
         call rzero(IZTAMR1FC,nt2)
         nt2 = LZ2*LZ2*2
         call rzero(IZAMR2CF,nt2)
         call rzero(IZTAMR2CF,nt2)
         call rzero(IZAMR2FC,nt2)
         call rzero(IZTAMR2FC,nt2)
      elseif(IFAXIS) then       ! axisymmetric case
         nt2 = LY1*LY1*2
         call rzero(IAAMR1CF,nt2)
         call rzero(IATAMR1CF,nt2)
         call rzero(IAAMR1FC,nt2)
         call rzero(IATAMR1FC,nt2)
         nt2 = LY2*LY2*2
         call rzero(IAAMR2CF,nt2)
         call rzero(IATAMR2CF,nt2)
         call rzero(IAAMR2FC,nt2)
         call rzero(IATAMR2FC,nt2)
      endif
      ! crs (vertex) interpolation
      nt2 = AMR_LXC*AMR_LXC*2
      call rzero(IXAMRCCF,nt2)
      call rzero(IXTAMRCCF,nt2)
      ! multiplicity arrays
      nt2 = LX1*LY1*LZ1
      call rone(IMAMR1,nt2)
      nt2 = LX2*LY2*LZ2
      call rone(IMAMR2,nt2)

      ! get interpolation operators
      ! Volume interpolation
      ! MESH M1
      ! coarse -> fine
      ! X
      ! negative
      do il=1,NX1
         tmpl(il) = 0.5*(ZGM1(il,1) -1.0)
      enddo
      call IGLLM (IXAMR1CF,IXTAMR1CF,ZGM1(1,1),tmpl,NX1,NX1,NX1,NX1)
      ! positive; we use symmetry
      do jl=1,NX1
         do il=1,NX1
            IXAMR1CF(NX1-il+1,NX1-jl+1,2)  = IXAMR1CF(il,jl,1)
            IXTAMR1CF(NX1-il+1,NX1-jl+1,2) = IXTAMR1CF(il,jl,1)
         enddo
      enddo

      ! Y
      ! negative
      do il=1,NY1
         tmpl(il) = 0.5*(ZGM1(il,2) -1.0)
      enddo
      call IGLLM (IYAMR1CF,IYTAMR1CF,ZGM1(1,2),tmpl,NY1,NY1,NY1,NY1)
      ! positive; we use symmetry
      do jl=1,NY1
         do il=1,NY1
            IYAMR1CF(NY1-il+1,NY1-jl+1,2)  = IYAMR1CF(il,jl,1)
            IYTAMR1CF(NY1-il+1,NY1-jl+1,2) = IYTAMR1CF(il,jl,1)
         enddo
      enddo

      if (IF3D) then
      ! Z
      ! negative
         do il=1,NZ1
            tmpl(il) = 0.5*(ZGM1(il,3) -1.0)
         enddo
         call IGLLM (IZAMR1CF,IZTAMR1CF,ZGM1(1,3),tmpl,NZ1,NZ1,NZ1,NZ1)
      ! positive; we use symmetry
         do jl=1,NZ1
            do il=1,NZ1
               IZAMR1CF(NZ1-il+1,NZ1-jl+1,2)  = IZAMR1CF(il,jl,1)
               IZTAMR1CF(NZ1-il+1,NZ1-jl+1,2) = IZTAMR1CF(il,jl,1)
            enddo
         enddo
      else
         IZAMR1CF(NZ1,NZ1,1) = 1.0
         IZTAMR1CF(NZ1,NZ1,1) = 1.0
         IZAMR1CF(NZ1,NZ1,2) = 1.0
         IZTAMR1CF(NZ1,NZ1,2) = 1.0
      endif

!     fine -> coarse
!     X
!     negative
      nt2 = NX1/2 + mod(NX1,2)
      do il=1,nt2
         tmpl(il) = 2.0*ZGM1(il,1) + 1.0
      enddo
      call IGLLM (IXAMR1FC,IXTAMR1FC,ZGM1(1,1),tmpl,NX1,nt2,NX1,NX1)
!     positive; we use symmetry
      do jl=1,NX1
         do il=1,nt2
            IXAMR1FC(NX1-il+1,NX1-jl+1,2)  = IXAMR1FC(il,jl,1)
            IXTAMR1FC(NX1-jl+1,NX1-il+1,2) = IXTAMR1FC(jl,il,1)
         enddo
      enddo

!     Y
!     negative
      nt2 = NY1/2 + mod(NY1,2)
      do il=1,nt2
         tmpl(il) = 2.0*ZGM1(il,2) + 1.0
      enddo
      call IGLLM (IYAMR1FC,IYTAMR1FC,ZGM1(1,2),tmpl,NY1,nt2,NY1,NY1)
!     positive; we use symmetry
      do jl=1,NY1
         do il=1,nt2
            IYAMR1FC(NY1-il+1,NY1-jl+1,2)  = IYAMR1FC(il,jl,1)
            IYTAMR1FC(NY1-jl+1,NY1-il+1,2) = IYTAMR1FC(jl,il,1)
         enddo
      enddo

      if (IF3D) then
!     Z
!     negative
         nt2 = NZ1/2 + mod(NZ1,2)
         do il=1,nt2
            tmpl(il) = 2.0*ZGM1(il,3) + 1.0
         enddo
         call IGLLM (IZAMR1FC,IZTAMR1FC,ZGM1(1,3),tmpl,NZ1,nt2,NZ1,NZ1)
!     positive; we use symmetry
         do jl=1,NZ1
            do il=1,nt2
               IZAMR1FC(NZ1-il+1,NZ1-jl+1,2)  = IZAMR1FC(il,jl,1)
               IZTAMR1FC(NZ1-jl+1,NZ1-il+1,2) = IZTAMR1FC(jl,il,1)
            enddo
         enddo
      else
         nt2 = NZ1/2 + mod(NZ1,2)
         IZAMR1FC(nt2,NZ1,1) = 1.0
         IZTAMR1FC(NZ1,nt2,1) = 1.0
         IZAMR1FC(nt2,NZ1,2) = 1.0
         IZTAMR1FC(NZ1,nt2,2) = 1.0
      endif

!     MESH M2
!     coarse -> fine
!     X
!     negative
      do il=1,NX2
         tmpl(il) = 0.5*(ZGM2(il,1) -1.0)
      enddo
      if (IFSPLIT) then         ! P-N-P_N
         call IGLLM (IXAMR2CF,IXTAMR2CF,ZGM2(1,1),tmpl,NX2,NX2,NX2,NX2)
      else                      ! P-N-P_N-2
         call IGLM (IXAMR2CF,IXTAMR2CF,ZGM2(1,1),tmpl,NX2,NX2,NX2,NX2)
      endif

!     positive; we use symmetry
      do jl=1,NX2
         do il=1,NX2
            IXAMR2CF(NX2-il+1,NX2-jl+1,2)  = IXAMR2CF(il,jl,1)
            IXTAMR2CF(NX2-il+1,NX2-jl+1,2) = IXTAMR2CF(il,jl,1)
         enddo
      enddo

!     Y
!     negative
      do il=1,NY2
         tmpl(il) = 0.5*(ZGM2(il,2) -1.0)
      enddo
      if (IFSPLIT) then         ! P-N-P_N
         call IGLLM (IYAMR2CF,IYTAMR2CF,ZGM2(1,2),tmpl,NY2,NY2,NY2,NY2)
      else                      ! P-N-P_N-2
         call IGLM (IYAMR2CF,IYTAMR2CF,ZGM2(1,2),tmpl,NY2,NY2,NY2,NY2)
      endif
!     positive; we use symmetry
      do jl=1,NY2
         do il=1,NY2
            IYAMR2CF(NY2-il+1,NY2-jl+1,2)  = IYAMR2CF(il,jl,1)
            IYTAMR2CF(NY2-il+1,NY2-jl+1,2) = IYTAMR2CF(il,jl,1)
         enddo
      enddo

      if (IF3D) then
!     Z
!     negative
         do il=1,NZ2
            tmpl(il) = 0.5*(ZGM2(il,3) -1.0)
         enddo
         if (IFSPLIT) then      ! P-N-P_N
            call IGLLM (IZAMR2CF,IZTAMR2CF,ZGM2(1,3),tmpl,
     $           NZ2,NZ2,NZ2,NZ2)
         else                   ! P-N-P_N-2
            call IGLM (IZAMR2CF,IZTAMR2CF,ZGM2(1,3),tmpl,
     $           NZ2,NZ2,NZ2,NZ2)
         endif
!     positive; we use symmetry
         do jl=1,NZ2
            do il=1,NZ2
               IZAMR2CF(NZ2-il+1,NZ2-jl+1,2)  = IZAMR2CF(il,jl,1)
               IZTAMR2CF(NZ2-il+1,NZ2-jl+1,2) = IZTAMR2CF(il,jl,1)
            enddo
         enddo
      else
         IZAMR2CF(NZ2,NZ2,1) = 1.0
         IZTAMR2CF(NZ2,NZ2,1) = 1.0
         IZAMR2CF(NZ2,NZ2,2) = 1.0
         IZTAMR2CF(NZ2,NZ2,2) = 1.0
      endif

!     fine -> coarse
!     X
!     negative
      nt2 = NX2/2 + mod(NX2,2)
      do il=1,nt2
         tmpl(il) = 2.0*ZGM2(il,1) + 1.0
      enddo
      if (IFSPLIT) then         ! P-N-P_N
         call IGLLM (IXAMR2FC,IXTAMR2FC,ZGM2(1,1),tmpl,NX2,nt2,NX2,NX2)
      else                      ! P-N-P_N-2
         call IGLM (IXAMR2FC,IXTAMR2FC,ZGM2(1,1),tmpl,NX2,nt2,NX2,NX2)
      endif
!     positive; we use symmetry
      do jl=1,NX2
         do il=1,nt2
            IXAMR2FC(NX2-il+1,NX2-jl+1,2)  = IXAMR2FC(il,jl,1)
            IXTAMR2FC(NX2-jl+1,NX2-il+1,2) = IXTAMR2FC(jl,il,1)
         enddo
      enddo

!     Y
!     negative
      nt2 = NY2/2 + mod(NY2,2)
      do il=1,nt2
         tmpl(il) = 2.0*ZGM2(il,2) + 1.0
      enddo
      if (IFSPLIT) then         ! P-N-P_N
         call IGLLM (IYAMR2FC,IYTAMR2FC,ZGM2(1,2),tmpl,NY2,nt2,NY2,NY2)
      else                      ! P-N-P_N-2
         call IGLM (IYAMR2FC,IYTAMR2FC,ZGM2(1,2),tmpl,NY2,nt2,NY2,NY2)
      endif
!     positive; we use symmetry
      do jl=1,NY2
         do il=1,nt2
            IYAMR2FC(NY2-il+1,NY2-jl+1,2)  = IYAMR2FC(il,jl,1)
            IYTAMR2FC(NY2-jl+1,NY2-il+1,2) = IYTAMR2FC(jl,il,1)
         enddo
      enddo

      if (IF3D) then
!     Z
!     negative
         nt2 = NZ2/2 + mod(NZ2,2)
         do il=1,nt2
            tmpl(il) = 2.0*ZGM2(il,3) + 1.0
         enddo
         if (IFSPLIT) then      ! P-N-P_N
            call IGLLM (IZAMR2FC,IZTAMR2FC,ZGM2(1,3),tmpl,
     $           NZ2,nt2,NZ2,NZ2)
         else                   ! P-N-P_N-2
            call IGLM (IZAMR2FC,IZTAMR2FC,ZGM2(1,3),tmpl,
     $           NZ2,nt2,NZ2,NZ2)
         endif
!     positive; we use symmetry
         do jl=1,NZ2
            do il=1,nt2
               IZAMR2FC(NZ2-il+1,NZ2-jl+1,2)  = IZAMR2FC(il,jl,1)
               IZTAMR2FC(NZ2-jl+1,NZ2-il+1,2) = IZTAMR2FC(jl,il,1)
            enddo
         enddo
      else
         nt2 = NZ2/2 + mod(NZ2,2)
         IZAMR2FC(nt2,NZ2,1) = 1.0
         IZTAMR2FC(NZ2,nt2,1) = 1.0
         IZAMR2FC(nt2,NZ2,2) = 1.0
         IZTAMR2FC(NZ2,nt2,2) = 1.0
      endif


!     special treatment of axisymmetric 2D case
      if (.not.IF3D.and.IFAXIS) then
!     copy current interpolation operators
!     mesh 1
         call copy(ICAMR1CF,IYAMR1CF,NY1*NY1*2)
         call copy(ICAMR1FC,IYAMR1FC,NY1*NY1*2)
         call copy(ICTAMR1CF,IYTAMR1CF,NY1*NY1*2)
         call copy(ICTAMR1FC,IYTAMR1FC,NY1*NY1*2)
!     mesh 2
         call copy(ICAMR2CF,IYAMR2CF,NY2*NY2*2)
         call copy(ICAMR2FC,IYAMR2FC,NY2*NY2*2)
         call copy(ICTAMR2CF,IYTAMR2CF,NY2*NY2*2)
         call copy(ICTAMR2FC,IYTAMR2FC,NY2*NY2*2)

!     get intrpolation operators
!     mesh 1
!     coarse -> fine
!     Y
!     negative
         do il=1,NY1
            tmpl(il) = 0.5*(ZAM1(il) -1.0)
         enddo
         call IGLJM (IAAMR1CF,IATAMR1CF,ZAM1,tmpl,NY1,NY1,NY1,NY1)
!     positive; we use symmetry
         do jl=1,NY1
            do il=1,NY1
               IAAMR1CF(NY1-il+1,NY1-jl+1,2)  = IAAMR1CF(il,jl,1)
               IATAMR1CF(NY1-il+1,NY1-jl+1,2) = IATAMR1CF(il,jl,1)
            enddo
         enddo

!     fine -> coarse
!     Y
!     negative
         nt2 = NY1/2 + mod(NY1,2)
         do il=1,nt2
            tmpl(il) = 2.0*ZAM1(il) + 1.0
         enddo
         call IGLJM (IAAMR1FC,IATAMR1FC,ZAM1,tmpl,NY1,nt2,NY1,NY1)
!     positive; we use symmetry
         do jl=1,NY1
            do il=1,nt2
               IAAMR1FC(NY1-il+1,NY1-jl+1,2)  = IAAMR1FC(il,jl,1)
               IATAMR1FC(NY1-jl+1,NY1-il+1,2) = IATAMR1FC(jl,il,1)
            enddo
         enddo

!     MESH M2
!     coarse -> fine
!     Y
!     negative
         do il=1,NY2
            tmpl(il) = 0.5*(ZAM2(il) -1.0)
         enddo
         if (IFSPLIT) then      ! P-N-P_N
            call IGLJM (IAAMR2CF,IATAMR2CF,ZAM2,tmpl,
     $           NY2,NY2,NY2,NY2)
         else                   ! P-N-P_N-2
            call IGJM (IAAMR2CF,IATAMR2CF,ZAM2,tmpl,
     $           NY2,NY2,NY2,NY2)
         endif
!     positive; we use symmetry
         do jl=1,NY2
            do il=1,NY2
               IAAMR2CF(NY2-il+1,NY2-jl+1,2)  = IAAMR2CF(il,jl,1)
               IATAMR2CF(NY2-il+1,NY2-jl+1,2) = IATAMR2CF(il,jl,1)
            enddo
         enddo

!     fine -> coarse
!     Y
!     negative
         nt2 = NY2/2 + mod(NY2,2)
         do il=1,nt2
            tmpl(il) = 2.0*ZAM2(il) + 1.0
         enddo
         if (IFSPLIT) then      ! P-N-P_N
            call IGLJM (IAAMR2FC,IATAMR2FC,ZAM2,tmpl,
     $           NY2,nt2,NY2,NY2)
         else                   ! P-N-P_N-2
            call IGJM (IAAMR2FC,IATAMR2FC,ZAM2,tmpl,
     $           NY2,nt2,NY2,NY2)
         endif
!     positive; we use symmetry
         do jl=1,NY2
            do il=1,nt2
               IAAMR2FC(NY2-il+1,NY2-jl+1,2)  = IAAMR2FC(il,jl,1)
               IATAMR2FC(NY2-jl+1,NY2-il+1,2) = IATAMR2FC(jl,il,1)
            enddo
         enddo

      endif                     ! axisymmetric

!     for crs interpolation (vertices only)
!     MESH CRS
      call zwgll(zgmc,wgmc,AMR_LXC)
!     coarse -> fine
!     X
!     negative
      do il=1,AMR_LXC
         tmpl(il) = 0.5*(zgmc(il) -1.0)
      enddo
      call IGLLM (IXAMRCCF,IXTAMRCCF,zgmc,tmpl,AMR_LXC,AMR_LXC,
     $            AMR_LXC,AMR_LXC)
!     positive; we use symmetry
      do jl=1,AMR_LXC
         do il=1,AMR_LXC
            IXAMRCCF(AMR_LXC-il+1,AMR_LXC-jl+1,2)  = IXAMRCCF(il,jl,1)
            IXTAMRCCF(AMR_LXC-il+1,AMR_LXC-jl+1,2) = IXTAMRCCF(il,jl,1)
         enddo
      enddo

!     for multiplicity
!     mesh 1
!     X
      if (mod(NX1,2).eq.1) then
         il = NX1/2 + 1
         do kl= 1, NZ1
            do jl=1, NY1
               IMAMR1(il,jl,kl) = IMAMR1(il,jl,kl) +1.0
            enddo
         enddo
      endif

!     Y
      if (mod(NY1,2).eq.1) then
         jl = NY1/2 + 1
         do kl= 1, NZ1
            do il=1, NX1
               IMAMR1(il,jl,kl) = IMAMR1(il,jl,kl) +1.0
            enddo
         enddo
         if (mod(NX1,2).eq.1) then
            il = NX1/2 + 1
            do kl= 1, NZ1
               IMAMR1(il,jl,kl) = IMAMR1(il,jl,kl) +1.0
            enddo
         endif
      endif

!     Z
      if (IF3D) then
         if (mod(NZ1,2).eq.1) then
            kl = NZ1/2 + 1
            do jl= 1, NY1
               do il=1, NX1
                  IMAMR1(il,jl,kl) = IMAMR1(il,jl,kl) +1.0
               enddo
            enddo
            if (mod(NX1,2).eq.1) then
               il = NX1/2 + 1
               do jl=1,NY1
                  IMAMR1(il,jl,kl) = IMAMR1(il,jl,kl) +1.0
               enddo
            endif
            if (mod(NY1,2).eq.1) then
               jl = NY1/2 + 1
               do il=1,NX1
                  IMAMR1(il,jl,kl) = IMAMR1(il,jl,kl) +1.0
               enddo
            endif
            if (mod(NX1,2).eq.1.and.mod(NY1,2).eq.1) then
               il = NX1/2 + 1
               jl = NY1/2 + 1
               IMAMR1(il,jl,kl) = IMAMR1(il,jl,kl) +1.0
            endif
         endif
      endif

!     calculate inverse
      nt2 = NX1*NY1*NZ1
      call invcol1(IMAMR1,nt2)

!     mesh 2
!     X
      if (mod(NX2,2).eq.1) then
         il = NX2/2 + 1
         do kl= 1, NZ2
            do jl=1, NY2
               IMAMR2(il,jl,kl) = IMAMR2(il,jl,kl) +1.0
            enddo
         enddo
      endif

!     Y
      if (mod(NY2,2).eq.1) then
         jl = NY2/2 + 1
         do kl= 1, NZ2
            do il=1, NX2
               IMAMR2(il,jl,kl) = IMAMR2(il,jl,kl) +1.0
            enddo
         enddo
         if (mod(NX2,2).eq.1) then
            il = NX2/2 + 1
            do kl= 1, NZ2
               IMAMR2(il,jl,kl) = IMAMR2(il,jl,kl) +1.0
            enddo
         endif
      endif

!     Z
      if (IF3D) then
         if (mod(NZ2,2).eq.1) then
            kl = NZ2/2 + 1
            do jl= 1, NY2
               do il=1, NX2
                  IMAMR2(il,jl,kl) = IMAMR2(il,jl,kl) +1.0
               enddo
            enddo
            if (mod(NX2,2).eq.1) then
               il = NX2/2 + 1
               do jl=1,NY2
                  IMAMR2(il,jl,kl) = IMAMR2(il,jl,kl) +1.0
               enddo
            endif
            if (mod(NY2,2).eq.1) then
               jl = NY2/2 + 1
               do il=1,NX2
                  IMAMR2(il,jl,kl) = IMAMR2(il,jl,kl) +1.0
               enddo
            endif
            if (mod(NX2,2).eq.1.and.mod(NY2,2).eq.1) then
               il = NX2/2 + 1
               jl = NY2/2 + 1
               IMAMR2(il,jl,kl) = IMAMR2(il,jl,kl) +1.0
            endif
         endif
      endif

!     calculate inverse
      nt2 = NX2*NY2*NZ2
      call invcol1(IMAMR2,nt2)

!     to get proper J-1 on faces and edges for fast diagonalisation method
!     I assume here LX1=LY1=LZ1, so only one array is needed
      call ftovecl(IMAMR1F,IMAMR1,1,NX1,NY1,NZ1)
      call etovec(IMAMR1E,1,IMAMR1,NX1,NY1,NZ1)

#ifdef DEBUG
!     testing
      call io_file_freeid(iunit, ierr)
      write(str,'(i2.2)') NID
      open(unit=iunit,file='ref_coeff.txt'//str)
      write(iunit,*) 'IXAMR1CF 1; m1'
      do il=1,nx1
         write(iunit,*) il,(IXAMR1CF(il,jl,1),jl=1,nx1)
      enddo
      write(iunit,*) 'IXAMR1CF 2'
      do il=1,nx1
         write(iunit,*) il,(IXAMR1CF(il,jl,2),jl=1,nx1)
      enddo
      write(iunit,*) 'IXAMR1FC 1; m1'
      do il=1,nx1
         write(iunit,*) il,(IXAMR1FC(il,jl,1),jl=1,nx1)
      enddo
      write(iunit,*) 'IXAMR1FC 2'
      do il=1,nx1
         write(iunit,*) il,(IXAMR1FC(il,jl,2),jl=1,nx1)
      enddo

      write(iunit,*) 'IXTAMR1CF 1; m1'
      do il=1,nx1
         write(iunit,*) il,(IXTAMR1CF(il,jl,1),jl=1,nx1)
      enddo
      write(iunit,*) 'IXTAMR1CF 2'
      do il=1,nx1
         write(iunit,*) il,(IXTAMR1CF(il,jl,2),jl=1,nx1)
      enddo
      write(iunit,*) 'IXTAMR1FC 1; m1'
      do il=1,nx1
         write(iunit,*) il,(IXTAMR1FC(il,jl,1),jl=1,nx1)
      enddo
      write(iunit,*) 'IXTAMR1FC 2'
      do il=1,nx1
         write(iunit,*) il,(IXTAMR1FC(il,jl,2),jl=1,nx1)
      enddo

      write(iunit,*) 'IYAMR1CF 1; m1'
      do il=1,nx1
         write(iunit,*) il,(IYAMR1CF(il,jl,1),jl=1,nx1)
      enddo
      write(iunit,*) 'IYAMR1CF 2'
      do il=1,nx1
         write(iunit,*) il,(IYAMR1CF(il,jl,2),jl=1,nx1)
      enddo
      write(iunit,*) 'IYAMR1FC 1; m1'
      do il=1,nx1
         write(iunit,*) il,(IYAMR1FC(il,jl,1),jl=1,nx1)
      enddo
      write(iunit,*) 'IYAMR1FC 2'
      do il=1,nx1
         write(iunit,*) il,(IYAMR1FC(il,jl,2),jl=1,nx1)
      enddo

      write(iunit,*) 'IYTAMR1CF 1; m1'
      do il=1,nx1
         write(iunit,*) il,(IYTAMR1CF(il,jl,1),jl=1,nx1)
      enddo
      write(iunit,*) 'IYTAMR1CF 2'
      do il=1,nx1
         write(iunit,*) il,(IYTAMR1CF(il,jl,2),jl=1,nx1)
      enddo
      write(iunit,*) 'IYTAMR1FC 1; m1'
      do il=1,nx1
         write(iunit,*) il,(IYTAMR1FC(il,jl,1),jl=1,nx1)
      enddo
      write(iunit,*) 'IYTAMR1FC 2'
      do il=1,nx1
         write(iunit,*) il,(IYTAMR1FC(il,jl,2),jl=1,nx1)
      enddo

!     mesh 2
      write(iunit,*) 'IXAMR2CF 1; m2'
      do il=1,nx2
         write(iunit,*) il,(IXAMR2CF(il,jl,1),jl=1,nx2)
      enddo
      write(iunit,*) 'IXAMR2CF 2'
      do il=1,nx2
         write(iunit,*) il,(IXAMR2CF(il,jl,2),jl=1,nx2)
      enddo
      write(iunit,*) 'IXAMR2FC 1; m2'
      do il=1,nx2
         write(iunit,*) il,(IXAMR2FC(il,jl,1),jl=1,nx2)
      enddo
      write(iunit,*) 'IXAMR2FC 2'
      do il=1,nx2
         write(iunit,*) il,(IXAMR2FC(il,jl,2),jl=1,nx2)
      enddo

      write(iunit,*) 'IXTAMR2CF 1; m2'
      do il=1,nx2
         write(iunit,*) il,(IXTAMR2CF(il,jl,1),jl=1,nx2)
      enddo
      write(iunit,*) 'IXTAMR2CF 2'
      do il=1,nx2
         write(iunit,*) il,(IXTAMR2CF(il,jl,2),jl=1,nx2)
      enddo
      write(iunit,*) 'IXTAMR2FC 1; m2'
      do il=1,nx2
         write(iunit,*) il,(IXTAMR2FC(il,jl,1),jl=1,nx2)
      enddo
      write(iunit,*) 'IXTAMR2FC 2'
      do il=1,nx2
         write(iunit,*) il,(IXTAMR2FC(il,jl,2),jl=1,nx2)
      enddo


      write(iunit,*) 'IYAMR2CF 1; m2'
      do il=1,nx2
         write(iunit,*) il,(IYAMR2CF(il,jl,1),jl=1,nx2)
      enddo
      write(iunit,*) 'IYAMR2CF 2'
      do il=1,nx2
         write(iunit,*) il,(IYAMR2CF(il,jl,2),jl=1,nx2)
      enddo
       write(iunit,*) 'IYTAMR2CF 1; m2'
      do il=1,nx2
         write(iunit,*) il,(IYTAMR2CF(il,jl,1),jl=1,nx2)
      enddo
      write(iunit,*) 'IYTAMR2CF 2'
      do il=1,nx2
         write(iunit,*) il,(IYTAMR2CF(il,jl,2),jl=1,nx2)
      enddo

      write(iunit,*) 'IXAMRCCF 1; crs'
      do il=1,amr_lxc
         write(iunit,*) il,(IXAMRCCF(il,jl,1),jl=1,amr_lxc)
      enddo
      write(iunit,*) 'IXAMRCCF 2; crs'
      do il=1,amr_lxc
         write(iunit,*) il,(IXAMRCCF(il,jl,2),jl=1,amr_lxc)
      enddo
      write(iunit,*) 'IXTAMRCCF 1; crs'
      do il=1,amr_lxc
         write(iunit,*) il,(IXTAMRCCF(il,jl,1),jl=1,amr_lxc)
      enddo
      write(iunit,*) 'IXTAMRCCF 2; crs'
      do il=1,amr_lxc
         write(iunit,*) il,(IXTAMRCCF(il,jl,2),jl=1,amr_lxc)
      enddo

      close(iunit)
#endif
      return
      end subroutine
!=======================================================================
!> @brief Get interpolation parameters for octal interpolation; hsmg solver
      subroutine amr_genwz_hsmg
      implicit none

      include 'SIZE'
      include 'INPUT'           ! IF3D, IFSPLIT
      include 'WZ'              ! ZGM(123)
      include 'HSMG'
      include 'AMR_REFINE_HSMG'

!     local variables
      integer level,il, jl, kl     !     loop indexes
      integer nt1, nzt1, nt2

!     tmp  array to calculate interpolation positions
      integer nmax
      parameter (nmax=84) ! like in speclib
      real tmpl(nmax)

#ifdef DEBUG
!     for testing
      character*2 str
      integer iunit, ierr
#endif
!-----------------------------------------------------------------------
!     zero arrays
!     Volume interpolation
      nt2 = LXM*LXM*2*LMGN
      call rzero(IXAMRHCF,nt2)
      call rzero(IXTAMRHCF,nt2)
      call rzero(IXAMRHFC,nt2)
      call rzero(IXTAMRHFC,nt2)

!     multiplicity arrays
      nt2 = LXM*LXM*LXM*LMGN
      call rone(IMAMRH,nt2)
      nt2 = LXM*LXM*LMGN
      call rone(IMAMRHF,nt2)

!     get information about multi-grid levels
      call hsmg_setup_mg_nx

!     we skip coarse-grid solver
      if (MG_LMAX.lt.3) return

!     get coordinates
      call hsmg_setup_semhat

!     get interpolation operators
!     I skip level 1 as for P_N-P_N-2 this is covered by coarse grid interpolation
      do level=2,MG_LMAX-1
!     coarse -> fine
!     X
!     negative
         nt1 = MG_NH(level)
         do il=1,nt1
            tmpl(il) = 0.5*(MG_ZH(il,level) -1.0)
         enddo
         call IGLLM (IXAMRHCF(1,1,level),IXTAMRHCF(1,1,level),
     $              MG_ZH(1,level),tmpl,nt1,nt1,nt1,nt1)
!     positive; we use symmetry
         do jl=1,nt1
            do il=1,nt1
               IXAMRHCF(nt1-il+1+(nt1-jl)*nt1,2,level)  =
     $              IXAMRHCF(il+(jl-1)*nt1,1,level)
               IXTAMRHCF(nt1-il+1+(nt1-jl)*nt1,2,level) =
     $              IXTAMRHCF(il+(jl-1)*nt1,1,level)
            enddo
         enddo

!     fine -> coarse
!     X
!     negative
         nt2 = nt1/2 + mod(nt1,2)
         do il=1,nt2
            tmpl(il) = 2.0*MG_ZH(il,level) + 1.0
         enddo
         call IGLLM (IXAMRHFC(1,1,level),IXTAMRHFC(1,1,level),
     $              MG_ZH(1,level),tmpl,nt1,nt2,nt1,nt1)
!     positive; we use symmetry
         do jl=1,nt1
            do il=1,nt2
               IXAMRHFC(nt1-il+1+(nt1-jl)*nt1,2,level)  =
     $              IXAMRHFC(il+(jl-1)*nt1,1,level)
               IXTAMRHFC(nt1-jl+1+(nt1-il)*nt1,2,level) =
     $              IXTAMRHFC(jl+(il-1)*nt1,1,level)
            enddo
         enddo
      enddo ! level

!     for multiplicity
!     I skip level 1 as for P_N-P_N-2 this is covered by coarse grid interpolation
      do level=2,MG_LMAX-1
         nt1 = MG_NH(level)+2 ! because we work on Schwarz
         if (IF3D) then
            nzt1 = nt1
         else
            nzt1 = 1
         endif
         if (mod(nt1,2).eq.1) then
!     X
            il = nt1/2 + 1
            do kl=1, nzt1
               do jl=1, nt1
                  IMAMRH(il+(jl-1 +(kl-1)*nt1)*nt1,level) =
     $                   IMAMRH(il+(jl-1 +(kl-1)*nt1)*nt1,level) +1.0
               enddo
            enddo
!     Y
            jl = nt1/2 + 1
            do kl=1, nzt1
               do il=1, nt1
                  IMAMRH(il+(jl-1 +(kl-1)*nt1)*nt1,level) =
     $                   IMAMRH(il+(jl-1 +(kl-1)*nt1)*nt1,level) +1.0
               enddo
            enddo
            il = nt1/2 + 1
            do kl=1, nzt1
               IMAMRH(il+(jl-1 +(kl-1)*nt1)*nt1,level) =
     $               IMAMRH(il+(jl-1 +(kl-1)*nt1)*nt1,level) +1.0
            enddo
!     Z
            if (IF3D) then
               kl = nt1/2 + 1
               do jl=1, nt1
                  do il=1, nt1
                     IMAMRH(il+(jl-1 +(kl-1)*nt1)*nt1,level) =
     $                   IMAMRH(il+(jl-1 +(kl-1)*nt1)*nt1,level) +1.0
                  enddo
               enddo
               il = nt1/2 + 1
               do jl=1, nt1
                  IMAMRH(il+(jl-1 +(kl-1)*nt1)*nt1,level) =
     $               IMAMRH(il+(jl-1 +(kl-1)*nt1)*nt1,level) +1.0
               enddo
               jl = nt1/2 + 1
               do il=1, nt1
                  IMAMRH(il+(jl-1 +(kl-1)*nt1)*nt1,level) =
     $               IMAMRH(il+(jl-1 +(kl-1)*nt1)*nt1,level) +1.0
               enddo
               il = nt1/2 + 1
               IMAMRH(il+(jl-1 +(kl-1)*nt1)*nt1,level) =
     $               IMAMRH(il+(jl-1 +(kl-1)*nt1)*nt1,level) +1.0
            endif

         endif

!     calculate inverse
         nt2 = nt1*nt1*nzt1
         call invcol1(IMAMRH(1,level),nt2)

!     to get proper J-1 on faces for fast diagonalisation method
!     I assume here LX1=LY1=LZ1, so only one array is needed
         call ftovecl(IMAMRHF(1,level),IMAMRH(1,level),1,nt1,nt1,nzt1)

      enddo  ! level

      return
      end subroutine
!=======================================================================
!> @brief Map single variable of coarse element to the fine one
!! @param[out]   vf      fine element vector
!! @param[in]    vc      coarse element vector
!! @param[in]    iel     element number
!! @param[in]    ch_pos  child position
!! @param[in]    limesh   mesh mark (velocity, pressure, mhd)
!! @param[in]    tmp     work array tmp(lnx,lny,lnz,2)
!! @param[in]    lnx     element size (x)
!! @param[in]    lny     element size (y)
!! @param[in]    lnz     element size (z)
      subroutine amr_mapcf(vf,vc,iel,ch_pos,limesh,tmp,lnx,lny,lnz)
      implicit none

      include 'SIZE'
      include 'INPUT'           ! IF3D, IFSPLIT
      include 'GEOM'            ! IFRZER
      include 'AMR'
      include 'AMR_REFINE'

      ! argument list
      integer lnx,lny,lnz
      real vf(lnx,lny,lnz), vc(lnx,lny,lnz)
      integer iel,ch_pos(3) ! ch_pos is child position in 3D
      integer limesh
      real tmp(lnx,lny,lnz,2) ! work array

      ! local variables
      integer nxy, nyz
      integer iz  ! loop index
!-----------------------------------------------------------------------
      nxy = lnx*lny
      nyz = lny*lnz

!      ! Use the appropriate derivative- and interpolation operator in
!      ! the y-direction (= radial direction if axisymmetric).

      if (limesh.eq.1) then ! velocity mesh

         if (IFAXIS) then
            if (IFRZER(iel)) then
               call copy(IYTAMR1CF,IATAMR1CF,lny*lny*2)
            else
               call copy(IYTAMR1CF,ICTAMR1CF,lny*lny*2)
            endif
         endif

         if (IF3D) then
            call mxm(IXAMR1CF(1,1,ch_pos(1)),lnx,vc,lnx,
     $           tmp(1,1,1,1),nyz)
            do iz = 1,lnz
               call mxm(tmp(1,1,iz,1),lnx,IYTAMR1CF(1,1,ch_pos(2)),
     $              lny,tmp(1,1,iz,2),lny)
            enddo
            call mxm(tmp(1,1,1,2),nxy,IZTAMR1CF(1,1,ch_pos(3)),lnz,
     $           vf,lnz)
         else
            call mxm(IXAMR1CF(1,1,ch_pos(1)),lnx,vc,lnx,tmp,nyz)
            call mxm(tmp,lnx,IYTAMR1CF(1,1,ch_pos(2)),lny,vf,lny)
         endif
      elseif (limesh.eq.2) then !pressure mesh

         if (IFAXIS) then
            if (IFRZER(iel)) then
               call copy(IYTAMR2CF,IATAMR2CF,lny*lny*2)
            else
               call copy(IYTAMR2CF,ICTAMR2CF,lny*lny*2)
            endif
         endif

         if (IF3D) then
            call mxm(IXAMR2CF(1,1,ch_pos(1)),lnx,vc,lnx,
     $           tmp(1,1,1,1),nyz)
            do iz = 1,lnz
               call mxm(tmp(1,1,iz,1),lnx,IYTAMR2CF(1,1,ch_pos(2)),
     $              lny,tmp(1,1,iz,2),lny)
            enddo
            call mxm(tmp(1,1,1,2),nxy,IZTAMR2CF(1,1,ch_pos(3)),lnz,
     $           vf,lnz)
         else
            call mxm(IXAMR2CF(1,1,ch_pos(1)),lnx,vc,lnx,tmp,nyz)
            call mxm(tmp,lnx,IYTAMR2CF(1,1,ch_pos(2)),lny,vf,lny)
         endif
      else  ! there should be as well place for mhd mesh
         call amr_abort('ERROR: amr_mapcf; wrong limesh')
      endif

      return
      end subroutine
!=======================================================================
!> @brief Map single variable of fine element to the coarse one
!! @param[out]   vc      coarse element vector
!! @param[in]    vf      fine element vector
!! @param[in]    iel     element number
!! @param[in]    ch_pos  child position
!! @param[in]    limesh   mesh mark (velocity, pressure, mhd)
!! @param[in]    tmp     work array tmp(lnx,lny,lnz,2)
!! @param[in]    lnx     element size (x)
!! @param[in]    lny     element size (y)
!! @param[in]    lnz     element size (z)
      subroutine amr_mapfc(vc,vf,iel,ch_pos,limesh,tmp,lnx,lny,lnz)
      implicit none

      include 'SIZE'
      include 'INPUT'           ! IF3D, IFSPLIT
      include 'GEOM'            ! IFRZER
      include 'AMR'
      include 'AMR_REFINE'

      ! argument list
      integer lnx,lny,lnz
      real vc(lnx,lny,lnz), vf(lnx,lny,lnz)
      integer iel,ch_pos(3)   !     ch_pos is child position in 3D
      integer limesh
      real tmp(lnx,lny,lnz,2) ! work array

      ! local variables
      integer nxy, nyz
      integer iz  ! loop index
!-----------------------------------------------------------------------
      nxy = lnx*lny
      nyz = lny*lnz

!      ! Use the appropriate derivative- and interpolation operator in
!      ! the y-direction (= radial direction if axisymmetric).

      if (limesh.eq.1) then ! velocity mesh

         if (IFAXIS) then
            if (IFRZER(iel)) then
               call copy(IYTAMR1FC,IATAMR1FC,lny*lny*2)
            else
               call copy(IYTAMR1FC,ICTAMR1FC,lny*lny*2)
            endif
         endif

         if (IF3D) then
            call mxm(IXAMR1FC(1,1,ch_pos(1)),lnx,vf,lnx,
     $           tmp(1,1,1,1),nyz)
            do iz = 1,lnz
               call mxm(tmp(1,1,iz,1),lnx,IYTAMR1FC(1,1,ch_pos(2)),
     $              lny,tmp(1,1,iz,2),lny)
            enddo
            call mxm(tmp(1,1,1,2),nxy,IZTAMR1FC(1,1,ch_pos(3)),lnz,
     $           vc,lnz)
         else
            call mxm(IXAMR1FC(1,1,ch_pos(1)),lnx,vf,lnx,tmp,nyz)
            call mxm(tmp,lnx,IYTAMR1FC(1,1,ch_pos(2)),lny,vc,lny)
         endif
      elseif (limesh.eq.2) then !pressure mesh

         if (IFAXIS) then
            if (IFRZER(iel)) then
               call copy(IYTAMR2FC,IATAMR2FC,lny*lny*2)
            else
               call copy(IYTAMR2FC,ICTAMR2FC,lny*lny*2)
            endif
         endif

         if (IF3D) then
            call mxm(IXAMR2FC(1,1,ch_pos(1)),lnx,vf,lnx,
     $           tmp(1,1,1,1),nyz)
            do iz = 1,lnz
               call mxm(tmp(1,1,iz,1),lnx,IYTAMR2FC(1,1,ch_pos(2)),
     $              lny,tmp(1,1,iz,2),lny)
            enddo
            call mxm(tmp(1,1,1,2),nxy,IZTAMR2FC(1,1,ch_pos(3)),lnz,
     $           vc,lnz)
         else
            call mxm(IXAMR2FC(1,1,ch_pos(1)),lnx,vf,lnx,tmp,nyz)
            call mxm(tmp,lnx,IYTAMR2FC(1,1,ch_pos(2)),lny,vc,lny)
         endif
      else  ! there should be as well place for mhd mesh
         call amr_abort('ERROR: amr_mapfc; wrong limesh')
      endif

      return
      end subroutine
!=======================================================================
!> @brief Perform single element refinement operation of given variable
!! @param[inout] vcf     refined vector
!! @param[in]    el_lst  element list
!! @param[in]    limesh   mesh mark (velocity, pressure, mhd)
!! @param[in]    tmp     work array tmp(lnx,lny,lnz,3)
!! @param[in]    lnx     element size (x)
!! @param[in]    lny     element size (y)
!! @param[in]    lnz     element size (z)
!! @param[in]    leln    element count
      subroutine amr_refine_vs(vcf,el_lst,limesh,tmp,
     $           lnx,lny,lnz,leln)
      implicit none

      include 'SIZE'
      include 'AMR'

      ! argument list
      integer lnx,lny,lnz,leln
      real vcf(lnx,lny,lnz,leln)
      integer el_lst(AMR_NVRT)
      integer limesh
      real tmp(lnx,lny,lnz,3) ! work array

      ! local variables
      integer ch_pos(3)   !     child position in 3D
      integer iel         !     element position
      integer il, nl
!-----------------------------------------------------------------------
      ! I assume el_lst(1) gives position of the coarse block
      ! and final ch_pos() = 1,1,1
      ! copy coarse element
      nl=lnx*lny*lnz
      call copy(tmp,vcf(1,1,1,el_lst(1)),nl)
      ! loop over all the children
      do il= 1,AMR_NVRT
         ! get child position
         iel = el_lst(il)
         ch_pos(3) = (il-1)/4 +1
         ch_pos(2) = mod((il-1)/2,2) +1
         ch_pos(1) = mod(il-1,2) +1

         ! refine
         call amr_mapcf(vcf(1,1,1,iel),tmp(1,1,1,1),iel,ch_pos,
     $                       limesh,tmp(1,1,1,2),lnx,lny,lnz)

      enddo

      return
      end subroutine
!=======================================================================
!> @brief Perform single element coarsening operation of given variable
!! @param[inout] vfc       coarsened vector
!! @param[in]    el_lst    element list
!! @param[in]    limesh   mesh mark (velocity, pressure, mhd)
!! @param[in]    tmp     work array tmp(lnx,lny,lnz,3)
!! @param[in]    lnx     element size (x)
!! @param[in]    lny     element size (y)
!! @param[in]    lnz     element size (z)
!! @param[in]    leln    element count
      subroutine amr_coarse_vs(vfc,el_lst,limesh,tmp,
     $           lnx,lny,lnz,leln)
      implicit none

      include 'SIZE'
      include 'AMR'
      include 'AMR_REFINE'

      ! argument list
      integer lnx,lny,lnz,leln
      real vfc(lnx,lny,lnz,leln)
      integer el_lst(AMR_NVRT)
      integer limesh
      real tmp(lnx,lny,lnz,3) ! work array

      ! local variables
      integer ch_pos(3)    !     child position in 3D
      integer iel          !     element position
      integer il, nl
!-----------------------------------------------------------------------
      ! I assume el_lst(1) gives position of the coarse block
      ! and inital ch_pos() = 1,1,1
      nl=lnx*lny*lnz
      ! loop over all the children
      do il= 1,AMR_NVRT
         ! get child position
         iel = el_lst(il)
         ch_pos(3) = (il-1)/4 +1
         ch_pos(2) = mod((il-1)/2,2) +1
         ch_pos(1) = mod(il-1,2) +1
         ! coarsen
         call amr_mapfc(tmp(1,1,1,1),vfc(1,1,1,iel),iel,ch_pos,
     $                       limesh,tmp(1,1,1,2),lnx,lny,lnz)
         ! sum contributions
         if (il.eq.1) then
            call copy(vfc(1,1,1,iel),tmp,nl)
         else
            call add2(vfc(1,1,1,el_lst(1)),tmp,nl)
         endif
      enddo
      ! take into account multiplicity of points
      if (limesh.eq.1) then ! velocity mesh
         call col2(vfc(1,1,1,el_lst(1)),IMAMR1,nl)
      elseif (limesh.eq.2) then !pressure mesh
         call col2(vfc(1,1,1,el_lst(1)),IMAMR2,nl)
      else  ! there should be as well place for mhd mesh
         call amr_abort('ERROR: amr_coarsen_vs; wrong limesh')
      endif

      return
      end subroutine
!=======================================================================
!> @brief Perform all refinement operations for given variable
!! @details Threre are two possible implementations. In following rourtine
!! we performe all block refinemnts for a single variable.
!! @param[inout] vcf     refined vector
!! @param[in]    limesh   mesh mark (velocity, pressure, mhd)
!! @param[in]    tmp     work array tmp(lnx,lny,lnz,3)
!! @param[in]    lnx     element size (x)
!! @param[in]    lny     element size (y)
!! @param[in]    lnz     element size (z)
!! @param[in]    leln    element count
      subroutine amr_refine_vm(vcf,limesh,tmp,lnx,lny,lnz,leln)
      implicit none

      include 'SIZE'
      include 'AMR'

      ! argument list
      integer lnx,lny,lnz,leln
      real vcf(lnx,lny,lnz,leln)
      integer limesh
      real tmp(lnx,lny,lnz,3) ! work array

      ! local variables
      integer el_lst(AMR_NVRT)   ! local element list for refinement
      integer il, jl ! loop index
!-----------------------------------------------------------------------
      do il=0,AMR_RFN_NR-1,AMR_NVRT
        do jl=1,AMR_NVRT
            el_lst(jl) = AMR_GLGL_RFN(3,il+jl)
        enddo
        call amr_refine_vs(vcf,el_lst,limesh,tmp,lnx,lny,lnz,leln)
      enddo

      return
      end subroutine
!=======================================================================
!> @brief Perform all coarsening operations for given variable
!! @details Threre are two possible implementations. In following rourtine
!! we performe all block coarsenings for a single variable.
!! @param[inout] vfc       coarsened vector
!! @param[in]    limesh   mesh mark (velocity, pressure, mhd)
!! @param[in]    tmp     work array tmp(lnx,lny,lnz,3)
!! @param[in]    lnx     element size (x)
!! @param[in]    lny     element size (y)
!! @param[in]    lnz     element size (z)
!! @param[in]    leln    element count
      subroutine amr_coarse_vm(vfc,limesh,tmp,
     $           lnx,lny,lnz,leln)
      implicit none

      include 'SIZE'
      include 'AMR'
      include 'AMR_REFINE'

      ! argument list
      integer lnx,lny,lnz,leln
      real vfc(lnx,lny,lnz,leln)
      integer limesh
      real tmp(lnx,lny,lnz,3) ! work array

      ! local variables
      integer el_lst(AMR_NVRT)   ! local element list for refinement
      integer il, jl ! loop index
!-----------------------------------------------------------------------
      do il=1,AMR_CRS_NR
        do jl=1,AMR_NVRT
            el_lst(jl) = AMR_GLGL_CRS(2,jl,il)
        enddo
        call amr_coarse_vs(vfc,el_lst,limesh,tmp,lnx,lny,lnz,leln)
      enddo

      return
      end subroutine
!=======================================================================
!> @brief Refinement, transfer, coarsening of single vector
!! @details This routine performs local refinement followed by data transfer
!!    and local coarsening of a single variable
!! @param[inout] var     refined vector
!! @param[inout] vrt     transfer mark array (work array)
!! @param[in]    tmp     work array tmp(lnx,lny,lnz,3)
!! @param[in]    limesh   mesh mark (velocity, pressure, mhd)
!! @param[in]    lnx     element size (x)
!! @param[in]    lny     element size (y)
!! @param[in]    lnz     element size (z)
!! @param[in]    lbuf    element count
      subroutine amr_refine_coarsen_var(var,vrt,tmp,
     $             limesh,lnx,lny,lnz,lbuf)
      implicit none

      include 'SIZE'
      include 'AMR'

      ! argument list
      integer lnx,lny,lnz,lbuf
      integer limesh
      real var(lnx,lny,lnz,lbuf)
      integer vrt(2,lbuf)
      real tmp(lnx,lny,lnz,3) ! work array

      ! local variables
      integer itmp
      real t1, t2

      ! functions
      real dnekclock
!-----------------------------------------------------------------------
      ! local refinement
      t1 = dnekclock()
      call amr_refine_vm(var,limesh,tmp,lnx,lny,lnz,lbuf)
      t2 = dnekclock()
      AMR_TCL = AMR_TCL + t2 - t1

      ! data redistribution
      itmp = lnx*lny*lnz
      call amr_vec_transfer(var,vrt,itmp,lbuf)
      t1 = dnekclock()
      AMR_TCC = AMR_TCC + t1 - t2

      ! local coarsening
      call amr_coarse_vm(var,limesh,tmp,lnx,lny,lnz,lbuf)
      t2 = dnekclock()
      AMR_TCL = AMR_TCL + t2 - t1

      return
      end subroutine
!=======================================================================
!> @brief Refinement, transfer, coarsening for all variables
!! @details This routine performs local refinement followed by data transfer
!!    and local coarsening of all necessary variables. First given variable
!!    is coppied to bigger array to avoid problem with element number
!!    exceeding lelt during refinement stage.
!! @remarks This routine uses global scratch space SCRNS, SCRUZ
      subroutine amr_refine_coarsen
      implicit none

      include 'SIZE'
      include 'INPUT'           ! IF3D
      include 'SOLN'
      include 'TSTEP'           ! NBDINP
      include 'GEOM'            ! IFGEOM
      include 'AMR'

      ! local variables
      integer il, jl ! loop index
      integer ntot, ntoto
      integer limesh  ! mesh mark (velocity, pressure, mhd)

      ! work arrays
      real tmpv(LX1,LY1,LZ1,3), tmpvp(LX2,LY2,LZ2,3),
     $     tmpb(LBX1,LBY1,LBZ1,3), tmpbp(LBX2,LBY2,LBZ2,3)
      ! big arrays to perform refinement
      integer lbuf    ! buffer size (element number)
      parameter (lbuf=LELT*AMR_NVRT)
      real vr(LX1*LY1*LZ1*lbuf)
      integer vi(2*lbuf)
      common /scrns/ vr
      common /scruz/ vi
!-----------------------------------------------------------------------
      ! mesh 1
      limesh = 1
      ! coordinates
      ntot  = lx1*ly1*lz1*NELT
      ntoto = lx1*ly1*lz1*AMR_NELT_O

      call copy(vr,XM1,ntoto)
      call amr_refine_coarsen_var(vr,vi,tmpv,limesh,
     $                                 LX1,LY1,LZ1,lbuf)
      call copy(XM1,vr,ntot)

      call copy(vr,YM1,ntoto)
      call amr_refine_coarsen_var(vr,vi,tmpv,limesh,
     $                                 LX1,LY1,LZ1,lbuf)
      call copy(YM1,vr,ntot)

      if(IF3D) then
         call copy(vr,ZM1,ntoto)
         call amr_refine_coarsen_var(vr,vi,tmpv,limesh,
     $                                 LX1,LY1,LZ1,lbuf)
         call copy(ZM1,vr,ntot)
      endif

      ! velocity
      ntot = lx1*ly1*lz1*NELV
      ntoto = lx1*ly1*lz1*AMR_NELV_O

      call copy(vr,VX,ntoto)
      call amr_refine_coarsen_var(vr,vi,tmpv,limesh,
     $                                 LX1,LY1,LZ1,lbuf)
      call copy(VX,vr,ntot)

      call copy(vr,VY,ntoto)
      call amr_refine_coarsen_var(vr,vi,tmpv,limesh,
     $                                 LX1,LY1,LZ1,lbuf)
      call copy(VY,vr,ntot)

      if(IF3D) then
         call copy(vr,VZ,ntoto)
         call amr_refine_coarsen_var(vr,vi,tmpv,limesh,
     $                                 LX1,LY1,LZ1,lbuf)
         call copy(VZ,vr,ntot)
      endif

      ! arrays for time integration AB[XYZ][12]
      if (IFTRAN) then
         call copy(vr,ABX1,ntoto)
         call amr_refine_coarsen_var(vr,vi,tmpv,limesh,
     $                                 LX1,LY1,LZ1,lbuf)
         call copy(ABX1,vr,ntot)

         call copy(vr,ABY1,ntoto)
         call amr_refine_coarsen_var(vr,vi,tmpv,limesh,
     $                                 LX1,LY1,LZ1,lbuf)
         call copy(ABY1,vr,ntot)

         if(IF3D) then
            call copy(vr,ABZ1,ntoto)
            call amr_refine_coarsen_var(vr,vi,tmpv,limesh,
     $                                 LX1,LY1,LZ1,lbuf)
            call copy(ABZ1,vr,ntot)
         endif

         call copy(vr,ABX2,ntoto)
         call amr_refine_coarsen_var(vr,vi,tmpv,limesh,
     $                                 LX1,LY1,LZ1,lbuf)
         call copy(ABX2,vr,ntot)

         call copy(vr,ABY2,ntoto)
         call amr_refine_coarsen_var(vr,vi,tmpv,limesh,
     $                                 LX1,LY1,LZ1,lbuf)
         call copy(ABY2,vr,ntot)

         if(IF3D) then
            call copy(vr,ABZ2,ntoto)
            call amr_refine_coarsen_var(vr,vi,tmpv,limesh,
     $                                 LX1,LY1,LZ1,lbuf)
            call copy(ABZ2,vr,ntot)
         endif

      endif

      ! lag arrays velocity
      if (IFFLOW.or.IFMHD) then
         do il =1, NBDINP-1

            call copy(vr,VXLAG(1,1,1,1,il),ntoto)
            call amr_refine_coarsen_var(vr,vi,tmpv,limesh,
     $                                 LX1,LY1,LZ1,lbuf)
            call copy(VXLAG(1,1,1,1,il),vr,ntot)

            call copy(vr,VYLAG(1,1,1,1,il),ntoto)
            call amr_refine_coarsen_var(vr,vi,tmpv,limesh,
     $                                 LX1,LY1,LZ1,lbuf)
            call copy(VYLAG(1,1,1,1,il),vr,ntot)

            if(IF3D) then
               call copy(vr,VZLAG(1,1,1,1,il),ntoto)
               call amr_refine_coarsen_var(vr,vi,tmpv,limesh,
     $                                 LX1,LY1,LZ1,lbuf)
               call copy(VZLAG(1,1,1,1,il),vr,ntot)
            endif
         enddo
      endif

      ! mesh 2
      limesh = 2

      ! pressure
      ntot = lx2*ly2*lz2*NELV
      ntoto = lx2*ly2*lz2*AMR_NELV_O

      call copy(vr,PR,ntoto)
      call amr_refine_coarsen_var(vr,vi,tmpvp,limesh,
     $                                 LX2,LY2,LZ2,lbuf)
      call copy(PR,vr,ntot)

      ! lag arrays pressure PRLAG
      if (NBD.ge.2) then
         call copy(vr,AMR_PRLAG(1,1,1,1,1),ntoto)
         call amr_refine_coarsen_var(vr,vi,tmpvp,limesh,
     $                                 LX2,LY2,LZ2,lbuf)
         call copy(AMR_PRLAG(1,1,1,1,1),vr,ntot)
         ! make pressuere lag arrays consistent
         call copy(PRLAG(1,1,1,1,1),vr,ntot)
      endif

      if (NBD.eq.3) then
         call copy(vr,AMR_PRLAG(1,1,1,1,2),ntoto)
         call amr_refine_coarsen_var(vr,vi,tmpvp,limesh,
     $                                 LX2,LY2,LZ2,lbuf)
         call copy(AMR_PRLAG(1,1,1,1,2),vr,ntot)
      endif

      ! user defined divergence USRDIV
      ntot = lx2*ly2*lz2*NELT
      ntoto = lx2*ly2*lz2*AMR_NELT_O

      call copy(vr,USRDIV,ntoto)
      call amr_refine_coarsen_var(vr,vi,tmpvp,limesh,
     $                                 LX2,LY2,LZ2,lbuf)
      call copy(USRDIV,vr,ntot)

      ! temperature and passive scalars
      if (IFHEAT) then
      ! mesh 1
         limesh = 1
         do il=2,NFIELD
            ! varialbe T; temperature nad passive scalars
            ! this should depened on mesh type (T vs V)
            ntot  = lx1*ly1*lz1*NELT
            ntoto = lx1*ly1*lz1*AMR_NELT_O

            call copy(vr,T(1,1,1,1,il-1),ntoto)
            call amr_refine_coarsen_var(vr,vi,tmpv,limesh,
     $                                 LX1,LY1,LZ1,lbuf)
            call copy(T(1,1,1,1,il-1),vr,ntot)

            ! arrays time integration VGRADT[12]
            call copy(vr,VGRADT1(1,1,1,1,il-1),ntoto)
            call amr_refine_coarsen_var(vr,vi,tmpv,limesh,
     $                                 LX1,LY1,LZ1,lbuf)
            call copy(VGRADT1(1,1,1,1,il-1),vr,ntot)

            call copy(vr,VGRADT2(1,1,1,1,il-1),ntoto)
            call amr_refine_coarsen_var(vr,vi,tmpv,limesh,
     $                                 LX1,LY1,LZ1,lbuf)
            call copy(VGRADT2(1,1,1,1,il-1),vr,ntot)

            ! lag arrays TLAG
            do jl =1, NBDINP-1
               call copy(vr,TLAG(1,1,1,1,jl,il-1),ntoto)
               call amr_refine_coarsen_var(vr,vi,tmpv,limesh,
     $                                 LX1,LY1,LZ1,lbuf)
               call copy(TLAG(1,1,1,1,jl,il-1),vr,ntot)
            enddo
         enddo                  ! passive scalar loop
      endif                     ! IFHEAT

      ! mhd
      ! NOT TESTED AND NOT SURE IT IS OK!!!!!!!!!!!!!!!!
      if (IFMHD) then
      ! mesh 3
         limesh = 3
         ntot  = lbx1*lby1*lbz1*NELV
         ntoto = lbx1*lby1*lbz1*AMR_NELV_O

         ! magnetic field
         call copy(vr,BX,ntoto)
         call amr_refine_coarsen_var(vr,vi,tmpb,limesh,
     $                                 LBX1,LBY1,LBZ1,lbuf)
         call copy(BX,vr,ntot)

         call copy(vr,BY,ntoto)
         call amr_refine_coarsen_var(vr,vi,tmpb,limesh,
     $                                 LBX1,LBY1,LBZ1,lbuf)
         call copy(BY,vr,ntot)

         if (IF3D) then
            call copy(vr,BZ,ntoto)
            call amr_refine_coarsen_var(vr,vi,tmpb,limesh,
     $                                 LBX1,LBY1,LBZ1,lbuf)
            call copy(BZ,vr,ntot)
         endif

         ! arrays for time integration BB[XYZ][12]
         call copy(vr,BBX1,ntoto)
         call amr_refine_coarsen_var(vr,vi,tmpb,limesh,
     $                                 LBX1,LBY1,LBZ1,lbuf)
         call copy(BBX1,vr,ntot)

         call copy(vr,BBY1,ntoto)
         call amr_refine_coarsen_var(vr,vi,tmpb,limesh,
     $                                 LBX1,LBY1,LBZ1,lbuf)
         call copy(BBY1,vr,ntot)

         if (IF3D) then
            call copy(vr,BBZ1,ntoto)
            call amr_refine_coarsen_var(vr,vi,tmpb,limesh,
     $                                 LBX1,LBY1,LBZ1,lbuf)
            call copy(BBZ1,vr,ntot)
         endif

         call copy(vr,BBX2,ntoto)
         call amr_refine_coarsen_var(vr,vi,tmpb,limesh,
     $                                 LBX1,LBY1,LBZ1,lbuf)
         call copy(BBX2,vr,ntot)

         call copy(vr,BBY2,ntoto)
         call amr_refine_coarsen_var(vr,vi,tmpb,limesh,
     $                                 LBX1,LBY1,LBZ1,lbuf)
         call copy(BBY2,vr,ntot)

         if (IF3D) then
            call copy(vr,BBZ2,ntoto)
            call amr_refine_coarsen_var(vr,vi,tmpb,limesh,
     $                                 LBX1,LBY1,LBZ1,lbuf)
            call copy(BBZ2,vr,ntot)
         endif

         ! lag arrays
         do il =1, NBDINP-1
            call copy(vr,BXLAG(1,il),ntoto)
            call amr_refine_coarsen_var(vr,vi,tmpb,limesh,
     $                                 LBX1,LBY1,LBZ1,lbuf)
            call copy(BXLAG(1,il),vr,ntot)

            call copy(vr,BYLAG(1,il),ntoto)
            call amr_refine_coarsen_var(vr,vi,tmpb,limesh,
     $                                 LBX1,LBY1,LBZ1,lbuf)
            call copy(BYLAG(1,il),vr,ntot)

            if (IF3D) then
               call copy(vr,BZLAG(1,il),ntoto)
               call amr_refine_coarsen_var(vr,vi,tmpb,limesh,
     $                                 LBX1,LBY1,LBZ1,lbuf)
               call copy(BZLAG(1,il),vr,ntot)
            endif
         enddo

         ! mesh 4
         limesh = 4
         ntot  = lbx2*lby2*lbz2*NELV
         ntoto = lbx2*lby2*lbz2*AMR_NELV_O

         call copy(vr,PM,ntoto)
         call amr_refine_coarsen_var(vr,vi,tmpbp,limesh,
     $                                 LBX2,LBY2,LBZ2,lbuf)
         call copy(PM,vr,ntot)

         ! lag arrays
         if (NBDINP.eq.3) then
            call copy(vr,PMLAG,ntoto)
            call amr_refine_coarsen_var(vr,vi,tmpbp,limesh,
     $                                 LBX2,LBY2,LBZ2,lbuf)
            call copy(PMLAG,vr,ntot)
         endif
      endif


      ! moving mesh; BM1LAG
      if (IFGEOM) then
         call amr_abort
     $        ('Error: amr_refine_local_v no moving mesh yet.')
      endif

      ! This could be place for perturbation, but interpolation would not
      ! give correct base flow structure, so I skip it for now

      return
      end subroutine
!=======================================================================
