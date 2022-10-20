!> @file amr_dssum.f
!! @ingroup nekamr
!! @brief Routines relater to nonconforming stiffness summation
!! @author Adam Peplinski
!! @date Mar 6, 2020
!=======================================================================
!> @brief Setup nonconforming interpolation arrays
      subroutine amr_set_noncon()
      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'TSTEP'
      include 'PARALLEL'
      include 'AMR'
      include 'AMR_NONCON'

      ! local variables
      integer ntot, nel
!-----------------------------------------------------------------------
      ! get interpolation matrices for direct stiffness summation
      call amr_set_jmat()

      ! multigrid solver; THIS SHOULD BE MOVED TO OTHER PLACE
      if (int(param(43)).eq.0) call amr_set_jmat_hsmg()

      ! here I set mjultiplication arrays fof Jt and Jinv
      if (ifheat) then
        ifield = 2
        nel = nelt
      else
        ifield = 1
        nel = nelv
      endif
      ntot = lx1*ly1*lz1*nel
      call rone(AMR_JTMLT,ntot)
      call amr_apply_jt(AMR_JTMLT,nx1,ny1,nz1,nel)
      call fgslib_gs_op(gsh_fld(ifield),AMR_JTMLT,1,1,0)
      call invcol1 (AMR_JTMLT,ntot)

      call rone(AMR_JIMLT,ntot)
      call amr_apply_ji(AMR_JIMLT,nx1,ny1,nz1,nel)
      call fgslib_gs_op(gsh_fld(ifield),AMR_JIMLT,1,1,0)
      call invcol1 (AMR_JIMLT,ntot)

      return
      end subroutine
!=======================================================================
!> @brief Setup nonconforming interpolation arrays for dssum
      subroutine amr_set_jmat()
      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'AMR'
      include 'AMR_REFINE'
      include 'AMR_TOPOL'
      include 'AMR_NONCON'

      ! local variables
      integer iel, ifc, im, il, jl
      integer lnx2
      integer lposx, lposy
!-----------------------------------------------------------------------
      ! reset counter
      im = 0

      if (IF3D) then
         do iel=1,nelt
            ! hanging element
            if (AMR_HNG_ELM(iel).eq.1) then
               ! face loop
               do ifc=1,AMR_NFCS
                  ! hanging face
                  lposx = AMR_HNG_FCS(ifc,iel)
                  if (lposx.ne.-1) then
                     ! mortar element pointer
                     im=im+1
                     amr_mrtr_fc(ifc,iel) = im
                     amr_nonc_fcs_f(im) = ifc
                     amr_nonc_fcs_e(im) = iel
                     ! face position
                     lposy = lposx/2 +1
                     lposx = mod(lposx,2) +1
                     ! face orientaation; I assume relative face aligmnment is correctly
                     ! treated in node numbering routine, so nothing special here
                     amr_ifjmt(im) = .FALSE.
                     ! copy operator
                     if (ifc.eq.1.or.ifc.eq.3) then ! xz face
                        ! GLL points; velocity temperature
                        ! X
                        lnx2 = NX1*NX1
                        call copy(amr_jm_fc(1,1,1,im),
     $                            IXAMR1CF(1,1,lposx),lnx2)
                        ! Z
                        lnx2 = NZ1*NZ1
                        call copy(amr_jm_fc(1,1,2,im),
     $                            IZTAMR1CF(1,1,lposy),lnx2)
                        ! inverse interpolation
                        ! X
                        lnx2 = NX1*NX1
                        call copy(amr_jm_fci(1,1,1,im),
     $                            IXAMR1FC(1,1,lposx),lnx2)
                        ! Z
                        lnx2 = NZ1*NZ1
                        call copy(amr_jm_fci(1,1,2,im),
     $                            IZTAMR1FC(1,1,lposy),lnx2)
                        ! GL points on GLL elements; presure preconditioner
                        if (IFSPLIT) then ! P-N-P_N
                           call amr_abort('Unsuported, amr_set_jmat')
                        else
                           ! X
                           lnx2 = NX1*NX1
                           call rzero(amr_jmlcf2(1,1,1,im),lnx2)
                           call rzero(amr_jmlfc2(1,1,1,im),lnx2)
                           do jl = 1,NX2
                              do il = 1,NX2
                                 amr_jmlcf2(il+1,jl+1,1,im) =
     $                                IXAMR2CF(il,jl,lposx)
                                 amr_jmlfc2(il+1,jl+1,1,im) =
     $                                IXAMR2FC(il,jl,lposx)
                              enddo
                           enddo
                           ! Z
                           lnx2 = NZ1*NZ1
                           call rzero(amr_jmlcf2(1,1,2,im),lnx2)
                           call rzero(amr_jmlfc2(1,1,2,im),lnx2)
                           do jl = 1,NZ2
                              do il = 1,NZ2
                                 amr_jmlcf2(il+1,jl+1,2,im) =
     $                                IZTAMR2CF(il,jl,lposy)
                                 amr_jmlfc2(il+1,jl+1,2,im) =
     $                                IZTAMR2FC(il,jl,lposy)
                              enddo
                           enddo

                        endif
                        ! crs grid (vertices only); presure preconditioner
                        lnx2 = AMR_LXC*AMR_LXC
                        ! X
                        call copy(amr_jmc_fc(1,1,1,im),
     $                            IXAMRCCF(1,1,lposx),lnx2)
                        ! Z
                        call copy(amr_jmc_fc(1,1,2,im),
     $                            IXTAMRCCF(1,1,lposy),lnx2)

                     elseif (ifc.eq.2.or.ifc.eq.4) then ! yz face
                        ! GLL points; velocity temperature
                        ! Y
                        lnx2 = NY1*NY1
                        call copy(amr_jm_fc(1,1,1,im),
     $                            IYAMR1CF(1,1,lposx),lnx2)
                        ! Z
                        lnx2 = NZ1*NZ1
                        call copy(amr_jm_fc(1,1,2,im),
     $                            IZTAMR1CF(1,1,lposy),lnx2)
                        ! inverse interpolation
                        ! Y
                        lnx2 = NY1*NY1
                        call copy(amr_jm_fci(1,1,1,im),
     $                            IYAMR1FC(1,1,lposx),lnx2)
                        ! Z
                        lnx2 = NZ1*NZ1
                        call copy(amr_jm_fci(1,1,2,im),
     $                            IZTAMR1FC(1,1,lposy),lnx2)
                        ! GL points on GLL elements; presure preconditioner
                        if (IFSPLIT) then ! P-N-P_N
                           call amr_abort('Unsuported, amr_set_jmat')
                        else
                           ! Y
                           lnx2 = NY1*NY1
                           call rzero(amr_jmlcf2(1,1,1,im),lnx2)
                           call rzero(amr_jmlfc2(1,1,1,im),lnx2)
                           do jl = 1,NY2
                              do il = 1,NY2
                                 amr_jmlcf2(il+1,jl+1,1,im) =
     $                                IYAMR2CF(il,jl,lposx)
                                 amr_jmlfc2(il+1,jl+1,1,im) =
     $                                IYAMR2FC(il,jl,lposx)
                              enddo
                           enddo

!     Z
                           lnx2 = NZ1*NZ1
                           call rzero(amr_jmlcf2(1,1,2,im),lnx2)
                           call rzero(amr_jmlfc2(1,1,2,im),lnx2)
                           do jl = 1,NZ2
                              do il = 1,NZ2
                                 amr_jmlcf2(il+1,jl+1,2,im) =
     $                                IZTAMR2CF(il,jl,lposy)
                                 amr_jmlfc2(il+1,jl+1,2,im) =
     $                                IZTAMR2FC(il,jl,lposy)
                              enddo
                           enddo

                        endif
                        ! crs grid (vertices only); presure preconditioner
                        lnx2 = AMR_LXC*AMR_LXC
                        ! Y
                        call copy(amr_jmc_fc(1,1,1,im),
     $                            IXAMRCCF(1,1,lposx),lnx2)
                        ! Z
                        call copy(amr_jmc_fc(1,1,2,im),
     $                            IXTAMRCCF(1,1,lposy),lnx2)

                     elseif (ifc.eq.5.or.ifc.eq.6) then ! xy face
                        ! GLL points; velocity temperature
                        ! X
                        lnx2 = NX1*NX1
                        call copy(amr_jm_fc(1,1,1,im),
     $                            IXAMR1CF(1,1,lposx),lnx2)
                        ! Y
                        lnx2 = NY1*NY1
                        call copy(amr_jm_fc(1,1,2,im),
     $                            IYTAMR1CF(1,1,lposy),lnx2)
                        ! inverse interpolation
                        ! X
                        lnx2 = NX1*NX1
                        call copy(amr_jm_fci(1,1,1,im),
     $                            IXAMR1FC(1,1,lposx),lnx2)
                        ! Y
                        lnx2 = NY1*NY1
                        call copy(amr_jm_fci(1,1,2,im),
     $                            IYTAMR1FC(1,1,lposy),lnx2)
                        ! GL points on GLL elements; presure preconditioner
                        if (IFSPLIT) then ! P-N-P_N
                           call amr_abort('Unsuported, amr_set_jmat')
                        else
                           ! X
                           lnx2 = NX1*NX1
                           call rzero(amr_jmlcf2(1,1,1,im),lnx2)
                           call rzero(amr_jmlfc2(1,1,1,im),lnx2)
                           do jl = 1,NX2
                              do il = 1,NX2
                                 amr_jmlcf2(il+1,jl+1,1,im) =
     $                                IXAMR2CF(il,jl,lposx)
                                 amr_jmlfc2(il+1,jl+1,1,im) =
     $                                IXAMR2FC(il,jl,lposx)
                              enddo
                           enddo

                           ! Y
                           lnx2 = NY1*NY1
                           call rzero(amr_jmlcf2(1,1,2,im),lnx2)
                           call rzero(amr_jmlfc2(1,1,2,im),lnx2)
                           do jl = 1,NY2
                              do il = 1,NY2
                                 amr_jmlcf2(il+1,jl+1,2,im) =
     $                                IYTAMR2CF(il,jl,lposy)
                                 amr_jmlfc2(il+1,jl+1,1,im) =
     $                                IYTAMR2FC(il,jl,lposy)
                              enddo
                           enddo

                        endif
                        ! crs grid (vertices only); presure preconditioner
                        lnx2 = AMR_LXC*AMR_LXC
                        ! X
                        call copy(amr_jmc_fc(1,1,1,im),
     $                            IXAMRCCF(1,1,lposx),lnx2)
                        ! Y
                        call copy(amr_jmc_fc(1,1,2,im),
     $                            IXTAMRCCF(1,1,lposy),lnx2)

                     endif
                  else
                     amr_mrtr_fc(ifc,iel)=0
                  endif
               enddo
            else
               do ifc=1,AMR_NFCS
                  amr_mrtr_fc(ifc,iel)=0
               enddo
            endif
         enddo
         ! save number of hanging faces
         amr_mrtr_fcs_m = im

         ! edge operators
         im = 0
         do iel=1,nelt
            ! hanging element
            if (AMR_HNG_ELM(iel).eq.1) then
               ! loop over edges
               do ifc=1,AMR_NEDG
                  ! hanging edge
                  ! I check only hanging edges inrelated to hanging faces
                  lposx = AMR_HNG_EDG(ifc,iel)
                  if (lposx.eq.0.or.lposx.eq.1) then
                     ! mortar element pointer
                     im=im+1
                     amr_mrtr_ed(ifc,iel) = im
                     amr_nonc_edg_f(im) = ifc
                     amr_nonc_edg_e(im) = iel
                     ! edge position
                     lposx = lposx + 1
                     ! like for faces I assume all orientations are consitent at this stage
                     ! (corrected by proper node numbering)
                     if(ifc.lt.5) then ! edges along x
                        lnx2 = NX1*NX1
                        call copy(amr_jm_ed(1,1,im),
     $                            IXAMR1CF(1,1,lposx),lnx2)
                        ! inverse interpolation
                        call copy(amr_jm_edi(1,1,im),
     $                            IXAMR1FC(1,1,lposx),lnx2)
                        ! crs grid (vertices only); presure preconditioner
                        lnx2 = AMR_LXC*AMR_LXC
                        call copy(amr_jmc_ed(1,1,im),
     $                            IXAMRCCF(1,1,lposx),lnx2)
                     else if(ifc.lt.9) then ! edges along y
                        lnx2 = NY1*NY1
                        call copy(amr_jm_ed(1,1,im),
     $                            IYAMR1CF(1,1,lposx),lnx2)
                        ! inverse interpolation
                        call copy(amr_jm_edi(1,1,im),
     $                            IYAMR1FC(1,1,lposx),lnx2)
                        ! crs grid (vertices only); presure preconditioner
                        lnx2 = AMR_LXC*AMR_LXC
                        call copy(amr_jmc_ed(1,1,im),
     $                            IXAMRCCF(1,1,lposx),lnx2)
                     else ! edges along z
                        lnx2 = NZ1*NZ1
                        call copy(amr_jm_ed(1,1,im),
     $                            IZAMR1CF(1,1,lposx),lnx2)
                        ! inverse interpolation
                        call copy(amr_jm_edi(1,1,im),
     $                            IZAMR1FC(1,1,lposx),lnx2)
                        ! crs grid (vertices only); presure preconditioner
                        lnx2 = AMR_LXC*AMR_LXC
                        call copy(amr_jmc_ed(1,1,im),
     $                            IXAMRCCF(1,1,lposx),lnx2)
                     endif
                  else
                     amr_mrtr_ed(ifc,iel) = 0
                  endif
               enddo
            else
               do ifc=1,12
                  amr_mrtr_ed(ifc,iel) = 0
               enddo
            endif
         enddo
         ! save number of hangign edges
         amr_mrtr_edg_m = im
      else
         do iel=1,nelt
            ! hanging element
            if (AMR_HNG_ELM(iel).eq.1) then
               ! face loop
               do ifc=1,AMR_NFCS
                  ! hanging face
                  lposx = AMR_HNG_FCS(ifc,iel)
                  if (lposx.ne.-1) then
                     ! mortar element pointer
                     im=im+1
                     amr_mrtr_fc(ifc,iel) = im
                     amr_nonc_fcs_f(im) = ifc
                     amr_nonc_fcs_e(im) = iel
                     ! face position
                     lposx = lposx +1
                     ! face orientaation; I assume relative face aligmnment is correctly
                     ! treated in node numbering routine, so nothing special here
                     amr_ifjmt(im) = .FALSE.
                     ! copy operator
                     if (ifc.eq.1.or.ifc.eq.3) then ! x face
                        ! GLL points; velocity temperature
                        ! X
                        lnx2 = NX1*NX1
                        call copy(amr_jm_fc(1,1,1,im),
     $                            IXAMR1CF(1,1,lposx),lnx2)
                        ! inverse interpolation
                        call copy(amr_jm_fci(1,1,1,im),
     $                            IXAMR1FC(1,1,lposx),lnx2)
                        ! GL points on GLL elements; presure preconditioner
                        if (IFSPLIT) then ! P-N-P_N
                           call amr_abort('Unsuported, amr_set_jmat')
                        else
                           ! X
                           lnx2 = NX1*NX1
                           call rzero(amr_jmlcf2(1,1,1,im),lnx2)
                           call rzero(amr_jmlfc2(1,1,1,im),lnx2)
                           do jl = 1,NX2
                              do il = 1,NX2
                                 amr_jmlcf2(il+1,jl+1,1,im) =
     $                                IXAMR2CF(il,jl,lposx)
                                 amr_jmlfc2(il+1,jl+1,1,im) =
     $                                IXAMR2FC(il,jl,lposx)
                              enddo
                           enddo

                        endif
                        ! crs grid (vertices only); presure preconditioner
                        lnx2 = AMR_LXC*AMR_LXC
                        ! X
                        call copy(amr_jmc_fc(1,1,1,im),
     $                            IXAMRCCF(1,1,lposx),lnx2)

                     elseif (ifc.eq.2.or.ifc.eq.4) then ! y face
                        ! GLL points; velocity temperature
                        ! Y
                        lnx2 = NY1*NY1
                        call copy(amr_jm_fc(1,1,1,im),
     $                            IYAMR1CF(1,1,lposx),lnx2)
                        ! inverse interpolation
                        call copy(amr_jm_fci(1,1,1,im),
     $                            IYAMR1FC(1,1,lposx),lnx2)
                        ! GL points on GLL elements; presure preconditioner
                        if (IFSPLIT) then ! P-N-P_N
                           call amr_abort('Unsuported, amr_set_jmat')
                        else
                           ! Y
                           lnx2 = NY1*NY1
                           call rzero(amr_jmlcf2(1,1,1,im),lnx2)
                           call rzero(amr_jmlfc2(1,1,1,im),lnx2)
                           do jl = 1,NY2
                              do il = 1,NY2
                                 amr_jmlcf2(il+1,jl+1,1,im) =
     $                                IYAMR2CF(il,jl,lposx)
                                 amr_jmlfc2(il+1,jl+1,1,im) =
     $                                IYAMR2FC(il,jl,lposx)
                              enddo
                           enddo

                        endif
                        ! crs grid (vertices only); presure preconditioner
                        lnx2 = AMR_LXC*AMR_LXC
                        ! Y
                        call copy(amr_jmc_fc(1,1,1,im),
     $                            IXAMRCCF(1,1,lposx),lnx2)

                     endif
                  else
                     amr_mrtr_fc(ifc,iel)=0
                  endif
               enddo
            else
               do ifc=1,AMR_NFCS
                  amr_mrtr_fc(ifc,iel)=0
               enddo
            endif
         enddo
         ! save number of hanging faces
         amr_mrtr_fcs_m = im
      endif

      return
      end subroutine
!=======================================================================
!> @brief Setup nonconforming interpolation arrays for multigrid solver
!     NOTICE!!! This routine depend on the variables set in set_jmat, so
!     it has to be called just after it
      subroutine amr_set_jmat_hsmg()
      implicit none

      include 'SIZE'
      include 'HSMG'
      include 'INPUT'
      include 'AMR'
      include 'AMR_TOPOL'
      include 'AMR_NONCON'
      include 'AMR_REFINE_HSMG'

      ! local variables
      integer iel, ifc, im, il, jl, kl
      integer lnx2, lnx
      integer lposx, lposy
      integer ierr
!-----------------------------------------------------------------------
      ! reset error mark
      ierr = 0

      if (IF3D) then
         ! loop over hanging faces
         do im=1, amr_mrtr_fcs_m
            ifc = amr_nonc_fcs_f(im)
            iel = amr_nonc_fcs_e(im)
            ! hanging face position
            lposx = AMR_HNG_FCS(ifc,iel)

            lposy = lposx/2 +1
            lposx = mod(lposx,2) +1
            ! copy operator
            ! I assume here (like in the rest of hsmg code) that LX1.eq.LY1.eq.LZ1,
            ! so I do not distinguis between different faces. Code taking into
            ! account different direction can be found in set_jmat
            do il=2,MG_LMAX-1 ! loop ove levels
               lnx = MG_NH(il)
               lnx2 = lnx*lnx
               call copy(amr_jm_fc_mg(1,1,im,il),
     $                   IXAMRHCF(1,lposx,il),lnx2)
               call copy(amr_jm_fc_mg(1,2,im,il),
     $                   IXTAMRHCF(1,lposy,il),lnx2)
            enddo
            ! for Schwarz operator
            do il=2,MG_LMAX-1 ! loop ove levels
               lnx = MG_NH(il)+2
               lnx2 = lnx*lnx
               call rzero(amr_jmlCF_mg(1,1,im,il),lnx2)
               call rzero(amr_jmlFC_mg(1,1,im,il),lnx2)
               call rzero(amr_jmlCF_mg(1,2,im,il),lnx2)
               call rzero(amr_jmlFC_mg(1,2,im,il),lnx2)
               do kl=1,MG_NH(il)
                  do jl=1,MG_NH(il)
                     amr_jmlCF_mg(jl+1+kl*lnx,1,im,il) =
     $                       IXAMRHCF(jl+(kl-1)*MG_NH(il),lposx,il)
                     amr_jmlFC_mg(jl+1+kl*lnx,1,im,il) =
     $                       IXAMRHFC(jl+(kl-1)*MG_NH(il),lposx,il)
                     amr_jmlCF_mg(jl+1+kl*lnx,2,im,il) =
     $                       IXTAMRHCF(jl+(kl-1)*MG_NH(il),lposy,il)
                     amr_jmlFC_mg(jl+1+kl*lnx,2,im,il) =
     $                       IXTAMRHFC(jl+(kl-1)*MG_NH(il),lposy,il)
                  enddo
               enddo
            enddo
         enddo ! loop over hanging faces

         ! edge operators
         ! loop over hanging edges
         do im=1, amr_mrtr_edg_m
            ifc = amr_nonc_edg_f(im)
            iel = amr_nonc_edg_e(im)
            ! hanging face position
            lposx = AMR_HNG_EDG(ifc,iel)

            lposx = lposx +1
            ! copy operator
            ! I assume here (like in the rest of hsmg code) that LX1.eq.LY1, so
            ! I do not distinguis between different faces. Code taking into
            ! account different direction can be found in set_jmat
            do il=2,MG_LMAX-1 ! loop ove levels
               lnx = MG_NH(il)
               lnx2 = lnx*lnx
               call copy(amr_jm_ed_mg(1,im,il),
     $                   IXAMRHCF(1,lposx,il),lnx2)
            enddo
         enddo ! loop over hanging edges

      else ! if3d
         !loop over hanging faces
         do im=1, amr_mrtr_fcs_m
            ifc = amr_nonc_fcs_f(im)
            iel = amr_nonc_fcs_e(im)
            ! hanging face position
            lposx = AMR_HNG_FCS(ifc,iel)

            lposx = lposx +1
            ! copy operator
            ! I assume here (like in the rest of hsmg code) that LX1.eq.LY1, so
            ! I do not distinguis between different faces. Code taking into
            ! account different direction can be found in set_jmat
            do il=2,MG_LMAX-1 ! loop ove levels
               lnx = MG_NH(il)
               lnx2 = lnx*lnx
               call copy(amr_jm_fc_mg(1,1,im,il),
     $                   IXAMRHCF(1,lposx,il),lnx2)
            enddo
            ! for Schwarz operator
            do il=2,MG_LMAX-1 ! loop ove levels
               lnx = MG_NH(il)+2
               lnx2 = lnx*lnx
               call rzero(amr_jmlCF_mg(1,1,im,il),lnx2)
               call rzero(amr_jmlFC_mg(1,1,im,il),lnx2)
               do kl=1,MG_NH(il)
                  do jl=1,MG_NH(il)
                     amr_jmlCF_mg(jl+1+kl*lnx,1,im,il) =
     $                       IXAMRHCF(jl+(kl-1)*MG_NH(il),lposx,il)
                     amr_jmlFC_mg(jl+1+kl*lnx,1,im,il) =
     $                       IXAMRHFC(jl+(kl-1)*MG_NH(il),lposx,il)
                  enddo
               enddo
            enddo
         enddo ! loop over hanging faces
      endif ! if3d

      return
      end subroutine
!=======================================================================
!> @brief Apply transposed interpolation operator (translating child data)
!! @param[inout]    ul           field array
!! @param[in]       nx, ny, nz   element size
!! @param[in]       nel          element number
      subroutine amr_apply_jt(ul,nx,ny,nz,nel)
      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'AMR'
      include 'AMR_TOPOL'
      include 'AMR_NONCON'

      ! argument list
      integer nx, ny, nz, nel
      real ul(nx,ny,nz,nel)

      ! local variables
      integer ie, iface, iedge    ! loop index
      integer im
      integer lface
      parameter (lface=lx1*lz1)
      real uin(lface,2*ldim),uout(lface)
      common /nonctmp/ uin, uout
!-----------------------------------------------------------------------
      ! face part
      do ie = 1 , nel
         if (AMR_HNG_ELM(ie).eq.1) then
           ! Note, we zero out u() on this face after extracting, for
           ! consistency reasons discovered during Jerry's thesis.
           ! Thus,  "ftovec_0" rather than ftovec().   (iface -- Ed notation)
           do iface = 1, AMR_NFCS
              im = amr_mrtr_fc(iface,ie)
              if (im.ne.0) then
                 call ftovec_0l(uin(1,iface),ul(1,1,1,ie),
     $                          iface,nx,ny,nz)
              endif
           enddo
           do iface=1, AMR_NFCS
              im = amr_mrtr_fc(iface,ie)
              if (im.ne.0) then
                 if (if3d) then
                   call matvec3t(uout,amr_jm_fc(1,1,1,im),uin(1,iface),
     $                           amr_ifjmt(im),nx,nx)
                 else
                   call matvect(uout,amr_jm_fc(1,1,1,im),uin(1,iface),
     $                          nx,nx)
                 endif
                 call vectof_addl(ul(1,1,1,ie),uout,iface,nx,ny,nz)
              endif
           enddo
         endif
      enddo

      ! edge part
      if (IF3D) then
        do ie = 1, nel
          if (AMR_HNG_ELM(ie).eq.1) then
            do iedge = 1, AMR_NEDG
                im = amr_mrtr_ed(iedge,ie)
                if (im.ne.0) then
                    call etovec(uin,iedge,ul(1,1,1,ie),nx,ny,nz)
                    call matvect (uout,amr_jm_ed(1,1,im),uin,nx,nx)
                    call vectoe(ul(1,1,1,ie),iedge,uout,nx,ny,nz)
                endif
            enddo
          endif
        enddo
      endif

      return
      end subroutine
!=======================================================================
!> @brief Apply interpolation operator (parent onto children)
!! @param[inout]    ul           field array
!! @param[in]       nx, ny, nz   element size
!! @param[in]       nel          element number
      subroutine amr_apply_j(ul,nx,ny,nz,nel)
      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'AMR'
      include 'AMR_TOPOL'
      include 'AMR_NONCON'

      ! argument list
      integer nx, ny, nz, nel
      real ul(nx,ny,nz,nel)

      ! local variables
      integer ie, iface, iedge    ! loop index
      integer im
      integer lface
      parameter (lface=lx1*lz1)
      real uin(lface,2*ldim),uout(lface)
      common /nonctmp/ uin, uout
!-----------------------------------------------------------------------
      ! face part
      do ie = 1, nel
         if (AMR_HNG_ELM(ie).eq.1) then
           do iface = 1, AMR_NFCS
              im = amr_mrtr_fc(iface,ie)
              if (im.ne.0) then
                 call ftovecl(uin(1,iface),ul(1,1,1,ie),iface,nx,ny,nz)
              endif
           enddo
           do iface=1, AMR_NFCS
              im = amr_mrtr_fc(iface,ie)
              if (im.ne.0) then
                 call matvec3(uout,amr_jm_fc(1,1,1,im),uin(1,iface),
     $                        amr_ifjmt(im),nx,nz)
                 call vectofl(ul(1,1,1,ie),uout,iface,nx,ny,nz)
              endif
           enddo
         endif
      enddo

      ! edge part
      if (IF3D) then
        do ie = 1, nel
          if (AMR_HNG_ELM(ie).eq.1) then
            do iedge = 1, AMR_NEDG
                im = amr_mrtr_ed(iedge,ie)
                if (im.ne.0) then
                    call etovec(uin,iedge,ul(1,1,1,ie),nx,ny,nz)
                    call mxm (amr_jm_ed(1,1,im),nx,uin,nx,uout,1)
                    call vectoe(ul(1,1,1,ie),iedge,uout,nx,ny,nz)
                endif
            enddo
          endif
        enddo
      endif

      return
      end subroutine
!=======================================================================
!> @brief Apply inverse interpolation operator (children onto parent)
!! @param[inout]    ul           field array
!! @param[in]       nx, ny, nz   element size
!! @param[in]       nel          element number
      subroutine amr_apply_ji(ul,nx,ny,nz,nel)
      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'AMR'
      include 'AMR_TOPOL'
      include 'AMR_REFINE'
      include 'AMR_NONCON'

      ! argument list
      integer nx, ny, nz, nel
      real ul(nx,ny,nz,nel)

      ! local variables
      integer ie, iface, iedge    ! loop index
      integer im
      integer lface, nface
      parameter (lface=lx1*lz1)
      real uin(lface,2*ldim),uout(lface)
      common /nonctmp/ uin, uout
!-----------------------------------------------------------------------
      ! face part
      nface = nx*nz
      do ie = 1, nel
         if (AMR_HNG_ELM(ie).eq.1) then
           do iface = 1, AMR_NFCS
              im = amr_mrtr_fc(iface,ie)
              if (im.ne.0) then
                 call ftovecl(uin(1,iface),ul(1,1,1,ie),iface,nx,ny,nz)
              endif
           enddo
           do iface=1, AMR_NFCS
              im = amr_mrtr_fc(iface,ie)
              if (im.ne.0) then
                 call matvec3(uout,amr_jm_fci(1,1,1,im),uin(1,iface),
     $                        amr_ifjmt(im),nx,nz)
                 call col2(uout,IMAMR1F,nface)
                 call vectofl(ul(1,1,1,ie),uout,iface,nx,ny,nz)
              endif
           enddo
         endif
      enddo

      ! edge part
      if (IF3D) then
        do ie = 1, nel
          if (AMR_HNG_ELM(ie).eq.1) then
            do iedge = 1, AMR_NEDG
                im = amr_mrtr_ed(iedge,ie)
                if (im.ne.0) then
                    call etovec(uin,iedge,ul(1,1,1,ie),nx,ny,nz)
                    call mxm (amr_jm_edi(1,1,im),nx,uin,nx,uout,1)
                    call col2(uout,IMAMR1E,nx)
                    call vectoe(ul(1,1,1,ie),iedge,uout,nx,ny,nz)
                endif
            enddo
          endif
        enddo
      endif

      return
      end subroutine
!=======================================================================
!> @brief Apply transposed interpolation operator for Schwarz
!! @param[inout]    ul           field array
!! @param[in]       jm           interpolation matrix
!! @param[in]       nx, ny, nz   element size
!! @param[in]       nel          element number
      subroutine amr_apply_jt_sch(ul,jm,nx,ny,nz,nel)
      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'AMR'
      include 'AMR_TOPOL'
      include 'AMR_NONCON'

      ! argument list
      integer nx, ny, nz, nel
      real ul(nx,ny,nz,nel)
      real jm(nx,nx,2,maxmor)

      ! local variables
      integer ie, iface   ! loop index
      integer im
      integer lface
      parameter (lface=lx1*lz1)
      real uin(lface,2*ldim),uout(lface)
      common /nonctmp/ uin, uout
!-----------------------------------------------------------------------
      ! face part
      do ie = 1 , nel
         if (AMR_HNG_ELM(ie).eq.1) then
           ! Note, we zero out u() on this face after extracting, for
           ! consistency reasons discovered during Jerry's thesis.
           ! Thus,  "ftovec_0" rather than ftovec().   (iface -- Ed notation)
           do iface = 1, AMR_NFCS
              im = amr_mrtr_fc(iface,ie)
              if (im.ne.0) then
                 call ftovec_0l(uin(1,iface),ul(1,1,1,ie),
     $                          iface,nx,ny,nz)
              endif
           enddo
           do iface=1, AMR_NFCS
              im = amr_mrtr_fc(iface,ie)
              if (im.ne.0) then
                 if (if3d) then
                   call matvec3t(uout,jm(1,1,1,im),uin(1,iface),
     $                           amr_ifjmt(im),nx,nx)
                 else
                   call matvect(uout,jm(1,1,1,im),uin(1,iface),
     $                          nx,nx)
                 endif
                 call vectof_addl(ul(1,1,1,ie),uout,iface,nx,ny,nz)
              endif
           enddo
         endif
      enddo

      ! we work on GL points so no edge part
      return
      end subroutine
!=======================================================================
!> @brief Apply interpolation operator for Schwarz
!! @param[inout]    ul           field array
!! @param[in]       jm           interpolation matrix
!! @param[in]       nx, ny, nz   element size
!! @param[in]       nel          element number
      subroutine amr_apply_j_sch(ul,jm,nx,ny,nz,nel)
      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'AMR'
      include 'AMR_TOPOL'
      include 'AMR_NONCON'

      ! argument list
      integer nx, ny, nz, nel
      real ul(nx,ny,nz,nel)
      real jm(nx,nx,2,maxmor)

      ! local variables
      integer ie, iface   ! loop index
      integer im
      integer lface
      parameter (lface=lx1*lz1)
      real uin(lface,2*ldim),uout(lface)
      common /nonctmp/ uin, uout
!-----------------------------------------------------------------------
      ! face part
      do ie = 1, nel
         if (AMR_HNG_ELM(ie).eq.1) then
           do iface = 1, AMR_NFCS
              im = amr_mrtr_fc(iface,ie)
              if (im.ne.0) then
                 call ftovecl(uin(1,iface),ul(1,1,1,ie),iface,nx,ny,nz)
              endif
           enddo
           do iface=1, AMR_NFCS
              im = amr_mrtr_fc(iface,ie)
              if (im.ne.0) then
                 call matvec3(uout,jm(1,1,1,im),uin(1,iface),
     $                        amr_ifjmt(im),nx,nz)
                 call vectofl(ul(1,1,1,ie),uout,iface,nx,ny,nz)
              endif
           enddo
         endif
      enddo

      ! we work on GL points so no edge part
      return
      end subroutine
!=======================================================================
!> @brief Apply inverse interpolation operator for Schwarz
!! @param[inout]    ul           field array
!! @param[in]       jm           interpolation matrix
!! @param[in]       nx, ny, nz   element size
!! @param[in]       nel          element number
      subroutine amr_apply_ji_sch(ul,jm,nx,ny,nz,nel)
      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'AMR'
      include 'AMR_TOPOL'
      include 'AMR_REFINE'
      include 'AMR_NONCON'

      ! argument list
      integer nx, ny, nz, nel
      real ul(nx,ny,nz,nel)
      real jm(nx,nx,2,maxmor)

      ! local variables
      integer ie, iface   ! loop index
      integer im
      integer lface, nface
      parameter (lface=lx1*lz1)
      real uin(lface,2*ldim),uout(lface)
      common /nonctmp/ uin, uout
!-----------------------------------------------------------------------
      ! face part
      nface = nx*nz
      do ie = 1, nel
         if (AMR_HNG_ELM(ie).eq.1) then
           do iface = 1, AMR_NFCS
              im = amr_mrtr_fc(iface,ie)
              if (im.ne.0) then
                 call ftovecl(uin(1,iface),ul(1,1,1,ie),iface,nx,ny,nz)
              endif
           enddo
           do iface=1, AMR_NFCS
              im = amr_mrtr_fc(iface,ie)
              if (im.ne.0) then
                 call matvec3(uout,jm(1,1,1,im),uin(1,iface),
     $                        amr_ifjmt(im),nx,nz)
                 call col2(uout,IMAMR1F,nface)
                 call vectofl(ul(1,1,1,ie),uout,iface,nx,ny,nz)
              endif
           enddo
         endif
      enddo

      ! we work on GL points so no edge part
      return
      end subroutine
!=======================================================================
!> @brief Apply transposed interpolation operator (single element)
!! @param[inout]    ul           field array
!! @param[in]       jmf          interpolation matrix face
!! @param[in]       jme          interpolation matrix edge
!! @param[in]       nx, ny, nz   element size
!! @param[in]       iel          element number
      subroutine amr_apply_jt_loc(ul,jmf,jme,nx,ny,nz,iel)
      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'AMR'
      include 'AMR_TOPOL'
      include 'AMR_NONCON'

      ! argument list
      integer nx, ny, nz, iel
      real ul(nx,ny,nz)
      real jmf(nx,nx,2,maxmor),jme(nx,nx,maxmor)

      ! local variables
      integer iedge, iface   ! loop index
      integer im
      integer lface
      parameter (lface=lx1*lz1)
      real uin(lface,2*ldim),uout(lface)
      common /nonctmp/ uin, uout
!-----------------------------------------------------------------------
      ! face part
      if (AMR_HNG_ELM(iel).eq.1) then
        ! Note, we zero out u() on this face after extracting, for
        ! consistency reasons discovered during Jerry's thesis.
        ! Thus,  "ftovec_0" rather than ftovec().   (iface -- Ed notation)
        do iface = 1, AMR_NFCS
           im = amr_mrtr_fc(iface,iel)
           if (im.ne.0) then
              call ftovec_0l(uin(1,iface),ul,iface,nx,ny,nz)
           endif
        enddo
        do iface=1, AMR_NFCS
           im = amr_mrtr_fc(iface,iel)
           if (im.ne.0) then
              if (if3d) then
                call matvec3t(uout,jmf(1,1,1,im),uin(1,iface),
     $                        amr_ifjmt(im),nx,nx)
              else
                call matvect(uout,jmf(1,1,1,im),uin(1,iface),nx,nx)
              endif
              call vectof_addl(ul,uout,iface,nx,ny,nz)
           endif
        enddo
      endif

      ! edge part
      if (IF3D) then
        if (AMR_HNG_ELM(iel).eq.1) then
          do iedge = 1, AMR_NEDG
              im = amr_mrtr_ed(iedge,iel)
              if (im.ne.0) then
                  call etovec(uin,iedge,ul,nx,ny,nz)
                  call matvect (uout,jme(1,1,im),uin,nx,nx)
                  call vectoe(ul,iedge,uout,nx,ny,nz)
              endif
          enddo
        endif
      endif

      return
      end subroutine
!=======================================================================
!> @brief Apply interpolation operator (single element)
!! @param[inout]    ul           field array
!! @param[in]       jmf          interpolation matrix face
!! @param[in]       jme          interpolation matrix edge
!! @param[in]       nx, ny, nz   element size
!! @param[in]       iel          element number
      subroutine amr_apply_j_loc(ul,jmf,jme,nx,ny,nz,iel)
      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'AMR'
      include 'AMR_TOPOL'
      include 'AMR_NONCON'

      ! argument list
      integer nx, ny, nz, iel
      real ul(nx,ny,nz)
      real jmf(nx,nx,2,maxmor),jme(nx,nx,maxmor)

      ! local variables
      integer iedge, iface   ! loop index
      integer im
      integer lface
      parameter (lface=lx1*lz1)
      real uin(lface,2*ldim),uout(lface)
      common /nonctmp/ uin, uout
!-----------------------------------------------------------------------
      ! face part
      if (AMR_HNG_ELM(iel).eq.1) then
        do iface = 1, AMR_NFCS
           im = amr_mrtr_fc(iface,iel)
           if (im.ne.0) then
              call ftovecl(uin(1,iface),ul,iface,nx,ny,nz)
           endif
        enddo
        do iface=1, AMR_NFCS
           im = amr_mrtr_fc(iface,iel)
           if (im.ne.0) then
              call matvec3(uout,jmf(1,1,1,im),uin(1,iface),
     $                     amr_ifjmt(im),nx,nz)
              call vectofl(ul,uout,iface,nx,ny,nz)
           endif
        enddo
      endif

      ! edge part
      if (IF3D) then
        if (AMR_HNG_ELM(iel).eq.1) then
          do iedge = 1, AMR_NEDG
            im = amr_mrtr_ed(iedge,iel)
            if (im.ne.0) then
               call etovec(uin,iedge,ul,nx,ny,nz)
               call mxm (jme(1,1,im),nx,uin,nx,uout,1)
               call vectoe(ul,iedge,uout,nx,ny,nz)
            endif
          enddo
        endif
      endif

      return
      end subroutine
!=======================================================================
!> @brief Zero field value at children faces
!! @param[inout]    ul           field array
!! @param[in]       nx, ny, nz   element size
!! @param[in]       nel          element number
      subroutine amr_zero_ch(ul,nx,ny,nz,nel)
      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'TSTEP'
      include 'AMR'
      include 'AMR_TOPOL'

      ! argument list
      integer nx,ny,nz,nel
      real ul(nx,ny,nz,nel)

      ! local variables
      integer el,fc,if_hng
      real zero
!-----------------------------------------------------------------------
      zero  = 0.0

      ! face
      do el=1,nelfld(IFIELD)
        if (AMR_HNG_ELM(el).ne.0) then
           do fc=1,AMR_NFCS
              if (AMR_HNG_FCS(fc,el).ne.-1) then
                 call facev(ul,el,fc,zero,nx,ny,nz)
              endif
           enddo
        endif
      enddo
      ! edge
      if (IF3D) then
        do el=1,nelfld(IFIELD)
          if (AMR_HNG_ELM(el).ne.0) then
             do fc=1,AMR_NEDG
                if_hng = AMR_HNG_EDG(fc,el)
                if (if_hng.eq.0.or.if_hng.eq.1) then
                  call ctoe(ul(1,1,1,el),fc,zero,nx,ny,nz)
                endif
             enddo
          endif
        enddo
      endif

      return
      end subroutine
!=======================================================================
!> @brief Fill nonconforming children with constant
!! @param[inout]    vec          array
!! @param[in]       nv           vector length
!! @param[in]       cnst         constant value
      subroutine amr_fill_ch(vec,nv,cnst)
      implicit none
      include 'SIZE'
      include 'AMR'
      include 'AMR_NONCON'

      ! argument list
      real vec(nv),cnst
      integer nv

      ! local variables
      integer il
!-----------------------------------------------------------------------
      do il=1,nv
         if (amr_qmlt(il).eq.0) vec(il)=cnst
      enddo

      return
      end subroutine
!=======================================================================
!> @brief Project u on children using mask
!! @param[inout]    u            field array
!! @param[in]       mask         mask array
!! @param[in]       nx, ny, nz   element size
!! @param[in]       nel          element number
      subroutine amr_h1_proj_msk(u,mask,nx,ny,nz,nel)
      implicit none

      include 'SIZE'
      include 'CTIMER'
      include 'INPUT'
      include 'PARALLEL'
      include 'TSTEP'
      include 'MVGEOM' ! wmult
      include 'SOLN'   ! vmult, tmult

      ! argument list
      integer nx,ny,nz,nel
      real u(nx,ny,nz,nel),mask(nx,ny,nz,nel)

      ! scratch arrays
      integer lface
      parameter (lface=lx1*ly1)
      real uin(lface,2*ldim),uout(lface)
      common /nonctmp/ uin, uout

      ! local variables
      integer ntot
      real timee
!-----------------------------------------------------------------------
      if(ifsync) call nekgsync()

#ifdef TIMER
      if (icalld.eq.0) then
         tdsmx=0.
         tdsmn=0.
      endif
      icalld=icalld+1
      etime1=dnekclock()
#endif

      ntot = nx*ny*nz*nel
      !                    ~  ~T
      ! Implement   :=   J Q  Q  Mu
      !
      !             T
      if (ifield.eq.0) then ! mesh
         call col2(u,wmult,ntot)
      elseif  (ifield.eq.1) then ! velocity
         call col2(u,vmult,ntot)
      elseif  (ifield.ge.2.and.ifield.le.nfield) then ! temperature and passive scalar
         call col2(u,tmult(1,1,1,1,ifield-1),ntot)
         call col2(u,tmult(1,1,1,1,ifield-1),ntot)
      endif
      !             ~ ~T
      ! This is the Q Q  part
      !
      call fgslib_gs_op(gsh_fld(ifield),u,1,1,0)
      call col2   (u,mask,ntot)
      !
      ! This is the J  part,  interpolating parent solution onto child
      !
      call apply_J(u,nx,ny,nz,nel)

#ifdef TIMER
      timee=(dnekclock()-etime1)
      tdsum=tdsum+timee
      ndsum=icalld
      tdsmx=max(timee,tdsmx)
      tdsmn=min(timee,tdsmn)
#endif
      return
      end subroutine
!=======================================================================
!> @brief Direct stiffness summation for Schwarz operator
!! @param[inout]    ul           field array
!! @param[in]       nx, ny, nz   element size
!! @param[in]       nel          element number
      subroutine amr_dssum_sch(ul,nx,ny,nz,nel)
      implicit none

      include 'SIZE'
      include 'CTIMER'
      include 'INPUT'
      include 'PARALLEL'
      include 'TSTEP'
      include 'AMR'
      include 'AMR_NONCON'

      ! argument list
      integer nx, ny, nz, nel
      real ul(nx,ny,nz,nel)

      ! local variables
      integer ifldt

      ! declare timer
      real timee
!-----------------------------------------------------------------------
      ifldt = ifield
      if (ifldt.eq.ifldmhd) ifldt = 1

      if (ifsync) call nekgsync()

#ifdef TIMER
      if (icalld.eq.0) then
         tdsmx=0.
         tdsmn=0.
      endif
      icalld=icalld+1
      etime1=dnekclock()
#endif
      !              T     -1
      ! This is the J  or J   part,  translating/interpolating child data to parent
      ! for Schwarz operations J^-1 is more effcient than J^T, but it gives
      ! not symmetric operator, so for PCG J^T must be used
      if (int(param(42)).eq.1) then
         call amr_apply_jt_sch(ul,amr_jmlCF2,nx,ny,nz,nel)
      else
         call amr_apply_ji_sch(ul,amr_jmlFC2,nx,ny,nz,nel)
      endif
      !             ~ ~T
      ! This is the Q Q  part
      call fgslib_gs_op(gsh_fld(ifldt),ul,1,1,0)  ! 1 ==> +

      ! This is the J  part,  interpolating parent solution onto child
      call amr_apply_j_sch(ul,amr_jmlCF2,nx,ny,nz,nel)

#ifdef TIMER
      timee=(dnekclock()-etime1)
      tdsum=tdsum+timee
      ndsum=icalld
      tdsmx=max(timee,tdsmx)
      tdsmn=min(timee,tdsmn)
#endif

      return
      end subroutine
!=======================================================================
!> @brief Direct stiffness summation on vector field with mask
!! @param[inout]    va,vb,vc     vector field
!! @param[in]       ma,mb,mc     mask arrays
!! @param[in]       nx, ny, nz   element size
!! @param[in]       nel          element number
      subroutine amr_opdssum_msk(va,vb,vc,ma,mb,mc,nx,ny,nz,nel)
      implicit none

      include 'SIZE'
      include 'INPUT'

      ! argumen list
      integer nx,ny,nz, nel
      real va(nx,ny,nz, nel),vb(nx,ny,nz, nel),vc(nx,ny,nz, nel)
      real ma(nx,ny,nz, nel),mb(nx,ny,nz, nel),mc(nx,ny,nz, nel)
!-----------------------------------------------------------------------
      if (ifcyclic) then
         call rotate_cyc(va,vb,vc,1)
         call dssum_msk(va,ma,nx,ny,nz)
         call dssum_msk(vb,mb,nx,ny,nz)
         if (IF3D) call dssum_msk(vc,mc,nx,ny,nz)
         call rotate_cyc  (va,vb,vc,0)
      else
         call dssum_msk(va,ma,nx,ny,nz)
         call dssum_msk(vb,mb,nx,ny,nz)
         if (IF3D) call dssum_msk(vc,mc,nx,ny,nz)
      endif

      return
      end subroutine
!=======================================================================
!> @brief Project vector field on children
!! @param[inout]    va,vb,vc     vector field
!! @param[in]       nx, ny, nz   element size
!! @param[in]       nel          element number
      subroutine amr_oph1_proj(va,vb,vc,nx,ny,nz,nel)
      implicit none

      include 'SIZE'
      include 'INPUT'

      ! argumen list
      integer nx,ny,nz, nel
      real va(nx,ny,nz, nel),vb(nx,ny,nz, nel),vc(nx,ny,nz, nel)
!-----------------------------------------------------------------------
      if (ifcyclic) then
         call rotate_cyc(va,vb,vc,1)
         call h1_proj(va,nx,ny,nz)
         call h1_proj(vb,nx,ny,nz)
         if (IF3D) call h1_proj(vc,nx,ny,nz)
         call rotate_cyc  (va,vb,vc,0)
      else
         call h1_proj(va,nx,ny,nz)
         call h1_proj(vb,nx,ny,nz)
         if (IF3D) call h1_proj(vc,nx,ny,nz)
      endif

      return
      end subroutine
!=======================================================================
