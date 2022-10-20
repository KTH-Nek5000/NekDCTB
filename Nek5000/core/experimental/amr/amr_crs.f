!> @file amr_crs.f
!! @ingroup nekamr
!! @brief Coarse grid solver related routines
!! @author Adam Peplinski
!! @date Nov 28, 2016
!=======================================================================
!> @brief Generate crs base function.
!! @details Generate coarse base function to calculate coarse preconditioner.
!! @note This routine assumes nx1=ny1 in 2D and nx1=ny1=nz1 in 3D
!! @todo Move it to the other file
      subroutine amr_gen_crs_basis()
      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'AMR'
      include 'AMR_CRS'

      ! locall variables
      ! 1D profiles
      real z01(lx1),z10(lx1),z1h(lx1),zh1(lx1),z0h(lx1),zh0(lx1)
      real zr(lx1),zs(lx1),zt(lx1),zrm(lx1),zsm(lx1),ztm(lx1)

      integer il, jl, kl, ml  ! loop index
!-----------------------------------------------------------------------
      ! get 1D profiles
      call zwgll(zr,zs,nx1)

      do il=1,nx1
         z10(il) = 0.5*(1.0-zr(il)) !   1-->0
         z01(il) = 0.5*(1.0+zr(il)) !   0-->1
         z0h(il) = 0.5*z01(il)      !   0-->0.5
         zh0(il) = 0.5*z10(il)      ! 0.5-->0
         z1h(il) = z10(il)+z0h(il)  !   1-->0.5
         zh1(il) = z01(il)+zh0(il)  ! 0.5-->1
      enddo

      ! depending on grid dimension
      if (IF3D) then
         ! loop over all verticies
         do il = 1, AMR_NVRT
            ! different base functions

            ! combinations of: 1-->0 (r) and 1-->0 (s) and 1-->0 (t)
            ! non-hanging node, no hanging faces sharing the node
            ! (any not refined element and refined concave corner)
            if (mod(il,2).eq.0) then
               call copy(zr,z01,nx1)
            else
               call copy(zr,z10,nx1)
            endif
            if (il.eq.3.or.il.eq.4.or.il.eq.7.or.il.eq.8) then
               call copy(zs,z01,nx1)
            else
               call copy(zs,z10,nx1)
            endif
            if (il.gt.4) then
               call copy(zt,z01,nx1)
            else
               call copy(zt,z10,nx1)
            endif

            do jl=1,nz1
               do kl=1,ny1
                  do ml=1,nx1
                     AMR_CRS_BASE(ml,kl,jl,1,il) = zr(ml)*zs(kl)*zt(jl)
                  enddo
               enddo
            enddo

            ! combinations of: 1-->0.5 and 1-->0.5 and 1-->0.5
            ! non-hanging node, hanging faces perpendicular to r and s and t
            ! (refined convex corner)
            if (mod(il,2).eq.0) then
               call copy(zr,zh1,nx1)
               call copy(zrm,z10,nx1)
            else
               call copy(zr,z1h,nx1)
               call copy(zrm,z01,nx1)
            endif
            if (il.eq.3.or.il.eq.4.or.il.eq.7.or.il.eq.8) then
               call copy(zs,zh1,nx1)
               call copy(zsm,z10,nx1)
            else
               call copy(zs,z1h,nx1)
               call copy(zsm,z01,nx1)
            endif
            if (il.gt.4) then
               call copy(zt,zh1,nx1)
               call copy(ztm,z10,nx1)
            else
               call copy(zt,z1h,nx1)
               call copy(ztm,z01,nx1)
            endif

            do jl=1,nz1
               do kl=1,ny1
                  do ml=1,nx1
                     AMR_CRS_BASE(ml,kl,jl,2,il) = zr(ml)*zs(kl)*zt(jl)
     $                 - 0.125*zrm(ml)*zsm(kl)*ztm(jl)
                  enddo
               enddo
            enddo

            ! combinations of: 1-->0.5 and 1-->0.5 and 1-->0.5
            ! non-hanging node, hanging faces perpendicular to r and s
            ! (refined convex edge)
            if (mod(il,2).eq.0) then
               call copy(zr,zh1,nx1)
               call copy(zrm,z10,nx1)
            else
               call copy(zr,z1h,nx1)
               call copy(zrm,z01,nx1)
            endif
            if (il.eq.3.or.il.eq.4.or.il.eq.7.or.il.eq.8) then
               call copy(zs,zh1,nx1)
               call copy(zsm,z10,nx1)
            else
               call copy(zs,z1h,nx1)
               call copy(zsm,z01,nx1)
            endif
            if (il.gt.4) then
               call copy(zt,zh1,nx1)
               call copy(ztm,zh1,nx1)
            else
               call copy(zt,z1h,nx1)
               call copy(ztm,z1h,nx1)
            endif

            do jl=1,nz1
               do kl=1,ny1
                  do ml=1,nx1
                     AMR_CRS_BASE(ml,kl,jl,3,il) = zr(ml)*zs(kl)*zt(jl)
     $                 - 0.25*zrm(ml)*zsm(kl)*ztm(jl)
                  enddo
               enddo
            enddo

            ! combinations of: 1-->0.5 and 1-->0.5 and 1-->0.5
            ! non-hanging node, hanging faces perpendicular to r and t
            ! (refined convex edge)
            if (mod(il,2).eq.0) then
               call copy(zr,zh1,nx1)
               call copy(zrm,z10,nx1)
            else
               call copy(zr,z1h,nx1)
               call copy(zrm,z01,nx1)
            endif
            if (il.eq.3.or.il.eq.4.or.il.eq.7.or.il.eq.8) then
               call copy(zs,zh1,nx1)
               call copy(zsm,zh1,nx1)
            else
               call copy(zs,z1h,nx1)
               call copy(zsm,z1h,nx1)
            endif
            if (il.gt.4) then
               call copy(zt,zh1,nx1)
               call copy(ztm,z10,nx1)
            else
               call copy(zt,z1h,nx1)
               call copy(ztm,z01,nx1)
            endif

            do jl=1,nz1
               do kl=1,ny1
                  do ml=1,nx1
                     AMR_CRS_BASE(ml,kl,jl,4,il) = zr(ml)*zs(kl)*zt(jl)
     $                 - 0.25*zrm(ml)*zsm(kl)*ztm(jl)
                  enddo
               enddo
            enddo

            ! combinations of: 1-->0.5 and 1-->0.5 and 1-->0.5
            ! non-hanging node, hanging faces perpendicular to s and t
            ! (refined convex edge)
            if (mod(il,2).eq.0) then
               call copy(zr,zh1,nx1)
               call copy(zrm,zh1,nx1)
            else
               call copy(zr,z1h,nx1)
               call copy(zrm,z1h,nx1)
            endif
            if (il.eq.3.or.il.eq.4.or.il.eq.7.or.il.eq.8) then
               call copy(zs,zh1,nx1)
               call copy(zsm,z10,nx1)
            else
               call copy(zs,z1h,nx1)
               call copy(zsm,z01,nx1)
            endif
            if (il.gt.4) then
               call copy(zt,zh1,nx1)
               call copy(ztm,z10,nx1)
            else
               call copy(zt,z1h,nx1)
               call copy(ztm,z01,nx1)
            endif

            do jl=1,nz1
               do kl=1,ny1
                  do ml=1,nx1
                     AMR_CRS_BASE(ml,kl,jl,5,il) = zr(ml)*zs(kl)*zt(jl)
     $                 - 0.25*zrm(ml)*zsm(kl)*ztm(jl)
                  enddo
               enddo
            enddo

            ! combinations of: 1-->0.5 and 1-->0.5 and 0.5-->0
            ! hanging node, hanging faces perpendicular to r and s
            ! (refined convex edge)
            if (mod(il,2).eq.0) then
               call copy(zr,zh1,nx1)
               call copy(zrm,z10,nx1)
            else
               call copy(zr,z1h,nx1)
               call copy(zrm,z01,nx1)
            endif
            if (il.eq.3.or.il.eq.4.or.il.eq.7.or.il.eq.8) then
               call copy(zs,zh1,nx1)
               call copy(zsm,z10,nx1)
            else
               call copy(zs,z1h,nx1)
               call copy(zsm,z01,nx1)
            endif
            if (il.gt.4) then
               call copy(zt,z0h,nx1)
               call copy(ztm,z0h,nx1)
            else
               call copy(zt,zh0,nx1)
               call copy(ztm,zh0,nx1)
            endif

            do jl=1,nz1
               do kl=1,ny1
                  do ml=1,nx1
                     AMR_CRS_BASE(ml,kl,jl,6,il) = zr(ml)*zs(kl)*zt(jl)
     $                 - 0.25*zrm(ml)*zsm(kl)*ztm(jl)
                  enddo
               enddo
            enddo

            ! combinations of: 1-->0.5 and 0.5-->0 and 1-->0.5
            ! hanging node, hanging faces perpendicular to r and t
            ! (refined convex edge)
            if (mod(il,2).eq.0) then
               call copy(zr,zh1,nx1)
               call copy(zrm,z10,nx1)
            else
               call copy(zr,z1h,nx1)
               call copy(zrm,z01,nx1)
            endif
            if (il.eq.3.or.il.eq.4.or.il.eq.7.or.il.eq.8) then
               call copy(zs,z0h,nx1)
               call copy(zsm,z0h,nx1)
            else
               call copy(zs,zh0,nx1)
               call copy(zsm,zh0,nx1)
            endif
            if (il.gt.4) then
               call copy(zt,zh1,nx1)
               call copy(ztm,z10,nx1)
            else
               call copy(zt,z1h,nx1)
               call copy(ztm,z01,nx1)
            endif

            do jl=1,nz1
               do kl=1,ny1
                  do ml=1,nx1
                     AMR_CRS_BASE(ml,kl,jl,7,il) = zr(ml)*zs(kl)*zt(jl)
     $                 - 0.25*zrm(ml)*zsm(kl)*ztm(jl)
                  enddo
               enddo
            enddo

            ! combinations of: 0.5-->0 and 1-->0.5 and 1-->0.5
            ! hanging node, hanging faces perpendicular to s and t
            ! (refined convex edge)
            if (mod(il,2).eq.0) then
               call copy(zr,z0h,nx1)
               call copy(zrm,z0h,nx1)
            else
               call copy(zr,zh0,nx1)
               call copy(zrm,zh0,nx1)
            endif
            if (il.eq.3.or.il.eq.4.or.il.eq.7.or.il.eq.8) then
               call copy(zs,zh1,nx1)
               call copy(zsm,z10,nx1)
            else
               call copy(zs,z1h,nx1)
               call copy(zsm,z01,nx1)
            endif
            if (il.gt.4) then
               call copy(zt,zh1,nx1)
               call copy(ztm,z10,nx1)
            else
               call copy(zt,z1h,nx1)
               call copy(ztm,z01,nx1)
            endif

            do jl=1,nz1
               do kl=1,ny1
                  do ml=1,nx1
                     AMR_CRS_BASE(ml,kl,jl,8,il) = zr(ml)*zs(kl)*zt(jl)
     $                 - 0.25*zrm(ml)*zsm(kl)*ztm(jl)
                  enddo
               enddo
            enddo

            ! combinations of: 1-->0 and 1-->0.5 and 1-->0.5
            ! non-hanging node, hanging face perpendicular to r
            ! (refined face)
            if (mod(il,2).eq.0) then
               call copy(zr,z01,nx1)
            else
               call copy(zr,z10,nx1)
            endif
            if (il.eq.3.or.il.eq.4.or.il.eq.7.or.il.eq.8) then
               call copy(zs,zh1,nx1)
            else
               call copy(zs,z1h,nx1)
            endif
            if (il.gt.4) then
               call copy(zt,zh1,nx1)
            else
               call copy(zt,z1h,nx1)
            endif

            do jl=1,nz1
               do kl=1,ny1
                  do ml=1,nx1
                     AMR_CRS_BASE(ml,kl,jl,9,il) = zr(ml)*zs(kl)*zt(jl)
                  enddo
               enddo
            enddo

            ! combinations of: 1-->0.5 and 1-->0 and 1-->0.5
            ! non-hanging node, hanging face perpendicular to s
            ! (refined face)
            if (mod(il,2).eq.0) then
               call copy(zr,zh1,nx1)
            else
               call copy(zr,z1h,nx1)
            endif
            if (il.eq.3.or.il.eq.4.or.il.eq.7.or.il.eq.8) then
               call copy(zs,z01,nx1)
            else
               call copy(zs,z10,nx1)
            endif
            if (il.gt.4) then
               call copy(zt,zh1,nx1)
            else
               call copy(zt,z1h,nx1)
            endif

            do jl=1,nz1
               do kl=1,ny1
                  do ml=1,nx1
                     AMR_CRS_BASE(ml,kl,jl,10,il) =zr(ml)*zs(kl)*zt(jl)
                  enddo
               enddo
            enddo

            ! combinations of: 1-->0.5 and 1-->0.5 and 1-->0
            ! non-hanging node, hanging face perpendicular to t
            ! (refined face)
            if (mod(il,2).eq.0) then
               call copy(zr,zh1,nx1)
            else
               call copy(zr,z1h,nx1)
            endif
            if (il.eq.3.or.il.eq.4.or.il.eq.7.or.il.eq.8) then
               call copy(zs,zh1,nx1)
            else
               call copy(zs,z1h,nx1)
            endif
            if (il.gt.4) then
               call copy(zt,z01,nx1)
            else
               call copy(zt,z10,nx1)
            endif

            do jl=1,nz1
               do kl=1,ny1
                  do ml=1,nx1
                     AMR_CRS_BASE(ml,kl,jl,11,il) =zr(ml)*zs(kl)*zt(jl)
                  enddo
               enddo
            enddo

            ! combinations of: 1-->0 and 1-->0.5 and 0.5-->0
            ! hanging node, hanging face perpendicular to r or t
            ! (refined face)
            if (mod(il,2).eq.0) then
               call copy(zr,z01,nx1)
            else
               call copy(zr,z10,nx1)
            endif
            if (il.eq.3.or.il.eq.4.or.il.eq.7.or.il.eq.8) then
               call copy(zs,zh1,nx1)
            else
               call copy(zs,z1h,nx1)
            endif
            if (il.gt.4) then
               call copy(zt,z0h,nx1)
            else
               call copy(zt,zh0,nx1)
            endif

            do jl=1,nz1
               do kl=1,ny1
                  do ml=1,nx1
                     AMR_CRS_BASE(ml,kl,jl,12,il) =zr(ml)*zs(kl)*zt(jl)
                  enddo
               enddo
            enddo

            ! combinations of: 1-->0.5 and 1-->0 and 0.5-->0
            ! hanging node, hanging face perpendicular to s or t
            ! (refined face)
            if (mod(il,2).eq.0) then
               call copy(zr,zh1,nx1)
            else
               call copy(zr,z1h,nx1)
            endif
            if (il.eq.3.or.il.eq.4.or.il.eq.7.or.il.eq.8) then
               call copy(zs,z01,nx1)
            else
               call copy(zs,z10,nx1)
            endif
            if (il.gt.4) then
               call copy(zt,z0h,nx1)
            else
               call copy(zt,zh0,nx1)
            endif

            do jl=1,nz1
               do kl=1,ny1
                  do ml=1,nx1
                     AMR_CRS_BASE(ml,kl,jl,13,il) =zr(ml)*zs(kl)*zt(jl)
                  enddo
               enddo
            enddo

            ! combinations of: 1-->0 and 0.5-->0 and 1-->0.5
            ! hanging node, hanging face perpendicular to r or s
            ! (refiend face)
            if (mod(il,2).eq.0) then
               call copy(zr,z01,nx1)
            else
               call copy(zr,z10,nx1)
            endif
            if (il.eq.3.or.il.eq.4.or.il.eq.7.or.il.eq.8) then
               call copy(zs,z0h,nx1)
            else
               call copy(zs,zh0,nx1)
            endif
            if (il.gt.4) then
               call copy(zt,zh1,nx1)
            else
               call copy(zt,z1h,nx1)
            endif

            do jl=1,nz1
               do kl=1,ny1
                  do ml=1,nx1
                     AMR_CRS_BASE(ml,kl,jl,14,il) =zr(ml)*zs(kl)*zt(jl)
                  enddo
               enddo
            enddo


            ! combinations of: 1-->0 and 0.5-->0 and 0.5-->0
            ! hanging node, hanging face perpendicular to r or s or t
            ! (refined face)
            if (mod(il,2).eq.0) then
               call copy(zr,z01,nx1)
            else
               call copy(zr,z10,nx1)
            endif
            if (il.eq.3.or.il.eq.4.or.il.eq.7.or.il.eq.8) then
               call copy(zs,z0h,nx1)
            else
               call copy(zs,zh0,nx1)
            endif
            if (il.gt.4) then
               call copy(zt,z0h,nx1)
            else
               call copy(zt,zh0,nx1)
            endif

            do jl=1,nz1
               do kl=1,ny1
                  do ml=1,nx1
                     AMR_CRS_BASE(ml,kl,jl,15,il) =zr(ml)*zs(kl)*zt(jl)
                  enddo
               enddo
            enddo

            ! combinations of: 1-->0 and 1-->0 and 1-->0.5
            ! non-hanging node, no hanging faces, hanging edge along t
            ! (refined concave edge)
            if (mod(il,2).eq.0) then
               call copy(zr,z01,nx1)
            else
               call copy(zr,z10,nx1)
            endif
            if (il.eq.3.or.il.eq.4.or.il.eq.7.or.il.eq.8) then
               call copy(zs,z01,nx1)
            else
               call copy(zs,z10,nx1)
            endif
            if (il.gt.4) then
               call copy(zt,zh1,nx1)
            else
               call copy(zt,z1h,nx1)
            endif

            do jl=1,nz1
               do kl=1,ny1
                  do ml=1,nx1
                     AMR_CRS_BASE(ml,kl,jl,16,il) = zr(ml)*zs(kl)*zt(jl)
                  enddo
               enddo
            enddo

            ! combinations of: 1-->0 and 1-->0.5 and 1-->0
            ! non-hanging node, no hanging faces, hanging edge along s
            ! (refined concave edge)
            if (mod(il,2).eq.0) then
               call copy(zr,z01,nx1)
            else
               call copy(zr,z10,nx1)
            endif
            if (il.eq.3.or.il.eq.4.or.il.eq.7.or.il.eq.8) then
               call copy(zs,zh1,nx1)
            else
               call copy(zs,z1h,nx1)
            endif
            if (il.gt.4) then
               call copy(zt,z01,nx1)
            else
               call copy(zt,z10,nx1)
            endif

            do jl=1,nz1
               do kl=1,ny1
                  do ml=1,nx1
                     AMR_CRS_BASE(ml,kl,jl,17,il) = zr(ml)*zs(kl)*zt(jl)
                  enddo
               enddo
            enddo

            ! combinations of: 1-->0.5 and 1-->0 and 1-->0
            ! non-hanging node, no hanging faces, hanging edge along r
            ! (refined concave edge)
            if (mod(il,2).eq.0) then
               call copy(zr,zh1,nx1)
            else
               call copy(zr,z1h,nx1)
            endif
            if (il.eq.3.or.il.eq.4.or.il.eq.7.or.il.eq.8) then
               call copy(zs,z01,nx1)
            else
               call copy(zs,z10,nx1)
            endif
            if (il.gt.4) then
               call copy(zt,z01,nx1)
            else
               call copy(zt,z10,nx1)
            endif

            do jl=1,nz1
               do kl=1,ny1
                  do ml=1,nx1
                     AMR_CRS_BASE(ml,kl,jl,18,il) = zr(ml)*zs(kl)*zt(jl)
                  enddo
               enddo
            enddo

            ! combinations of: 1-->0 and 1-->0 and 0.5-->0
            ! hanging node, no hanging faces, hanging edge along t
            ! (refined concave edge)
            if (mod(il,2).eq.0) then
               call copy(zr,z01,nx1)
            else
               call copy(zr,z10,nx1)
            endif
            if (il.eq.3.or.il.eq.4.or.il.eq.7.or.il.eq.8) then
               call copy(zs,z01,nx1)
            else
               call copy(zs,z10,nx1)
            endif
            if (il.gt.4) then
               call copy(zt,z0h,nx1)
            else
               call copy(zt,zh0,nx1)
            endif

            do jl=1,nz1
               do kl=1,ny1
                  do ml=1,nx1
                     AMR_CRS_BASE(ml,kl,jl,19,il) = zr(ml)*zs(kl)*zt(jl)
                  enddo
               enddo
            enddo

            ! combinations of: 1-->0 and 0.5-->0 and 1-->0
            ! hanging node, no hanging faces, hanging edge along s
            ! (refined concave edge)
            if (mod(il,2).eq.0) then
               call copy(zr,z01,nx1)
            else
               call copy(zr,z10,nx1)
            endif
            if (il.eq.3.or.il.eq.4.or.il.eq.7.or.il.eq.8) then
               call copy(zs,z0h,nx1)
            else
               call copy(zs,zh0,nx1)
            endif
            if (il.gt.4) then
               call copy(zt,z01,nx1)
            else
               call copy(zt,z10,nx1)
            endif

            do jl=1,nz1
               do kl=1,ny1
                  do ml=1,nx1
                     AMR_CRS_BASE(ml,kl,jl,20,il) = zr(ml)*zs(kl)*zt(jl)
                  enddo
               enddo
            enddo

            ! combinations of: 0.5-->0 and 1-->0 and 1-->0
            ! hanging node, no hanging faces, hanging edge along r
            ! (refined concave edge)
            if (mod(il,2).eq.0) then
               call copy(zr,z0h,nx1)
            else
               call copy(zr,zh0,nx1)
            endif
            if (il.eq.3.or.il.eq.4.or.il.eq.7.or.il.eq.8) then
               call copy(zs,z01,nx1)
            else
               call copy(zs,z10,nx1)
            endif
            if (il.gt.4) then
               call copy(zt,z01,nx1)
            else
               call copy(zt,z10,nx1)
            endif

            do jl=1,nz1
               do kl=1,ny1
                  do ml=1,nx1
                     AMR_CRS_BASE(ml,kl,jl,21,il) = zr(ml)*zs(kl)*zt(jl)
                  enddo
               enddo
            enddo

         enddo  ! vertex loop
      else ! IF3D
         ! loop over all verticies
         do il = 1, AMR_NVRT
            ! different base functions

            ! combinations of: 1-->0 (r) and 1-->0 (s)
            ! non-hanging node, no hanging faces sharing the node
            if (mod(il,2).eq.0) then
               call copy(zr,z01,nx1)
            else
               call copy(zr,z10,nx1)
            endif
            if (il.ge.3) then
               call copy(zs,z01,nx1)
            else
               call copy(zs,z10,nx1)
            endif

            do jl=1,ny1
               do kl=1,nx1
                  AMR_CRS_BASE(kl,jl,1,1,il) = zr(kl)*zs(jl)
               enddo
            enddo

            ! combinations of: 1-->0.5 and 1-->0.5
            ! non-hanging node, hanging faces perpendicular to r and s
            if (mod(il,2).eq.0) then
               call copy(zr,zh1,nx1)
               call copy(zrm,zh0,nx1)
            else
               call copy(zr,z1h,nx1)
               call copy(zrm,z0h,nx1)
            endif
            if (il.ge.3) then
               call copy(zs,zh1,nx1)
               call copy(zsm,zh0,nx1)
            else
               call copy(zs,z1h,nx1)
               call copy(zsm,z0h,nx1)
            endif

            do jl=1,ny1
               do kl=1,nx1
                  AMR_CRS_BASE(kl,jl,1,2,il) = zr(kl)*zs(jl)
     $                             - zrm(kl)*zsm(jl)
               enddo
            enddo

            ! combinations of: 1-->0.5 and 1-->0
            ! non-hanging node, hanging face perpendicular to s
            if (mod(il,2).eq.0) then
               call copy(zr,zh1,nx1)
            else
               call copy(zr,z1h,nx1)
            endif
            if (il.ge.3) then
               call copy(zs,z01,nx1)
            else
               call copy(zs,z10,nx1)
            endif

            do jl=1,ny1
               do kl=1,nx1
                  AMR_CRS_BASE(kl,jl,1,3,il) = zr(kl)*zs(jl)
               enddo
            enddo

            ! combinations of: 1-->0 and 1-->0.5
            ! non-hanging node, hanging face perpendicular to r
            if (mod(il,2).eq.0) then
               call copy(zr,z01,nx1)
            else
               call copy(zr,z10,nx1)
            endif
            if (il.ge.3) then
               call copy(zs,zh1,nx1)
            else
               call copy(zs,z1h,nx1)
            endif

            do jl=1,ny1
               do kl=1,nx1
                  AMR_CRS_BASE(kl,jl,1,4,il) = zr(kl)*zs(jl)
               enddo
            enddo

            ! combinations of: 0.5-->0 and 1-->0
            ! hanging node, hanging face perpendicular to r or s
            do jl=1,ny1
               do kl=1,nx1
                  AMR_CRS_BASE(kl,jl,1,5,il) = 0.5*
     $                                  AMR_CRS_BASE(kl,jl,1,1,il)
               enddo
            enddo

         enddo  ! vertex loop
      endif  ! IF3D

      return
      end subroutine
!=======================================================================
!> @brief Generate vertex to crs base mapping
!! @todo Move it to the other file
      subroutine amr_gen_crs_map()
      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'TOPOL'
      include 'AMR'
      include 'AMR_TOPOL'
      include 'AMR_CRS'

!     locall variables
      integer ie, iv, ieg  ! loop index
      integer ierr    ! error mark
      integer ir, is, it  ! hanging mark along r, s, t
      integer ivc   ! vertex number touching hanging face
      integer itmp  ! dummy variable

!     vertex to face mapping; symmetric notation
      integer ivface(3,8)
      save ivface
      data ivface /1,3,5,  2,3,5,  1,4,5,  2,4,5,
     $             1,3,6,  2,3,6,  1,4,6,  2,4,6/

!     face vertex to element verex mapping
      integer ivftove (0:3,6)
      save ivftove
      data ivftove / 1,3,5,7,  2,4,6,8,  1,2,5,6,
     $               3,4,7,8,  1,2,3,4,  5,6,7,8/

!     inverse of ICEDG to store edge orientation; 3D runs only
!     vertex-vertex pair => rst orinetation
      integer icedg_inv(8,8)
!     vertex to edge mapping including rst orientation
      integer ivetoedg(3,8)

!     functions
      integer iglmax
!-----------------------------------------------------------------------
!     for 3D invert ICEDG to get edge orientation
!     this should be done only once
      if (IF3D) then
         call ifill(icedg_inv,-1,64)
         do ie = 1,12 ! loop over edges
            icedg_inv(ICEDG(1,ie),ICEDG(2,ie)) = ICEDG(3,ie)
            icedg_inv(ICEDG(2,ie),ICEDG(1,ie)) = ICEDG(3,ie)
            ivetoedg(ICEDG(3,ie),ICEDG(1,ie)) = ie
            ivetoedg(ICEDG(3,ie),ICEDG(2,ie)) = ie
         enddo
      endif

!     reset mapping
      ie = AMR_NVRT*NELT
      call ione(AMR_CRS_MAP,ie)

      ierr = 0
!     depending on grid dimension
      if (IF3D) then
!     loop ove elements
         do ie = 1,nelt
            if (AMR_HNG_ELM(ie).gt.0) then ! any face or edge hanging?
!     loop over element verticies
               do iv = 1, AMR_NVRT
!     check hanging faces; remember neks doesn't use symmetric notation
                  ir = AMR_HNG_FCS(EFACE(ivface(1,iv)),ie)
                  is = AMR_HNG_FCS(EFACE(ivface(2,iv)),ie)
                  it = AMR_HNG_FCS(EFACE(ivface(3,iv)),ie)
                  ivc = -1
                  if (AMR_HNG_VRT(iv,ie).eq.2) then ! edge hanging node
                     select case (ir)
                         case(-1)! non-hanging r
                            select case(is)
                               case(-1) ! non-hanging s
                                  select case(it)
                                     case(-1) ! non-hanging t
!     look for concave edges
                                       ir = AMR_HNG_EDG
     $                                        (ivetoedg(1,iv),ie)
                                       is = AMR_HNG_EDG
     $                                        (ivetoedg(2,iv),ie)
                                       it = AMR_HNG_EDG
     $                                        (ivetoedg(3,iv),ie)

                                       select case(ir)
                                         case (-1)  ! not hanging
                                           select case (is)
                                             case (-1)  ! not hanging
                                               select case (it)
                                                 case (-1)  ! not hanging
                                                  ierr = 1
                                                 case (0:1)
                                                  AMR_CRS_MAP(iv,ie)=19
                                                 case default
                                                   ierr = 1
                                               end select
                                             case (0:1)
                                               AMR_CRS_MAP(iv,ie) =20
                                             case default
                                               ierr = 1
                                           end select
                                         case (0:1)
                                           AMR_CRS_MAP(iv,ie) =21
                                         case default
                                           ierr = 1
                                       end select

!     possibly hanging edge of the corner element; nothing to do
                                     case(0:3)  ! hanging t
!     find the edge direction
!     get vertex number touching hanging face
                                        ivc = ivftove(it,ivface(3,iv))
!     get edge orientation
                                        ivc = icedg_inv(iv,ivc)
                                        select case(ivc)
                                           case (1) ! edge along r
                                              AMR_CRS_MAP(iv,ie) = 12
                                           case (2) ! edge along s
                                              AMR_CRS_MAP(iv,ie) = 13
                                           case default
                                              ierr = 1
                                        end select
                                     case default
                                        ierr = 1
                                  end select
                               case(0:3) ! hanging s
                                  select case(it)
                                     case(-1) ! non-hanging t
!     find the edge direction
!     get vertex number touching hanging face
                                        ivc = ivftove(is,ivface(2,iv))
!     get edge orientation
                                        ivc = icedg_inv(iv,ivc)
                                        select case(ivc)
                                           case (1) ! edge along r
                                              AMR_CRS_MAP(iv,ie) = 14
                                           case (3) ! edge along t
                                              AMR_CRS_MAP(iv,ie) = 13
                                           case default
                                              ierr = 1
                                        end select
                                     case(0:3)  ! hanging t
                                        AMR_CRS_MAP(iv,ie) = 8
                                     case default
                                        ierr = 1
                                  end select
                               case default
                                  ierr = 1
                            end select
                         case(0:3) ! hanging r
                            select case(is)
                               case(-1) ! non-hanging s
                                  select case(it)
                                     case(-1) ! non-hanging t
!     find the edge direction
!     get vertex number touching hanging face
                                        ivc = ivftove(ir,ivface(1,iv))
!     get edge orientation
                                        ivc = icedg_inv(iv,ivc)
                                        select case(ivc)
                                           case (2) ! edge along s
                                              AMR_CRS_MAP(iv,ie) = 14
                                           case (3) ! edge along t
                                              AMR_CRS_MAP(iv,ie) = 12
                                           case default
                                              ierr = 1
                                        end select
                                     case(0:3)  ! hanging t
                                        AMR_CRS_MAP(iv,ie) = 7
                                     case default
                                        ierr = 1
                                  end select
                               case(0:3) ! hanging s
                                  select case(it)
                                     case(-1) ! non-hanging t
                                        AMR_CRS_MAP(iv,ie) = 6
                                     case(0:3)  ! hanging t
                                        ierr = 1
                                     case default
                                        ierr = 1
                                  end select
                               case default
                                  ierr=1
                            end select
                         case default
                            ierr = 1
                     end select
                  else if (AMR_HNG_VRT(iv,ie).eq.1) then ! face hanging node
                     AMR_CRS_MAP(iv,ie) = 15
                  else if (AMR_HNG_VRT(iv,ie).eq.0) then ! non-hangign node
                     select case (ir)
                         case(-1)! non-hanging r
                            select case(is)
                               case(-1) ! non-hanging s
                                  select case(it)
                                     case(-1) ! non-hanging t
!     look for concave edges
                                       ir = AMR_HNG_EDG
     $                                        (ivetoedg(1,iv),ie)
                                       is = AMR_HNG_EDG
     $                                        (ivetoedg(2,iv),ie)
                                       it = AMR_HNG_EDG
     $                                        (ivetoedg(3,iv),ie)

                                       select case(ir)
                                         case (-1)  ! not hanging
                                           select case (is)
                                             case (-1)  ! not hanging
                                               select case (it)
                                                 case (-1)  ! not hanging
                                                 ! concave vertex or conforming mesh
                                                  AMR_CRS_MAP(iv,ie)=1
                                                 case (0:1)
                                                  AMR_CRS_MAP(iv,ie)=16
                                                 case default
                                                   ierr = 1
                                               end select
                                             case (0:1)
                                               AMR_CRS_MAP(iv,ie)=17
                                             case default
                                               ierr = 1
                                           end select
                                         case (0:1)
                                           AMR_CRS_MAP(iv,ie)=18
                                         case default
                                           ierr = 1
                                       end select

                                     case(0:3)  ! hanging t
                                        AMR_CRS_MAP(iv,ie) = 11
                                     case default
                                        ierr = 1
                                  end select
                               case(0:3) ! hanging s
                                  select case(it)
                                     case(-1) ! non-hanging t
                                        AMR_CRS_MAP(iv,ie) = 10
                                     case(0:3)  ! hanging t
                                        AMR_CRS_MAP(iv,ie) = 5
                                     case default
                                        ierr = 1
                                  end select
                               case default
                                  ierr = 1
                            end select
                         case(0:3) ! hanging r
                            select case(is)
                               case(-1) ! non-hanging s
                                  select case(it)
                                     case(-1) ! non-hanging t
                                        AMR_CRS_MAP(iv,ie) = 9
                                     case(0:3)  ! hanging t
                                        AMR_CRS_MAP(iv,ie) = 4
                                     case default
                                        ierr = 1
                                  end select
                               case(0:3) ! hanging s
                                  select case(it)
                                     case(-1) ! non-hanging t
                                        AMR_CRS_MAP(iv,ie) = 3
                                     case(0:3)  ! hanging t
                                        AMR_CRS_MAP(iv,ie) = 2
                                     case default
                                        ierr = 1
                                  end select
                               case default
                                  ierr=1
                            end select
                         case default
                            ierr = 1
                     end select
                  else
                     ierr = 1
                  endif  ! AMR_HNG_VRT
               enddo  ! vertex loop
            endif
         enddo  ! element loop
      else  ! IF3D
!     loop ove elements
         do ie = 1,nelt
            if (AMR_HNG_ELM(ie).eq.1) then ! any face or edge hanging?
!     loop over all verticies
               do iv = 1, AMR_NVRT
!     check hanging faces; remember neks doesn't use symmetric notation
                  ir = AMR_HNG_FCS(EFACE(ivface(1,iv)),ie)
                  is = AMR_HNG_FCS(EFACE(ivface(2,iv)),ie)
                  if (AMR_HNG_VRT(iv,ie).eq.1) then ! hanging node
                     if(ir.eq.-1.neqv.is.eq.-1) then
                       AMR_CRS_MAP(iv,ie) = 5
                     else
                        ierr = 1
                     endif
                  else ! non-hangign node
                     select case (ir)
                         case(-1)! non-hanging
                            select case(is)
                               case(-1) ! non-hanging
                                  AMR_CRS_MAP(iv,ie) = 1
                               case(0:1) ! hanging
                                  AMR_CRS_MAP(iv,ie) = 3
                               case default
                                  ierr = 1
                            end select
                         case(0:1) ! hanging
                            select case(is)
                               case(-1) ! non-hanging
                                  AMR_CRS_MAP(iv,ie) = 4
                               case(0:1) ! hanging
                                  AMR_CRS_MAP(iv,ie) = 2
                               case default
                                  ierr=1
                            end select
                         case default
                            ierr = 1
                     end select
                  endif  ! AMR_HNG_VRT
               enddo  ! vertex loop
            endif
         enddo  ! element loop
      endif

!     error check
      ierr = iglmax(ierr,1)
      if (ierr.eq.1) then
         call amr_log
     $     (AMR_LP_ERR,'Error:amr_gen_crs_map wrong hanging mark')
         call exitt
      endif
      return
      end
!=======================================================================
!> @brief Get crs base function.
!! @details Get coarse base function to calculate coarse preconditioner.
!! @param[out]   bcrs      crs base function part
!! @param[in]    iel       element number
!! @param[in]    ivr       vertex number
!! @todo Move it to the other file
      subroutine amr_get_crs_basis(bcrs,iel,ivr)
      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'AMR'
      include 'AMR_CRS'

!     argument list
      real bcrs(nx1,ny1,nz1)
      integer iel, ivr

!     locall variables
      integer il, jl, kl  ! loop index
!-----------------------------------------------------------------------
!     copy base function
      do kl=1,nz1
         do jl=1,ny1
            do il=1,nx1
               bcrs(il,jl,kl) =
     $                 AMR_CRS_BASE(il,jl,kl,AMR_CRS_MAP(ivr,iel),ivr)
            enddo
         enddo
      enddo

      return
      end
!=======================================================================
