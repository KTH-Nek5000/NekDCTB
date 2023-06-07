!> @file amr_trns.f
!! @ingroup nekamr
!! @brief Set of global transfer routines to redistribute elements.
!! @author Adam Peplinski
!! @date Jun 05, 2016
!=======================================================================
!> @brief Redistribute tree data between processors
!! @param[inout]  igrp     elelment greoup mark
!! @param[inout]  lvl      elelment level
!! @param[inout]  crvl     external curved surface mark
!! @param[inout]  bcl      boundary condition mark
!! @remarks This routine uses global scratch space SCRUZ
      subroutine amr_tree_transfer(igrpl,lvl,crvl,bcl)
      implicit none

      include 'SIZE'
      include 'PARALLEL'
      include 'INPUT'
      include 'TOPOL'
      include 'AMR'

      ! argument list
      integer igrpl(LELT),lvl(LELT),crvl(AMR_NFCS,LELT)
      integer bcl(AMR_NFCS,LELT)

      ! local variables
      ! integer array size
      integer isw
      parameter (isw=2+2+2*AMR_NFCS)

      ! transfer arrays
      integer vi(isw,LELT)
      integer*8 vl(1)           ! required by crystal rauter
      real vr(1)                ! required by crystal rauter
      common /scruz/ vi

      integer lnelt,itmp
      integer eg,il,jl,kl       ! loop index
      integer key(1)            ! required by crystal rauter; for sorting

      ! temporary storage for face renumbering
      integer bct(AMR_NFCS), crvt(AMR_NFCS)

#ifdef DEBUG
      ! for testing; take into accout that debugging can be run on less than 100 processes
      integer iunit, ierr
      character(len=2) str
#endif
!-----------------------------------------------------------------------
      call amr_log(AMR_LP_PRD,'Redistributing forest data.')

      ! take number of local p4est elements
      lnelt = AMR_NELT

      ! single send
      ! pack integer array
      do eg=1,lnelt
         ! global element number
         itmp = AMR_nelit + eg
         vi(1,eg) = itmp
         ! processor id
         vi(2,eg) = AMR_PN_NID(eg)

         ! element group
         itmp = 3
         vi(itmp,eg) = igrpl(eg)

         ! element level
         itmp = itmp+1
         vi(itmp,eg) = lvl(eg)

         ! B.C.
         do il=1,AMR_NFCS
            itmp = itmp+1
            vi(itmp,eg) = bcl(il,eg)
         enddo

         ! curvature
         do il=1,AMR_NFCS
            itmp = itmp+1
            vi(itmp,eg) = crvl(il,eg)
         enddo
      enddo

      ! transfer arrays
      call fgslib_crystal_tuple_transfer
     $     (cr_h,lnelt,LELT,vi,isw,vl,0,vr,0,2)

      ! test local element number
      if (lnelt.ne.NELT) then
         call amr_abort('Error: tree_transfer; lnelt /= nelt')
      endif

      ! sort elements acording to their global number
      key(1) = 1
      call fgslib_crystal_tuple_sort
     $     (cr_h,lnelt,vi,isw,vl,0,vr,0,key,1)

      ! put data back
      ! integer arrays
      do eg=1,NELT

!#ifdef DEBUG
         ! check of back communication arrays; just for testing
         ! global element number
         if (vi(1,eg).ne.LGLEL(eg)) then
            call amr_abort
     $           ('Error: tree_transfer; back comm. err.; el')
         endif
         ! origin process id
         if (vi(2,eg).ne.AMR_NP_NID(eg)) then
            call amr_abort
     $           ('Error: tree_transfer; back comm. err.; nid')
         endif
!#endif
         ! element group
         itmp = 3
         igrpl(eg) = vi(itmp,eg)

         ! element level
         itmp = itmp+1
         lvl(eg) = vi(itmp,eg)

         ! B.C.
         do il=1,AMR_NFCS
            itmp = itmp+1
            bcl(il,eg) = vi(itmp,eg)
         enddo

         ! curvature
         do il=1,AMR_NFCS
            itmp = itmp+1
            crvl(il,eg) = vi(itmp,eg)
         enddo
      enddo

      ! change ordering of faces
      ! curvature, BC
      do eg=1,NELT           !  SWAP TO PREPROCESSOR NOTATION
            call icopy(crvt,crvl(1,eg),AMR_NFCS)
            call icopy(bct,bcl(1,eg),AMR_NFCS)
            do il=1,AMR_NFCS
               crvl(il,eg) = crvt(EFACE1(il))
               bcl(il,eg) = bct(EFACE1(il))
            enddo
      enddo

#ifdef DEBUG
      ! testing
      write(str,'(i2.2)') nid
      call io_file_freeid(iunit, ierr)
      open(unit=iunit,file='tree_dat.txt'//str)
      write(iunit,*) 'Global info', nid, nelgt, nelgv, nelt, nelv
      do il=1,lnelt
         write(iunit,*) 'element', il, LGLEL(il),vi(1,il)
         write(iunit,*) 'group, level', igrpl(il),lvl(il)
         write(iunit,*) 'curvature, bc'
         do jl=1,AMR_NFCS
            write(iunit,*) jl,crvl(jl,il),bcl(jl,il)
         enddo
      enddo
      close(iunit)
#endif
      return
      end subroutine
!=======================================================================
!> @brief Global communication for error estimator mark
!! @remarks This routine uses global scratch space SCRUZ
      subroutine amr_mark_transfer()
      implicit none

      include 'SIZE'
      include 'PARALLEL'
      include 'INPUT'
      include 'TOPOL'
      include 'AMR'

      ! local variables
      ! integer and real array sizes
      integer isw
      parameter (isw=(1+2))

      ! transfer arrays
      integer vi(isw,LELT)
      integer*8 vl(1)       ! required by crystal rauter
      real vr(1)            ! required by crystal rauter
      common /scruz/ vi

      integer lnelt,eg
      integer key(1)            ! required by crystal rauter; for sorting
!-----------------------------------------------------------------------
      ! redistribute integer refinement mark back to p4est element distribution
      call amr_log(AMR_LP_PRD,'Redistributing MARK array.')

      ! take number of local nek elements
      lnelt = NELT

      ! single send
      ! pack integer array
      do eg=1,lnelt
        ! global element number
        vi(1,eg) = LGLEL(eg)
        ! processor id
        vi(2,eg) = AMR_NP_NID(eg)
        ! mark array
        vi(3,eg) = AMR_MARK(eg)
      enddo

      ! transfer arrays
      call fgslib_crystal_tuple_transfer
     $     (cr_h,lnelt,LELT,vi,isw,vl,0,vr,0,2)

      ! test local element number
      if (lnelt.ne.AMR_NELT) then
         call amr_abort('Error: mark_transfer; lnelt /= nelt')
      endif

      ! sort elements acording to their global number
      key(1) = 1
      call fgslib_crystal_tuple_sort(cr_h,lnelt,vi,isw,vl,0,vr,0,key,1)

      ! put data back
      do eg=1,AMR_NELT

!#ifdef DEBUG
        ! check of back communication arrays; just for testing
        ! global element number
        if (vi(1,eg).ne.(AMR_NELIT+eg)) then
           call amr_abort
     $          ('Error: mark_transfer; back comm. err.; el')
        endif
        ! origin processor id
        if (vi(2,eg).ne.AMR_PN_NID(eg)) then
           call amr_abort
     $          ('Error: mark_transfer; back comm. err.; nid')
        endif
!#endif

         ! mark array
         AMR_MARK(eg) = vi(3,eg)
      enddo

      return
      end subroutine
!=======================================================================
!> @brief Global communication for refinement history; MAP part
!! @remarks This routine uses global scratch space SCRUZ
      subroutine amr_hst_map_transfer()
      implicit none

      include 'SIZE'
      include 'PARALLEL'
      include 'INPUT'
      include 'TOPOL'
      include 'AMR'

      ! local variables
      ! integer and real array sizes
      integer isw
      parameter (isw=(3+1))

      ! transfer arrays
      integer vi(isw,LELT)
      integer*8 vl(1)        ! required by crystal rauter
      real vr(1)             ! required by crystal rauter
      common /scruz/ vi

      integer lnelt
      integer eg,jl             ! loop index
      integer key(1)            ! required by crystal rauter; for sorting

#ifdef DEBUG
      ! for testing
      integer iunit,ierr
      character(len=2) str1, str2
      ! call number
      integer icalled
      save icalled
      data icalled /0/
#endif
!-----------------------------------------------------------------------
      call amr_log(AMR_LP_PRD,'Redistributing history MAP.')

      ! take number of local p4est elements
      lnelt = AMR_NELT

      ! count array entries
      jl = 0

      ! single send
      ! pack integer array
      do eg=1,lnelt
        ! if the element was unchanged
        if (AMR_GLGL_MAP(1,eg).ne.0) then
            ! count array entries
            jl=jl+1

            ! old global element number
            vi(1,jl) = AMR_GLGL_MAP(1,eg)
            ! old processor id
            vi(2,jl) = AMR_GLGL_MAP(3,eg)
            ! old local element number
            vi(3,jl) = AMR_GLGL_MAP(2,eg)

            ! current global element number
            vi(4,jl) = eg + AMR_NELIT
        endif
      enddo

      ! test local element number
      if (jl.ne.AMR_MAP_NR) then
         call amr_abort('Error: map_transfer; wrong jl')
      endif
      ! correct number of array entries
      lnelt = AMR_MAP_NR

      ! reset position array; whole array as those values will be used
      do eg=1,AMR_NVRT*LELT
         AMR_GLGL_MAP(1,eg) = 0
         AMR_GLGL_MAP(2,eg) = 0
         AMR_GLGL_MAP(3,eg) = -1
      enddo

      ! transfer arrays
      call fgslib_crystal_tuple_transfer
     $     (cr_h,lnelt,LELT,vi,isw,vl,0,vr,0,2)

      if(lnelt.gt.0) then
        ! test local element number
        if (lnelt.gt.AMR_NELT_O) then
           call amr_abort('Error: map_transfer;wrong lnelt')
        endif

        ! sort elements acording to their old global number
        key(1) = 1
        call fgslib_crystal_tuple_sort
     $   (cr_h,lnelt,vi,isw,vl,0,vr,0,key,1)

        ! put data back
        do eg=1,lnelt
            ! map array
            AMR_GLGL_MAP(1,vi(3,eg)) = vi(4,eg)
            !AMR_GLGL_MAP(2,vi(3,eg)) = GLLEL(vi(4,eg))
            AMR_GLGL_MAP(3,vi(3,eg)) = GLLNID(vi(4,eg))
        enddo
      endif

      ! save current number of "valid" elements
      AMR_MAP_NR = lnelt

#ifdef DEBUG
      ! testing
      ! count call number
      icalled = icalled+1
      write(str1,'(i2.2)') NID
      write(str2,'(i2.2)') icalled
      call io_file_freeid(iunit, ierr)
      open(unit=iunit,file='MAP_dat.txt'//str1//'i'//str2)
      write(iunit,*) 'Global info', nid, nelgt, nelgv, nelt, nelv
      write(iunit,*) 'Global info_old', AMR_NELGT_O, AMR_NELT_O,
     $                AMR_NELV_O
      do eg=1,AMR_NELT_O
         write(iunit,*) eg, (AMR_GLGL_MAP(jl,eg),jl=1,3)
      enddo
      close(iunit)
#endif

      return
      end subroutine
!=======================================================================
!> @brief Global communication for refinement history; RFN part
!! @remarks This routine uses global scratch space SCRNS, SCRUZ
      subroutine amr_hst_rfn_transfer()
      implicit none

      include 'SIZE'
      include 'PARALLEL'
      include 'INPUT'
      include 'TOPOL'
      include 'AMR'

      ! local variables
      ! integer and real array sizes
      integer isw,rsw
      parameter (isw=(3+2))
      parameter (rsw=1)
      integer lbuf    ! buffer size (element number)
      parameter (lbuf=LELT*AMR_NVRT)

      ! transfer arrays
      real vr(rsw)
      integer vi(isw,lbuf)
      integer*8 vl              ! required by crystal rauter
      common /scrns/ vr
      common /scruz/ vi

      integer lnelt
      integer eg,itmp1,itmp2,il,jl
      integer key(1)            ! required by crystal rauter; for sorting

#ifdef DEBUG
      ! for testing
      integer iunit,ierr
      character(len=2) str
      ! call number
      integer icalled, icalled_w
      save icalled
      data icalled /0/
      parameter (icalled_w=6)
#endif
!-----------------------------------------------------------------------
      call amr_log(AMR_LP_PRD,'Redistributing history RFN.')

      ! take number of local refined p4est elements
      lnelt = AMR_RFN_NR

      ! single send
      ! pack integer array
      do eg=1,lnelt
        ! old global parent element number
        vi(1,eg) = AMR_GLGL_RFN(2,eg)
        ! old parent element processor id
        vi(2,eg) = AMR_GLGL_RFN(4,eg)
        ! old local parent element number
        vi(3,eg) = AMR_GLGL_RFN(3,eg)

        ! current global element number
        vi(4,eg) = AMR_GLGL_RFN(1,eg)
        ! child position
        vi(5,eg) = AMR_GLGL_RFN(5,eg)

      enddo

      ! transfer arrays
      call fgslib_crystal_tuple_transfer
     $     (cr_h,lnelt,lbuf,vi,isw,vl,0,vr,0,2)


      if (lnelt.gt.0) then

        ! check arrays size
        if (lnelt.gt.lbuf) then
           call amr_abort
     $          ('Error: rfn_transfer; too many refined blocks 1')
        endif

        ! test local element number; it must be multiplication of AMR_NVRT
        if (mod(lnelt,AMR_NVRT).ne.0) then
           call amr_abort('Error: rfn_transfer; wrong lnelt')
        endif

        ! sort elements acording to their global number
        key(1) = 1
        call fgslib_crystal_tuple_sort
     $       (cr_h,lnelt,vi,isw,vl,0,vr,0,key,1)

        ! put data back
        ! count new allocated elements
        il = AMR_NELT_O
        do eg=1,lnelt,AMR_NVRT
            itmp1 = il
            do jl=0,AMR_NVRT-1
                ! which child
                itmp2 = vi(5,eg+jl)
                ! new global element number
                AMR_GLGL_RFN(1,eg+itmp2)  = vi(4,eg+jl)
                ! local parent position sanity check
                if (vi(3,eg).ne.vi(3,eg+jl)) call amr_abort
     $                  ('Error: rfn_transfer; wrong parent number')
                AMR_GLGL_RFN(2,eg+itmp2)  = vi(3,eg)
                ! local child position
                if (itmp2.eq.0) then
                    AMR_GLGL_RFN(3,eg+itmp2)  = vi(3,eg)
                else
                    il=il+1
                    AMR_GLGL_RFN(3,eg+itmp2)  = itmp1+itmp2
                endif
            enddo
        enddo

        ! check arrays size
        itmp1=il
        if (il.gt.lbuf) call amr_abort
     $          ('Error: rfn_transfer; too many refined blocks 2')

        ! fill map array checking if everything is correct
        do eg=1,lnelt
            il = AMR_GLGL_RFN(3,eg)
            if (AMR_GLGL_MAP(1,il).eq.0) then
                AMR_GLGL_MAP(1,il) = AMR_GLGL_RFN(1,eg)
                AMR_GLGL_MAP(3,il) = GLLNID(AMR_GLGL_RFN(1,eg))
            else
                call amr_abort
     $               ('Error: rfn_transfer; index allready used')
            endif
        enddo
      else
        itmp1=AMR_NELT_O
      endif
      ! save current number of "valid" elements
      AMR_RFN_NR = lnelt
      AMR_RFN_NR_S = itmp1

#ifdef DEBUG
      ! testing
      ! count call number
      icalled = icalled+1
      if (icalled.eq.icalled_w) then
      write(str,'(i2.2)') nid
      call io_file_freeid(iunit, ierr)
      ! RFN
      open(unit=iunit,file='RFN_dat.txt'//str)
      write(iunit,*) 'Global info', nid, AMR_RFN_NR
      do eg=1,AMR_RFN_NR
         write(iunit,*) eg, (AMR_GLGL_RFN(il,eg),il=1,3)
      enddo
      close(iunit)
      ! modified MAP
      open(unit=iunit,file='MAP2_dat.txt'//str)
      write(iunit,*) 'Global info', nid, nelgt, nelgv, nelt, nelv,
     $                AMR_RFN_NR_S
      write(iunit,*) 'Global info_old', AMR_NELGT_O, AMR_NELT_O,
     $                AMR_NELV_O
      do eg=1,AMR_RFN_NR_S
         write(iunit,*) eg, (AMR_GLGL_MAP(il,eg),il=1,3)
      enddo
      close(iunit)
      endif
#endif

      return
      end subroutine
!=======================================================================
!> @brief Global communication for refinement history; CRS part
!! @remarks This routine uses global scratch space SCRNS, SCRUZ, SCREV
      subroutine amr_hst_crs_transfer
      implicit none

      include 'SIZE'
      include 'PARALLEL'
      include 'INPUT'
      include 'TOPOL'
      include 'AMR'

      ! local variables
      ! integer and real array sizes
      integer isw,rsw
      parameter (isw=(3+2))
      parameter (rsw=1)
      integer lbuf    ! buffer size (element number)
      parameter (lbuf=LELT*AMR_NVRT)

      ! transfer arrays
      real vr(rsw)
      integer vi(isw,lbuf)
      integer*8 vl              ! required by crystal rauter
      common /scrns/ vr
      common /scruz/ vi

      ! tmp array
      integer vi_tmp(lbuf)
      common /screv/ vi_tmp

      integer lnelt

      integer eg,itmp,itmp2,itmp3,il,jl, kl
      integer key(1)            ! required by crystal rauter; for sorting

      ! functions
      integer igl_running_sum

#ifdef DEBUG
      ! for testing
      integer iunit,ierr
      character(len=2) str
      ! call number
      integer icalled, icalled_w
      save icalled
      data icalled /0/
      parameter (icalled_w=6)
#endif
!-----------------------------------------------------------------------
      call amr_log(AMR_LP_PRD,'Redistributing history CRS.')

      ! first get tmp global numbering of all elements to be removed
      ! count number of elements on lower number processors
      itmp  = igl_running_sum(AMR_CRS_NR)

      if (AMR_CRS_NR.gt.0) then
        ! calculate tmp global element number
        ! starting point
        itmp = NELGT+(itmp - AMR_CRS_NR)*(AMR_NVRT-1)
        do il = 1, AMR_CRS_NR
            do jl=2,AMR_NVRT
                itmp = itmp + 1
                AMR_GLGL_CRS(1,jl,il) = itmp
            enddo
        enddo

        ! make local copy of gllnid; this can contain communication and
        ! will be used twice
        do il = 1, AMR_CRS_NR
           vi_tmp(il) = GLLNID(AMR_GLGL_CRS(1,1,il))
        enddo

        ! fill transfer array,
        ! this one will be send to processors owning the elements for coarsening
        ! count array entries
        kl=0
        do il = 1, AMR_CRS_NR
            do jl=1,AMR_NVRT
                kl=kl+1
                ! old global child element number
                vi(1,kl) = AMR_GLGL_CRS(2,jl,il)
                ! old child element processor id
                vi(2,kl) = AMR_GLGL_CRS(4,jl,il)
                ! old local child element number
                vi(3,kl) = AMR_GLGL_CRS(3,jl,il)

                ! new element global number
                vi(4,kl) = AMR_GLGL_CRS(1,jl,il)
                ! new parent element processor id
                vi(5,kl) = vi_tmp(il)
            enddo
        enddo
        lnelt = kl
      else
        lnelt = 0
      endif

      ! transfer arrays
      call fgslib_crystal_tuple_transfer
     $     (cr_h,lnelt,lbuf,vi,isw,vl,0,vr,0,2)

      ! local number of elements to be send for coarsening
      AMR_CRS_NR_S = lnelt

      if (lnelt.gt.0) then
        ! sort elements acording to their global number
        key(1) = 1
        call fgslib_crystal_tuple_sort
     $       (cr_h,lnelt,vi,isw,vl,0,vr,0,key,1)

        ! fill MAP arrays
        do il=1,lnelt
          itmp = vi(3,il)
          if (AMR_GLGL_MAP(1,itmp).eq.0) then
            AMR_GLGL_MAP(1,itmp) = vi(4,il)
            AMR_GLGL_MAP(3,itmp) = vi(5,il)
          else
            call amr_abort('Error: crs_transfer; index allready used')
          endif
        enddo
      endif

      ! this one will be send to processors owning the elements resulting from coarsening
      if (AMR_CRS_NR.gt.0) then
        ! fill transfer array
        ! count array entries
        kl=0
        do il = 1, AMR_CRS_NR
            itmp = AMR_GLGL_CRS(1,1,il)
            do jl=1,AMR_NVRT
                kl=kl+1
                ! new parent element global number
                vi(1,kl) = itmp
                ! new parent element processor id
                vi(2,kl) = vi_tmp(il)
                ! new child element global number
                vi(3,kl) = AMR_GLGL_CRS(1,jl,il)
                ! child position
                vi(4,kl) = jl
            enddo
        enddo
        lnelt = kl
      else
        lnelt = 0
      endif

      ! transfer arrays
      call fgslib_crystal_tuple_transfer
     $     (cr_h,lnelt,lbuf,vi,isw,vl,0,vr,0,2)

      if (lnelt.gt.0) then

        ! check arrays size
        if (lnelt.gt.lbuf) call amr_abort
     $          ('Error: crs_transfer; too many coarsening blocks 1')

        ! test local element number; it must be multiplication of AMR_NVRT
        if (mod(lnelt,AMR_NVRT).ne.0) then
           call amr_abort('Error: crs_transfer; wrong lnelt')
        endif

        ! sort elements acording to their global number
        key(1) = 1
        call fgslib_crystal_tuple_sort
     $   (cr_h,lnelt,vi,isw,vl,0,vr,0,key,1)

        ! put data back
        ! count new allocated elements
        il = NELT
        do eg=1,lnelt/AMR_NVRT
            kl=(eg-1)*AMR_NVRT
            itmp = vi(1,kl+1)
            itmp2 = il
            do jl=1,AMR_NVRT
                ! sanity check; global parent position
                if (itmp.ne.vi(1,kl+jl)) call amr_abort
     $                  ('Error: crs_transfer; wrong parent number')

                ! which child
                itmp3 = vi(4,kl+jl)
                ! new global element number
                AMR_GLGL_CRS(1,itmp3,eg)  = vi(3,kl+jl)

                ! local child position
                if (itmp3.eq.1) then
                    AMR_GLGL_CRS(2,itmp3,eg)  = GLLEL(itmp)
                else
                    il=il+1
                    AMR_GLGL_CRS(2,itmp3,eg)  = itmp2 + itmp3 - 1
                endif
            enddo
        enddo

        ! check arrays size
        itmp=il
        if (itmp.gt.lbuf) then
           call amr_abort
     $          ('Error: crs_transfer; too many coarsening blocks 2')
        endif

        lnelt = lnelt/AMR_NVRT

      endif
      ! save current number of "valid" elements
      AMR_CRS_NR = lnelt

      ! check if there are still some zeros in MAP
      do eg=1,AMR_RFN_NR_S
        if (AMR_GLGL_MAP(1,eg).eq.0) then
           call amr_abort('Error: crs_transfer; inconsistent MAP')
        endif
      enddo

#ifdef DEBUG
      ! testing
      ! count call number
      icalled = icalled+1
      if (icalled.eq.icalled_w) then
      write(str,'(i2.2)') nid
      call io_file_freeid(iunit, ierr)
      ! CRS
      open(unit=iunit,file='CRS_dat.txt'//str)
      write(iunit,*) 'Global info', nid, AMR_CRS_NR
      do eg=1,AMR_CRS_NR
        do il=1,AMR_NVRT
            write(iunit,*) eg,il,(AMR_GLGL_CRS(jl,il,eg),jl=1,2)
        enddo
      enddo
      close(iunit)
      ! modified MAP
      open(unit=iunit,file='MAP3_dat.txt'//str)
      write(iunit,*) 'Global info', nid, nelgt, nelgv, nelt, nelv,
     $                AMR_RFN_NR_S
      write(iunit,*) 'Global info_old', AMR_NELGT_O, AMR_NELT_O,
     $                AMR_NELV_O
      do eg=1,AMR_RFN_NR_S
         write(iunit,*) eg, (AMR_GLGL_MAP(jl,eg),jl=1,3)
      enddo
      close(iunit)
      endif
#endif

      return
      end subroutine
!=======================================================================
!> @brief Redistribute single variable after refinement
!! @param[inout] vr      redistributed vector
!! @param[inout] vi      transfer mark array (work array)
!! @param[in]    rsw     array size
!! @param[inout] lbuff   array size
      subroutine amr_vec_transfer(vr, vi, rsw, lbuff)
      implicit none

      include 'SIZE'
      include 'PARALLEL'
      include 'INPUT'
      include 'TOPOL'
      include 'AMR'

      ! argument list
      integer rsw, lbuff
      real vr(rsw,lbuff)
      integer vi(2,lbuff)

      ! local variables
      ! integer array size
      integer isw
      parameter (isw=2)
      integer*8 vl              ! required by crystal rauter
      integer lnelt,eg
      integer key(1)            ! required by crystal rauter; for sorting
!-----------------------------------------------------------------------
      ! take number of local p4est elements
      lnelt = AMR_RFN_NR_S

      ! single send
      ! pack integer array
      do eg=1,lnelt
         ! global element number
         vi(1,eg) = AMR_GLGL_MAP(1,eg)
         ! processor id
         vi(2,eg) = AMR_GLGL_MAP(3,eg)
      enddo

      ! min aray size
      eg = lbuff
      ! transfer arrays
      call fgslib_crystal_tuple_transfer
     $     (cr_h,lnelt,eg,vi,isw,vl,0,vr,rsw,2)

      ! test local element number
      if (lnelt.ne.(NELT+AMR_CRS_NR*(AMR_NVRT-1))) then
         call amr_abort('Error: v1_transfer; lnelt /= nelt')
      endif

      ! sort elements acording to their global number
      key(1) = 1
      call fgslib_crystal_tuple_sort
     $     (cr_h,lnelt,vi,isw,vl,0,vr,rsw,key,1)

      return
      end subroutine
!=======================================================================
!> @brief Redistribute node data.
!! @details Redistribute global node numberring and hanging node mark.
!! @param[inout] lnodes  redistributed vector
!! @param[in]    vnode   number of nodes
!! @param[inout] lnelt   local number of elements
!! @remarks This routine uses global scratch space SCRUZ
      subroutine amr_node_transfer(lnodes,vnode,lnelt)
      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'
      include 'TOPOL'
      include 'AMR'
      include 'AMR_TOPOL'

      ! argument list
      integer vnode ! number of nodes; counted vertices, faces, edges
      integer lnelt ! local number of elements; p4est count
      integer*8 lnodes(vnode*LELT) ! global numberring of nodes

      ! local variables
      integer el, fl ! loop index
      integer itmp, itmpv(AMR_NFCS) ! dummy variables

      ! transfer arrays
      integer vi_size
      parameter (vi_size = 2+1+AMR_NVRT+AMR_NFCS+AMR_NEDG)
      real vr                   ! required by c.r.
      integer vi(vi_size*LELT)        ! processor and element information
      common /scruz/ vi
      integer key(1)            ! required by crystal rauter; for sorting
!-----------------------------------------------------------------------
      ! data transfer
      ! single send
      ! pack processor, element and hang arrays info
      do el=1,lnelt
        ! global element number
        itmp = AMR_NELIT + el
        vi(1+(el-1)*vi_size) = itmp
        ! processor id
        vi(2+(el-1)*vi_size) = AMR_PN_NID(el)
        ! hang arrays
        itmp = 3
        ! element
        vi(itmp+(el-1)*vi_size) = AMR_HNG_ELM(el)
        itmp = itmp+1
        ! vertex
        do fl=1,AMR_NVRT
            vi(itmp+(el-1)*vi_size) = AMR_HNG_VRT(fl,el)
            itmp = itmp+1
        enddo
        ! face
        do fl=1,AMR_NFCS
            vi(itmp+(el-1)*vi_size) = AMR_HNG_FCS(fl,el)
            itmp = itmp+1
        enddo
        ! edge
#if N_DIM == 3
        do fl = 1, AMR_NEDG
            vi(itmp+(el-1)*vi_size) = AMR_HNG_EDG(fl,el)
            itmp = itmp+1
        enddo
#endif
      enddo

      ! transfer arrays
      call fgslib_crystal_tuple_transfer
     $     (cr_h,lnelt,LELT,vi,vi_size,lnodes,vnode,vr,0,2)

      ! test local element number
      if (lnelt.ne.NELT) call amr_abort
     $     ('Error: amr_topol_transfer; lnelt /= nelt')

      ! sort elements acording to their global number
      key(1) = 1
      call fgslib_crystal_tuple_sort
     $     (cr_h,lnelt,vi,vi_size,lnodes,vnode,vr,0,key,1)

      ! move data back
      ! hang arrays
      do el=1,lnelt
        itmp = 3
        ! element
        AMR_HNG_ELM(el) = vi(itmp+(el-1)*vi_size)
        itmp = itmp+1
        ! vertex
        do fl = 1, AMR_NVRT
            AMR_HNG_VRT(fl,el) = vi(itmp+(el-1)*vi_size)
            itmp = itmp+1
        enddo
        ! face
        do fl = 1, AMR_NFCS
            AMR_HNG_FCS(fl,el) = vi(itmp+(el-1)*vi_size)
            itmp = itmp+1
        enddo
        ! edge
#if N_DIM == 3
        do fl = 1, AMR_NEDG
           AMR_HNG_EDG(fl,el) = vi(itmp+(el-1)*vi_size)
           itmp = itmp+1
        enddo
#endif
      enddo

      ! change ordering of faces; due to inconsistent ordering of bc and cbc arrays
      ! AMR_HNG_FCS
      do el=1,NELT           !  swap to preprocessor notation
        call icopy(itmpv,AMR_HNG_FCS(1,el),AMR_NFCS)
        do fl=1,AMR_NFCS
            AMR_HNG_FCS(fl,el) = itmpv(EFACE1(fl))
        enddo
      enddo

      return
      end subroutine
!=======================================================================
!> @brief Redistribute aligment data.
!! @details Redistribute orientation of faces and edges.
!! @param[inout] lnelt   local number of elements
!! @remarks This routine uses global scratch space SCRUZ
      subroutine amr_algn_transfer(lnelt)
      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'
      include 'TOPOL'
      include 'AMR'
      include 'AMR_TOPOL'

      ! argument list
      integer lnelt ! local number of elements; p4est count


      ! local variables
      integer el, fl ! loop index
      integer itmp, itmpv(AMR_NFCS) ! dummy variables

      ! transfer arrays
      integer*8 vl              ! required by crystal rauter
      integer vi_size
      parameter (vi_size = 2+AMR_NFCS+AMR_NEDG)
      real vr                   ! required by c.r.
      integer vi(vi_size*LELT)        ! processor and element information
      common /scruz/ vi
      integer key(1)            ! required by crystal rauter; for sorting
!-----------------------------------------------------------------------
      ! data transfer
      ! single send
      ! pack processor, element and hang arrays info
      do el=1,lnelt
        ! global element number
        itmp = AMR_NELIT + el
        vi(1+(el-1)*vi_size) = itmp
        ! processor id
        vi(2+(el-1)*vi_size) = AMR_PN_NID(el)
        ! alignment arrays
        itmp = 3
        ! face
        do fl=1,AMR_NFCS
            vi(itmp+(el-1)*vi_size) = AMR_ALGN_FCS(fl,el)
            itmp = itmp+1
        enddo
        ! edge
#if N_DIM == 3
        do fl = 1, AMR_NEDG
            vi(itmp+(el-1)*vi_size) = AMR_ALGN_EDG(fl,el)
            itmp = itmp+1
        enddo
#endif
      enddo

      ! transfer arrays
      call fgslib_crystal_tuple_transfer
     $     (cr_h,lnelt,LELT,vi,vi_size,vl,0,vr,0,2)

      ! test local element number
      if (lnelt.ne.NELT) call amr_abort
     $     ('Error: amr_algn_transfer; lnelt /= nelt')

      ! sort elements acording to their global number
      key(1) = 1
      call fgslib_crystal_tuple_sort
     $     (cr_h,lnelt,vi,vi_size,vl,0,vr,0,key,1)

      ! move data back
      ! alignment arrays
      do el=1,lnelt
        itmp = 3
        ! face
        do fl = 1, AMR_NFCS
            AMR_ALGN_FCS(fl,el) = vi(itmp+(el-1)*vi_size)
            itmp = itmp+1
        enddo
        ! edge
#if N_DIM == 3
        do fl = 1, AMR_NEDG
           AMR_ALGN_EDG(fl,el) = vi(itmp+(el-1)*vi_size)
           itmp = itmp+1
        enddo
#endif
      enddo

      return
      end subroutine
!=======================================================================
!> @brief Redistribute family data.
!! @param[inout] lnelt   local number of elements
!! @remarks This routine uses global scratch space SCRUZ
      subroutine amr_fml_transfer(lnelt)
      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'
      include 'TOPOL'
      include 'AMR'
      include 'AMR_TOPOL'

      ! argument list
      integer lnelt ! local number of elements; p4est count

      ! local variables
      integer el, fl ! loop index
      integer itmp, itmpv(AMR_NFCS) ! dummy variables

      ! transfer arrays
      integer*8 vl              ! required by crystal rauter
      integer vi_size
      parameter (vi_size = 2+2)
      real vr                   ! required by c.r.
      integer vi(vi_size*LELT)        ! processor and element information
      common /scruz/ vi
      integer key(1)            ! required by crystal rauter; for sorting
!-----------------------------------------------------------------------
      ! data transfer
      ! single send
      ! pack processor, element and hang arrays info
      do el=1,lnelt
        ! global element number
        itmp = AMR_NELIT + el
        vi(1+(el-1)*vi_size) = itmp
        ! processor id
        vi(2+(el-1)*vi_size) = AMR_PN_NID(el)
        ! family arrays
        itmp = 3
        do fl=1,2
            vi(itmp+(el-1)*vi_size) = AMR_FML_MARK(fl,el)
            itmp = itmp+1
        enddo
      enddo

      ! transfer arrays
      call fgslib_crystal_tuple_transfer
     $     (cr_h,lnelt,LELT,vi,vi_size,vl,0,vr,0,2)

      ! test local element number
      if (lnelt.ne.NELT) call amr_abort
     $     ('Error: amr_fml_transfer; lnelt /= nelt')

      ! sort elements acording to their global number
      key(1) = 1
      call fgslib_crystal_tuple_sort
     $     (cr_h,lnelt,vi,vi_size,vl,0,vr,0,key,1)

      ! move data back
      ! family arrays
      do el=1,lnelt
        itmp = 3
        do fl = 1, 2
            AMR_FML_MARK(fl,el) = vi(itmp+(el-1)*vi_size)
            itmp = itmp+1
        enddo
      enddo

      return
      end subroutine
!=======================================================================
!> @brief Redistribute graph data.
!! @param[inout] lnelt   local number of elements
!! @param[out]   gnelt   global element number list
!! @param[inout] grph    graph
!! @param[inout] goff    graph offset
!! @param[inout] ewgt    edge weight
!! @param[inout] elprt           element partitioning
!! @param[in]    lgnel   element list size
!! @param[in]    lgrph   graph array size
!! @param[in]    lgoff   graph offset array size
!! @remarks This routine uses global scratch space SCRUZ
      subroutine amr_graph_transfer(lnelt,gnelt,grph,goff,ewgt,
     $           elprt,lgnel,lgrph,lgoff)
      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'
      include 'TOPOL'
      include 'AMR'

      ! argument list
      integer lnelt, lgnel, lgrph, lgoff
      integer gnelt(lgnel), grph(lgrph), goff(lgoff), ewgt(lgrph),
     $        elprt(lgoff)

      ! local variables
      integer il, jl ! loop index
      integer itmp1, itmp2 ! dummy variables

      ! transfer arrays
      integer*8 vl              ! required by crystal rauter
      integer vi_s, vi_b
      parameter (vi_s = 2+2, vi_b = LELT*LX1*LY1*LZ1)
      real vr                   ! required by c.r.
      integer vi(vi_s*vi_b) ! processor and element information
      common /scruz/ vi
      integer key(1)            ! required by crystal rauter; for sorting
!-----------------------------------------------------------------------
      ! pack processor, element and graph arrays info
      itmp2 = 0
      do il=1,lnelt
         itmp1 = AMR_nelit + il
         do jl=goff(il)+1,goff(il+1)
            ! global element number
            vi(itmp2*vi_s + 1) = itmp1
            ! processor id
            vi(itmp2*vi_s + 2) = elprt(il)
            ! graph edge
            vi(itmp2*vi_s + 3) = grph(jl)
            ! edge weight
            vi(itmp2*vi_s + 4) = ewgt(jl)
            itmp2 = itmp2 + 1
        enddo
      enddo

      ! transfer arrays
      call fgslib_crystal_tuple_transfer
     $     (cr_h,itmp2,vi_b,vi,vi_s,vl,0,vr,0,2)

      ! check array sizes
      if (itmp2.gt.vi_b) then
         call amr_abort
     $        ('Error: amr_graph_transfer; buffer too small')
      endif

      ! sort elements acording to their global number
      key(1) = 1
      call fgslib_crystal_tuple_sort
     $     (cr_h,itmp2,vi,vi_s,vl,0,vr,0,key,1)

      ! move data back
      jl = 1
      gnelt(jl) = vi(1)
      goff(jl) = 0
      do il=1,itmp2
        grph(il) = vi((il-1)*vi_s + 3)
        ewgt(il) = vi((il-1)*vi_s + 4)
        if(gnelt(jl).ne.vi((il-1)*vi_s + 1)) then
           jl = jl + 1
           gnelt(jl) = vi((il-1)*vi_s + 1)
           goff(jl)  = il-1
        endif
      enddo
      goff(jl+1) = itmp2
      lnelt = jl

      return
      end subroutine
!=======================================================================
!> @brief Redistribute graph node's weights and coordinates.
!! @param[inout] lnelt   local number of elements
!! @param[inout] nwgt    node weight
!! @param[inout] ncrd    node coordinates
!! @param[in]    elprt   element partitioning
!! @param[in]    lgoff   array size
!! @remarks This routine uses global scratch space SCRUZ, SCRNS
      subroutine amr_grnd_transfer(lnelt,nwgt,ncrd,elprt,lgoff)
      implicit none

      include 'SIZE'
      include 'INPUT'
      include 'PARALLEL'
      include 'TOPOL'
      include 'AMR'
!      include 'AMR_TOPOL'

      ! argument list
      integer lnelt, lgoff
      integer nwgt(lgoff), elprt(lgoff)
      real ncrd(LDIM*lgoff)

      ! local variables
      integer il, jl ! loop index

      ! transfer arrays
      integer*8 vl              ! required by crystal rauter
      integer vi_s, vi_b
      parameter (vi_s = 2+1, vi_b = AMR_NVRT*LELT+1)
      real vr(LDIM*vi_b)    ! node coordinates
      integer vi(vi_s*vi_b) ! processor and element information
      common /scrns/ vr
      common /scruz/ vi
      integer key(1)            ! required by crystal rauter; for sorting
!-----------------------------------------------------------------------
      ! pack data, element and graph arrays info
      do il=1,lnelt
            ! global element number
            vi((il-1)*vi_s + 1) = AMR_nelit + il
            ! processor id
            vi((il-1)*vi_s + 2) = elprt(il)
            ! node weight
            vi((il-1)*vi_s + 3) = nwgt(il)
            ! coordinates
            do jl = 1, LDIM
               vr((il-1)*LDIM + jl) = ncrd((il-1)*LDIM + jl)
            enddo
      enddo

      ! transfer arrays
      call fgslib_crystal_tuple_transfer
     $     (cr_h,lnelt,vi_b,vi,vi_s,vl,0,vr,LDIM,2)

      ! check array sizes
      if (lnelt.gt.vi_b) then
         call amr_abort
     $        ('Error: amr_grnd_transfer; buffer too small')
      endif


      ! sort elements acording to their global number
      key(1) = 1
      call fgslib_crystal_tuple_sort
     $     (cr_h,lnelt,vi,vi_s,vl,0,vr,LDIM,key,1)

      ! unpack data
      do il=1,lnelt
        nwgt(il)  = vi((il-1)*vi_s + 3)
        do jl = 1, LDIM
           ncrd((il-1)*LDIM + jl) = vr((il-1)*LDIM + jl)
        enddo
      enddo

      return
      end subroutine
!=======================================================================
!> @brief Transfer single integer between p4est and nek5000 distributions
!! @param[inout] lnelt   local number of elements
!! @param[inout] gnelt   global element number
!! @param[inout] elmrk   element mark (partitioning/local number)
!! @param[inout] eldst   element destination
!! @param[in]    lgoff   array size
!! @remarks This routine uses global scratch space SCRUZ
      subroutine amr_int_transfer(lnelt,gnelt,elmrk,eldst,lgoff)
      implicit none

      include 'SIZE'
      include 'PARALLEL'
      include 'AMR'

      ! argument list
      integer lnelt, lgoff
      integer gnelt(lgoff), elmrk(lgoff), eldst(lgoff)

      ! local variables
      integer il, jl ! loop index
      integer itmp   ! dummy variable

      ! transfer arrays
      integer*8 vl              ! required by crystal rauter
      integer vi_s, vi_b
      parameter (vi_s = 2+1, vi_b = LELT)
      real vr(1)    ! node coordinates
      integer vi(vi_s*vi_b) ! processor and element information
      common /scruz/ vi
      integer key(1)            ! required by crystal rauter; for sorting
!-----------------------------------------------------------------------
      ! pack data, element and graph arrays info
      do il=1,lnelt
            ! global element number
            vi((il-1)*vi_s + 1) = gnelt(il)
            ! destination process id
            vi((il-1)*vi_s + 2) = eldst(il)
            ! element mark
            vi((il-1)*vi_s + 3) = elmrk(il)
      enddo

      ! transfer arrays
      call fgslib_crystal_tuple_transfer
     $     (cr_h,lnelt,vi_b,vi,vi_s,vl,0,vr,0,2)

      ! check array sizes
      if (lnelt.gt.vi_b) then
         call amr_abort
     $    ('Error: amr_int_transfer; buffer too small')
      endif

      ! test local element number
      if (lnelt.ne.AMR_NELT) call amr_abort
     $     ('Error: amr_int_transfer; lnelt /= AMR_NELT')

      ! sort elements acording to their global number
      key(1) = 1
      call fgslib_crystal_tuple_sort
     $     (cr_h,lnelt,vi,vi_s,vl,0,vr,0,key,1)

      ! unpack data
      do il=1,lnelt
        gnelt(il) = vi((il-1)*vi_s + 1)
        eldst(il) = vi((il-1)*vi_s + 2)
        elmrk(il) = vi((il-1)*vi_s + 3)
      enddo

      return
      end subroutine
!=======================================================================
