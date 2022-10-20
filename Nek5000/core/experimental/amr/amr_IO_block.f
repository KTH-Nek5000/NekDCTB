!> @file amr_IO_block.f
!! @ingroup nekamr
!! @brief block data to initialise common block for I/O routines for nekamr
!! @details Following Nek5000 standard I keep block data in seaprate file
!! @author Adam Peplinski
!! @date Mar 7, 2016
!! @note this file is identical with Nek_framework/io/io_tools/io_tools_block.f
!=======================================================================
      block data io_common_init
!     keeep track of max iunit generated
      integer io_iunit_min, io_iunit_max
      common /io_iunit/ io_iunit_min, io_iunit_max
      data io_iunit_min/200/
      data io_iunit_max/200/
      end
