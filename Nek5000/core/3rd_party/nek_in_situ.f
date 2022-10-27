c-----------------------------------------------------------------------
      subroutine in_situ_init()
#ifdef VISIT
      call visit_init()
#elif CATALYST
      call catalyst_init()
#elif ADIOS2
      call adios2_init()
#elif ASYNC_CATALYST
      call async_catalyst_init()
#endif
      end
c-----------------------------------------------------------------------
      subroutine in_situ_check()
#ifdef VISIT
      call visit_check()
#elif CATALYST
      call catalyst_process()
#elif ADIOS2
      call adios2_write()
#elif ASYNC_CATALYST
      call async_catalyst_write()
#endif
      end
c-----------------------------------------------------------------------
      subroutine in_situ_end()
#ifdef VISIT
      call visit_end()
#elif CATALYST
      call catalyst_end()
#elif ADIOS2
      call adios2_end()
#elif ASYNC_CATALYST
      call async_catalyst_end()
#endif
      end
c-----------------------------------------------------------------------

