$(OBJDIR)/frame.o           :${COMPRESS_SRC}/framework_files/frame.f;              $(F77) -c $(FL2) -I./ $< -o $@
$(OBJDIR)/mntrlog_block.o   :${COMPRESS_SRC}/framework_files/mntrlog_block.f;      $(F77) -c $(FL2) -I./ $< -o $@
$(OBJDIR)/mntrlog.o         :${COMPRESS_SRC}/framework_files/mntrlog.f;            $(F77) -c $(FL2) -I./ $< -o $@
$(OBJDIR)/mntrtmr_block.o   :${COMPRESS_SRC}/framework_files/mntrtmr_block.f;      $(F77) -c $(FL2) -I./ $< -o $@
$(OBJDIR)/mntrtmr.o         :${COMPRESS_SRC}/framework_files/mntrtmr.f;            $(F77) -c $(FL2) -I./ $< -o $@
$(OBJDIR)/rprm_block.o      :${COMPRESS_SRC}/framework_files/rprm_block.f;         $(F77) -c $(FL2) -I./ $< -o $@
$(OBJDIR)/rprm.o            :${COMPRESS_SRC}/framework_files/rprm.f;               $(F77) -c $(FL2) -I./ $< -o $@
$(OBJDIR)/io_tools_block.o  :${COMPRESS_SRC}/framework_files/io_tools_block.f;     $(F77) -c $(FL2) -I./ $< -o $@
$(OBJDIR)/io_tools.o        :${COMPRESS_SRC}/framework_files/io_tools.f;           $(F77) -c $(FL2) -I./ $< -o $@
$(OBJDIR)/sperri_trunc.o        :${COMPRESS_SRC}/framework_files/sperri_trunc.f;           $(F77) -c $(FL2) -I./ $< -o $@
$(OBJDIR)/io_trunc.o        :${COMPRESS_SRC}/io_trunc.f;                           $(F77) -c $(FL2) -I./ $< -o $@
$(OBJDIR)/test_io_trunc.o        :${COMPRESS_SRC}/test_io_trunc.f;                           $(F77) -c $(FL2) -I./ $< -o $@

