$(OBJDIR)/io_stream.o        :${STREAM_SRC}/io_stream.f;                           $(F77) -c $(FL2) -I./ $< -o $@
