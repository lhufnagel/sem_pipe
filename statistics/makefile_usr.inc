$(OBJDIR)/uparam.o               :src-usr/uparam.f;                            $(F77) -c $(FL2) $< -o $@
$(OBJDIR)/IO_tools.o             :src-usr/IO_tools.f;                          $(F77) -c $(FL2) $< -o $@
$(OBJDIR)/avg_read.o             :src-usr/avg_read.f;                          $(F77) -c $(FL2) $< -o $@
