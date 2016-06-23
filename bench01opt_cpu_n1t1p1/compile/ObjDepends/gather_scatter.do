gather_scatter.o ObjDepends/gather_scatter.do: gather_scatter.f90
gather_scatter.o: kinds_mod.o
gather_scatter.o: communicate.o
gather_scatter.o: constants.o
gather_scatter.o: blocks.o
gather_scatter.o: distribution.o
gather_scatter.o: domain.o
gather_scatter.o: domain_size.o
gather_scatter.o: exit_mod.o
