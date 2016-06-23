domain.o ObjDepends/domain.do: domain.f90
domain.o: kinds_mod.o
domain.o: constants.o
domain.o: communicate.o
domain.o: broadcast.o
domain.o: blocks.o
domain.o: distribution.o
domain.o: exit_mod.o
domain.o: io_types.o
domain.o: boundary.o
domain.o: domain_size.o
