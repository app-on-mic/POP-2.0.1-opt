io_binary.o ObjDepends/io_binary.do: io_binary.f90
io_binary.o: kinds_mod.o
io_binary.o: domain_size.o
io_binary.o: domain.o
io_binary.o: constants.o
io_binary.o: boundary.o
io_binary.o: communicate.o
io_binary.o: broadcast.o
io_binary.o: gather_scatter.o
io_binary.o: exit_mod.o
io_binary.o: io_types.o
