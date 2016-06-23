solvers.o ObjDepends/solvers.do: solvers.f90
solvers.o: kinds_mod.o
solvers.o: blocks.o
solvers.o: distribution.o
solvers.o: domain.o
solvers.o: constants.o
solvers.o: boundary.o
solvers.o: global_reductions.o
solvers.o: gather_scatter.o
solvers.o: broadcast.o
solvers.o: grid.o
solvers.o: io.o
solvers.o: time_management.o
solvers.o: exit_mod.o
