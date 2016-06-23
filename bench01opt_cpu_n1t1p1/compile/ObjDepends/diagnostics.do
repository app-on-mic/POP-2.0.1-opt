diagnostics.o ObjDepends/diagnostics.do: diagnostics.f90
diagnostics.o: domain.o
diagnostics.o: constants.o
diagnostics.o: prognostic.o
diagnostics.o: time_management.o
diagnostics.o: io.o
diagnostics.o: broadcast.o
diagnostics.o: global_reductions.o
diagnostics.o: grid.o
diagnostics.o: solvers.o
diagnostics.o: forcing.o
diagnostics.o: timers.o
diagnostics.o: exit_mod.o
