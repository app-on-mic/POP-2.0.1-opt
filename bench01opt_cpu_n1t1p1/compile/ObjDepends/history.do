history.o ObjDepends/history.do: history.f90
history.o: kinds_mod.o
history.o: domain.o
history.o: constants.o
history.o: prognostic.o
history.o: grid.o
history.o: io.o
history.o: broadcast.o
history.o: time_management.o
history.o: forcing.o
history.o: forcing_shf.o
history.o: exit_mod.o
