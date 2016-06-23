vmix_rich.o ObjDepends/vmix_rich.do: vmix_rich.f90
vmix_rich.o: kinds_mod.o
vmix_rich.o: blocks.o
vmix_rich.o: distribution.o
vmix_rich.o: domain.o
vmix_rich.o: constants.o
vmix_rich.o: grid.o
vmix_rich.o: broadcast.o
vmix_rich.o: io.o
vmix_rich.o: state_mod.o
vmix_rich.o: time_management.o
vmix_rich.o: exit_mod.o
