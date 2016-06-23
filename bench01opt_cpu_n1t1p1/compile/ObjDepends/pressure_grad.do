pressure_grad.o ObjDepends/pressure_grad.do: pressure_grad.f90
pressure_grad.o: kinds_mod.o
pressure_grad.o: blocks.o
pressure_grad.o: constants.o
pressure_grad.o: operators.o
pressure_grad.o: grid.o
pressure_grad.o: broadcast.o
pressure_grad.o: communicate.o
pressure_grad.o: io_types.o
pressure_grad.o: state_mod.o
pressure_grad.o: time_management.o
pressure_grad.o: exit_mod.o
