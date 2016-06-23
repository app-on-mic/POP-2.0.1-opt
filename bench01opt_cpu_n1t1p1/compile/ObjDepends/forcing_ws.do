forcing_ws.o ObjDepends/forcing_ws.do: forcing_ws.f90
forcing_ws.o: kinds_mod.o
forcing_ws.o: constants.o
forcing_ws.o: blocks.o
forcing_ws.o: distribution.o
forcing_ws.o: domain.o
forcing_ws.o: boundary.o
forcing_ws.o: io.o
forcing_ws.o: forcing_tools.o
forcing_ws.o: time_management.o
forcing_ws.o: grid.o
forcing_ws.o: exit_mod.o
