barotropic.o ObjDepends/barotropic.do: barotropic.f90
barotropic.o: kinds_mod.o
barotropic.o: blocks.o
barotropic.o: domain.o
barotropic.o: constants.o
barotropic.o: prognostic.o
barotropic.o: boundary.o
barotropic.o: solvers.o
barotropic.o: operators.o
barotropic.o: grid.o
barotropic.o: time_management.o
barotropic.o: global_reductions.o
barotropic.o: forcing.o
barotropic.o: forcing_ap.o
barotropic.o: tavg.o
