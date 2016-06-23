#!/bin/bash
[ ! -f ${NETCDF_CPUINC}/netcdf.h ] && echo "netcdf.h(CPU) not found" && exit -1
[ ! -f ${NETCDF_CPULIB}/libnetcdf.a ] && echo "libnetcdf.a(CPU) not found" && exit -1
rm -rf pop
pushd compile
DEBUG=
#DEBUG="-traceback -check all -debug all"

#mpiicc  -DPOSIX -I${NETCDF_CPUINC} -mcmodel=medium -ipo -O3 -xMIC-AVX512 -fp-model precise -ftz -g -c printtrace.c
mpiicc  -DPOSIX -I${NETCDF_CPUINC} -mcmodel=medium -ipo -O3 -xMIC-AVX512 -fp-model precise -ftz -g -c fix_64.c
mpiifort -mcmodel=medium -ipo -I${NETCDF_CPUINC} -IObjDepends -O3 -xMIC-AVX512 -fp-model fast=2 -assume byterecl -ftz -qopt-report=5 -qopt-report-phase=loop,vec,openmp -ipo  -free -g ${DEBUG} -c kinds_mod.f90
mpiifort -mcmodel=medium -ipo -I${NETCDF_CPUINC} -IObjDepends -O3 -xMIC-AVX512 -fp-model fast=2 -assume byterecl -ftz -qopt-report=5 -qopt-report-phase=loop,vec,openmp -ipo  -free -g ${DEBUG} -c domain_size.f90
mpiifort -mcmodel=medium -ipo -I${NETCDF_CPUINC} -IObjDepends -O3 -xMIC-AVX512 -fp-model fast=2 -assume byterecl -ftz -qopt-report=5 -qopt-report-phase=loop,vec,openmp -ipo  -free -g ${DEBUG} -c constants.f90
mpiifort -mcmodel=medium -ipo -I${NETCDF_CPUINC} -IObjDepends -O3 -xMIC-AVX512 -fp-model fast=2 -assume byterecl -ftz -qopt-report=5 -qopt-report-phase=loop,vec,openmp -ipo  -free -g ${DEBUG} -c communicate.f90
mpiifort -mcmodel=medium -ipo -I${NETCDF_CPUINC} -IObjDepends -O3 -xMIC-AVX512 -fp-model fast=2 -assume byterecl -ftz -qopt-report=5 -qopt-report-phase=loop,vec,openmp -ipo  -free -g ${DEBUG} -c exit_mod.f90
mpiifort -mcmodel=medium -ipo -I${NETCDF_CPUINC} -IObjDepends -O3 -xMIC-AVX512 -fp-model fast=2 -assume byterecl -ftz -qopt-report=5 -qopt-report-phase=loop,vec,openmp -ipo  -free -g ${DEBUG} -c blocks.f90
mpiifort -mcmodel=medium -ipo -I${NETCDF_CPUINC} -IObjDepends -O3 -xMIC-AVX512 -fp-model fast=2 -assume byterecl -ftz -qopt-report=5 -qopt-report-phase=loop,vec,openmp -ipo  -free -g ${DEBUG} -c distribution.f90
mpiifort -mcmodel=medium -ipo -I${NETCDF_CPUINC} -IObjDepends -O3 -xMIC-AVX512 -fp-model fast=2 -assume byterecl -ftz -qopt-report=5 -qopt-report-phase=loop,vec,openmp -ipo  -free -g ${DEBUG} -c broadcast.f90
mpiifort -mcmodel=medium -ipo -I${NETCDF_CPUINC} -IObjDepends -O3 -xMIC-AVX512 -fp-model fast=2 -assume byterecl -ftz -qopt-report=5 -qopt-report-phase=loop,vec,openmp -ipo  -free -g ${DEBUG} -c io_types.f90
mpiifort -cpp -mcmodel=medium -ipo -I${NETCDF_CPUINC} -IObjDepends -O3 -xMIC-AVX512 -fp-model fast=2 -assume byterecl -ftz -qopt-report=5 -qopt-report-phase=loop,vec,openmp -ipo  -free -g ${DEBUG} -c boundary.f90
mpiifort -mcmodel=medium -ipo -I${NETCDF_CPUINC} -IObjDepends -O3 -xMIC-AVX512 -fp-model fast=2 -assume byterecl -ftz -qopt-report=5 -qopt-report-phase=loop,vec,openmp -ipo  -free -g ${DEBUG} -c domain.f90
mpiifort -mcmodel=medium -ipo -I${NETCDF_CPUINC} -IObjDepends -O3 -xMIC-AVX512 -fp-model fast=2 -assume byterecl -ftz -qopt-report=5 -qopt-report-phase=loop,vec,openmp -ipo  -free -g ${DEBUG} -c gather_scatter.f90
mpiifort -mcmodel=medium -ipo -I${NETCDF_CPUINC} -IObjDepends -O3 -xMIC-AVX512 -fp-model fast=2 -assume byterecl -ftz -qopt-report=5 -qopt-report-phase=loop,vec,openmp -ipo  -free -g ${DEBUG} -c io_netcdf.f90
mpiifort -mcmodel=medium -ipo -I${NETCDF_CPUINC} -IObjDepends -O3 -xMIC-AVX512 -fp-model fast=2 -assume byterecl -ftz -qopt-report=5 -qopt-report-phase=loop,vec,openmp -ipo  -free -g ${DEBUG} -c io_binary.f90
mpiifort -mcmodel=medium -ipo -I${NETCDF_CPUINC} -IObjDepends -O3 -xMIC-AVX512 -fp-model fast=2 -assume byterecl -ftz -qopt-report=5 -qopt-report-phase=loop,vec,openmp -ipo  -free -g ${DEBUG} -c io.f90
mpiifort -mcmodel=medium -ipo -I${NETCDF_CPUINC} -IObjDepends -O3 -xMIC-AVX512 -fp-model fast=2 -assume byterecl -ftz -qopt-report=5 -qopt-report-phase=loop,vec,openmp -ipo  -free -g ${DEBUG} -c global_reductions.f90
mpiifort -mcmodel=medium -ipo -I${NETCDF_CPUINC} -IObjDepends -O3 -xMIC-AVX512 -fp-model fast=2 -assume byterecl -ftz -qopt-report=5 -qopt-report-phase=loop,vec,openmp -ipo  -free -g ${DEBUG} -c grid.f90
mpiifort -mcmodel=medium -ipo -I${NETCDF_CPUINC} -IObjDepends -O3 -xMIC-AVX512 -fp-model fast=2 -assume byterecl -ftz -qopt-report=5 -qopt-report-phase=loop,vec,openmp -ipo  -free -g ${DEBUG} -c prognostic.f90
mpiifort -mcmodel=medium -ipo -I${NETCDF_CPUINC} -IObjDepends -O3 -xMIC-AVX512 -fp-model fast=2 -assume byterecl -ftz -qopt-report=5 -qopt-report-phase=loop,vec,openmp -ipo  -free -g ${DEBUG} -c time_management.f90
mpiifort -mcmodel=medium -ipo -I${NETCDF_CPUINC} -IObjDepends -O3 -xMIC-AVX512 -fp-model fast=2 -assume byterecl -ftz -qopt-report=5 -qopt-report-phase=loop,vec,openmp -ipo  -free -g ${DEBUG} -c solvers.f90
mpiifort -mcmodel=medium -ipo -I${NETCDF_CPUINC} -IObjDepends -O3 -xMIC-AVX512 -fp-model fast=2 -assume byterecl -ftz -qopt-report=5 -qopt-report-phase=loop,vec,openmp -ipo  -free -g ${DEBUG} -c forcing_tools.f90
mpiifort -mcmodel=medium -ipo -I${NETCDF_CPUINC} -IObjDepends -O3 -xMIC-AVX512 -fp-model fast=2 -assume byterecl -ftz -qopt-report=5 -qopt-report-phase=loop,vec,openmp -ipo  -free -g ${DEBUG} -c forcing_ws.f90
mpiifort -mcmodel=medium -ipo -I${NETCDF_CPUINC} -IObjDepends -O3 -xMIC-AVX512 -fp-model fast=2 -assume byterecl -ftz -qopt-report=5 -qopt-report-phase=loop,vec,openmp -ipo  -free -g ${DEBUG} -c forcing_shf.f90
mpiifort -mcmodel=medium -ipo -I${NETCDF_CPUINC} -IObjDepends -O3 -xMIC-AVX512 -fp-model fast=2 -assume byterecl -ftz -qopt-report=5 -qopt-report-phase=loop,vec,openmp -ipo  -free -g ${DEBUG} -c forcing_sfwf.f90
mpiifort -mcmodel=medium -ipo -I${NETCDF_CPUINC} -IObjDepends -O3 -xMIC-AVX512 -fp-model fast=2 -assume byterecl -ftz -qopt-report=5 -qopt-report-phase=loop,vec,openmp -ipo  -free -g ${DEBUG} -c forcing_pt_interior.f90
mpiifort -mcmodel=medium -ipo -I${NETCDF_CPUINC} -IObjDepends -O3 -xMIC-AVX512 -fp-model fast=2 -assume byterecl -ftz -qopt-report=5 -qopt-report-phase=loop,vec,openmp -ipo  -free -g ${DEBUG} -c forcing_s_interior.f90
mpiifort -mcmodel=medium -ipo -I${NETCDF_CPUINC} -IObjDepends -O3 -xMIC-AVX512 -fp-model fast=2 -assume byterecl -ftz -qopt-report=5 -qopt-report-phase=loop,vec,openmp -ipo  -free -g ${DEBUG} -c forcing_ap.f90
mpiifort -mcmodel=medium -ipo -I${NETCDF_CPUINC} -IObjDepends -O3 -xMIC-AVX512 -fp-model fast=2 -assume byterecl -ftz -qopt-report=5 -qopt-report-phase=loop,vec,openmp -ipo  -free -g ${DEBUG} -c ice.f90
mpiifort -mcmodel=medium -ipo -I${NETCDF_CPUINC} -IObjDepends -O3 -xMIC-AVX512 -fp-model fast=2 -assume byterecl -ftz -qopt-report=5 -qopt-report-phase=loop,vec,openmp -ipo  -free -g ${DEBUG} -c timers.f90
mpiifort -mcmodel=medium -ipo -I${NETCDF_CPUINC} -IObjDepends -O3 -xMIC-AVX512 -fp-model fast=2 -assume byterecl -ftz -qopt-report=5 -qopt-report-phase=loop,vec,openmp -ipo  -free -g ${DEBUG} -c forcing_coupled.f90
mpiifort -mcmodel=medium -ipo -I${NETCDF_CPUINC} -IObjDepends -O3 -xMIC-AVX512 -fp-model fast=2 -assume byterecl -ftz -qopt-report=5 -qopt-report-phase=loop,vec,openmp -ipo  -free -g ${DEBUG} -c tavg.f90
mpiifort -mcmodel=medium -ipo -I${NETCDF_CPUINC} -IObjDepends -O3 -xMIC-AVX512 -fp-model fast=2 -assume byterecl -ftz -qopt-report=5 -qopt-report-phase=loop,vec,openmp -ipo  -free -g ${DEBUG} -c forcing.f90
mpiifort -mcmodel=medium -ipo -I${NETCDF_CPUINC} -IObjDepends -O3 -xMIC-AVX512 -fp-model fast=2 -assume byterecl -ftz -qopt-report=5 -qopt-report-phase=loop,vec,openmp -ipo  -free -g ${DEBUG} -c diagnostics.f90
mpiifort -mcmodel=medium -ipo -I${NETCDF_CPUINC} -IObjDepends -O3 -xMIC-AVX512 -fp-model fast=2 -assume byterecl -ftz -qopt-report=5 -qopt-report-phase=loop,vec,openmp -ipo  -free -g ${DEBUG} -c state_mod.f90
mpiifort -mcmodel=medium -ipo -I${NETCDF_CPUINC} -IObjDepends -O3 -xMIC-AVX512 -fp-model fast=2 -assume byterecl -ftz -qopt-report=5 -qopt-report-phase=loop,vec,openmp -ipo  -free -g ${DEBUG} -c operators.f90
mpiifort -mcmodel=medium -ipo -I${NETCDF_CPUINC} -IObjDepends -O3 -xMIC-AVX512 -fp-model fast=2 -assume byterecl -ftz -qopt-report=5 -qopt-report-phase=loop,vec,openmp -ipo  -free -g ${DEBUG} -c vmix_const.f90
mpiifort -mcmodel=medium -ipo -I${NETCDF_CPUINC} -IObjDepends -O3 -xMIC-AVX512 -fp-model fast=2 -assume byterecl -ftz -qopt-report=5 -qopt-report-phase=loop,vec,openmp -ipo  -free -g ${DEBUG} -c vmix_rich.f90
mpiifort -mcmodel=medium -ipo -I${NETCDF_CPUINC} -IObjDepends -O3 -xMIC-AVX512 -fp-model fast=2 -assume byterecl -ftz -qopt-report=5 -qopt-report-phase=loop,vec,openmp -ipo  -free -g ${DEBUG} -c vmix_kpp.f90
mpiifort -mcmodel=medium -ipo -I${NETCDF_CPUINC} -IObjDepends -O3 -xMIC-AVX512 -fp-model fast=2 -assume byterecl -ftz -qopt-report=5 -qopt-report-phase=loop,vec,openmp -ipo  -free -g ${DEBUG} -c vertical_mix.f90
mpiifort -mcmodel=medium -ipo -I${NETCDF_CPUINC} -IObjDepends -O3 -xMIC-AVX512 -fp-model fast=2 -assume byterecl -ftz -qopt-report=5 -qopt-report-phase=loop,vec,openmp -ipo  -free -g ${DEBUG} -c hmix_gm.f90
mpiifort -mcmodel=medium -ipo -I${NETCDF_CPUINC} -IObjDepends -O3 -xMIC-AVX512 -fp-model fast=2 -assume byterecl -ftz -qopt-report=5 -qopt-report-phase=loop,vec,openmp -ipo  -free -g ${DEBUG} -c advection.f90
mpiifort -mcmodel=medium -ipo -I${NETCDF_CPUINC} -IObjDepends -O3 -xMIC-AVX512 -fp-model fast=2 -assume byterecl -ftz -qopt-report=5 -qopt-report-phase=loop,vec,openmp -ipo  -free -g ${DEBUG} -c pressure_grad.f90
mpiifort -mcmodel=medium -ipo -I${NETCDF_CPUINC} -IObjDepends -O3 -xMIC-AVX512 -fp-model fast=2 -assume byterecl -ftz -qopt-report=5 -qopt-report-phase=loop,vec,openmp -ipo  -free -g ${DEBUG} -c topostress.f90
mpiifort -mcmodel=medium -ipo -I${NETCDF_CPUINC} -IObjDepends -O3 -xMIC-AVX512 -fp-model fast=2 -assume byterecl -ftz -qopt-report=5 -qopt-report-phase=loop,vec,openmp -ipo  -free -g ${DEBUG} -c hmix_del2.f90
mpiifort -mcmodel=medium -ipo -I${NETCDF_CPUINC} -IObjDepends -O3 -xMIC-AVX512 -fp-model fast=2 -assume byterecl -ftz -qopt-report=5 -qopt-report-phase=loop,vec,openmp -ipo  -free -g ${DEBUG} -c hmix_del4.f90
mpiifort -mcmodel=medium -ipo -I${NETCDF_CPUINC} -IObjDepends -O3 -xMIC-AVX512 -fp-model fast=2 -assume byterecl -ftz -qopt-report=5 -qopt-report-phase=loop,vec,openmp -ipo  -free -g ${DEBUG} -c hmix_aniso.f90
mpiifort -mcmodel=medium -ipo -I${NETCDF_CPUINC} -IObjDepends -O3 -xMIC-AVX512 -fp-model fast=2 -assume byterecl -ftz -qopt-report=5 -qopt-report-phase=loop,vec,openmp -ipo  -free -g ${DEBUG} -c horizontal_mix.f90
mpiifort -mcmodel=medium -ipo -I${NETCDF_CPUINC} -IObjDepends -O3 -xMIC-AVX512 -fp-model fast=2 -assume byterecl -ftz -qopt-report=5 -qopt-report-phase=loop,vec,openmp -ipo  -free -g ${DEBUG} -c baroclinic.f90
mpiifort -mcmodel=medium -ipo -I${NETCDF_CPUINC} -IObjDepends -O3 -xMIC-AVX512 -fp-model fast=2 -assume byterecl -ftz -qopt-report=5 -qopt-report-phase=loop,vec,openmp -ipo  -free -g ${DEBUG} -c barotropic.f90
mpiifort -mcmodel=medium -ipo -I${NETCDF_CPUINC} -IObjDepends -O3 -xMIC-AVX512 -fp-model fast=2 -assume byterecl -ftz -qopt-report=5 -qopt-report-phase=loop,vec,openmp -ipo  -free -g ${DEBUG} -c current_meters.f90
mpiifort -mcmodel=medium -ipo -I${NETCDF_CPUINC} -IObjDepends -O3 -xMIC-AVX512 -fp-model fast=2 -assume byterecl -ftz -qopt-report=5 -qopt-report-phase=loop,vec,openmp -ipo  -free -g ${DEBUG} -c drifters.f90
mpiifort -mcmodel=medium -ipo -I${NETCDF_CPUINC} -IObjDepends -O3 -xMIC-AVX512 -fp-model fast=2 -assume byterecl -ftz -qopt-report=5 -qopt-report-phase=loop,vec,openmp -ipo  -free -g ${DEBUG} -c history.f90
mpiifort -mcmodel=medium -ipo -I${NETCDF_CPUINC} -IObjDepends -O3 -xMIC-AVX512 -fp-model fast=2 -assume byterecl -ftz -qopt-report=5 -qopt-report-phase=loop,vec,openmp -ipo  -free -g ${DEBUG} -c hydro_sections.f90
mpiifort -mcmodel=medium -ipo -I${NETCDF_CPUINC} -IObjDepends -O3 -xMIC-AVX512 -fp-model fast=2 -assume byterecl -ftz -qopt-report=5 -qopt-report-phase=loop,vec,openmp -ipo  -free -g ${DEBUG} -c surface_hgt.f90
mpiifort -mcmodel=medium -ipo -I${NETCDF_CPUINC} -IObjDepends -O3 -xMIC-AVX512 -fp-model fast=2 -assume byterecl -ftz -qopt-report=5 -qopt-report-phase=loop,vec,openmp -ipo  -free -g ${DEBUG} -c xdisplay.f90
mpiifort -mcmodel=medium -ipo -I${NETCDF_CPUINC} -IObjDepends -O3 -xMIC-AVX512 -fp-model fast=2 -assume byterecl -ftz -qopt-report=5 -qopt-report-phase=loop,vec,openmp -ipo  -free -g ${DEBUG} -c restart.f90
mpiifort -mcmodel=medium -ipo -I${NETCDF_CPUINC} -IObjDepends -O3 -xMIC-AVX512 -fp-model fast=2 -assume byterecl -ftz -qopt-report=5 -qopt-report-phase=loop,vec,openmp -ipo  -free -g ${DEBUG} -c movie.f90
mpiifort -mcmodel=medium -ipo -I${NETCDF_CPUINC} -IObjDepends -O3 -xMIC-AVX512 -fp-model fast=2 -assume byterecl -ftz -qopt-report=5 -qopt-report-phase=loop,vec,openmp -ipo  -free -g ${DEBUG} -c output.f90
mpiifort -mcmodel=medium -ipo -I${NETCDF_CPUINC} -IObjDepends -O3 -xMIC-AVX512 -fp-model fast=2 -assume byterecl -ftz -qopt-report=5 -qopt-report-phase=loop,vec,openmp -ipo  -free -g ${DEBUG} -c initial.f90
mpiifort -mcmodel=medium -ipo -I${NETCDF_CPUINC} -IObjDepends -O3 -xMIC-AVX512 -fp-model fast=2 -assume byterecl -ftz -qopt-report=5 -qopt-report-phase=loop,vec,openmp -ipo  -free -g ${DEBUG} -c step_mod.f90
mpiifort -mcmodel=medium -ipo -I${NETCDF_CPUINC} -IObjDepends -O3 -xMIC-AVX512 -fp-model fast=2 -assume byterecl -ftz -qopt-report=5 -qopt-report-phase=loop,vec,openmp -ipo  -free -g ${DEBUG} -c POP.f90

popd
# pop_9-off
mpiifort -o pop -mcmodel=medium -ipo compile/fix_64.o compile/advection.o compile/baroclinic.o compile/barotropic.o compile/blocks.o compile/constants.o compile/current_meters.o compile/diagnostics.o compile/distribution.o compile/domain.o compile/drifters.o compile/exit_mod.o compile/forcing_ap.o compile/forcing_coupled.o compile/forcing.o compile/forcing_pt_interior.o compile/forcing_sfwf.o compile/forcing_shf.o compile/forcing_s_interior.o compile/forcing_tools.o compile/forcing_ws.o compile/grid.o compile/history.o compile/hmix_aniso.o compile/hmix_del2.o compile/hmix_del4.o compile/hmix_gm.o compile/horizontal_mix.o compile/hydro_sections.o compile/ice.o compile/initial.o compile/io_binary.o compile/io.o compile/io_netcdf.o compile/io_types.o compile/kinds_mod.o compile/movie.o compile/operators.o compile/output.o compile/POP.o compile/pressure_grad.o compile/prognostic.o compile/restart.o compile/solvers.o compile/state_mod.o compile/step_mod.o compile/surface_hgt.o compile/tavg.o compile/time_management.o compile/timers.o compile/topostress.o compile/vertical_mix.o compile/vmix_const.o compile/vmix_kpp.o compile/vmix_rich.o compile/xdisplay.o compile/domain_size.o compile/boundary.o compile/broadcast.o compile/communicate.o compile/gather_scatter.o compile/global_reductions.o  -L${NETCDF_CPULIB} -lnetcdf   

