Recipe: Building and Running POP for Intel® Xeon Phi™ Processors
Author: Lin, Junmin (Intel)

I. Overview
This article provides a recipe for how to obtain, compile, and run an
optimized version of POP (Parallel Ocean Program) 2.0.1 with a  bench01 (0.1
degree high resolution) workload on Intel® Xeon® processors and Intel® Xeon
Phi™ processors.
The source for this version of POP 2.0.1 as well as the bench01 workload can
be obtained by contacting Prof. Zhenya Song at songroy@fio.org.cn.

II. Introduction
POP (Parallel Ocean Program) is an ocean circulation model that solves the
three-dimensional primitive equations for ocean dynamics. It consists of the
baroclinic and barotropic solvers that respectively solve 3D equations
explicitly and 2D surface pressure implicitly. 
POP was developed by LANL USA and is widely used in ocean and climate
research. The official web site of POP is located at
http://oceans11.lanl.gov/trac/POP. It is also incorporated into FIO-ESM (First
Institute of Oceanography-Earth System Model) as the ocean component.
This version of POP is optimized for the performance on both Intel Xeon
processors and Intel Xeon Phi processors. Optimizations in this package
include:
•	The addition of two improved block distribution algorithms for better
load balancing.
•	Vectorization in KPP (K-Profile Parameterization) vertical mixing
scheme. 
•	Use of CA-PCG (Communication-Avoiding Preconditioned Conjugate
Gradient) in the barotropic solver to reduce global communication.
•	Fine-grained optimization to improve data locality. 
•	Compiler options tuning. 

III. Preliminaries
1.	To build this package, install Intel® MPI Library 5.0 or higher and
Intel® Parallel Studio XE 2016 or higher on your host system (2016 update1 has
a compile-time error; you’ll need to upgrade to the 2016 update2 version).
2.	Please contact Prof. Zhenya Song at songroy@fio.org.cn to get the
optimized POP source and bench01 workload packages. Please specify that you’d
like the version used for the Intel Recipes.
3.	Set up the Intel MPI Library and Intel® Fortran Compiler environments:
> source
> /opt/intel/compilers_and_libraries_<version>/linux/mpi/bin64/mpivars.sh
> source
> /opt/intel/compilers_and_libraries_<version>/linux/bin/compilervars.sh
> intel64
4.	Set up the other environmental variables.
export I_MPI_FABRICS=shm:tcp
export I_MPI_PIN=enable
ulimit –s unlimited
export OMP_STACKSIZE=1000M
5.	In order to run POP on an Intel Xeon Phi processor, reboot the system
with the quadrant cluster mode and cache memory mode via BIOS settings. Please
refer to Intel® Xeon Phi™ x200 Processor - Memory Modes and Cluster Modes:
Configuration and Use Cases for more details on memory configuration.

IV. Build POP for the Intel® Xeon® processor
1.	Unpack the source code to any directory of /home/<user> 
> tar xvfj POP_opt_20160512.tar.bz2
This creates the POP_build directory.
2.	Build the executables for the Intel Xeon processor with specified rank
size.
> cd /home/<user>/POP_build
> ./batchbuild_hsw.sh -m cpu -a hsw_opt -t hsw_opt_bench01 -p "28,36"
This builds executables for the Intel Xeon processor that run with 28 and 36
MPI processes. The executables are located at the path of
/home/<user>/POP_build/bench01opt_cpu_n1t1p1, with the names of
pop.hsw_opt_bench01.p28 and pop.hsw_opt_bench01.p36, respectively.

V. Build POP for the Intel Xeon Phi processor
1.	Build the executables for the Intel Xeon processor with specified rank
size.
> cd /home/<user>/POP_build
> ./batchbuild_knl.sh -m cpu -a knl_opt -t knl_opt_bench01 -p "62,108"
This builds the executables for the Intel Xeon processor that run with 62 and
108 MPI processes. The executables are located at the path of
/home/<user>/POP_build/bench01opt_cpu_n1t1p1, with the names of
pop.knl_opt_bench01.p62 and pop.knl_opt_bench01.p108, respectively.

VI. Run POP with the bench01 workload on the Intel Xeon processor
1.	Copy the bench01 workload package to the same directory where the
source package is located, and unpack it.
	> cd /home/<user>
	> tar xvfj POP_bench01_bin.tar.bz2
This unpacks the bench01 input data set to the path of
/home/<user>/POP_build/bench01
2.	Check the soft links of topography.0.1 and grid.0.1 files in the
/home/<user>/POP_build/bench01opt_cpu_n1t1p1. Make sure they are valid links
to ../bench01/topography.0.1 and ../bench01/grid.0.1
If not, correct them.
	> cd /home/<user>/POP_build/bench01opt_cpu_n1t1p1
	> ln –s ../bench01/topography.0.1 .
	> ln –s ../bench01/grid.0.1 .
3.	Run POP with the bench01 workload on the Intel Xeon processor.
	> cd /home/<user>/POP_build
	> ./batchrun_hsw.sh -t hsw_opt_bench01 -p "28,36"
This runs POP with 28 MPI ranks and 36 MPI ranks in turn. 
4.	Check the performance. In addition to the screen output, there are
also .log files generated in the POP_build/bench01opt_cpu_n1t1p1 directory. A
successful run of POP shows the timers as follows.
Timer number 10 =      316.27 seconds
  Timer stats (node): min =        312.53 seconds
                      max =        316.27 seconds
                      mean=        315.26 seconds
…
---------------------------------------------------------------------
POP exiting...
 Successful completion of POP run

---------------------------------------------------------------------
The timer 10, 12, and 13 are the total time, baroclinic time, and barotropic
time, respectively.

VII. Run POP with bench01 workload on the Intel Xeon Phi processor
Steps 1 and 2 are the same as steps 1 and 2 in section VI. If you have already
done them in section VI, you can skip them here.
1.	Run POP with bench01 workload on the Intel Xeon Phi processor.
	> cd /home/<user>/POP_build
	> ./batchrun_knl.sh -t knl_opt_bench01 -p "62,108"
This will run POP with 62 MPI ranks and 108 MPI ranks in turn. 
2.	Check the performance, as described in section VI, step 4.

VIII. Performance gain
For the bench01 workload, the following graph shows the speedup achieved from
the Intel Xeon Phi processor, compared to the Intel Xeon processor. As you can
see, we get:
•	Up to 1.29x faster with the Intel® Xeon Phi™ processor 7210 compared
to the two-socket Intel® Xeon® processor E5-2697 v4. 
•	Up to 1.41x faster with the Intel® Xeon Phi™ processor 7250 compared
to the two-socket Intel Xeon processor E5-2697 v4. 
 
Comments on performance improvement on Intel Xeon Phi:
•	POP with the bench01 workload is memory bandwidth bound, and therefore
POP benefits from MCDRAM greatly.
•	POP has good parallel scalability with MPI ranks, and benefits from
more cores. However, since POP is memory bandwidth bound and is sensitive to
the average L2 cache size per rank, 2 ranks per core outperform 4 ranks per
core. In addition, block-wise data decomposition determines that more than 108
ranks could gain no better load balance at the cost of increased memory
traffic. And OpenMP threading could also increase the demand on memory
bandwidth as much as MPI parallelization. Therefore, the best performance on
Intel® Xeon Phi™ 7250 is achieved with 108 MPI ranks even though there are as
many as 272 threads available.
•	POP is well vectorized, since most of the loops are well structured.
However, due to its memory bound nature, the performance improvement from the
added register size available with AVX512 is limited, which is also observed
with AVX2 on Intel Xeon platform.

Testing platform configuration:
Intel Xeon processor E5-2697 v4: Dual-socket Intel Xeon processor E5-2697 v4,
2.3 GHz, 18 cores/socket, 36 cores, 72 threads (HT and Turbo ON), DDR4 128 GB,
2400 MHz, Oracle Linux* Server release 6.7.
Intel Xeon Phi processor 7210 (64 cores): Intel Xeon Phi processor 7210, 64
cores, 256 threads, 1300 MHz core freq. (HT and Turbo ON), 1600 MHz uncore
freq., MCDRAM 16 GB 6.4 GT/s, BIOS 10D28, DDR4 96 GB, 2133 MHz, Red Hat 7.2,
quad cluster mode, MCDRAM cache memory mode.
Intel Xeon Phi processor 7250 (68 cores): Intel Xeon Phi processor 7250, 68
cores, 272 threads, 1400 MHz core freq. (HT and Turbo ON), 1700 MHz uncore
freq., MCDRAM 16 GB 7.2 GT/s, BIOS 10D28, DDR4 96 GB, 2400 MHz, Red Hat 7.2,
quad cluster mode, MCDRAM cache memory mode.

