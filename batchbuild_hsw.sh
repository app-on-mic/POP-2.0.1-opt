#!/bin/bash

# Argument = [-b|--build 0|1=default] [-m|--mode cpu|symm|both=default] [-a|--arch hsw=default|knl] [-t|--batch "avx|hsw|knl|default=<null>binary ext"] [-p|--proclist "num1,num2,...,numN|default=28,36"] [-r|--run 0|1=default]
usage()
{
  echo $0 '[-b|--build 0|1=default] [-m|--mode cpu|symm|both=default] [-a|--arch hsw=default|knl] [-t|--batch "avx|hsw|knl|default=<null>binary ext"] [-p|--proclist "num1,num2,...,numN|default=28,36"] [-r|--run 0|1=default]'
}
ARG_BUILD=1
ARG_MODE=cpu
ARG_ARCH=hsw
ARG_BATCH=unknown_ext
ARG_PROCLIST=28,36
ARG_RUN=1

while [[ $# -ge 1 ]]
do
key="$1"

case $key in
    -b|--build)
    ARG_BUILD="$2"
    shift # past argument
    ;;
    -m|--mode)
    ARG_MODE="$2"
    shift # past argument
    ;;
    -a|--arch)
    ARG_ARCH="$2"
    shift # past argument
    ;;
    -t|--batch)
    ARG_BATCH="$2"
    shift # past argument
    ;;
    -p|--proclist)
    ARG_PROCLIST="$2"
    shift # past argument
    ;;
    -r|--run)
    ARG_RUN="$2"
    shift # past argument
    ;;
    *)
    usage && exit -1
            # unknown option
    ;;
esac
shift # past argument or value
done
# turn off running if in batch build mode
[[ "$ARG_BATCH" != "unknown_ext" ]] && ((ARG_RUN=0))
case $ARG_ARCH in
hsw_opt ) buildcpucmd=mk_hsw_opt.sh;;
* ) buildcpucmd=mk_hsw_opt.sh;;
esac

#echo "ARG_BUILD=" $ARG_BUILD
#echo "ARG_MODE=" $ARG_MODE
#echo "ARG_ARCH=" $ARG_ARCH
#echo "ARG_BATCH=" $ARG_BATCH
#echo "ARG_PROCLIST=" $ARG_PROCLIST
#echo "ARG_RUN=" $ARG_RUN
#exit 0

PROCLIST=`echo $ARG_PROCLIST | sed -e "s/[,;]/ /g"`
NETCDF_CPUPATH=`pwd`/local
NETCDF_MICPATH=`pwd`/local_knl
SRCRD_SFX=_cpu_n0t1p1
############################################
# check netcdf path
read -p "please input the path to NETCDF(CPU)[${NETCDF_CPUPATH}]" -t 30 a
[ -z $a ] || NETCDF_CPUPATH=$a
echo -n "checking netcdf.h..."
find ${NETCDF_CPUPATH} -maxdepth 2 -name netcdf.h -print
[ $? -ne 0 ] && echo "please check the NETCDF(CPU) PATH" && exit -1;
echo -n "checking libnetcdf.a..."
find ${NETCDF_CPUPATH} -maxdepth 2 -name libnetcdf.a -print
[ $? -ne 0 ] && echo "please check the NETCDF(CPU) PATH" && exit -1; 
read -p "please input the path to NETCDF(MIC)[${NETCDF_MICPATH}]" -t 30 a
[ -z $a ] || NETCDF_MICPATH=$a
echo -n "checking netcdf.h..."
find ${NETCDF_MICPATH} -maxdepth 2 -name netcdf.h -print
[ $? -ne 0 ] && echo "please check the NETCDF(MIC) PATH" && exit -1;
echo -n "checking libnetcdf.a..."
find ${NETCDF_MICPATH} -maxdepth 2 -name libnetcdf.a -print
[ $? -ne 0 ] && echo "please check the NETCDF(MIC) PATH" && exit -1; 

echo "netcdf check successfullly."
export NETCDF_CPUINC=${NETCDF_CPUPATH}/include
export NETCDF_CPULIB=${NETCDF_CPUPATH}/lib
export NETCDF_MICINC=${NETCDF_MICPATH}/include
export NETCDF_MICLIB=${NETCDF_MICPATH}/lib

############################################
# CPU_N: CPU MPI rank# per node
# MIC_N: MIC MPI rank# per node
# SYM_N: MPI rank# per node(CPU+MIC,calculated)
# S: source dir of test case(e.g.,bench01d), block_sizes are defined in domain_size.F90
# P: node number
############################################

# blocksize=180x100
function base_distrb()
{
  local result
  case $1 in
  28 ) result=20;;
  36 ) result=16;;
  56 ) result=9;;
  64 ) result=9;;
  66 ) result=9;;
  68 ) result=9;;
  84 ) result=6;;
  112 ) result=6;;
  140 ) result=4;;
  168 ) result=4;;
  196 ) result=4;;
  224 ) result=4;;
  252 ) result=4;;
  280 ) result=2;;
  * ) result=9;;
  esac
  echo $result
}

# blocksize=180x100(default)
function opt_distrb()
{
  local result
  ((result=(420-1)/$1 + 1))
  echo $result
}

# decide block#/rank by workload
# $1: workload; $2: rank#
function choose_distrb()
{
  local result
  case $1 in
  "01base" ) result=`base_distrb $2`;;
  "01opt" ) result=`opt_distrb $2`;;
  * ) result=9;;
  esac
  echo $result
}

for CPU_N in ${PROCLIST}; do
for S in 01opt ; do # 01base 01dov 02dov 
for P in 1; do
((NP=CPU_N*P))
CASE=bench${S}_cpu_n1t1p1

if [[ "$ARG_MODE" == "cpu" || "$ARG_MODE" == "both" ]]; then
#########################################
# build & run pure CPU
#########################################
if [ "$ARG_BUILD" == "1" ]; then
pushd ${CASE}
SYM_N=$CPU_N
((MIC_N=SYM_N-CPU_N))
#########################################
# configure constants in source code
# 1. configure block distribution parameters
sed -i '/nblock_mic_pp = .* ,&! /{
    s/= [0-9]*/= '"1"'/g
}
/nproc_cpu_pn = .* ,&! /{
    s/= [0-9]*/= '"$CPU_N"'/g
}
/nproc_mic_pn = .* ,&! /{
    s/= [0-9]*/= '"$MIC_N"'/g
}
/nproc_pn = .* ,&! /{
    s/= [0-9]*/= '"$SYM_N"'/g
}
/nnode = .* ! /{
    s/= [0-9]*/= '"$P"'/g
}' compile/distribution.f90
# get block_size_x/y from domain_size.F90 under current dir, which was copied from bench01d directory
BLOCK_SIZE_X=`cat domain_size.F90 | sed -n '/ block_size_x = .*,/{s/^.*= \([0-9]*\).*$/\1/p}' | head -n 1`
BLOCK_SIZE_Y=`cat domain_size.F90 | sed -n '/ block_size_y = .*/{s/^.*= \([0-9]*\).*$/\1/p}' | head -n 1`
NX_BKFACTOR=1
NY_BKFACTOR=1
# 2. configure block creation parameter
sed -i '/ nx_bkfactor = .* ,/{
    s/= [0-9]*/= '"$NX_BKFACTOR"'/g
}
/ ny_bkfactor = .*/{
    s/= [0-9]*/= '"$NY_BKFACTOR"'/g
}' compile/blocks.f90
sed -i '/ block_size_x = .*,/{
    s/= [0-9]*/= '"$BLOCK_SIZE_X"'/g
}
/ block_size_y = .*/{
    s/= [0-9]*/= '"$BLOCK_SIZE_Y"'/g
}' compile/domain_size.f90
MAXBLKS=`choose_distrb ${S} ${NP}`
sed -i '/max_blocks_.* = /{
    s/= [0-9]*/= '"$MAXBLKS"'/g
}' compile/domain_size.f90
cat compile/domain_size.f90
# 3. configure capcg parameter
CACG_SSTEP=8
CACG_MATLEVEL=1
sed -i '/ cacg_sstep = .*/{
    s/= [0-9]*/= '"$CACG_SSTEP"'/g
}
/ cacg_matlevel = .*/{
    s/= [0-9]*/= '"$CACG_MATLEVEL"'/g
}' compile/blocks.f90
#########################################
# build the CPU executable
cp ../${buildcpucmd} ./
./${buildcpucmd} |& tee buildpurecpu.log
[ ! -f pop ] && echo "build pop (pure CPU) failed." && exit -1
if [ "ARG_BATCH" != "unknown_ext" ]; then
cp -p pop pop.${ARG_BATCH}.p${NP}
fi
popd
else
pushd ${CASE}
[ ! -f pop.${ARG_ARCH}.p${NP} ] && echo "pop.cpu (pure CPU) not exist, need rebuild." && exit -1
cp -p pop.${ARG_ARCH}.p${NP} pop
popd
fi
if [ "$ARG_RUN" == "1" ]; then
#########################################
# setup pop_in for run
pushd ${CASE}
sed -i 's/nprocs_clinic = .*/nprocs_clinic = '"$NP"'/g;s/nprocs_tropic = .*/nprocs_tropic = '"$NP"'/g; ' pop_in
popd
./e3_cpu_bdw.sh bench${S} ${CPU_N} 1 0 ${P} ${ARG_ARCH}
#./e3_cpu.sh bench${S} ${CPU_N} 1 0 ${P}
fi
fi  # ARG_MODE

done
done
done

