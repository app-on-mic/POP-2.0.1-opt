#!/bin/bash

# Argument = [-t|--batch "avx|hsw|knl|default=<null>binary ext"] [-p|--proclist "num1,num2,...,numN|default=28,36"]
usage()
{
  echo $0 '[-t|--batch "avx|hsw|knl|default=<null>binary ext"] [-p|--proclist "num1,num2,...,numN|default=28,36"]'
}
ARG_BATCH=hsw_opt_bench01
ARG_PROCLIST=28,36

while [[ $# -ge 1 ]]
do
key="$1"

case $key in
    -t|--batch)
    ARG_BATCH="$2"
    shift # past argument
    ;;
    -p|--proclist)
    ARG_PROCLIST="$2"
    shift # past argument
    ;;
    *)
    usage && exit -1
            # unknown option
    ;;
esac
shift # past argument or value
done

#echo "ARG_BATCH=" $ARG_BATCH
#echo "ARG_PROCLIST=" $ARG_PROCLIST
#exit 0

EXT=$ARG_BATCH
PROCLIST=`echo $ARG_PROCLIST | sed -e "s/[,;]/ /g"`

############################################################
# check core# and generate exclude list string
# input: $1: hostname
# output: exclude list string
############################################################
function gen_excl_list()
{
MPI_CORE=`ssh $1 lscpu | grep "Core(s) per socket" | awk '{print $NF}'`
MPI_SOCKET=`ssh $1 lscpu | grep "Socket(s):" | awk '{print $NF}'`
MPI_HT=`ssh $1 lscpu | grep "Thread(s) per core:" | awk '{print $NF}'`
EXCL_LIST=
((EXC_SINGLE=(MPI_CORE*MPI_SOCKET-CPU_N*CPU_T)/MPI_SOCKET))
if ((EXC_SINGLE > 0)); then
  for ((i=0;i<MPI_HT;i++));do
    for ((j=0;j<MPI_SOCKET;j++));do
      [[ ! -z "$EXCL_LIST" ]] && EXCL_LIST=${EXCL_LIST}","
      ((STOP=(i*MPI_SOCKET+j+1)*MPI_CORE-1))
      ((START=STOP-EXC_SINGLE+1))
      EXCL_LIST=${EXCL_LIST}"$START-$STOP"
    done
  done
fi
  echo $EXCL_LIST
}

ulimit -s unlimited
export OMP_STACKSIZE=1000M

RUNDIR=bench01opt_cpu_n1t1p1
pushd $RUNDIR
cp -p pop_in.hsw pop_in

for NP in $PROCLIST; do
sed -i 's/nprocs_clinic = .*/nprocs_clinic = '"$NP"'/g;s/nprocs_tropic = .*/nprocs_tropic = '"$NP"'/g; ' pop_in

((CPU_N=NP))
((CPU_T=1))
EXCL_LIST=`gen_excl_list "localhost"`

mpirun -n ${NP} -env OMP_NUM_THREADS=1 -env I_MPI_PIN_PROCESSOR_LIST=allcores:map=bunch -env I_MPI_PIN_PROCESSOR_EXCLUDE_LIST=${EXLC_LIST} ./pop.${EXT}.p${NP} |& tee pop_${EXT}_p${NP}.log

done

popd
