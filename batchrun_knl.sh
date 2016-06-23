#!/bin/bash

# Argument = [-t|--batch "avx|hsw|knl|default=<null>binary ext"] [-p|--proclist "num1,num2,...,numN|default=28,36"]
usage()
{
  echo $0 '[-t|--batch "avx|hsw|knl|default=<null>binary ext"] [-p|--proclist "num1,num2,...,numN|default=28,36"]'
}
ARG_BATCH=knl_opt_bench01
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

ulimit -s unlimited
export OMP_STACKSIZE=1000M

RUNDIR=bench01opt_cpu_n1t1p1
pushd $RUNDIR
cp -p pop_in.knl pop_in

for NP in $PROCLIST; do
sed -i 's/nprocs_clinic = .*/nprocs_clinic = '"$NP"'/g;s/nprocs_tropic = .*/nprocs_tropic = '"$NP"'/g; ' pop_in

((DM=256/NP))
mpirun -n ${NP} -env I_MPI_PIN_DOMAIN=${DM} -env OMP_NUM_THREADS=1 ./pop.${EXT}.p${NP} |& tee pop_${EXT}_p${NP}.log

done

popd
