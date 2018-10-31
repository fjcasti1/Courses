#!/usr/bin/env bash
delay="${1:-'0'}"
job_com="${2:-main}"

timeLimit="0-00:21"
echo "Delay = ${delay} hours"
echo "Time Limit = ${timeLimit}"

OVERHEAD_THREADS=28
OVERHEAD_CPUS=14
OVERHEAD_JOBS=14
MKL_CPUS=1

trapped() {
  echo 'TRAPPED -- QUITTING'
  exit 70
}

#my_job() {

#}

#export -f my_job

trap "trapped" 1 2 3 4 5 6 7 8 

sbatch --comment="Sweep ${job_prefix} ${job_com}" << EOF
#!/bin/bash
#SBATCH -p parallel
#SBATCH -t ${timeLimit}
#SBATCH --begin=now+${delay}hour
#SBATCH --nodes=1-1
#SBATCH --ntasks=1
#SBATCH --mincpus=$OVERHEAD_THREADS
#SBATCH --mail-type ALL
#SBATCH --mail-user fjcasti1@asu.edu
#SBATCH -o "testOMP.out"
#SBATCH -e "testOMP.err"

module load intel/2018x

export MKL_NUM_THREADS=$MKL_CPUS
ulimit -s unlimited

export OMP_STACKSIZE=16M
echo "Jobs:"
export OMP_NUM_THREADS=4
time ./test_henon.exe >> out4.log 
export OMP_NUM_THREADS=2
time ./test_henon.exe >> out2.log 
export OMP_NUM_THREADS=1
time ./test_henon.exe >> out1.log 
EOF
#pcmd=$HOME/.local/bin/parallel
#\$pcmd -v -j $OVERHEAD_JOBS --col-sep='\s+' my_job :::: $job_rec 
