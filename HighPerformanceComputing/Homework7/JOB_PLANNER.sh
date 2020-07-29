#!/usr/bin/env bash
delay="${1:-'0'}"
job_com="${2:-main}"

timeLimit="0-00:10"
echo "Delay = ${delay} hours"
echo "Time Limit = ${timeLimit}"

OVERHEAD_THREADS=28
OVERHEAD_CPUS=14
OVERHEAD_JOBS=14
MKL_CPUS=1
outfile0=results
Nprocs=16

rm -rf $outfile
for (( N=2; N<=${Nprocs}; N=N+2))
do
  if (( N < 10 )); then
    outfile="${outfile0}0$N.res"
  else
    outfile="${outfile0}$N.res"
  fi
  echo "===========================" >> ${outfile}
  echo "With $N Processing Elements" >> ${outfile}
  echo "===========================" >> ${outfile}

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

module load intel-mpi/2018x

export MKL_NUM_THREADS=$MKL_CPUS
ulimit -s unlimited

export OMP_STACKSIZE=16M
(time mpirun -n $N ./mpihenon.exe henon.txt) &>> ${outfile} 
EOF
done

while [[ $(squeue -u $USER | wc -l) -gt "1" ]]
do
  sleep 3s
done
for filename in ./*.res;
do
  cat $filename >> results.res
  rm $filename
done
#pcmd=$HOME/.local/bin/parallel
#\$pcmd -v -j $OVERHEAD_JOBS --col-sep='\s+' my_job :::: $job_rec 
