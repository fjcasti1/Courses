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
Nprocs=10
outfile=results.res

echo "=========================" >> ${outfile}
echo "With $Nprocs Processing Element" >> ${outfile}
echo "=========================" >> ${outfile}

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
echo "$Nprocs"
(time mpirun -n $Nprocs ./mpihenon.exe henon.txt) &>> ${outfile} 
EOF
