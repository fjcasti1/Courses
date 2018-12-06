#!/usr/bin/env bash
delay="${1:-'0'}"
job_com="${2:-main}"

timeLimit="0-00:06"
echo "Delay = ${delay} hours"
echo "Time Limit = ${timeLimit}"

OVERHEAD_THREADS=28
job_prefix=Coarray
outfile0=resultsCoarray
Nprocs=16
outfile="${outfile0}${Nprocs}.res"

rm -rf $outfile

echo "===========================" >> ${outfile}
echo "With ${Nprocs} Processing Elements" >> ${outfile}
echo "===========================" >> ${outfile}

sbatch --comment="${job_prefix} ${job_com}" << EOF
#!/bin/bash
#SBATCH -p parallel
#SBATCH -t ${timeLimit}
#SBATCH --begin=now+${delay}hour
#SBATCH --nodes=1-1
#SBATCH --ntasks=1
#SBATCH --mincpus=$OVERHEAD_THREADS
#SBATCH --mail-type ALL
#SBATCH --mail-user fjcasti1@asu.edu
#SBATCH -o "testCoarray.out"
#SBATCH -e "testCoarray.err"

module load intel-mpi/2018x
export FOR_COARRAY_NUM_IMAGES=${Nprocs}

ulimit -s unlimited

(time mpirun -n 1 jacobicoarray.exe) &>> ${outfile} 
EOF
