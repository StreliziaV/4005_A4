#!/bin/bash
#SBATCH --job-name=your_job_name # Job name
#SBATCH --nodes=1                    # Run all processes on a single node	
#SBATCH --ntasks=1                   # number of processes = 1 
#SBATCH --cpus-per-task=20      # Number of CPU cores allocated to each process
#SBATCH --partition=Project            # Partition name: Project or Debug (Debug is default)
#SBATCH --output openmp.out         ## filename of the output

# cd /nfsmnt/119010369/4005_A4/
./openmpg 1000 2000 4
# ./openmp 1000 100 20
# ./openmp 1000 100 40
# ./openmp 1000 100 80
# ./openmp 1000 100 120
# ./openmp 1000 100 200