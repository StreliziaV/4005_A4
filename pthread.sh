#!/bin/bash
#SBATCH --job-name=your_job_name # Job name
#SBATCH --nodes=1                    # Run all processes on a single node	
#SBATCH --ntasks=1                   # number of processes = 1 
#SBATCH --cpus-per-task=20      # Number of CPU cores allocated to each process
#SBATCH --partition=Project            # Partition name: Project or Debug (Debug is default)
#SBATCH --output pthread.out         ## filename of the output

# cd /nfsmnt/119010369/4005_A4/
./pthreadg 1000 10000 4
# ./pthread 1000 100 20
# ./pthread 1000 100 40
# ./pthread 1000 100 80
# ./pthread 1000 100 120
# ./pthread 1000 100 200