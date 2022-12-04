#!/bin/bash
#SBATCH --job-name=mpip # Job name
#SBATCH --nodes=1                    # Run all processes on a single node	
#SBATCH --ntasks=20                   # number of processes = 20
#SBATCH --cpus-per-task=1      # Number of CPU cores allocated to each process (please use 1 here, in comparison with pthread)
#SBATCH --partition=Project            # Partition name: Project or Debug (Debug is default)
#SBATCH --output mpip.out         ## filename of the output

# cd /nfsmnt/119010355/CSC4005_2022Fall_Demo/project4_template/
# mpirun -np num_process ./mpig size num_thread max_iteration
mpirun -np 4 ./mpipg 1000 2 1000
# mpirun -np 20 ./mpi 1000 100
# mpirun -np 40 ./mpi 1000 100