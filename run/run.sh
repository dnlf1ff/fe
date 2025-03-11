#!/bin/bash
#SBATCH --job-name=chgt5ot      # Job name
#SBATCH --output=calc.%j.log      # Output file (%j will be replaced with job ID)
#SBATCH --error=calc.%j.log        # Error file (%j will be replaced with job ID)
#SBATCH --nodes=1                                # Number of nodes
#SBATCH --ntasks=8                             # Total number of tasks (MPI processes)
#SBATCH --time=12:00:00                          # Time limit hrs:min:sec
#SBATCH --partition=gpu		  	 # Partition name (update based on your cluster)
#SBATCH --nodelist=n008

echo "SLURM_NTASKS: $SLURM_NTASKS"

if [ -z "$SLURM_NTASKS" ] || [ "$SLURM_NTASKS" -le 0 ]; then
	echo "Error: SLURM_NTASKS is not set or is less than or equal to 0"
	exit 1
fi

python post.py > std.x
