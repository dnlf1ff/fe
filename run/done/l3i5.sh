#!/bin/bash
#SBATCH --job-name=l3i5   
#SBATCH --output=l3i5.%j.log   
#SBATCH --error=l3i5.%j.log  
#SBATCH --nodes=1      
#SBATCH --ntasks=1
#SBATCH --time=12:00:00    
#SBATCH --partition=gpu	
#SBATCH --nodelist=n008
#SBATCH --gres=gpu:0

echo "SLURM_NTASKS: $SLURM_NTASKS"

if [ -z "$SLURM_NTASKS" ] || [ "$SLURM_NTASKS" -le 0 ]; then
	echo "Error: SLURM_NTASKS is not set or is less than or equal to 0"
	exit 1
fi


python elastic.py chgTot_l3i5 >> l3i5.x 
