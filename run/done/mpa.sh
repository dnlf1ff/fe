#!/bin/bash
#SBATCH --job-name=mpa   
#SBATCH --output=mpa.%j.log   
#SBATCH --error=mpa.%j.log  
#SBATCH --nodes=1      
#SBATCH --ntasks=1
#SBATCH --time=12:00:00    
#SBATCH --partition=gpu2	
#SBATCH --nodelist=n014
#SBATCH --gres=gpu:0

echo "SLURM_NTASKS: $SLURM_NTASKS"

if [ -z "$SLURM_NTASKS" ] || [ "$SLURM_NTASKS" -le 0 ]; then
	echo "Error: SLURM_NTASKS is not set or is less than or equal to 0"
	exit 1
fi

export FLASH=1

python elastic.py mpa >> mpa.x 
