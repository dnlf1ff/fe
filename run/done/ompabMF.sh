#!/bin/bash
#SBATCH --job-name=ompabMF   
#SBATCH --output=ompabMF.%j.log   
#SBATCH --error=ompabMF.%j.log  
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
python elastic.py ompa_bMF_fp3_felsc_ft1 >> ompabMF.x
