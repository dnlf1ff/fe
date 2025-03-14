#!/bin/bash
#SBATCH --job-name=ompaMF   
#SBATCH --output=ompaMF.%j.log   
#SBATCH --error=ompaMF.%j.log  
#SBATCH --nodes=1      
#SBATCH --ntasks=4
#SBATCH --time=12:00:00    
#SBATCH --partition=gpu2080	
#SBATCH --nodelist=Lex100
#SBATCH --gres=gpu:1

module load cuda/12.1.1

echo "SLURM_NTASKS: $SLURM_NTASKS"

if [ -z "$SLURM_NTASKS" ] || [ "$SLURM_NTASKS" -le 0 ]; then
	echo "Error: SLURM_NTASKS is not set or is less than or equal to 0"
	exit 1
fi

export FLASH=1

python elastic.py ompa_MF_bMF_scFT2_mpa_FT >> ompaMF.x
