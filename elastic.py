import os, sys
from glob import glob
import numpy as np

from embark.util import freeLoader
from embark.conf import taskCONF
from embark.mlp import sevennRunner

from ase.units import GPa
from ase.io import read, write
from ase.optimize import FIRE2
from ase.filters import StrainFilter

from matscipy.elasticity import measure_triclinic_elastic_constants, fit_elastic_constants, full_3x3x3x3_to_Voigt_6x6
from modulus_util import *

cwd=os.getcwd()
conf=taskCONF(cwd)

fe=read('poscar/Fe.poscar', format='vasp')

def strain_opt(runner, mlp):
    fe_bulk=fe.copy()
    fe_bulk=fe_bulk.repeat((4,4,4))
    fe_bulk.calc=runner.calculator
    jar=conf.load_pickle(f'{mlp}/{mlp}+all',task=False)

    pre_opt=runner.calc(fe_bulk)
    jar['elastic'][f'{mlp}-pre']=pre_opt
    sf=StrainFilter(fe_bulk)
    fire=FIRE2(sf,logfile=os.path.join(os.environ['LOG'],mlp,f'{mlp}.strain.log'),trajectory=os.path.join(os.environ['TRAJ'],mlp,f'{mlp}.strain.traj'))
    fire.run(fmax=1e-3)

    post_opt=runner.calc(fe_bulk)
    jar['elastic'][f'{mlp}-post']=post_opt

    write(f'poscar/{mlp}/{mlp}.strain.poscar',fe_bulk,format='vasp')
    conf.save_pickle(jar, f'{mlp}/{mlp}+elastic',task=False)

def elastics(runner, mlp):
    fe_bulk=read(f'poscar/{mlp}/{mlp}.strain.poscar', format='vasp')
    fe_bulk.calc=runner.calculator

    C_least_squares, _ = fit_elastic_constants(fe_bulk, verbose=False, logfile=os.path.join(os.environ['LOG'],mlp,f'{mlp}.fit.elastic.log'),optimizer=FIRE2, steps=100, fmax=1e-3,delta=0.005,N_steps=5)
    jar=conf.load_pickle(f'{mlp}/{mlp}+elastic',task=False)
    post_opt=runner.calc(fe_bulk)
    jar['elastic'][f'{mlp}-least_squares']=post_opt
#    jar['elastic']['task_fr']=['C11', 'C12', 'C13', 'C14', 'C15', 'C16','C22', 'C23', 'C24', 'C25', 'C26', 'C33', 'C34', 'C35', 'C36','C44', 'C45', 'C46', 'C55', 'C56', 'C66', 'bulk_voigt', 'bulk_reuss', 'bulk_vrh','shear_voigt', 'shear_reuss', 'shear_vrh']

    elastic_df={}
    for i in range(6):
        for j in range(6):
            elastic_df[f'C{i+1}{j+1}'] = C_least_squares[i][j]/GPa
    elastic_df['bulk_voigt'] = get_voigt_bulk_modulus(C_least_squares)/GPa
    try:
        elastic_df['bulk_reuss'] = get_reuss_bulk_modulus(np.linalg.inv(C_least_squares/GPa))
        elastic_df['bulk_vrh'] = get_vrh_bulk_modulus(C_least_squares)/GPa
    except:
        elastic_df['bulk_reuss'] = np.nan
        elastic_df['bulk_vrh'] = np.nan
    elastic_df['shear_voigt'] = get_voigt_shear_modulus(C_least_squares)/GPa
    try:
        elastic_df['shear_reuss'] = get_reuss_shear_modulus(np.linalg.inv(C_least_squares/GPa))
        elastic_df['shear_vrh'] = get_vrh_shear_modulus(C_least_squares)/GPa
    except:
        elastic_df['shear_reuss']=np.nan
        elastic_df['shear_vrh']=np.nan
    jar['elastic'][f'{mlp}-least_squares']=elastic_df
    conf.save_pickle(jar, f'{mlp}/{mlp}+elastic',task=False)
    conf.save_pickle(jar, f'{mlp}/{mlp}+all',task=False)

if __name__ == '__main__':
    mlp=sys.argv[1]
    mlp_list=[mlp]
    for mlp in mlp_list:
        print(f'{mlp} start')
        runner=sevennRunner(mlp)
        strain_opt(runner, mlp)
        elastics(runner, mlp)
