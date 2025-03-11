import os, sys
from embark.util import freeLoader, make_dir
from embark.conf import designatedLoader, taskCONF
from embark.tool import aseTool
from embark.mlp import sevennRunner
import gc
from ase.io import read, write
from ase.build import bulk, make_supercell, surface
from ase.constraints import FixAtoms, FixSymmetry
from ase.optimize import FIRE2
from ase.filters import UnitCellFilter, StrainFilter
import pandas as pd
import numpy as np

cwd = os.getcwd()
conf=taskCONF(dir=cwd)

fe=read("poscar/Fe.poscar")
fe_sc=fe.repeat([4,4,4])
write("poscar/Fe_sc.poscar", fe_sc, format='vasp')

def make_df(mlp):
    df={}
    df['bulk']={}
    task_bulk=['nion','lengths','angles','potential_energy','force','stress','volume']
    df['bulk']['task']=task_bulk
    df['surface']={}
    task_surface=['nion','A_100','A_110','A_111','E_100','E_110','E_111']
    df['surface']['task']=task_surface
    df['vacancy']={}
    task_vac=['nion','E_vac']
    df['vacancy']['task']=task_vac
    df['elastic']={}
    task_elastic=['nion','B','C_prime','C11','C12','C44']
    df['elastic']['task']=task_elastic
    conf.save_pickle(df, f'{mlp}/{mlp}', task=False)
    return df

def _bulk(atoms,runner, mlp, df):
    print('relaxation for 4 by 4 by 4 cell')
    fe_bulk=atoms.copy()
    fe_bulk.calc=runner.calculator

    pre_opt=runner.calc(fe_bulk)
    df['bulk'][f'{mlp}-pre']=pre_opt

    fe_bulk.set_constraint(FixSymmetry(fe_bulk))
    usc=UnitCellFilter(fe_bulk)
    fire=FIRE2(usc,logfile=os.path.join(os.environ['LOG'],f'{mlp}/{mlp}.log'),trajectory=os.path.join(os.environ['TRAJ'],f'{mlp}/{mlp}.traj'))
    fire.run(fmax=1e-3)

    post_opt=runner.calc(fe_bulk)
    write(f'poscar/{mlp}/Fe_{mlp}.poscar', fe_bulk,format='vasp')
    df['bulk'][f'{mlp}-post']=post_opt

    conf.save_pickle(df, f'{mlp}/{mlp}', task=False)
    del fe_bulk
    gc.collect()
    print('relaxation for 4 by 4 by 4 cell -- done')

def vacancy(runner, mlp, df):
    print('letzzs get our vacancy formation E')
    fe_defect=read(f'poscar/{mlp}/Fe_{mlp}.poscar')
    del fe_defect[0]
    fe_defect.calc=runner.calculator
    usc=UnitCellFilter(fe_defect)
    fire=FIRE2(usc, logfile=os.path.join(os.environ['LOG'],f'{mlp}/{mlp}.vac.log'), trajectory=os.path.join(os.environ['TRAJ'],f'{mlp}/{mlp}.vac.traj'))
    fire.run(fmax=1e-3)

    post_opt=runner.calc(fe_defect)
    df['vacancy'][f'{mlp}-post']=post_opt
    conf.save_pickle(df, f'{mlp}/{mlp}+vac', task=False)

    print('relaxation for vacancy -- done')
    del fe_defect
    gc.collect()


def surface(runner, mlp, df):
    a0=conf.load_pickle(f'{mlp}/{mlp}', task=False)['bulk'][f'{mlp}-post'][1][0]
    for hkl in ['100','110','111']:
        slab=runner.surface_energy('Fe',hkl=hkl, a=a0)
        slab.calc=runner.calculator
        write(f'poscar/Fe_{hkl}.poscar', slab, format='vasp')
        pre_opt=runner.calc(slab)
        df['surface'][f'{mlp}-{hkl}-pre']=pre_opt
        conf.save_pickle(df, f'{mlp}/{mlp}+surface', task=False) 

        usc=UnitCellFilter(slab)
        fire=FIRE2(usc,logfile=os.path.join(os.environ['LOG'],f'{mlp}/{mlp}.{hkl}.log'),trajectory=os.path.join(os.environ['TRAJ'],f'{mlp}/{mlp}.{hkl}.traj'))
        fire.run(fmax=1e-3)
        post_opt=runner.calc(slab)
        df['surface'][f'{mlp}-{hkl}-post']=post_opt
        write(f'poscar/{mlp}/Fe_{hkl}_{mlp}.poscar',slab,format='vasp')
        conf.save_pickle(df,f'{mlp}/{mlp}+surface', task=False)

        del slab
        gc.collect()
        print(f'surface energy for {hkl} -- done')
    print('surface energies -- done')
    del a0
    gc.collect()

def surface_bulk(runner, mlp):
    a0=conf.load_pickle(f'{mlp}/{mlp}', task=False)['bulk'][f'{mlp}-post'][1][0]
    df=conf.load_pickle(f'{mlp}/{mlp}+surface', task=False)
    fe=bulk('Fe','bcc',a=a0)
    fe=fe.repeat((4,4,20))
    write('poscar/Fe_surface_bulk.poscar',fe,format='vasp')

    fe.calc=runner.calculator
    pre_opt=runner.calc(fe)
    df['surface'][f'{mlp}-bulk-pre']=pre_opt
    conf.save_pickle(df, f'{mlp}/{mlp}+all', task=False)

    fe.set_constraint(FixSymmetry(fe))
    usc=UnitCellFilter(fe)
    fire=FIRE2(usc,logfile=os.path.join(os.environ['LOG'],f'{mlp}/{mlp}.surface.bulk.log'),trajectory=os.path.join(os.environ['TRAJ'],f'{mlp}/{mlp}.surface.bulk.traj'))
    fire.run(fmax=1e-3)
    post_opt=runner.calc(fe)
    df['surface'][f'{mlp}-bulk-post']=post_opt

    write(f'poscar/{mlp}/Fe_surface_bulk_{mlp}.poscar',fe,format='vasp')
    conf.save_pickle(df,f'{mlp}/{mlp}+all', task=False)

    del a0, fe
    gc.collect()


def makedirs(mlp):
    dirs=['log','traj','jar','poscar']
    for d in dirs:
        make_dir(os.path.join(d,mlp),return_path=False)

def main():
    mlp_list=[sys.argv[1]]
    for mlp in mlp_list:
        makedirs(mlp)
        print('mlp', mlp)
        df=make_df(mlp)
        runner=sevennRunner(mlp)
        _bulk(fe_sc, runner, mlp, df)
    #    strains(runner, mlp, df)
        vacancy(runner, mlp, df)
        surface(runner, mlp, df)
        surface_bulk(runner, mlp)
        conf.save_pickle(df, f'{mlp}/{mlp}+all', task=False)
        del runner, df
        gc.collect()

if __name__ == '__main__':
    main()
