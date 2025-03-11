import numpy as np
from ase.io import read
from embark.util import freeLoader, UNITs
from embark.conf import taskCONF
from embark.tool import aseTool
#TODO: mkdir pre/calc/post

conf=taskCONF('.')
mlps=conf.mlps
fL=freeLoader
asetool=aseTool

def surface_energy(mlp, jar, hkl):
    atoms=read(f'poscar/{mlp}/Fe_{hkl}_{mlp}.poscar', format='vasp')
    A_surface=asetool.surface_area(atoms)
    n_surface=jar['surface'][f'{mlp}-{hkl}-post'][0]
    E_surface=jar['surface'][f'{mlp}-{hkl}-post'][3]
    n_bulk=jar['surface'][f'{mlp}-bulk-post'][0]
    E_bulk=jar['surface'][f'{mlp}-bulk-post'][3]
    if hkl=='100':
        return (E_surface-0.5*E_bulk)/A_surface*2, A_surface
    else:
        return (E_surface-E_bulk)/A_surface*2, A_surface

def vacancy_formation_energy(mlp, jar):
    E_defect=jar['vacancy'][f'{mlp}-post'][3]
    E_bulk=jar['bulk'][f'{mlp}-post'][3]
    n=jar['bulk'][f'{mlp}-post'][0]
    return E_defect-(n-1)*E_bulk/n

def make_df():
    df=pd.DataFrame()
    tasks=['a','E_vac','E_100','A_100','E_110','A_110','E_111','A_111','B','C_prime','C_44']
    df['task']=tasks
    df['unit']=['A','eV','J/m^2','A^2','J/m^2','A^2','J/m^2','A^2','GPa','GPa','GPa']
    conf.save_csv(df, 'pureFe')
    return df

def get_elastic(mlp, jar):
    B=jar['elastic'][f'{mlp}-least_squares']['bulk_voigt']
    C44=jar['elastic'][f'{mlp}-least_squares']['C44']
    C11=jar['elastic'][f'{mlp}-least_squares']['C11']
    C_12=jar['elastic'][f'{mlp}-least_squares']['C12']
    C_prime=(C11-C12)/2
    return B, C_prime, C44

def main():
    mlps=conf.mlps
    df=make_df()
    data=[]
    for mlp in mlps:
        jar=conf.load_pickle(f'{mlp}/{mlp}+all')
        data.append(df['a'][mlp]=jar['bulk'][f'{mlp}-post'][1][0])
        data.append(vacancy_formation_energy(mlp, jar))
        B, C_prime, C44=get_elastic(mlp, jar)
        data.append(B)
        data.append(C_prime)
        data.append(C44)
        for hkl in ['100','110','111']:
            E_surf, A_surf=surface_energy(mlp, jar, hkl)
            data.append(E_surf*16.0217662)
            data.append(A_surf)
    df[mlp]=data
    conf.save_csv(df, 'pureFe')
if __name__=='__main__':
    main()
