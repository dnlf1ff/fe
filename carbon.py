from embark.util import freeLoader
from embark.conf import taskCONF
from embark.mlp import sevennRunner

from ase import Atoms
from ase.io import read, write
from ase.build import bulk
from ase.lattice.cubic import BodyCenteredCubic
conf=taskCONF('.')

def create_sc(a, mlp)
    fe_unit=BodyCenteredCubic(directions=[1,0,0],[0,1,0],[0,0,1]], size=(1,1,1), symbol='Fe', pbc=True, latticeconstant=a)
    fe_sc=fe_unit.repeat((4,4,4))
    write(f'poscar/{mlp}_sc.poscar',sc, format='vasp')
    return fe_sc

def a1(fe,a,mlp):
    fe_a1 = fe.copy()
    del fe_a1[1]
    fe_a1.append(Atoms('C', positions=[(a/2,a/2,0)]))
    write(f'poscar/{mlp}_a1.poscar',fe_a1, format='vasp')
    return fe_a1

def a2(fe,a,mlp):
    fe_a2 = fe.copy()
    del fe_a2[1]
    fe_a2.append(Atoms('C', positions=[(0,0,a/2)]))
    write(f'poscar/{mlp}_a2.poscar',fe_a2,format='vasp')
    return fe_a2

#TODO: b, c
def b1(fe, a, mlp):
    fe_b1 = fe.copy()
    del fe_b1[1]
    del fe_b1[3]
    fe_b1.append(Atoms('C', positions=[(a/2,a/2,a)])))
    write(f'poscar/{mlp}_b1.poscar',fe_b1,format='vasp')
    return fe_b1

def b2(fe, a, mlp):
    fe_b2=fe.copy()
    del fe_b2[1]
    del fe_b2[3]
    fe_b2.append(Atoms('C',positions=[(0,0,a/2)]))
    write(f'poscar/{mlp}_b2.poscar',fe_b2,format='vasp')
    return fe_b2

def b3(fe, a, mlp):
    fe_b3=fe.copy()
    del fe_b3[1]
    del fe_b3[3]
    fe_b3.append(Atoms('C',positions=[(a/2,a/2,0)])
    write(f'poscar/{mlp}_b3.poscar', fe_b3, format='vasp')
    return fe_b3

def c1(fe, a, mlp):
    fe_c1=fe.copy()
    del fe_c1[0]
    del fe_c1[1]
    fe_c1.append(Atoms('C',positions=[(a/2,a/2,a)]))
    write(f'poscar/{mlp}_c1.poscar',fe_c1, format='vasp')
    return fe_c1


def c2(fe, a, mlp):
    fe_c2=fe.copy()
    del fe_c2[0]
    del fe_c2[1]
    fe_c1.append(Atoms('C',positions=[(a/2,a/2,0)]))
    write(f'poscar/{mlp}_c2.poscar',fe_c2 format='vasp')
    return fe_c2


def d(fe, a, mlp):
    fe_d = fe.copy()
    fe_d.append(Atoms('C', positions=[(a/2,a/2,0)]))
    fe_d.append(Atoms('C', positions=[(a/2,a/2,a)]))
    write(f'poscar/{mlp}_d.poscar',fe_d, format='vasp')
    return fe_d


def e(fe,a,mlp):
    fe_e = fe.copy()
    fe_e.append(Atoms('C', positions=[(a/2,a/2,a)]))
    fe_e.append(Atoms('C', positions=[(a/2,a,a/2)]))
    write(f'poscar/{mlp}_e.poscar',fe_e, format='vasp')
    return fe_e

def f(fe, a, mlp):
    fe_f=fe.copy()
    fe.append(Atoms('C',positions=[(0,a/2,0)]))
    fe.append(Atoms('C',positions=[(a,3*a/2,a)]))
    write(f'poscar/{mlp}_f.poscar',fe_f,format='vasp')

def i(fe,a,mlp):
    fe_i = fe.copy()
    fe_i.append(Atoms('C', positions=[(a/2,a/2,a)]))
    fe_i.append(Atoms('C', positions=[(0,0,a/2)]))
    write(f'poscar/{mlp}_i.poscar',fe_i,format='vasp')
    return fe_i

def g(fe,a,mlp):
    fe_g = fe.copy()
    fe_g.append(Atoms('C', positions=[(0,0,a/2)]))
    fe_g.append(Atoms('C', positions=[(0,a,a/2)]))
    write(f'poscar/{mlp}_g.poscar',fe_g,format='vasp')
    return fe_g

def h(fe, a, mlp):
    fe_h=fe.copy()
    fe_h.append(Atoms('C',positions=[(a/2,a/2,a)]))
    fe_h.append(Atoms('C',positions=[(a/2,a,a/2)]))
    write(f'poscar/{mlp}_h.poscar', fe_h, format='vasp')

def j(fe,a,mlp):
    fe_j = fe.copy()
    del fe_j[1]
    fe_j.append(Atoms('C', positions=[(a/2,a/2,0)]))
    fe_j.append(Atoms('C', positions=[(a/2,a/2,a)]))
    write(f'poscar/{mlp}_j.poscar',fe_j,format='vasp')
    return fe_j

def k(fe,a,mlp):
    fe_k = fe.copy()
    del fe_k[1]
    fe_k.append(Atoms('C', positions=[(a/2,a/2,0)]))
    fe_k.append(Atoms('C', positions=[(0,a/2,a/2)]))
    write(f'poscar/{mlp}_k.poscar',fe_k,format='vasp')
    return fe_k

def l(fe,a,mlp):
    fe_l = fe.copy()
    del fe_l[1]
    fe_l.append(Atoms('C', positions=[(a/2,a/2,a/4)]))
    fe_l.append(Atoms('C', positions=[(a/2,a/2,a*3/4)]))
    write(f'poscar/{mlp}_l.poscar',fe_l,format='vasp')
    return fe_l

def m(fe,a,mlp):
    fe_m = fe.copy()
    del fe_m[1]
    fe_m.append(Atoms('C', positions=[(a/2,a/4,a*3/4)]))
    fe_m.append(Atoms('C', positions=[(a/2,a*3/4,a/4)]))
    write(f'poscar/{mlp}_m.poscar',fe_m,format='vasp')
    return fe_m

def vac_1(fe,mlp):
    fe_vac1=fe.copy()
    del fe_vac1[1]
    write(f'poscar/{mlp}_vac1.posar',fe_vac1, format='vasp')
    return fe_vac1

def vac_13(fe,mlp)
    fe_vac13=fe.copy()
    del fe_vac13[1]
    del ve_vac13[3]
    write(f'poscsar/{mlp}_vac13.poscar',fe_vac13,format='vasp')

    return vac_13


def


def relax(atoms, runner):
    atoms.calc=runner.calculator
    atoms.set_constraints(FixSymmetry(atoms))
    ucf=UnitCellFilter(atoms)
    fire=FIRE2(ucf, logfile=f'log/{mlp}/C.{mlp}.log',traj=f'traj/{mlp}/C.{mlp}.traj')
    write(f'structure/{mlp}/C_{mlp}_opt.poscar',atoms,format='vasp')
    post_opt=runner.calc(atoms)
    return post_opt



def run(mlp,runner):
    jar=conf.load_pickle(f'jar/{mlp}/{mlp}+all',task=False)
    a=jar['bulk'][f'{mlp}-post'][1][0]
    
    fe=create_sc(a,mlp) 
    
    fe_a1=a1(fe,a,mlp)



    fe_a2=a2(fe,a,mlp)
    fe_b1=b1(fe,a,mlp)
    fe_b2=b2(fe,a,mlp)
    fe_b3=b3(fe,a,mlp)
    fe_c1=c1(fe,a,mlp)
    fe_c2=c2(fe,a,ml)
    fe_d=d(fe,a,mlp)
    fe_e=e(fe,a,mlp)
    fe_f=f(fe,a,mlp)
    fe_g=g(fe,a,mlp)
    fe_h=h(fe,a,mlp)
    fe_i=i(fe,a,mlp)
    fe_j=j(fe,a,mlp)
    fe_k=k(fe,a,mlp)
    fe_l=l(fe,a,mlp)
    fe_m=m(fe,a,mlp)

def main():
    jar=conf.load_pickle(f'jar/{mlp}/{mlp}+all',task=False)
    a0=jar['bulk'][f'{mlp}-post'][1][0]
    write('fe',fe,format='vasp')
    runner=sevennRunner(mlp)

if __name__ == '__main__':
    main()
