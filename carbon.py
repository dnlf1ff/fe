from embark.util import freeLoader
from embark.conf import taskCONF
from embark.mlp import sevennRunner

from ase import Atoms
from ase.io import read, write
from ase.build import bulk

conf=taskCONF('.')

def vacency(fe):
    fe_def = fe.copy()
    del fe_def[0]
    return fe_def

def interstitial(fe_def):
    fe_int = fe_def.copy()
    a=fe_int.cell.lengths()[0]
    fe_int.append(Atoms('C', positions=[(a/2,a/2,a/2)]))
    return fe_int

def a1(fe):
    fe_a1 = fe.copy()
    del fe_a1[0]
    a=fe.cell.lengths()[0]
    fe_a1.append(Atoms('C', positions=[(a/2,0,0)]))
    write('poscar/Fe_a1.poscar',fe_a1)
    return fe_a1

def a2(fe):
    fe_a2 = fe.copy()
    del fe_a2[0]
    a=fe.cell.lengths()[0]
    fe_a2.append(Atoms('C', positions=[(0,a/2,a/2)]))
    write('poscar/Fe_a2.poscar',fe_a2)
    return fe_a2

#TODO: b, c
def b1(fe):
    fe_b1 = fe.copy()
    del fe_b1[0]
    a=fe.cell.lengths()[0]
    fe_b1.append(Atoms('C', positions=[(0,a/2,0)]))
    write('poscar/Fe_b1.poscar',fe_b1)
    return fe_b1


def d(fe):
    fe_d = fe.copy()
    a=fe.cell.lengths()[0]
    fe_d.append(Atoms('C', positions=[(a/2,a/2,0)]))
    fe_d.append(Atoms('C', positions=[(a/2,a/2,a)]))
    write('poscar/Fe_d.poscar',fe_d)
    return fe_d

def h(fe):
    fe_h = fe.copy()
    a=fe.cell.lengths()[0]
    fe_h.append(Atoms('C', positions=[(a,a/2,a/2)]))
    fe_h.append(Atoms('C', positions=[(a/2,a/2,a)]))
    write('poscar/Fe_h.poscar',fe_h)
    return fe_h

def e(fe):
    fe_e = fe.copy()
    a=fe.cell.lengths()[0]
    fe_e.append(Atoms('C', positions=[(a/2,a/2,a)]))
    fe_e.append(Atoms('C', positions=[(0,a/2,0)]))
    write('poscar/Fe_e.poscar',fe_e)
    return fe_e

#TODO f

def i(fe):
    fe_i = fe.copy()
    a=fe.cell.lengths()[0]
    fe_i.append(Atoms('C', positions=[(a/2,a/2,a)]))
    fe_i.append(Atoms('C', positions=[(0,0,a/2)]))
    write('poscar/Fe_i.poscar',fe_i)
    return fe_i

def g(fe):
    fe_g = fe.copy()
    a=fe.cell.lengths()[0]
    fe_g.append(Atoms('C', positions=[(0,0,a/2)]))
    fe_g.append(Atoms('C', positions=[(a,0,a/2)]))
    write('poscar/Fe_g.poscar',fe_g)
    return fe_g

def j(fe):
    fe_j = fe.copy()
    del fe_j[0]
    a=fe.cell.lengths()[0]
    fe_j.append(Atoms('C', positions=[(a/2,a/2,0)]))
    fe_j.append(Atoms('C', positions=[(a/2,a/2,a)]))
    write('poscar/Fe_j.poscar',fe_j)
    return fe_j

def k(fe):
    fe_k = fe.copy()
    a=fe.cell.lengths()[0]
    del fe_k[0]
    fe_k.append(Atoms('C', positions=[(a/2,a/2,0)]))
    fe_k.append(Atoms('C', positions=[(0,a/2,a/2)]))
    write('poscar/Fe_k.poscar',fe_k)
    return fe_k

def l(fe):
    fe_l = fe.copy()
    a=fe.cell.lengths()[0]
    del fe_l[0]
    fe_l.append(Atoms('C', positions=[(a/2,a/2,a/4)]))
    fe_l.append(Atoms('C', positions=[(a/2,a/2,a*3/4)]))
    write('poscar/Fe_l.poscar',fe_l)
    return fe_l

def m(fe):
    fe_m = fe.copy()
    a=fe.cell.lengths()[0]
    del fe_m[0]
    fe_m.append(Atoms('C', positions=[(a/2,a/4,a*3/4)]))
    fe_m.append(Atoms('C', positions=[(a/2,a*3/4,a/4)]))
    write('poscar/Fe_m.poscar',fe_m)
    return fe_m

def main():
    jar=conf.load_pickle(f'jar/{mlp}/{mlp}+all',task=False)
    a0=jar['bulk'][f'{mlp}-post'][1][0]
    fe=bulk('Fe','bcc',a=a0).repeat((4,4,4))
    write('fe',fe,format='vasp')
    runner=sevennRunner(mlp)

if __name__ == '__main__':
    main()
