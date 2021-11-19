#!/usr/bin/env python

# Create input .cell file for the NEB

import numpy as np
import ase.io
import io
import clipboard
from ase.build import bulk, make_supercell
from ase.calculators.castep import Castep
from morse import MorsePotential
from ase.optimize import BFGS
from soprano.collection.generate import linspaceGen
from soprano.properties.linkage import LinkageList
from readts import TSFile

def read_forces(out, n_atoms):
    forces = []
    fix_cart = []
    fix = []
    while True:
        line = out.readline()
        fields = line.split()
        if len(fields) == 7:
            break
    for n in range(n_atoms):
        consd = np.array([0, 0, 0])
        fxyz = [0, 0, 0]
        for (i, force_component) in enumerate(fields[-4:-1]):
            if force_component.count("(cons'd)") > 0:
                consd[i] = 1
            fxyz[i] = float(force_component.replace(
                "(cons'd)", ''))
        if consd.all():
            fix.append(n)
        elif consd.any():
            fix_cart.append(FixCartesian(n, consd))
        forces.append(fxyz)
        line = out.readline()
        fields = line.split()
        
    return np.array(forces)


def read_neb_data(out):
    image, dof, tangent, real, real_proj, e_band, total = np.loadtxt(out, unpack=True)
    image = image.astype(int)
    dof = dof.astype(int)
    return image, dof, tangent, real, real_proj, e_band, total

def paste_forces(n_atoms):
    text = clipboard.paste()
    return read_forces(io.StringIO(text), n_atoms)

def paste_neb():
    text = clipboard.paste()
    return read_neb_data(io.StringIO(text))

A = 4.0
DE = 1.0   # eV
nnr = 2.55 # Angstrom, nearest neighbour distance

def makeMorseCalc():
    return MorsePotential(epsilon=DE, rho0=A, r0=nnr, apply_cutoff=False)

try:

    f_reac = io.read('CuVac-107-reac.cell')
    f_prod = io.read('CuVac-107-prod.cell')

except IOError:

    cu = bulk('Cu', a=nnr*2**0.5, cubic=True)

    cu3x3 = make_supercell(cu, np.eye(3)*3)

    # Remove the most central atom to create a vacancy
    C = np.ones(3)*cu3x3.get_cell()[0,0]*0.5
    pos = cu3x3.positions
    i1vac = np.argmin(np.linalg.norm(pos-C, axis=1))
    p1vac = pos[i1vac]

    cuvac = cu3x3.copy()
    del(cuvac[i1vac])

    # Now the one closest to it
    i2vac = np.argmin(np.linalg.norm(cuvac.positions-p1vac, axis=1))

    cuvac_reac = cuvac.copy()
    cuvac_prod = cuvac_reac.copy()
    cuvac_prod[i2vac].position = p1vac

    # Now optimise with Morse Potential
    # Save an unoptimized copy
    ase.io.write('CuVac-107.cell', cuvac_reac)

    cuvac_reac.calc = makeMorseCalc()
    cuvac_prod.calc = makeMorseCalc()

    dyn = BFGS(cuvac_reac)
    f_reac = cuvac_reac.get_forces()
    dyn.run(fmax=1e-3)

    print('Reactants converged at energy {0} eV'.format(cuvac_reac.get_potential_energy()))

    dyn = BFGS(cuvac_prod)
    f_prod = cuvac_prod.get_forces()
    dyn.run(fmax=1e-3)

    print('Products converged at energy {0} eV'.format(cuvac_prod.get_potential_energy()))
    
    ase.io.write('CuVac-107-reac.cell', cuvac_reac)
    ase.io.write('CuVac-107-prod.cell', cuvac_prod)



# Interpolate
lgen = linspaceGen(cuvac_reac, cuvac_prod, 3, True)
images = [im for im in lgen]
cuvac_intm = images[1] # for testing, use same intermediate and reactant 

# Now go for CASTEP
cuvac_reac.calc = Castep(directory='cuvac-neb-morse', label='CuVacMorse')
cuvac._rename_existing_dir = False
cuvac_reac.calc.cell.positions_abs_product = cuvac_prod
cuvac_reac.calc.cell.positions_abs_intermediate = cuvac_intm
cuvac_reac.calc.cell.fix_com = False
cuvac_reac.calc.cell.species_pot = [('Cu', '../Cu_C19_LDA_OTF.usp')]

# Convert to units
DE_kcal_mol = DE*23.060922344650095
# Derived constant
k = 2*A**2/nnr**2*DE_kcal_mol

n_steps = 10

cuvac_reac.calc.param.task = 'transitionstatesearch'
cuvac_reac.calc.param.iprint = 2
cuvac_reac.calc.param.tssearch_method = 'neb'
# cuvac_reac.calc.param.tssearch_neb_method = 'fire'
cuvac_reac.calc.param.tssearch_neb_max_iter = n_steps
cuvac_reac.calc.param.tssearch_max_path_points = 1
#cuvac_reac.calc.param.tssearch_neb_climbing = True
cuvac_reac.calc.param.tssearch_force_tol = 1e-3
cuvac_reac.calc.param.devel_code = """PP=T
pp: MORS=T MORS_CUT={CR} MORS_R={R} MORS_K={K} MORS_D={D} :endpp
""".format(CR=2.7*nnr, R=nnr, K=k, D=DE_kcal_mol)

cuvac_reac.calc.prepare_input_files(cuvac_reac)

import os, subprocess

for outfile in ['CuVacMorse.castep', 'CuVacMorse.ts']:
    outfile = os.path.join('cuvac-neb-morse', outfile)
    if os.path.exists(outfile):
        os.os.unlink(outfile)

castep = subprocess.Popen(['/Users/u1470235/gits/castep/obj/darwin_arm64_gfortran10--serial/castep.serial', 'CuVacMorse'],
                            cwd='cuvac-neb-morse')
castep.communicate()

# finally, run with ASE NEB using ode12r
from ase.neb import NEB, NEBOptimizer

neb = NEB(images, method='spline')

if os.path.exists('cuvac-neb-morse/CuVacMorse.ts'):
    tsf = TSFile('CuVacMorse', path='cuvac-neb-morse')
    rea = tsf.blocks['REA']
    tst = tsf.blocks['TST']
    pro = tsf.blocks['PRO']
    for image, src in zip(neb.images, [rea, tst, pro]):
        image.set_positions(src.get_positions(1)[0,:,:])

for image in neb.images:
    image.calc = makeMorseCalc()

f0_real = [ img.get_forces() for img in neb.images ]
f0_neb = neb.get_forces()

print(f'Number of images: {len(neb.images)}')

from ase.optimize import FIRE

opt = NEBOptimizer(neb, verbose=2)

f_neb = []
def store_forces():
    f_neb.append(neb.get_forces())

opt.attach(store_forces)

opt.run(fmax=1e-3, steps=n_steps)

f1_real = [ img.get_forces() for img in neb.images ]
f1_neb = neb.get_forces()
