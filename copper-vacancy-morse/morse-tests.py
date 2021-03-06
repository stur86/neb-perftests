#!/usr/bin/env python

# Create input .cell file for the NEB

import sys
import time
import numpy as np
import threading
import subprocess as sp
import argparse as ap
from ase import io
from ase.build import bulk, make_supercell
from ase.calculators.castep import Castep
from morse import MorsePotential
from ase.calculators.socketio import SocketIOCalculator
from ase.neb import NEB, SingleCalculatorNEB, NEBOptimizer
from ase.optimize import BFGS
from ase.optimize.ode import ODE12r
from soprano.collection.generate import linspaceGen
from soprano.properties.linkage import LinkageList
from readts import TSFile

A = 4.0
DE = 1.0   # eV
nnr = 2.55 # Angstrom, nearest neighbour distance

def makeMorseCalc():
    return MorsePotential(epsilon=DE, rho0=A, r0=nnr, apply_cutoff=False)

parser = ap.ArgumentParser()
parser.add_argument('--castep-cmd', '-cc', type=str, default='castep.mpi')
parser.add_argument('--force-recalc', '-f', action='store_true', default=False)
parser.add_argument('--no-castep', '-noc', action='store_true', default=False)
parser.add_argument('--no-ase', '-noa', action='store_true', default=False)
parser.add_argument('--no-socket', '-nos', action='store_true', default=False)

args = parser.parse_args()


try:

    cuvac_reac = io.read('CuVac-107-reac.cell')
    cuvac_prod = io.read('CuVac-107-prod.cell')
    has_opt_files = True

except IOError:

    has_opt_files = False


if not has_opt_files or args.force_recalc:

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
    io.write('CuVac-107.cell', cuvac_reac)

    cuvac_reac.calc = makeMorseCalc()
    cuvac_prod.calc = makeMorseCalc()

    dyn = BFGS(cuvac_reac)
    dyn.run(fmax=1e-3)

    print('Reactants converged at energy {0} eV'.format(cuvac_reac.get_potential_energy()))

    dyn = BFGS(cuvac_prod)
    dyn.run(fmax=1e-3)

    print('Products converged at energy {0} eV'.format(cuvac_prod.get_potential_energy()))
    
    io.write('CuVac-107-reac.cell', cuvac_reac)
    io.write('CuVac-107-prod.cell', cuvac_prod)

# Interpolate
lgen = linspaceGen(cuvac_reac, cuvac_prod, 3, True)
images = [im for im in lgen]
cuvac_intm = images[1] # for testing, use same intermediate and reactant 

# Now go for CASTEP
cuvac_castep = cuvac_reac.copy()
cuvac_castep.calc = Castep(directory='cuvac-neb-morse', label='CuVacMorse', castep_command=args.castep_cmd)
cuvac_castep.calc.cell.positions_abs_product = cuvac_prod
cuvac_castep.calc.cell.positions_abs_intermediate = cuvac_intm
cuvac_castep.calc.cell.fix_com = False
cuvac_castep.calc.cell.species_pot = [('Cu', '../Cu_C19_LDA_OTF.usp')]

# Convert to units
DE_kcal_mol = DE*23.060922344650095
# Derived constant
k = 2*A**2/nnr**2*DE_kcal_mol

n_path_points = 5
n_steps = 10

cuvac_castep.calc.param.task = 'transitionstatesearch'
cuvac_castep.calc.param.cut_off_energy = 100
cuvac_castep.calc.param.iprint = 2
cuvac_castep.calc.param.tssearch_method = 'neb'
cuvac_castep.calc.param.tssearch_neb_max_iter = n_steps
cuvac_castep.calc.param.tssearch_max_path_points = n_path_points
#cuvac_castep.calc.param.tssearch_neb_climbing = True
cuvac_castep.calc.param.tssearch_force_tol = 1e-3
cuvac_castep.calc.param.devel_code = """PP=T
pp: MORS=T MORS_CUT={CR} MORS_R={R} MORS_K={K} MORS_D={D} :endpp
""".format(CR=2.7*nnr, R=nnr, K=k, D=DE_kcal_mol)

if not args.no_castep:
    print('1:  RUNNING CASTEP TEST')
    cuvac_castep.get_potential_energy()

    # Read the TS file
    cuvacts = TSFile('CuVacMorse', 'cuvac-neb-morse')
    cuvacts_rea = cuvacts.blocks['REA']
    cuvacts_tst = cuvacts.blocks['TST']
    cuvacts_pro = cuvacts.blocks['PRO']
    cuvacts_imgs = cuvacts_rea[1] + cuvacts_tst[cuvacts_tst.last_index] + cuvacts_pro[1]

    for i, b in enumerate(cuvacts_imgs):
        # Save the atoms
        io.write('outputs/CuVac-CASTEP-{0}.cell'.format(i), b.atoms)
    # And the energies
    E_castep = [img.atoms.get_potential_energy() for img in cuvacts_imgs]
    np.savetxt('outputs/CuVac-CASTEP-E.dat', E_castep)


# Save the castep calculator, prepare for the ASE NEB
cast_calc = cuvac_castep.calc
cuvac_reac.calc = makeMorseCalc()

images = [cuvac_reac.copy()]
images += [cuvac_reac.copy() for i in range(n_path_points)]
images += [cuvac_prod.copy()]
neb = NEB(images)
# Interpolate linearly the potisions of the three middle images:
neb.interpolate()

# Set calculators:
for image in images:
    image.calc = makeMorseCalc()

# Optimize:
if not args.no_ase:
    print('2:  RUNNING ASE NEB TEST')
    optimizer = NEBOptimizer(neb, logfile='ase-neb.log')
    optimizer.run(fmax=0.04)

    # And write the files
    for i, a in enumerate(images):
        # Save the atoms
        io.write('outputs/CuVac-ASE-{0}.cell'.format(i), a)
    # And the energies
    E_ase = [a.get_potential_energy() for a in images]
    np.savetxt('outputs/CuVac-ASE-E.dat', E_ase)

# Now prepare for the socket
if not args.no_socket:
    cast_calc.param.task = 'socketdriver'
    cast_calc._directory = 'cuvac-socket-morse'
    # We keep the cutoff tiny to avoid initialising a big basis set
    

    def makeSocketCalc():
        return SocketIOCalculator(port=3141)

    cast_calc.prepare_input_files(cuvac_reac)

    images = [cuvac_reac.copy()]
    images += [cuvac_reac.copy() for i in range(n_path_points)]
    images += [cuvac_prod.copy()]
    neb = NEB(images, allow_shared_calculator=True)
    # Interpolate linearly the potisions of the three middle images:
    neb.interpolate()

    # Set calculators:
    scalc = makeSocketCalc()
    for image in images[1:-1]:
        image.calc = scalc

    # Optimize:
    print('3:  RUNNING ASE+CASTEP SOCKET NEB TEST')

    # Run in a parallel thread
    def optimization():
        optimizer = NEBOptimizer(neb, logfile='socket-neb.log')
        optimizer.run(fmax=0.04)

    optthr = threading.Thread(target=optimization)
    optthr.start()

    # And start CASTEP too
    time.sleep(1)
    proc = sp.Popen([args.castep_cmd, 'CuVacMorse'], cwd='./cuvac-socket-morse', stdout=sp.DEVNULL)

    optthr.join()
    print('Thread joined')
    proc.kill()
    print('CASTEP process killed')

    # And write the files
    for i, a in enumerate(images):
        # Change back the calculator (for the energy)
        a.calc = makeMorseCalc()
        # Save the atoms
        io.write('outputs/CuVac-SOCKET-{0}.cell'.format(i), a)
    # And the energies
    E_skt = [a.get_potential_energy() for a in images]
    np.savetxt('outputs/CuVac-SOCKET-E.dat', E_skt)


