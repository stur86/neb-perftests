#!/usr/bin/env python

# Create input .cell file for the NEB

import numpy as np
from ase import io
from ase.build import bulk, make_supercell
from ase.calculators.castep import Castep
from ase.calculators.morse import MorsePotential
from ase.optimize import BFGS
from soprano.collection.generate import linspaceGen
from soprano.properties.linkage import LinkageList

A = 4.0
DE = 1.0   # eV
nnr = 2.55 # Angstrom, nearest neighbour distance

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
def makeMorseCalc():
    return MorsePotential(epsilon=DE, rho0=A, r0=nnr)

cuvac_reac.calc = makeMorseCalc()
cuvac_prod.calc = makeMorseCalc()

dyn = BFGS(cuvac_reac)
dyn.run(fmax=1e-3)

dyn = BFGS(cuvac_prod)
dyn.run(fmax=1e-3)

# Interpolate
lgen = linspaceGen(cuvac_reac, cuvac_prod, 3, True)
images = [im for im in lgen]
cuvac_intm = images[1]

# run with ASE NEB using ode12r
from ase.neb import NEB, NEBOptimizer

neb = NEB(images, method='spline')
for image in neb.images:
    image.calc = makeMorseCalc()

opt = NEBOptimizer(neb)
opt.run(fmax=1e-3)

# Now go for CASTEP
cuvac_reac.calc = Castep(directory='cuvac-neb-morse', label='CuVacMorse')
cuvac_reac.calc.cell.positions_abs_product = cuvac_prod
cuvac_reac.calc.cell.positions_abs_intermediate = cuvac_intm
cuvac_reac.calc.cell.fix_com = False

# Morse potential

# Convert to units
DE_kcal_mol = DE*23.060922344650095
# Derived constant
k = 2*A**2/nnr**2*DE_kcal_mol

cuvac_reac.calc.param.task = 'transitionstatesearch'
cuvac_reac.calc.param.tssearch_method = 'neb'
cuvac_reac.calc.param.tssearch_neb_max_iter = 100
cuvac_reac.calc.param.tssearch_max_path_points = 5
#cuvac_reac.calc.param.tssearch_neb_climbing = True
cuvac_reac.calc.param.tssearch_force_tol = 1e-3
cuvac_reac.calc.param.devel_code = """PP=T
pp: MORS=T MORS_CUT={CR} MORS_R={R} MORS_K={K} MORS_D={D} :endpp
""".format(CR=2.2*nnr, R=nnr, K=k, D=DE_kcal_mol)

cuvac_reac.calc.prepare_input_files(cuvac_reac)