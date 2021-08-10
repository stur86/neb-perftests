import os
import numpy as np
from ase.build import bulk
from ase.optimize import BFGS
from ase.neb import NEB
from ase.io import read, Trajectory
from ase.geometry.geometry import get_distances
from ase.calculators.emt import EMT
from ase.calculators.castep import Castep
from ase.constraints import ExpCellFilter


def calc():
    return EMT()

def make_geomopt(a, name, force=False):

    from ase.io.castep import write_param

    calc = Castep(directory='cuvac-geomopt-{0}'.format(name),
                  label=name)
    if os.path.isdir(calc._directory) and not force:
        print('WARNING: path ' + calc._directory + ' exists, skipping')
        return
    calc.param.task = 'geometryoptimisation'
    calc.param.xc_functional = 'pbe'
    calc.param.max_scf_cycles = 100
    calc.param.cut_off_energy = 500.0
    calc.param.calculate_stress = True
    calc.param.write_cell_structure = True
    calc.cell.fix_all_cell = True
    calc.cell.kpoint_mp_grid = "4 4 4"

    a = a.copy()
    a.calc = calc

    calc.prepare_input_files(a, True)
    # Overwrite cell (we need fractional coordinates)
    a.write(calc._abs_path('{0}.cell'.format(name)), positions_frac=True)
    



N_intermediate = 3
N_cell = 2
initial = bulk('Cu', cubic=True)
initial *= N_cell

initial.calc = calc()

print("Starting bulk lattice parameter: {0} Ang".format(initial.get_cell()[0,0]))

# Relax the bulk
ecf = ExpCellFilter(initial, hydrostatic_strain=True)
qn = BFGS(ecf)
traj = Trajectory('CuVacBulk.traj', 'w', initial)
qn.attach(traj)
qn.run(fmax=0.05)

print("Relaxed bulk lattice parameter: {0} Ang".format(initial.get_cell()[0,0]))

initial.write('CuBulk.cell', positions_frac=True)

# place vacancy near centre of cell
D, D_len = get_distances(np.diag(initial.cell) / 2,
                         initial.positions,
                         initial.cell, initial.pbc)
vac_index = D_len.argmin()
vac_pos = initial.positions[vac_index]
del initial[vac_index]

# identify two opposing nearest neighbours of the vacancy
D, D_len = get_distances(vac_pos,
                         initial.positions,
                         initial.cell, initial.pbc)
D = D[0, :]
D_len = D_len[0, :]

nn_mask = np.abs(D_len - D_len.min()) < 1e-8
i1 = nn_mask.nonzero()[0][0]
i2 = ((D + D[i1])**2).sum(axis=1).argmin()

print(f'vac_index={vac_index} i1={i1} i2={i2} '
      f'distance={initial.get_distance(i1, i2, mic=True)}')

final = initial.copy()
final.positions[i1] = vac_pos

initial.calc = calc()
final.calc = calc()

qn = BFGS(initial)
qn.run(fmax=1e-3)
initial.write('initial.xyz')
make_geomopt(initial, 'initial')

qn = BFGS(final)
qn.run(fmax=1e-3)
final.write('final.xyz')
make_geomopt(final, 'final')

images = [initial]
for image in range(N_intermediate):
    image = initial.copy()
    image.calc = calc()
    images.append(image)
images.append(final)

neb = NEB(images, k=0.1)
neb.interpolate()

opt = BFGS(neb, trajectory='neb.traj')
opt.run(fmax=1e-2)

# Save the results
path = 'cuvac-ase-traj'
try:
    os.mkdir(path)
except FileExistsError:
    pass

for i, img in enumerate(images):
    img.write(os.path.join(path, 'cuvac-ase-{0}.cell'.format(i)))

# Energies?
E = [img.get_potential_energy() for img in images]
np.savetxt(os.path.join(path, 'cuvac-ase-E.dat'),
           np.array([np.arange(len(images)),
                     np.array(E)]).T)
