import os
import numpy as np
from ase.build import bulk
from ase.optimize import BFGS
from ase.neb import NEB
from ase.io import read
from ase.geometry.geometry import get_distances
from ase.calculators.emt import EMT
from ase.calculators.castep import Castep


def calc():
    return EMT()


def make_geomopt(a, name, force=False):
    calc = Castep(directory='cuvac-geomopt-{0}'.format(name),
                  label=name)
    if os.path.isdir(calc._directory) and not force:
        print('WARNING: path ' + calc._directory + ' exists, skipping')
        return
    calc.param.task = 'geometryoptimisation'
    calc.param.calculate_stress = True
    calc.param.write_cell_structure = True
    calc.cell.cell_constraints = "1 1 1\n0 0 0"
    calc.cell.kpoint_mp_grid = "3 3 3"
    a = a.copy()
    a.calc = calc

    calc.prepare_input_files(a, True)


N_intermediate = 3
N_cell = 2
initial = bulk('Cu', cubic=True)
initial *= N_cell

initial.write('CuBulk.cell')

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
