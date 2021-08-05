import numpy as np
from ase.build import bulk
from ase.optimize import BFGS
from ase.neb import NEB
from ase.io import read
from ase.geometry.geometry import get_distances
from ase.calculators.emt import EMT
   
def calc():
    return EMT()

N_intermediate = 3
N_cell = 2
initial = bulk('Cu', cubic=True)
initial *= N_cell

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

qn = BFGS(final)
qn.run(fmax=1e-3)
final.write('final.xyz')
    
# initial = read('initial.xyz')
# final = read('final.xyz')

images = [initial]
for image in range(N_intermediate):
    image = initial.copy()
    image.rattle()
    image.calc = calc()
    images.append(image)
images.append(final)

neb = NEB(images, k=0.1)
neb.interpolate()

for i, image in enumerate(neb.images):
    image.write(f'neb_{i:03d}.cell')

opt = BFGS(neb, trajectory='neb.traj')
opt.run(fmax=1e-2)
