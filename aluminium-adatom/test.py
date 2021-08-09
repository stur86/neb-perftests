import os
from ase.build import fcc100, add_adsorbate
from ase.constraints import FixAtoms
from ase.calculators.emt import EMT
from ase.optimize import QuasiNewton, BFGS
from ase.io import read
from ase.neb import NEB

if not os.path.exists('initial.traj') or not os.path.exists('final.traj'):
    # 2x2-Al(001) surface with 3 layers and an
    # Au atom adsorbed in a hollow site:
    slab = fcc100('Al', size=(2, 2, 3))
    add_adsorbate(slab, 'Au', 1.7, 'hollow')
    slab.center(axis=2, vacuum=4.0)

    slab.write('AlAu.cell')

    # Make sure the structure is correct:
    # view(slab)

    # Fix second and third layers:
    mask = [atom.tag > 1 for atom in slab]
    # print(mask)
    slab.set_constraint(FixAtoms(mask=mask))

    # Use EMT potential:
    slab.calc = EMT()

    # Initial state:
    qn = QuasiNewton(slab, trajectory='initial.traj')
    qn.run(fmax=0.05)

    # Final state:
    slab[-1].x += slab.get_cell()[0, 0] / 2
    qn = QuasiNewton(slab, trajectory='final.traj')
    qn.run(fmax=0.05)

initial = read('initial.traj')
final = read('final.traj')

constraint = FixAtoms(mask=[atom.tag > 1 for atom in initial])

images = [initial]
for i in range(3):
    image = initial.copy()
    image.calc = EMT()
    image.set_constraint(constraint)
    images.append(image)

images.append(final)

neb = NEB(images)
neb.interpolate()

for i, image in enumerate(neb.images):
    image.write(f'neb_{i:03d}.cell')

qn = BFGS(neb, trajectory='neb.traj')
qn.run(fmax=0.05)