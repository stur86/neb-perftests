#!/usr/bin/env python

# Create input .cell file for the NEB

from ase import io
from ase.calculators.castep import Castep
from soprano.collection.generate import linspaceGen

a1 = io.read('cuvac-geomopt-initial/initial-out.cell')
a2 = io.read('cuvac-geomopt-final/final-out.cell')

calc = Castep()
a1.calc = calc

calc.cell.kpoint_mp_grid = [4,4,4]
calc.cell.fix_com = False
calc.cell.positions_abs_product = a2

lgen = linspaceGen(a1, a2, 3, True)
images = [im for im in lgen]

# neb = NEB([a1, a1.copy(), a2])
# neb.interpolate()

calc.cell.positions_abs_intermediate = images[1]

for i, im in enumerate(images):
    io.write('CuVacNebIm{0}.cell'.format(i), im)

a1.write('cuvac-neb/CuVac.cell')