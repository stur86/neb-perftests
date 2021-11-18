#!/usr/bin/env python

# Create input .cell file for the NEB

from ase import io
from ase.calculators.castep import Castep
from soprano.collection.generate import linspaceGen

a1 = io.read('cuvac-geomopt-initial-fast/initial-out.cell')
a2 = io.read('cuvac-geomopt-final-fast/final-out.cell')

calc = Castep()
a1.calc = calc

calc.cell.kpoint_mp_grid = [2,2,2]
calc.cell.fix_com = False
calc.cell.positions_abs_product = a2

lgen = linspaceGen(a1, a2, 3, True)
images = [im for im in lgen]

calc.cell.positions_abs_intermediate = images[1]

a1.write('cuvac-neb-fast/CuVac.cell')

