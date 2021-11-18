#!/usr/bin/env python

import os
from ase import io
from ase.calculators.castep import Castep

# Number of images to test
n_images = [1, 2, 4, 8, 16]
n_cores_img = 8

base = io.read('CuVac.cell')
base.calc.merge_param('CuVac.param')
base.calc._seed = 'CuVac'

with open('SUBMIT.template.slurm') as f:
    slurm_template = f.read()

for im in n_images:
    folder = 'farm-{}'.format(im)
    base.calc.param.tssearch_max_path_points = im
    base.calc.param.num_farms = im
    base.calc._directory = folder
    base.calc.prepare_input_files(base, force_write=True)

    # Now write the slurm files
    with open(os.path.join(folder, 'SUBMIT.slurm'), 'w') as f:
        f.write(slurm_template.format(jobname=folder, seedname=base.calc._seed, 
                                      ntasks=im*n_cores_img))
