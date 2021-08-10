#!/usr/bin/env python

import argparse as ap
import subprocess as sp

parser = ap.ArgumentParser()
parser.add_argument('seedname', type=str, default=None,
                    help="""CASTEP seedname.""")
parser.add_argument('-n', type=int, default=16,
                    help="""Number of cores""")
parser.add_argument('-t', type=str, default='02:00:00',
                    help="""Runtime""")
parser.add_argument('-cp', type=str, default='castep.mpi',
                    help="""CASTEP path""")
args = parser.parse_args()

template = """#!/bin/bash
#SBATCH -p scarf
#SBATCH -n {n}
#SBATCH -t {t}
#SBATCH -o %J.log
#SBATCH -e %J.err

mpirun -n {n} {cp} {seedname}"""

script = template.format(seedname=args.seedname, t=args.t, n=args.n,
                         cp=args.cp)

# Submit
proc = sp.Popen(['sbatch'], stdin=sp.PIPE,
                stdout=sp.PIPE, stderr=sp.PIPE)
stdout, stderr = proc.communicate(script)

if stderr.strip() == '':
    print(stdout)
else:
    raise RuntimeError(stderr)
