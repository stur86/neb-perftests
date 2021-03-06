#!/usr/bin/env python

import os
import json
import argparse as ap
import subprocess as sp

pwd = os.path.split(os.path.abspath(__file__))[0]
defaults = json.load(open(os.path.join(pwd, 'csub_default.json')))

parser = ap.ArgumentParser()
parser.add_argument('seedname', type=str, default=None,
                    help="""CASTEP seedname.""")
parser.add_argument('-n', type=int, default=defaults['n'],
                    help="""Number of cores""")
parser.add_argument('-t', type=str, default=defaults['t'],
                    help="""Runtime""")
parser.add_argument('-cp', type=str, default=defaults['cp'],
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
stdout, stderr = proc.communicate(script.encode('utf-8'))

if stderr.strip().decode('utf-8') == '':
    print(stdout)
else:
    raise RuntimeError(stderr)
