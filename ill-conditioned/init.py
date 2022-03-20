# init.py
import signac
import numpy as np
import itertools


def grid(gridspec):
    for values in itertools.product(*gridspec.values()):
        yield dict(zip(gridspec.keys(), values))

import platform
plt = platform.system()

workspace = "./workspace/"

project = signac.init_project('gcmma', workspace=workspace)

statepoint_grid = {
        'case' : ["compliance", "mechanism"],
        'quad_degree' : np.array([0.0]),
        'p_norm_coeff' : np.array([0.0]),
        'penalty_param' : np.array([0.0])
        }

for sp in grid(statepoint_grid):
    print('Initializing job', sp)
    project.open_job(sp).init()

statepoint_grid = {
        'case' : ["stress"],
        'quad_degree' : np.array([0, 4, 20]),
        'p_norm_coeff' : np.array([8]),
        'penalty_param' : np.array([10.0])
        }

for sp in grid(statepoint_grid):
    print('Initializing job', sp)
    project.open_job(sp).init()
