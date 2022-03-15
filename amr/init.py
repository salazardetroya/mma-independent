# init.py
import signac
import numpy as np
from seurat import grid

import platform
plt = platform.system()

workspace = "./"
if plt == "Linux":
    workspace = "/p/lscratchh/miguel/topopt_gcmma/amr/"
elif plt == "Darwin":
    workspace = "./workspace/"

project = signac.init_project('gcmma', workspace=workspace)

statepoint_grid = {
        'case': ["compliance"],
        'coarsen_refine' : ["refine", "coarsen"],
        'ref_strategy' : ["A", "B"],
        'movlim' : np.array([0.1, 0.2])
        }

for sp in grid(statepoint_grid):
    print('Initializing job', sp)
    project.open_job(sp).init()

statepoint_grid = {
        'case': ["thermal_flow"],
        'coarsen_refine' : ["refine", "coarsen"],
        'ref_strategy' : ["A", "B"],
        'movlim' : np.array([0.1, 0.2])
        }

for sp in grid(statepoint_grid):
    print('Initializing job', sp)
    project.open_job(sp).init()
