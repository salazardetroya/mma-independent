# init.py
import signac
import numpy as np
import itertools
from seurat import grid

import platform

plt = platform.system()

workspace = "./"
if plt == "Linux":
    workspace = "/p/lscratchh/miguel/topopt_gcmma/uniform/"
elif plt == "Darwin":
    workspace = "./workspace/"

project = signac.init_project("gcmma", workspace=workspace)


statepoint_grid = {
    'nref': [1, 2, 3],
    'inner_product' : ["L2", "euclidean"],
    'eps_approach' : ['sigma', 'lb']
}

for sp in grid(statepoint_grid):
    print('Initializing job', sp)
    project.open_job(sp).init()
