Code accompanying the journal article "Another source of mesh dependence in topology optimization".
It requires [Firedrake](https://github.com/firedrakeproject/firedrake), [pyMMAopt](https://github.com/LLNL/pyMMAopt) and [Gmsh](https://gmsh.info/).

Generate the compliance mesh with Gmsh:
```
gmsh -2 beam_amr.geo
```
and similarly for the other `.geo` files.
[Signac](https://signac.io/) is necessary to use the below commands and run all the jobs.
Each script can be run separately without Signac.
```
python3 init.py
python3 project.py run
```
To visualize the results, [ParaView](https://www.paraview.org/) and [matplotlib](https://matplotlib.org/) are required.
