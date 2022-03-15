from firedrake import *
from firedrake_adjoint import *

from pyMMAopt import MMASolver
import argparse
import glob
import re
import os
import signal
from copy import copy

parser = argparse.ArgumentParser(description="Compliance problem with MMA")
parser.add_argument(
    "--nref",
    action="store",
    dest="nref",
    type=int,
    help="Number of mesh refinements",
    default=0,
)
parser.add_argument(
    "--inner_product",
    action="store",
    dest="inner_product",
    type=str,
    help="Inner product, euclidean or L2",
    default="L2",
)
parser.add_argument(
    "--output_dir",
    action="store",
    dest="output_dir",
    type=str,
    help="Directory for all the output",
    default="./",
)
args, unknown = parser.parse_known_args()
nref = args.nref
inner_product = args.inner_product
output_dir = args.output_dir

assert inner_product == "L2" or inner_product == "euclidean"
print(f"inner product is: {inner_product}")

m = RectangleMesh(20 * (2 ** nref), 10 * (2 ** nref), 2.0, 1.0, quadrilateral=True)
layers = 10 * (2 ** nref)
mesh = ExtrudedMesh(m, layers, layer_height=1.0 / layers)

NEUMANN = 2
DIRICHLET = 1

V = VectorFunctionSpace(mesh, "CG", 1)
u, v = TrialFunction(V), TestFunction(V)
print(f"# DOFS: {V.dim()}")

# Elasticity parameters
E, nu = 1e2, 0.3
mu, lmbda = Constant(E / (2 * (1 + nu))), Constant(E * nu / ((1 + nu) * (1 - 2 * nu)))

# Helmholtz solver
RHO = FunctionSpace(mesh, "DG", 0)
init_design = 0.15
rho = Function(RHO)

with stop_annotating():
    rho.interpolate(Constant(init_design))

af, b = TrialFunction(RHO), TestFunction(RHO)

filter_radius = Constant(0.00005)
x, y, z = SpatialCoordinate(mesh)
with stop_annotating():
    x_ = interpolate(x, RHO)
    y_ = interpolate(y, RHO)
    z_ = interpolate(z, RHO)
Delta_h = sqrt(jump(x_) ** 2 + jump(y_) ** 2 + jump(z_) ** 2)
aH = filter_radius * jump(af) / Delta_h * jump(b) * (dS_h + dS_v) + af * b * dx
LH = rho * b * dx

rhof = Function(RHO)
solver_params = {
    "ksp_type": "preonly",
    "pc_type": "lu",
    "pc_factor_mat_solver_type": "mumps",
    "mat_mumps_icntl_14": 200,
    "mat_mumps_icntl_24": 1,
}
solver_parameters_amg = {
    "ksp_type": "cg",
    "ksp_converged_reason": None,
    "ksp_rtol": 1e-7,
    "pc_type": "hypre",
    "pc_hypre_type": "boomeramg",
    "pc_hypre_boomeramg_max_iter": 100,
    "pc_hypre_boomeramg_coarsen_type": "HMIS",
    "pc_hypre_boomeramg_agg_nl": 1,
    "pc_hypre_boomeramg_strong_threshold": 0.25,
    "pc_hypre_boomeramg_interp_type": "ext+i",
    "pc_hypre_boomeramg_P_max": 4,
    "pc_hypre_boomeramg_relax_type_all": "sequential-Gauss-Seidel",
    "pc_hypre_boomeramg_grid_sweeps_all": 1,
    "pc_hypre_boomeramg_truncfactor": 0.3,
    "pc_hypre_boomeramg_max_levels": 8,
}
solve(aH == LH, rhof, solver_parameters=solver_parameters_amg)
rhofControl = Control(rhof)

eps = Constant(1e-4)
p = Constant(3.0)


def simp(rho):
    return eps + (Constant(1.0) - eps) * rho ** p


def epsilon(v):
    return sym(nabla_grad(v))


def sigma(v):
    return 2.0 * mu * epsilon(v) + lmbda * tr(epsilon(v)) * Identity(3)


a = inner(simp(rhof) * sigma(u), epsilon(v)) * dx
load = as_vector(
    [
        0.0,
        -1.0 * conditional(And(y > 0.4, And(y < 0.6, And(z > 0.4, z < 0.6))), 1, 0),
        0.0,
    ]
)
L = inner(load, v) * ds_v(NEUMANN)

u_sol = Function(V)


bcs = DirichletBC(V, Constant((0.0, 0.0, 0.0)), DIRICHLET)
# create rigid body modes
x, y, z = SpatialCoordinate(mesh)
b0 = Function(V)
b1 = Function(V)
b2 = Function(V)
b3 = Function(V)
b4 = Function(V)
b5 = Function(V)
with stop_annotating():
    b0.interpolate(Constant([1, 0, 0]))
    b1.interpolate(Constant([0, 1, 0]))
    b2.interpolate(Constant([0, 0, 1]))
    b3.interpolate(as_vector([-y, x, 0.0]))
    b4.interpolate(as_vector([-x, 0.0, z]))
    b5.interpolate(as_vector([0.0, -z, y]))
nullmodes = VectorSpaceBasis([b0, b1, b2, b3, b4, b5])
# Make sure they're orthonormal.
nullmodes.orthonormalize()

solver_parameters_gamg = {
    "ksp_type": "cg",
    "ksp_max_it": 300,
    "pc_type": "gamg",
    "mat_type": "aij",
    "ksp_converged_reason": None,
}
solve(
    a == L,
    u_sol,
    bcs=bcs,
    solver_parameters=solver_parameters_gamg,
    near_nullspace=nullmodes,
)
c = Control(rho)
J = assemble(Constant(1e3) * inner(u_sol, load) * ds_v(NEUMANN))
Vol = assemble(rhof * dx)
VolControl = Control(Vol)


with stop_annotating():
    Vlimit = assemble(Constant(1.0) * dx(domain=mesh)) * init_design

rho_viz_f = Function(RHO, name="rho")
plot_file = f"{output_dir}/design_{nref}_{inner_product}.pvd"
print(plot_file)
controls_f = File(plot_file, target_continuity=H1, mode="a")


import itertools

global_counter = itertools.count()

plot_freq = 10 if args.nref in [0, 1, 2, 3] else 100

def deriv_cb(j, dj, rho):
    with stop_annotating():
        iter = next(global_counter)
        if iter % plot_freq == 0:
            rho_viz_f.assign(rhofControl.tape_value())
            controls_f.write(rho_viz_f)


Jhat = ReducedFunctional(J, c, derivative_cb_post=deriv_cb)
Volhat = ReducedFunctional(Vol, c)

class VolumeConstraint(InequalityConstraint):
    def __init__(self, Vhat, Vlimit, VolControl):
        self.Vhat = Vhat
        self.Vlimit = float(Vlimit)
        self.VolControl = VolControl
        self.tmpvec = Function(RHO)

    def function(self, m):
        # Compute the integral of the control over the domain
        integral = self.VolControl.tape_value()
        with stop_annotating():
            value = -integral / self.Vlimit + 1.0
        return [value]

    def jacobian(self, m):
        with stop_annotating():
            gradients = self.Vhat.derivative()
            with gradients.dat.vec as v:
                v.scale(-1.0 / self.Vlimit)
        return [gradients]

    def output_workspace(self):
        return [0.0]

    def length(self):
        """Return the number of components in the constraint vector (here, one)."""
        return 1


lb = 0.0
ub = 1.0
problem = MinimizationProblem(
    Jhat, bounds=(lb, ub), constraints=[VolumeConstraint(Volhat, Vlimit, VolControl)]
)

parameters_mma = {
    "move": 0.2,
    "maximum_iterations": 1000,
    "m": 1,
    "IP": 0,
    "tol": 1e-6,
    "accepted_tol": 1e-4,
    "norm": inner_product,
    "gcmma": True,
    "output_dir": output_dir
}

loop = 0

checkpoints = glob.glob(f"{output_dir}/checkpoint*")
checkpoints_sorted = sorted(checkpoints,
                    key=lambda L: list(map(int, re.findall(r'iter_(\d+)\.h5', L)))[0])

if checkpoints_sorted:
    last_file = checkpoints_sorted[-1]
    loop = int(re.findall(r'iter_(\d+)\.h5', last_file)[0])
    parameters_mma["restart_file"] = last_file

solver = MMASolver(problem, parameters=parameters_mma)

results = solver.solve(loop=loop)
final_design = f"{output_dir}/final_design.pvd"
print(plot_file)
controls_f = File(final_design)
with stop_annotating():
    rho_viz_f.assign(rhofControl.tape_value())
    controls_f.write(rho_viz_f)

with open(f"{output_dir}/finished_{nref}_{inner_product}.txt", "w") as f:
    f.write("Done")
