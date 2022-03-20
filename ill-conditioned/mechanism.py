from firedrake import *
from firedrake_adjoint import *
from penalization import ramp

from pyMMAopt import MMASolver
import argparse

parser = argparse.ArgumentParser(description="Compliance mmechanism problem with MMA")
parser.add_argument(
    "--nref",
    action="store",
    dest="nref",
    type=int,
    help="Number of mesh refinements",
    default=0,
)
parser.add_argument(
    "--uniform",
    action="store",
    dest="uniform",
    type=int,
    help="Use uniform mesh",
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
args = parser.parse_args()
nref = args.nref
uniform = args.uniform
inner_product = args.inner_product
output_dir = args.output_dir

assert inner_product == "L2" or inner_product == "euclidean"
assert uniform == 0 or uniform == 1
print(f"inner product is: {inner_product}")

if uniform == 1:
    mesh = Mesh("./mechanism_uniform.msh")
else:
    mesh = Mesh("./mechanism_amr.msh")

if nref > 0:
    mh = MeshHierarchy(mesh, nref)
    mesh = mh[-1]
elif nref < 0:
    raise RuntimeError("Non valid mesh argument")

V = VectorFunctionSpace(mesh, "CG", 1)
u, v = TrialFunction(V), TestFunction(V)
print(f"# DOFS: {V.dim()}")

# Elasticity parameters
E, nu = 1.0, 0.3
mu, lmbda = Constant(E / (2 * (1 + nu))), Constant(E * nu / ((1 + nu) * (1 - 2 * nu)))

# Helmholtz solver
RHO = FunctionSpace(mesh, "DG", 0)
RHOF = FunctionSpace(mesh, "CG", 1)
with stop_annotating():
    rho = interpolate(Constant(0.1), RHO)
af, b = TrialFunction(RHOF), TestFunction(RHOF)

filter_radius = Constant(0.8)
aH = filter_radius * inner(grad(af), grad(b)) * dx + af * b * dx
LH = rho * b * dx

rhof = Function(RHOF)
solver_params = {
    "ksp_type": "preonly",
    "pc_type": "lu",
    "pc_factor_mat_solver_type": "mumps",
    "mat_mumps_icntl_14": 200,
    "mat_mumps_icntl_24": 1,
}
solve(aH == LH, rhof, solver_parameters=solver_params)
rhofControl = Control(rhof)


def epsilon(v):
    return sym(nabla_grad(v))


def sigma(v):
    return 2.0 * mu * epsilon(v) + lmbda * tr(epsilon(v)) * Identity(2)


DIRICHLET = 1
INLET = 2
OUTLET = 3
ROLLERS = 4
k_in = Constant(1.0 / 3.0)
k_out = Constant(1.0e-3 / 3.0)
eps = Constant(1e-5)
ramp_rho = ramp(rhof, ramp_p=20.0, val_0=1e-5)

p = Constant(3.0)


def simp(rho):
    return eps + (Constant(1.0) - eps) * rho ** p


a = (
    inner(simp(rhof) * sigma(u), epsilon(v)) * dx
    + inner(k_in * u[0], v[0]) * ds(INLET)
    + inner(k_out * u[0], v[0]) * ds(OUTLET)
)


load = Constant((10.0, 0.0))
L = inner(load, v) * ds(INLET)

u_sol = Function(V)


bcs = [
    DirichletBC(V, Constant((0.0, 0.0)), DIRICHLET),
    DirichletBC(V.sub(1), Constant(0.0), ROLLERS),
]

solve(a == L, u_sol, bcs=bcs, solver_parameters=solver_params)
u_control = Control(u_sol)
c = Control(rho)
n = FacetNormal(mesh)
J = assemble(inner(u_sol, n) * ds(OUTLET))
Vol = assemble(rhof * dx)
VolControl = Control(Vol)


with stop_annotating():
    Vlimit = assemble(Constant(1.0) * dx(domain=mesh)) * 0.3

rho_viz_f = Function(RHOF, name="rho filtered")
plot_file = f"{output_dir}/design_{uniform}_{inner_product}.pvd"
controls_f = File(plot_file)

import itertools

global_counter1 = itertools.count()


def deriv_cb(j, dj, rho):
    iter = next(global_counter1)
    if iter % 10 == 0:
        with stop_annotating():
            rho_viz_f.assign(rhofControl.tape_value())
            controls_f.write(rho_viz_f)


Jhat = ReducedFunctional(J, c, derivative_cb_post=deriv_cb)
Volhat = ReducedFunctional(Vol, c)


class VolumeConstraint(InequalityConstraint):
    def __init__(self, Vhat, Vlimit, VolControl):
        self.Vhat = Vhat
        self.Vlimit = float(Vlimit)
        self.VolControl = VolControl

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


lb = 1e-5
ub = 1.0
problem = MinimizationProblem(
    Jhat, bounds=(lb, ub), constraints=[VolumeConstraint(Volhat, Vlimit, VolControl)]
)

parameters_mma = {
    "move": 0.2,
    "maximum_iterations": 200,
    "m": 1,
    "IP": 0,
    "tol": 1e-6,
    "accepted_tol": 1e-4,
    "norm": inner_product,
    "gcmma": True,
}
solver = MMASolver(problem, parameters=parameters_mma)

rho_opt = solver.solve()

with open(f"{output_dir}/finished_{uniform}_{inner_product}.txt", "w") as f:
    f.write("Done")
