from firedrake import *
from firedrake_adjoint import *

from pyMMAopt import MMASolver
import argparse

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
print(f"mesh is : {uniform == 1}")

if uniform == 1:
    mesh = Mesh("./beam_uniform.msh")
else:
    mesh = Mesh("./beam_amr.msh")

if nref > 0:
    mh = MeshHierarchy(mesh, nref)
    mesh = mh[-1]
elif nref < 0:
    raise RuntimeError("Non valid mesh argument")

V = VectorFunctionSpace(mesh, "CG", 1)
u, v = TrialFunction(V), TestFunction(V)
print(f"# DOFS: {V.dim()}")

# Elasticity parameters
E, nu = 1e0, 0.3
mu, lmbda = Constant(E / (2 * (1 + nu))), Constant(E * nu / ((1 + nu) * (1 - 2 * nu)))

# Helmholtz solver
RHO = FunctionSpace(mesh, "DG", 0)
rho = interpolate(Constant(0.1), RHO)
af, b = TrialFunction(RHO), TestFunction(RHO)

filter_radius = Constant(0.2)
x, y = SpatialCoordinate(mesh)
x_ = interpolate(x, RHO)
y_ = interpolate(y, RHO)
Delta_h = sqrt(jump(x_) ** 2 + jump(y_) ** 2)
aH = filter_radius * jump(af) / Delta_h * jump(b) * dS + af * b * dx
LH = rho * b * dx

rhof = Function(RHO)
solver_params = {
    "ksp_type": "preonly",
    "pc_type": "lu",
    "pc_factor_mat_solver_type": "mumps",
    "mat_mumps_icntl_14": 200,
    "mat_mumps_icntl_24": 1,
}
solve(aH == LH, rhof, solver_parameters=solver_params)
rhofControl = Control(rhof)

eps = Constant(1e-5)
p = Constant(3.0)


def simp(rho):
    return eps + (Constant(1.0) - eps) * rho ** p


def epsilon(v):
    return sym(nabla_grad(v))


def sigma(v):
    return 2.0 * mu * epsilon(v) + lmbda * tr(epsilon(v)) * Identity(2)


DIRICHLET = 3
NEUMANN = 4

a = inner(simp(rhof) * sigma(u), epsilon(v)) * dx
load = Constant((0.0, -1.0))
L = inner(load, v) * ds(NEUMANN)

u_sol = Function(V)


bcs = DirichletBC(V, Constant((0.0, 0.0)), DIRICHLET)

solve(a == L, u_sol, bcs=bcs, solver_parameters=solver_params)
c = Control(rho)
J = assemble(Constant(1e-4) * inner(u_sol, load) * ds(NEUMANN))
Vol = assemble(rhof * dx)
VolControl = Control(Vol)


with stop_annotating():
    Vlimit = assemble(Constant(1.0) * dx(domain=mesh)) * 0.3

rho_viz_f = Function(RHO, name="rho")
plot_file = f"{output_dir}/design_{uniform}_{inner_product}.pvd"
print(plot_file)
controls_f = File(plot_file)


def deriv_cb(j, dj, rho):
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
