from firedrake import *
from firedrake_adjoint import *
from pyadjoint.placeholder import Placeholder

from pyMMAopt import MMASolver
import argparse

parser = argparse.ArgumentParser(description="Stress constrained topology optimization problem with MMA")
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
parser.add_argument(
    "--quad_degree",
    action="store",
    dest="quad_degree",
    type=int,
    help="Integration degree for the ramp function",
    default=8,
)
parser.add_argument(
    "--p_norm_coeff",
    action="store",
    dest="p_norm_coeff",
    type=int,
    help="ramp function approximation coefficient",
    default=8,
)
parser.add_argument(
    "--penalty_param",
    action="store",
    dest="penalty_param",
    type=float,
    help="Penalty parameter for the constraint",
    default=10,
)
args = parser.parse_args()
nref = args.nref
uniform = args.uniform
inner_product = args.inner_product
output_dir = args.output_dir
p_norm_number = args.p_norm_coeff
quad_degree = args.quad_degree
penalty_param = args.penalty_param

from firedrake import PETSc

print = lambda x: PETSc.Sys.Print(x, comm=COMM_SELF)
assert inner_product == "L2" or inner_product == "euclidean"
print(f"inner product is: {inner_product}")

if uniform == 1:
    mesh = Mesh("./lbracket_uniform.msh")
else:
    mesh = Mesh("./lbracket_amr.msh")

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
with stop_annotating():
    rho = interpolate(Constant(0.5), RHO)
af, b = TrialFunction(RHO), TestFunction(RHO)

filter_radius = Constant(1.2)
x, y = SpatialCoordinate(mesh)
with stop_annotating():
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


def epsilon(v):
    return sym(nabla_grad(v))


def sigma(v):
    return 2.0 * mu * epsilon(v) + lmbda * tr(epsilon(v)) * Identity(2)


DIRICHLET = 1
NEUMANN = 2

eps = Constant(1e-5)
p = Constant(3.0)


def simp(rho):
    return (eps + rho - eps * rho) ** (3.0)


a = inner(simp(rhof) * sigma(u), epsilon(v)) * dx(degree=2)
load = Constant((0.0, -0.5))
L = inner(load, v) * ds(NEUMANN)

u_sol = Function(V)


bcs = DirichletBC(V, Constant((0.0, 0.0)), DIRICHLET)

solve(a == L, u_sol, bcs=bcs, solver_parameters=solver_params)
c = Control(rho)
u_sol_control = Control(u_sol)


def vonmises2D(stress_tensor):
    return sqrt(
        stress_tensor[0, 0] * stress_tensor[0, 0]
        + stress_tensor[1, 1] * stress_tensor[1, 1]
        - stress_tensor[0, 0] * stress_tensor[1, 1]
        + 3.0 * stress_tensor[1, 0] * stress_tensor[0, 1]
    )


def stress_penalization(a):
    return (eps + a - eps * a) ** (0.5)


sigma_y = 1.5
p_norm_coeff = Constant(10.0)
Placeholder(p_norm_coeff)


def ramp_function(vonmises):
    return (1.0 + (vonmises ** 2 / sigma_y ** 2) ** p_norm_coeff) ** (
        1.0 / p_norm_coeff
    ) - 1.0


penalty_param = Constant(penalty_param)
J = assemble(
    rhof * dx
    + penalty_param
    * ramp_function(stress_penalization(rhof) * vonmises2D(sigma(u_sol)))
    * dx(degree=quad_degree)
)


Vol = assemble(rho * dx(1))
VolControl = Control(Vol)


with stop_annotating():
    Vlimit = assemble(Constant(1.0) * dx(1, domain=mesh)) * 1e-8

rho_viz_f = Function(RHO, name="rho")
plot_file = f"{output_dir}/design_{uniform}_{inner_product}.pvd"
controls_f = File(plot_file)
vonmises_plot_file = File(
    f"{output_dir}/vonmises_stress_ref_{nref}_uniform_{uniform}_inner_product_{inner_product}.pvd"
)
vonmises2D_viz = Function(RHO)


import itertools

global_counter1 = itertools.count()


def deriv_cb(j, dj, rho):
    iter = next(global_counter1)
    if iter % 10 == 0:
        with stop_annotating():
            print(f"PLOTTING")
            rho_viz_f.assign(rhofControl.tape_value())
            controls_f.write(rho_viz_f)
            vonmises2D_viz.assign(
                project(
                    stress_penalization(rho_viz_f)
                    * vonmises2D(sigma(u_sol_control.tape_value())),
                    RHO,
                    form_compiler_parameters={"degree": 0},
                    solver_parameters={"ksp_type": "preonly", "pc_type": "jacobi"},
                )
            )
            vonmises_plot_file.write(vonmises2D_viz)


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
        print(f"Volume: {integral}, Vlimit: {self.Vlimit}")
        with stop_annotating():
            value = -integral + self.Vlimit
        return [value]

    def jacobian(self, m):
        with stop_annotating():
            gradients = self.Vhat.derivative()
            with gradients.dat.vec as v:
                v.scale(-1.0)
        return [gradients]

    def output_workspace(self):
        return [0.0]

    def length(self):
        """Return the number of components in the constraint vector (here, one)."""
        return 1


lb = 1e-5
ub = 1.0
ramp_p_arr = [p_norm_number, 20.0]
max_iter_arr = [300, 150]
for ramp_p_val, max_iter in zip(ramp_p_arr, max_iter_arr):
    with stop_annotating():
        p_norm_coeff.assign(ramp_p_val)
    problem = MinimizationProblem(
        Jhat,
        bounds=(lb, ub),
        constraints=[VolumeConstraint(Volhat, Vlimit, VolControl)],
    )

    parameters_mma = {
        "move": 0.2,
        "maximum_iterations": max_iter,
        "m": 1,
        "IP": 0,
        "tol": 1e-6,
        "accepted_tol": 1e-4,
        "norm": inner_product,
        "gcmma": True,
    }

    solver = MMASolver(problem, parameters=parameters_mma)

    rho_opt = solver.solve()
    with stop_annotating():
        rho.assign(rho_opt)

with open(f"{output_dir}/finished_{uniform}_{inner_product}.txt", "w") as f:
    f.write("Done")
