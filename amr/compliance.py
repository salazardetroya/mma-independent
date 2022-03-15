import firedrake as fd
from firedrake import PETSc
from firedrake import inner, jump, sqrt, ds, dx, dS, sym, nabla_grad, tr, Identity, avg
from ufl.algebra import Abs
import firedrake_adjoint as fda
from pyMMAopt import MMASolver
import argparse
import numpy as np
from meshpy import triangle
from amr_tools import refine_mesh, coarsen_mesh, triangle_to_firedrake, pde_filter

try:
    from firedrake.cython import dmcommon
except ImportError:
    from firedrake.cython import dmplex as dmcommon

def print(x):
    return PETSc.Sys.Print(x)

def create_mesh_info(boundary_conditions):
    Lx, Ly = 100, 40
    points = [
        (0, 0),
        (Lx, 0),
        (Lx, Ly * 4.0/10.0),
        (Lx, Ly * 6.0/10.0),
        (Lx, Ly),
        (0, Ly)
    ]

    facets = [(i, (i + 1) % len(points)) for i in range(len(points))]
    markers = [0, 0, boundary_conditions["NEUMANN"], 0, 0, boundary_conditions["DIRICHLET"]]

    mesh_info = triangle.MeshInfo()
    mesh_info.set_points(points)
    mesh_info.set_facets(facets, facet_markers=markers)

    return mesh_info


def simp(rho, eps, p):
    return eps + (fd.Constant(1.0) - eps) * rho ** p


def epsilon(v):
    return sym(nabla_grad(v))


def sigma(v, mu, lmbda):
    return 2.0 * mu * epsilon(v) + lmbda * tr(epsilon(v)) * Identity(2)


def forward(V, rhof, boundary_conditions, solver_parameters=None):

    # Elasticity parameters
    E, nu = 1e0, 0.3
    mu, lmbda = fd.Constant(E / (2 * (1 + nu))), fd.Constant(E * nu / ((1 + nu) * (1 - 2 * nu)))

    eps = fd.Constant(1e-5)
    p = fd.Constant(3.0)

    u, v = fd.TrialFunction(V), fd.TestFunction(V)
    a = inner(simp(rhof, eps, p) * sigma(u, mu, lmbda), epsilon(v)) * dx
    L = inner(boundary_conditions["load"], v) * ds(boundary_conditions["NEUMANN"])

    u_sol = fd.Function(V)

    bcs = fd.DirichletBC(V, fd.Constant((0.0, 0.0)), boundary_conditions["DIRICHLET"])
    fd.solve(a == L, u_sol, bcs=bcs, solver_parameters=solver_parameters)
    fd.File("u_sol.pvd").write(u_sol)

    return u_sol


def opti_problem(c, u_sol, rhof, boundary_conditions, plot_history=None):

    with fda.stop_annotating():
        Vlimit = fd.assemble(fd.Constant(1.0) * dx(domain=u_sol.ufl_domain())) * 0.3

    lb = 1e-5
    ub = 1.0

    J = fd.assemble(fd.Constant(1e-4) * inner(u_sol, boundary_conditions["load"])
                    * ds(boundary_conditions["NEUMANN"]))
    Jhat = fda.ReducedFunctional(J, c, derivative_cb_post=plot_history)
    Vol = fd.assemble(rhof * dx)
    Volhat = fda.ReducedFunctional(Vol, c)
    VolControl = fda.Control(Vol)
    RHO = rhof.function_space()
    return fda.MinimizationProblem(
        Jhat, bounds=(lb, ub), constraints=[VolumeConstraint(RHO,
                                            Volhat, Vlimit, VolControl)]
    )


class VolumeConstraint(fda.InequalityConstraint):

    def __init__(self, RHO, Vhat, Vlimit, VolControl):
        self.Vhat = Vhat
        self.Vlimit = float(Vlimit)
        self.VolControl = VolControl
        self.tmpvec = fd.Function(RHO)

    def function(self, m):
        # Compute the integral of the control over the domain
        integral = self.VolControl.tape_value()
        with fda.stop_annotating():
            value = -integral / self.Vlimit + 1.0
        return [value]

    def jacobian(self, m):
        with fda.stop_annotating():
            gradients = self.Vhat.derivative()
            with gradients.dat.vec as v:
                v.scale(-1.0 / self.Vlimit)
        return [gradients]

    def output_workspace(self):
        return [0.0]

    def length(self):
        """Return the number of components in the constraint vector (here, one)."""
        return 1


def compliance_optimization(ref_strategy, inner_product, min_area, coarsen_refine, max_area, output_dir):

    assert inner_product == "L2" or inner_product == "euclidean"
    assert ref_strategy == "A" or ref_strategy == "B"
    print(f"inner product is: {inner_product}")

    boundary_conditions = {
        "DIRICHLET": 3,
        "NEUMANN": 4,
        "load": fd.Constant((0.0, -1.0))
    }

    mesh_info = create_mesh_info(boundary_conditions)
    triangle_mesh = triangle.build(mesh_info, max_volume=0.2)
    triangle_mesh.element_volumes.setup()
    mesh = triangle_to_firedrake(triangle_mesh)

    V = fd.VectorFunctionSpace(mesh, "CG", 1)
    print(f"# DOFS: {V.dim()}")

    RHO = fd.FunctionSpace(mesh, "DG", 0)
    with fda.stop_annotating():
        rho = fd.interpolate(fd.Constant(0.1), RHO)
        xold1_func = fd.Function(RHO).interpolate(rho)
        xold2_func = fd.Function(RHO).interpolate(rho)
        low_func = fd.Function(RHO).interpolate(rho)
        upp_func = fd.Function(RHO).interpolate(rho)

    if coarsen_refine == "coarsen":
        if ref_strategy == "A":
            max_iter_arr = np.array([100, 150, 600])
        else:
            max_iter_arr = np.array([150, 200, 600])
    elif coarsen_refine == "refine":
        if ref_strategy == "A":
            max_iter_arr = np.array([10, 80, 600])
        else:
            max_iter_arr = np.array([100, 150, 600])
    loop = 1

    solver_parameters = {
        "ksp_type": "preonly",
        "pc_type": "lu",
        "pc_factor_mat_solver_type": "mumps",
        "mat_mumps_icntl_14": 200,
        "mat_mumps_icntl_24": 1,
    }

    import itertools
    global_counter1 = itertools.count()
    for ref_iter, max_iter in enumerate(max_iter_arr):
        rhof = pde_filter(RHO, rho, solver_parameters=solver_parameters)
        u_sol = forward(V, rhof, boundary_conditions, solver_parameters=solver_parameters)
        rhofControl = fda.Control(rhof)
        uControl = fda.Control(u_sol)

        rho_viz_f = fd.Function(RHO, name="rho")
        plot_file = f"{output_dir}/design_{ref_strategy}_{inner_product}_{coarsen_refine}_{ref_iter}.pvd"
        controls_f = fd.File(plot_file)

        def deriv_cb(j, dj, rho):
            iter = next(global_counter1)
            if iter % 20 == 0:
                with fda.stop_annotating():
                    rho_viz_f.assign(rhofControl.tape_value())
                    controls_f.write(rho_viz_f)

        c = fda.Control(rho)
        problem = opti_problem(c, u_sol, rhof, boundary_conditions, plot_history=deriv_cb)

        parameters_mma = {
            "move": 0.1,
            "maximum_iterations": max_iter,
            "m": 1,
            "IP": 0,
            "tol": 1e-6,
            "accepted_tol": 1e-4,
            "norm": inner_product,
            "gcmma": True,
        }
        solver = MMASolver(problem, parameters=parameters_mma)

        results = solver.solve(xold1_func=xold1_func, xold2_func=xold2_func,
                               low_func=low_func, upp_func=upp_func, loop=loop)

        if loop >= max_iter_arr[-1]:
            break

        rho_opt = results["control"]
        xold1_func = results["xold1"]
        xold2_func = results["xold2"]
        low_func = results["low"]
        upp_func = results["upp"]
        loop = results["loop"]

        rhof_opt = rhofControl.tape_value()
        b = fd.TestFunction(RHO)
        error_indicator = fd.Function(RHO)

        with fda.stop_annotating():
            error_indicator_arr = (fd.assemble(Abs(jump(rhof_opt)) * avg(b) * dS).dat.data /
                                fd.assemble(avg(b) * dS).dat.data)
            with error_indicator.dat.vec as err_vec:
                err_vec[:] = error_indicator_arr

            areas = fd.project(fd.CellVolume(mesh), RHO)

            if coarsen_refine == "refine":
                triangle_mesh = refine_mesh(triangle_mesh, error_indicator, areas,
                                            min_area=min_area)
            else:
                triangle_mesh = coarsen_mesh(mesh_info, error_indicator, areas)

            mesh = triangle_to_firedrake(triangle_mesh)

        V = fd.VectorFunctionSpace(mesh, "CG", 1)
        print(f"# DOFS: {V.dim()}")

        RHO = fd.FunctionSpace(mesh, "DG", 0)
        with fda.stop_annotating():
            rho = fd.project(rho_opt, RHO)
            xold1_func = fd.project(xold1_func, RHO)
            xold2_func = fd.project(xold2_func, RHO)
            low_func = fd.project(low_func, RHO)
            upp_func = fd.project(upp_func, RHO)
        fda.get_working_tape().clear_tape()

    with open(f"{output_dir}/finished_{ref_strategy}_{inner_product}_{coarsen_refine}.txt", "w") as f:
        f.write("Done")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Compliance problem with MMA")
    parser.add_argument(
        "--ref_strategy",
        action="store",
        dest="ref_strategy",
        type=str,
        help="Strategy A or B for the AMR",
        default="A",
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
        "--min_area",
        action="store",
        dest="min_area",
        type=float,
        help="Minimum area an area can have to limit an explosion of DOFS",
        default=0.01,
    )
    parser.add_argument(
        "--coarsen_refine",
        action="store",
        dest="coarsen_refine",
        type=str,
        help="Decide whether to coarsen or refine",
        default="refine",
    )
    parser.add_argument(
        "--max_area",
        action="store",
        dest="max_area",
        type=float,
        help="Minimum area an area can have to limit an explosion of DOFS",
        default=1.0,
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
    ref_strategy = args.ref_strategy
    inner_product = args.inner_product
    coarsen_refine = args.coarsen_refine
    min_area = args.min_area
    max_area = args.max_area
    output_dir = args.output_dir
    compliance_optimization(ref_strategy, inner_product, min_area, coarsen_refine, max_area, output_dir)
