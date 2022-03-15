import argparse
from tkinter import W
import firedrake as fd
from firedrake import inner, dot, grad, div, dx, ds, dS, jump, sin, sqrt, pi, avg
import firedrake_adjoint as fda
from pyadjoint.placeholder import Placeholder
from pyMMAopt import MMASolver, ReducedInequality
from penalization import ramp
from meshpy import triangle
from amr_tools import pde_filter, triangle_to_firedrake, refine_mesh, coarsen_mesh
import numpy as np
from ufl.algebra import Abs

DOMAIN = 0
INMOUTH = 1
OUTMOUTHS = 2
INLET_WIDTH = 0.1
Y_INLET = 0.0

def print(x):
    return fd.PETSc.Sys.Print(x)

def create_triangle_mesh_info(boundary_conditions):
    points = [
        (-0.1, 0),
        (-0.1, -0.1),
        (0, -0.1),
        (0.0, -0.5),
        (1.0, -0.5),
        (1.0, -0.1),
        (1.1, -0.1),
        (1.1, 0.0),
    ]

    facets = [(i, (i + 1) % len(points)) for i in range(len(points))]
    WALLS = boundary_conditions["WALLS"]
    markers = [boundary_conditions["INLET"],
               WALLS,
               WALLS,
               WALLS,
               WALLS,
               WALLS,
               boundary_conditions["OUTLET"],
               WALLS,
              ]

    mesh_info = triangle.MeshInfo()
    mesh_info.set_points(points)
    mesh_info.set_facets(facets, facet_markers=markers)

    return mesh_info


def alpha(rho, Da, ramp_p):
    return fd.Constant(1.0 / Da) * ramp(rho, ramp_p=ramp_p)  # inv_ramp(rho, ramp_p=ramp_p)


def flow_problem(W, rhof, Re, Da, ramp_p, boundary_conditions, solver_parameters=None):

    v, q = fd.TestFunctions(W)

    up1 = fd.Function(W)
    u1, p1 = fd.split(up1)
    F1 = (1.0 / Re * inner(grad(u1), grad(v)) * dx
            + inner(dot(grad(u1), u1), v) * dx
            - p1 * div(v) * dx
            + div(u1) * q * dx
            + alpha(rhof, Da, ramp_p) * inner(u1, v) * dx
        )

    u_inflow = 1.0
    # Dirichelt boundary conditions
    _, y = fd.SpatialCoordinate(W.ufl_domain())
    inflow1 = fd.as_vector(
        [
            -u_inflow * sin(((y - (Y_INLET)) * pi) / INLET_WIDTH),
            0.0,
        ]
    )
    noslip = fd.Constant((0.0, 0.0))
    bcs1_1 = fd.DirichletBC(W.sub(0), noslip, boundary_conditions["WALLS"])
    bcs1_2 = fd.DirichletBC(W.sub(0), inflow1, boundary_conditions["INLET"])
    bcs1 = [bcs1_1, bcs1_2]

    problem = fd.NonlinearVariationalProblem(F1, up1, bcs=bcs1)
    solver = fd.NonlinearVariationalSolver(problem, solver_parameters=solver_parameters)
    solver.solve()

    return up1


def thermal_problem(T, rhof, u_flow, Pe, beta, boundary_conditions, solver_parameters=None):

    def generation(t):
        return rhof * beta * (fd.Constant(1.0) - t)

    tau = fd.TestFunction(T)
    t = fd.Function(T)
    mesh = T.ufl_domain()
    n = fd.FacetNormal(mesh)

    F = (inner(u_flow, grad(t)) * tau + fd.Constant(1.0 / Pe) * inner(grad(t), grad(tau))) * dx\
        - inner(fd.Constant(1.0 / Pe) * grad(t), n) * tau * ds(boundary_conditions["OUTLET"])\
        - generation(t) * tau * dx

    def L_T(t):
        return dot(u_flow, grad(t)) - fd.Constant(1.0 / Pe) * div(grad(t)) + rhof * beta * t

    def res(t):
        return dot(u_flow, grad(t)) - fd.Constant(1.0 / Pe) * div(grad(t)) - generation(t)

    R_U = res(t)
    beta_gls = 0.9
    h = fd.CellSize(mesh)
    tau_gls = beta_gls * (
        (4.0 * dot(u_flow, u_flow) / h ** 2)
         + 9.0 * (4.0 * fd.Constant(1.0 / Pe) / h ** 2) ** 2
         + (rhof * beta) ** 2
    ) ** (-0.5)
    degree = 4

    theta_U = L_T(tau)
    F_T = (F + tau_gls * inner(R_U, theta_U) * dx(degree=degree))

    bc1 = fd.DirichletBC(T, fd.Constant(0.0), boundary_conditions["INLET"])
    bcs = [bc1]
    problem_T = fd.NonlinearVariationalProblem(F_T, t, bcs=bcs)
    solver_T = fd.NonlinearVariationalSolver(
                    problem_T, solver_parameters=solver_parameters
                    )
    solver_T.solve()
    t.rename("Temperature")
    return t


def cost_functional(c, rhof, t, beta, deriv_cb=None):
    scale = 1e3
    J = fd.assemble(-rhof * fd.Constant(scale) * beta * (1 - t) * dx)
    return fda.ReducedFunctional(J, c, derivative_cb_post=deriv_cb)


def pressure_constraint(c, p, boundary_conditions, Pdrop):
    G = (fd.assemble(p * ds(boundary_conditions["INLET"]))
        - fd.assemble(p * ds(boundary_conditions["OUTLET"])))
    Gcontrol = fda.Control(G)
    Ghat = fda.ReducedFunctional(G, c)
    Glimit = Pdrop

    return Ghat, Gcontrol, Glimit


def navier_stokes_flow():  # sourcery no-metrics

    parser = argparse.ArgumentParser(description="")
    parser.add_argument(
        "--nu",
        dest="nu",
        type=float,
        action="store",
        default=1.0,
        help="Kinematic viscosity",
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
        "--output-dir",
        dest="output_dir",
        type=str,
        action="store",
        default="./",
        help="Directory where to save the result",
    )
    parser.add_argument(
        "--ref_strategy",
        action="store",
        dest="ref_strategy",
        type=str,
        help="Strategy A or B for the AMR",
        default="A",
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
        "--movlim",
        action="store",
        dest="movlim",
        type=float,
        help="MMA move limit",
        default=0.1,
    )
    parser.add_argument(
        "--output_dir",
        action="store",
        dest="output_dir",
        type=str,
        help="Directory for all the output",
        default="./",
    )

    opts = parser.parse_args()
    inner_product = opts.inner_product
    coarsen_refine = opts.coarsen_refine
    output_dir = opts.output_dir
    ref_strategy = opts.ref_strategy
    min_area = 6e-6
    max_area = 0.01

    solver_parameters_direct = {
        "snes_rtol": 1e-5,
        #"snes_atol": 1e-5,
        #"ksp_atol": 1e-5,
        #"ksp_rtol": 1e-5,
        "mat_type": "aij",
        "pc_type": "lu",
        "pc_factor_mat_solver_type": "mumps",
        #"ksp_converged_reason": None,
        #"ksp_type": "gmres",
        #"ksp_gmres_restart": 100,
        ##"pc_hypre_euclid_level": 3,
        #"ksp_pc_side": "right",
        #"pc_type": "hypre",
        #"pc_hypre_type": "euclid",
    }

    boundary_conditions = {
        "INLET": 1,
        "OUTLET": 2,
        "WALLS": 3
    }
    mesh_info = create_triangle_mesh_info(boundary_conditions)
    triangle_mesh = triangle.build(mesh_info, max_volume=0.00005) # 0.00005
    triangle_mesh.element_volumes.setup()

    mesh = triangle_to_firedrake(triangle_mesh)

    RHO = fd.FunctionSpace(mesh, "DG", 0)
    V = fd.VectorFunctionSpace(mesh, "CG", 2)
    Q = fd.FunctionSpace(mesh, "CG", 1)
    W = V * Q
    T = fd.FunctionSpace(mesh, 'CG', 1)
    print(f"Initial DOFs: {W.dim()}")
    ramp_p = fd.Constant(20.0)
    Re_val = 1.0
    Re = fd.Constant(Re_val)
    Pe = 1e4
    beta = fd.Constant(0.01)
    Pdrop = 70.0
    Da = 1e-6

    with fda.stop_annotating():
        rho = fd.interpolate(fd.Constant(0.0), RHO)
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

    import itertools
    global_counter1 = itertools.count()

    for ref_iter, max_iter in enumerate(max_iter_arr):
        rhof = pde_filter(RHO, rho, filter_radius=fd.Constant(1e-5),
                          solver_parameters=solver_parameters_direct)
        rhof_control = fda.Control(rhof)
        up1 = flow_problem(W, rhof, Re, Da, ramp_p, boundary_conditions, solver_parameters=solver_parameters_direct)
        u1, p1 = fd.split(up1)
        t = thermal_problem(T, rhof, u1, Pe, beta, boundary_conditions, solver_parameters=solver_parameters_direct)

        up_control = fda.Control(up1)
        t_control = fda.Control(t)

        c = fda.Control(rho)

        plot_file = f"{output_dir}/design_{ref_strategy}_{inner_product}_{coarsen_refine}_{ref_iter}.pvd"
        controls_f = fd.File(plot_file)
        vel_pvd = fd.File("velocity.pvd")
        t_pvd = fd.File("temperature.pvd")
        rho_viz_f = fd.Function(RHO, name="rho")
        t_viz = fd.Function(T, name="temperature")

        def deriv_cb(j, dj, rho):
            with fda.stop_annotating():
                iter = next(global_counter1)
                if iter % 20 == 0:
                    rho_viz_f.assign(rhof_control.tape_value())
                    #u_plot, p_plot = up_control.tape_value().split()
                    #u_plot.rename("Velocity")
                    #p_plot.rename("Pressure")
                    #t_viz.assign(t_control.tape_value())

                    #t_pvd.write(t_viz)
                    #vel_pvd.write(u_plot, p_plot)
                    controls_f.write(rho_viz_f)

        Jhat = cost_functional(c, rhof, t, beta, deriv_cb=deriv_cb)

        with fda.stop_annotating():
            n = fd.FacetNormal(u1.ufl_domain())
            inflow_flow = -1.0 * fd.assemble(inner(u1, n) * ds(boundary_conditions["INLET"]))

        Ghat, Gcontrol, Glimit = pressure_constraint(c, p1, boundary_conditions, Pdrop)

        # Bound constraints
        lb = 0.0
        ub = 1.0
        # Solve the optimisation problem with q = 0.01

        problem = fda.MinimizationProblem(
            Jhat,
            bounds=(lb, ub),
            constraints=[
                ReducedInequality(Ghat, Glimit, Gcontrol),
            ],
        )

        parameters_mma = {
            "move": opts.movlim,
            "maximum_iterations": max_iter,
            "m": 1,
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

        rhof_opt = rhof_control.tape_value()
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
                triangle_mesh = coarsen_mesh(mesh_info, error_indicator, areas, expand_factor=2.01)

            mesh = triangle_to_firedrake(triangle_mesh)

        V = fd.VectorFunctionSpace(mesh, "CG", 2)
        Q = fd.FunctionSpace(mesh, "CG", 1)
        W = V * Q
        print(f"# DOFS: {W.dim()}")
        T = fd.FunctionSpace(mesh, 'CG', 1)

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
    navier_stokes_flow()
