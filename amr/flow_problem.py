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
INLET_WIDTH = 0.2
Y_INLET = 0.0


def create_triangle_mesh_info(boundary_conditions):
    points = [
        (-0.2, 0),
        (0, 0),
        (1, 0),
        (1.2, 0),
        (1.2, 0.2),
        (1., 0.2),
        (1., 0.52),
        (1.2, 0.52),
        (1.2, 0.72),
        (1.0, 0.72),
        (1.0, 1.0),
        (1.2, 1.0),
        (1.2, 1.2),
        (1.0, 1.2),
        (0.0, 1.2),
        (0.0, 0.2),
        (-0.2, 0.2),
    ]

    facets = [(i, (i + 1) % len(points)) for i in range(len(points))]
    WALLS = boundary_conditions["WALLS"]
    markers = [WALLS,
               WALLS,
               WALLS,
               boundary_conditions["OUTLET1"],
               WALLS,
               WALLS,
               WALLS,
               boundary_conditions["OUTLET2"],
               WALLS,
               WALLS,
               WALLS,
               boundary_conditions["OUTLET3"],
               WALLS,
               WALLS,
               WALLS,
               WALLS,
               boundary_conditions["INLET1"],
              ]

    mesh_info = triangle.MeshInfo()
    mesh_info.set_points(points)
    mesh_info.set_facets(facets, facet_markers=markers)

    return mesh_info

def alpha(rho, Da, ramp_p):
    return fd.Constant(1.0 / Da) * ramp(rho, ramp_p=ramp_p)  # inv_ramp(rho, ramp_p=ramp_p)

def forward(W, rhof, Re, Da, ramp_p, boundary_conditions, solver_parameters=None):

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
            u_inflow * sin(((y - (Y_INLET)) * pi) / INLET_WIDTH),
            0.0,
        ]
    )
    noslip = fd.Constant((0.0, 0.0))
    bcs1_1 = fd.DirichletBC(W.sub(0), noslip, boundary_conditions["WALLS"])
    bcs1_2 = fd.DirichletBC(W.sub(0), inflow1, boundary_conditions["INLET1"])
    bcs1 = [bcs1_1, bcs1_2]

    problem = fd.NonlinearVariationalProblem(F1, up1, bcs=bcs1)
    solver = fd.NonlinearVariationalSolver(problem, solver_parameters=solver_parameters)
    solver.solve()

    return up1


def cost_functional(c, u1, rho, Re, Da, ramp_p, deriv_cb=None):
    J = fd.assemble(fd.Constant(1.0) / Re * inner(grad(u1), grad(u1)) * dx + alpha(rho, Da, ramp_p) * inner(u1, u1) * dx)
    return fda.ReducedFunctional(J, c, derivative_cb_post=deriv_cb)


def flow_constraint_functional(c, u1, boundary_conditions, inflow_flow):
    all_flow_outlets = []
    all_flow_outlets_ctrl = []
    n = fd.FacetNormal(u1.ufl_domain())
    for OUTLET in ["OUTLET1", "OUTLET2", "OUTLET3"]:
        outlet_flow = fd.assemble(inner(u1, n) * ds(boundary_conditions[OUTLET]))
        all_flow_outlets.append(outlet_flow)
        all_flow_outlets_ctrl.append(fda.Control(outlet_flow))

    deviation = (
        (
            1.0
            / 3.0
            * sum(
                (outlet_flow - inflow_flow / 3.0) ** 2
                for outlet_flow in all_flow_outlets
            )
        )
    ) ** 0.5

    Kcontrol = fda.Control(deviation)
    Khat = fda.ReducedFunctional(deviation, c)
    Klimit = inflow_flow * 0.001

    return Khat, Kcontrol, Klimit, all_flow_outlets_ctrl


def volume_constraint_functional(c, rhof):
    G = fd.assemble((fd.Constant(1.0) - rhof) * dx)
    Gcontrol = fda.Control(G)
    Ghat = fda.ReducedFunctional(G, c)
    with fda.stop_annotating():
        volume_domain = fd.assemble(fd.Constant(1.0) * dx(domain=rhof.ufl_domain()))
    Glimit = volume_domain / 3.5

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
    min_area = 2e-5
    max_area = 0.01

    solver_parameters_direct = {
        "snes_rtol": 1e-5,
        "ksp_type": "gmres",
        "mat_type": "aij",
        "pc_type": "lu",
        "pc_factor_mat_solver_type": "mumps",
    }

    boundary_conditions = {
        "INLET1": 1,
        "OUTLET1": 2,
        "OUTLET2": 3,
        "OUTLET3": 4,
        "WALLS": 5
    }
    mesh_info = create_triangle_mesh_info(boundary_conditions)
    triangle_mesh = triangle.build(mesh_info, max_volume=0.0002)
    triangle_mesh.element_volumes.setup()

    mesh = triangle_to_firedrake(triangle_mesh)

    RHO = fd.FunctionSpace(mesh, "DG", 0)
    V = fd.VectorFunctionSpace(mesh, "CG", 2)
    Q = fd.FunctionSpace(mesh, "CG", 1)
    W = V * Q
    print(f"Initial DOFs: {W.dim()}")
    ramp_p = fd.Constant(20.0)
    Re_val = 4.0
    Re = fd.Constant(Re_val)
    Da = 1e-5

    with fda.stop_annotating():
        rho = fd.interpolate(fd.Constant(0.5), RHO)
        xold1_func = fd.Function(RHO).interpolate(rho)
        xold2_func = fd.Function(RHO).interpolate(rho)
        low_func = fd.Function(RHO).interpolate(rho)
        upp_func = fd.Function(RHO).interpolate(rho)

    if coarsen_refine == "coarsen":
        if ref_strategy == "A":
            max_iter_arr = np.array([5, 10, 20, 100, 400])
        else:
            max_iter_arr = np.array([150, 200, 400])
    elif coarsen_refine == "refine":
        if ref_strategy == "A":
            max_iter_arr = np.array([5, 10, 20, 100, 400])
        else:
            max_iter_arr = np.array([100, 150, 400])
    loop = 1

    for ref_iter, max_iter in enumerate(max_iter_arr):
        rhof = pde_filter(RHO, rho, filter_radius=fd.Constant(1e-5),
                          solver_parameters=solver_parameters_direct)
        rhof_control = fda.Control(rhof)
        up1 = forward(W, rhof, Re, Da, ramp_p, boundary_conditions, solver_parameters=solver_parameters_direct)
        u1, p1 = fd.split(up1)

        up_control = fda.Control(up1)

        c = fda.Control(rho)

        plot_file = f"{output_dir}/design_{ref_strategy}_{inner_product}_{coarsen_refine}_{ref_iter}.pvd"
        controls_f = fd.File(plot_file)
        rho_viz_f = fd.Function(RHO, name="rho")

        import itertools

        global_counter1 = itertools.count()
        def deriv_cb(j, dj, rho):
            iter = next(global_counter1)
            if iter % 30 == 0:
                with fda.stop_annotating():
                    rho_viz_f.assign(rhof_control.tape_value())
                    u_plot, p_plot = up_control.tape_value().split()
                    u_plot.rename("Velocity")
                    p_plot.rename("Pressure")
                    controls_f.write(rho_viz_f)
                    print(f"{''.join([f'flow {i}: {val.tape_value()} ' for i, val in enumerate(all_flow_outlets_ctrl)])}")
                    print(f"inflow: {inflow_flow}")

        Jhat = cost_functional(c, u1, rhof, Re, Da, ramp_p, deriv_cb=deriv_cb)

        with fda.stop_annotating():
            n = fd.FacetNormal(u1.ufl_domain())
            inflow_flow = -1.0 * fd.assemble(inner(u1, n) * ds(boundary_conditions["INLET1"]))

        Khat, Kcontrol, Klimit, all_flow_outlets_ctrl = flow_constraint_functional(c, u1, boundary_conditions, inflow_flow)
        Ghat, Gcontrol, Glimit = volume_constraint_functional(c, rhof)

        # Bound constraints
        lb = 0.0
        ub = 1.0
        # Solve the optimisation problem with q = 0.01

        problem = fda.MinimizationProblem(
            Jhat,
            bounds=(lb, ub),
            constraints=[
                ReducedInequality(Ghat, Glimit, Gcontrol),
                ReducedInequality(Khat, Klimit, Kcontrol),
            ],
        )

        parameters_mma = {
            "move": 0.2,
            "maximum_iterations": max_iter,
            "m": 2,
            "norm": inner_product,
            "gcmma": False,
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

        with fda.stop_annotating():
            error_indicator = fd.assemble(Abs(jump(rhof_opt)) * avg(b) * dS)

            areas = fd.project(fd.CellVolume(mesh), RHO)
            fd.File("areas_errors.pvd").write(areas, error_indicator)

            if coarsen_refine == "refine":
                triangle_mesh = refine_mesh(triangle_mesh, error_indicator, areas,
                                            min_area=min_area)
            else:
                triangle_mesh = coarsen_mesh(mesh_info, error_indicator, areas)

            mesh = triangle_to_firedrake(triangle_mesh)

        V = fd.VectorFunctionSpace(mesh, "CG", 2)
        Q = fd.FunctionSpace(mesh, "CG", 1)
        W = V * Q
        print(f"# DOFS: {W.dim()}")

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
