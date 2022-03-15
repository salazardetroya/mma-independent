from meshpy import triangle
import numpy as np
import firedrake as fd
from firedrake import sqrt, jump, dS, dx
try:
    from firedrake.cython import dmcommon
except ImportError:
    from firedrake.cython import dmplex as dmcommon


def pde_filter(RHO, rho, filter_radius=fd.Constant(0.2), solver_parameters=None):
    mesh = RHO.ufl_domain()
    x, y = fd.SpatialCoordinate(mesh)
    x_ = fd.interpolate(x, RHO)
    y_ = fd.interpolate(y, RHO)
    Delta_h = sqrt(jump(x_) ** 2 + jump(y_) ** 2)
    af, b = fd.TrialFunction(RHO), fd.TestFunction(RHO)
    aH = filter_radius * jump(af) / Delta_h * jump(b) * dS + af * b * dx
    LH = rho * b * dx

    rhof = fd.Function(RHO)
    fd.solve(aH == LH, rhof, solver_parameters=solver_parameters)

    return rhof


def refine_mesh(triangle_mesh, error_indicator, areas, min_area=0.01):
    shrink = 100
    exponent = 2
    areas_array = areas.dat.data_ro
    max_err = error_indicator.dat.data_ro.max()
    for index, err in enumerate(error_indicator.dat.data_ro[:]):
        shrink_factor = shrink * (err / max_err)**exponent
        area = areas_array[index]
        triangle_mesh.element_volumes[index] = max(area / (1 + shrink_factor),
                                                   min_area)

    triangle_mesh = triangle.refine(triangle_mesh)
    triangle_mesh.element_volumes.setup()
    return triangle_mesh


def coarsen_mesh(mesh_info, error_indicator, areas, expand_factor=10.0):
    err_50 = np.percentile(error_indicator.dat.data_ro, 50)
    max_area = areas.dat.data_ro.max()
    min_area = areas.dat.data_ro.min()
    print(f"Min area: {min_area}")

    def coarsen_func(vertices, area):
        vert_origin, vert_destination, vert_apex = vertices
        bary_x = (vert_origin.x + vert_destination.x + vert_apex.x) / 3
        bary_y = (vert_origin.y + vert_destination.y + vert_apex.y) / 3

        error_at_bary = error_indicator.at((bary_x, bary_y))
        shrink_factor = (error_at_bary / err_50)

        if shrink_factor > 1.0:
            return area >= min_area * 2.0 # Prevent too small areas
        else:
            return False

    triangle_mesh = triangle.build(mesh_info, refinement_func=coarsen_func, max_volume=max_area * expand_factor)
    triangle_mesh.element_volumes.setup()
    return triangle_mesh


def triangle_to_firedrake(triangle_mesh, comm=fd.COMM_WORLD):

    elements, points = triangle_mesh.elements, triangle_mesh.points
    plex = fd.mesh._from_cell_list(2, elements, points, comm)

    markers = {
        tuple(sorted((v1, v2))): triangle_mesh.facet_markers[index]
        for index, (v1, v2) in enumerate(triangle_mesh.facets)
    }
    plex.createLabel(dmcommon.FACE_SETS_LABEL)
    plex.markBoundaryFaces("boundary_faces")
    boundary_faces = plex.getStratumIS("boundary_faces", 1).getIndices()
    offset = plex.getDepthStratum(0)[0]
    for face in boundary_faces:
        vertices = tuple(sorted([v - offset for v in plex.getCone(face)]))
        marker = markers[vertices]
        plex.setLabelValue(dmcommon.FACE_SETS_LABEL, face, marker)

    return fd.Mesh(plex, reorder=False)
