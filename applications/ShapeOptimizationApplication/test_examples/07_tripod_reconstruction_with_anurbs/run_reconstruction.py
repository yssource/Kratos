# importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.ShapeOptimizationApplication as KratosShape
import KratosMultiphysics.MeshingApplication as KratosMeshingApp
import KratosMultiphysics.StructuralMechanicsApplication as KratosCSM

# check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

# Import ANurbs library
import ANurbs as an
import numpy as np
import numpy.linalg as la
import scipy as sp

# Additional imports
import time
import os

# ======================================================================================================================================
# Parameters
# ======================================================================================================================================

# Parameters
parameters = KratosMultiphysics.Parameters("""
{
    "inpute_parameters":
    {
        "cad_filename"                  : "tripod.iga",
        "fem_filename"                  : "tripod.mdpa",
        "fe_refinement_level"           : 0,
        "update_variable_name"          : "SHAPE_CHANGE"
    },
    "output_parameters":
    {
        "output_folder"            : "01_Results",
        "result_geometry_filename" : "reconstructed_geometry.iga"
    }
}""")

# ======================================================================================================================================
# Helper functions
# ======================================================================================================================================

def OutputFEData(model_part, fem_output_filename, nodal_variables):
    from gid_output import GiDOutput
    nodal_results=nodal_variables
    gauss_points_results=[]
    VolumeOutput = True
    GiDPostMode = "Binary"
    GiDWriteMeshFlag = False
    GiDWriteConditionsFlag = True
    GiDMultiFileFlag = "Single"

    gig_io = GiDOutput(fem_output_filename, VolumeOutput, GiDPostMode, GiDMultiFileFlag, GiDWriteMeshFlag, GiDWriteConditionsFlag)
    gig_io.initialize_results(model_part)
    gig_io.write_results(1, model_part, nodal_results, gauss_points_results)
    gig_io.finalize_results()

def is_point_on_line(point, line_a, line_b, tolerance):
    pa = line_a - point
    ab = line_b - line_a

    ab2 = ab @ ab

    if ab2 == 0:
        if pa @ pa < tolerance:
            return True
        else:
            return False

    cross = pa[0] * ab[1] - pa[1] * ab[0]

    if cross**2 / ab2 < tolerance**2:
        ap_dot_ab = -pa @ ab
        d = np.sign(ap_dot_ab) * ap_dot_ab**2 / ab2

        if d < -tolerance**2:
            return False

        if d > ab2 + tolerance**2 + 2 * ab2 * tolerance:
            return False

        return True

    return False

def contains(point, polygon):
    inside = False

    i = 0
    j = len(polygon) - 1

    while i < len(polygon):
        U0 = polygon[i]
        U1 = polygon[j]

        if is_point_on_line(point, U0, U1, 0.01):
            return True

        if point[1] < U1[1]:
            if U0[1] <= point[1]:
                lhs = (point[1] - U0[1])*(U1[0] - U0[0])
                rhs = (point[0] - U0[0])*(U1[1] - U0[1])
                if lhs > rhs:
                    inside = not inside
        elif point[1] < U0[1]:
            lhs = (point[1] - U0[1]) * (U1[0] - U0[0])
            rhs = (point[0] - U0[0]) * (U1[1] - U0[1])
            if (lhs < rhs):
               inside = not inside

        j = i
        i += 1

    return inside

# ======================================================================================================================================
# Preprocessing
# ======================================================================================================================================

print("\n\n========================================================================================================")
print("> Start reconstruction...")
print("========================================================================================================")

# Measure time
start_time = time.time()

print("\n> Starting preprocessing...")

# Create results folder
output_folder = parameters["output_parameters"]["output_folder"].GetString()
if not os.path.exists( output_folder ):
    os.makedirs( output_folder )

# Read FE data
fem_input_filename = parameters["inpute_parameters"]["fem_filename"].GetString()
update_variable_name = parameters["inpute_parameters"]["update_variable_name"].GetString()
update_variable = KratosMultiphysics.KratosGlobals.GetVariable(update_variable_name)

model = KratosMultiphysics.Model()
fe_model_part = model.CreateModelPart("reconstruction_part")
fe_model_part.AddNodalSolutionStepVariable(update_variable)
fe_model_part.AddNodalSolutionStepVariable(KratosShape.CONTROL_POINT_UPDATE)

model_part_io = KratosMultiphysics.ModelPartIO(fem_input_filename[:-5])
model_part_io.ReadModelPart(fe_model_part)

# Refine if specified
prop_id = 1
prop = fe_model_part.Properties[prop_id]
mat = KratosCSM.LinearElasticPlaneStress2DLaw()
prop.SetValue(KratosMultiphysics.CONSTITUTIVE_LAW, mat.Clone())

refinement_level = parameters["inpute_parameters"]["fe_refinement_level"].GetInt()
for refinement_level in range(0,refinement_level):

    number_of_avg_elems = 10
    number_of_avg_nodes = 10
    nodal_neighbour_search = KratosMultiphysics.FindNodalNeighboursProcess(fe_model_part, number_of_avg_elems, number_of_avg_nodes)
    neighbour_calculator = KratosMultiphysics.FindElementalNeighboursProcess(fe_model_part,2,10)
    nodal_neighbour_search.Execute()
    neighbour_calculator.Execute()

    for elem in fe_model_part.Elements:
        elem.SetValue(KratosMultiphysics.SPLIT_ELEMENT,True)

    refine_on_reference = False
    interpolate_internal_variables = True
    Refine = KratosMeshingApp.LocalRefineTriangleMesh(fe_model_part)
    Refine.LocalRefineMesh(refine_on_reference, interpolate_internal_variables)

# Output FE-data with projected points
output_folder = parameters["output_parameters"]["output_folder"].GetString()
fem_output_filename = "fe_model_used_for_reconstruction"
fem_output_filename_with_path = os.path.join(output_folder,fem_output_filename)
nodal_variables = [parameters["inpute_parameters"]["update_variable_name"].GetString()]
OutputFEData(fe_model_part, fem_output_filename_with_path, nodal_variables)

# Read CAD data
cad_filename = parameters["inpute_parameters"]["cad_filename"].GetString()
model = an.Model.open(cad_filename)

print("> Preprocessing finished.")

# ======================================================================================================================================
# Reconstruction
# ======================================================================================================================================

# Prepare working variables
print("> Starting to prepare working variables....")
start_time_subprocess = time.time()

num_control_points = 0
for face_itr, face_i in enumerate(model.of_type('BrepFace')):
    surface_geometry = face_i.surface_geometry_3d().geometry
    num_control_points += surface_geometry.NbPolesU * surface_geometry.NbPolesV

num_fe_nodes = fe_model_part.NumberOfNodes()
node_coords = np.zeros((num_fe_nodes,3))
node_updates = np.zeros((num_fe_nodes,3))
node_coords_updated = np.zeros((num_fe_nodes,3))
for node_itr, node in enumerate(fe_model_part.Nodes):

    node_coords_i = np.array([node.X0, node.Y0, node.Z0])
    node_coords[node_itr,:] = node_coords_i

    node_update_i = np.array(node.GetSolutionStepValue(update_variable))
    node_updates[node_itr,:] = node_update_i

    node_coords_updated[node_itr,:] = node_coords_i+node_update_i

print("> Finished preparing working variables in " ,round( time.time()-start_time_subprocess, 3 ), " s.")

# Identify points pairs
print("> Starting identification of point pairs....")
start_time_subprocess = time.time()

point_pairs = []
for node_i in fe_model_part.Nodes:
    point_pairs.append([])

sorted_point_pairs = []
for face_itr, face_i in enumerate(model.of_type('BrepFace')):
    sorted_point_pairs.append([])

tessellation = an.CurveTessellation2D()

for face_itr, face_i in enumerate(model.of_type('BrepFace')):

    print(face_itr)

    surface = face_i.surface_3d()
    projection = an.PointOnSurfaceProjection3D(surface)

    boundary_polygon = []
    for trim in face_i.trims():
        tessellation.Compute(trim.curve_2d(), 0.01)
        for i in range(tessellation.NbPoints):
            boundary_polygon.append(tessellation.Point(i))

    # Bounding box and scaling as the bounding box is based on a simple pointwise discretization of surface
    min_x, min_y, min_z, max_x, max_y, max_z = projection.BoundingBox
    bounding_box_tolerance = 1
    min_x -= bounding_box_tolerance
    min_y -= bounding_box_tolerance
    min_z -= bounding_box_tolerance
    max_x += bounding_box_tolerance
    max_y += bounding_box_tolerance
    max_z += bounding_box_tolerance

    for node_itr, node_i in enumerate(fe_model_part.Nodes):

        node_coords_i = node_coords[node_itr]
        node_update_i = node_updates[node_itr]

        # Points outside bounding box are not considered
        if node_coords_i[0] < min_x or max_x < node_coords_i[0]:
            continue
        if node_coords_i[1] < min_y or max_y < node_coords_i[1]:
            continue
        if node_coords_i[2] < min_z or max_z < node_coords_i[2]:
            continue

        projection.Compute(Point=node_coords_i)
        projected_point_i = projection.Point
        projected_point_uv = np.array([projection.ParameterU, projection.ParameterV])

        is_inside = contains(projected_point_uv, boundary_polygon)

        if is_inside:
            point_pairs[node_itr].append([projected_point_i, projected_point_uv, face_itr])
            sorted_point_pairs[face_itr].append([node_itr, projected_point_i, projected_point_uv])

# Check results
# i = 0
# for itr, entry in enumerate(point_pairs):
#     if len(entry) >= 1:
#         for (x, y, z), _, face_itr in entry:
#             model.add({
#                 'Key': f'Point3D<{i}>',
#                 'Type': 'Point3D',
#                 'Location': [x, y, z],
#                 'Layer': f'PointPairs{face_itr}',
#             })
#             i += 1
#     if len(entry) == 0:
#         print("> WARNING: Missing point pair for point: ", itr)

# for face_itr, entry in enumerate(sorted_point_pairs):
#         for node_itr, (x, y, z), _ in entry:
#             model.add({
#                 'Key': f'Point3D<{i}>',
#                 'Type': 'Point3D',
#                 'Location': [x, y, z],
#                 'Layer': f'PointPairs{face_itr}',
#             })
#             i += 1

# Some additional output
# # Output FE-data with projected points
# for node_itr, node_i in enumerate(fe_model_part.Nodes):
#     node_coords_i = node_coords[node_itr]
#     projected_point_i = point_pairs[node_itr][0]
#     distance = projected_point_i - node_coords_i
#     node_i.SetSolutionStepValue(CONTROL_POINT_UPDATE,[distance[0],distance[1],distance[2]])

# output_folder = parameters["output_parameters"]["output_folder"].GetString()
# fem_output_filename = "fe_model_used_for_reconstruction"
# fem_output_filename_with_path = os.path.join(output_folder,fem_output_filename)
# nodal_variables = [parameters["inpute_parameters"]["update_variable_name"].GetString(), "CONTROL_POINT_UPDATE"]
# OutputFEData(fe_model_part, fem_output_filename_with_path, nodal_variables)

print("> Finished identification of point pairs in " ,round( time.time()-start_time_subprocess, 3 ), " s.")

# Assembly
print("> Starting assembly....")
start_time_subprocess = time.time()

dof_ids ={}
dofs = []

def get_dof_id(dof):
    if dof not in dof_ids:
        dof_ids[dof] = len(dof_ids)
        dofs.append(dof)
    return dof_ids[dof]

def get_dof(index):
    return dofs[index]

eq_ids ={}
eqs = []

def get_eq_id(eq):
    if eq not in eq_ids:
        eq_ids[eq] = len(eq_ids)
        eqs.append(eq)
    return eq_ids[eq]

def get_eq(index):
    return eqs[index]

# Assign dofs and count number of relevant control points
for face_itr, face_i in enumerate(model.of_type('BrepFace')):

    surface_geometry = face_i.surface_geometry_3d().geometry
    shape_function = an.SurfaceShapeEvaluator(DegreeU=surface_geometry.DegreeU, DegreeV=surface_geometry.DegreeV, Order=0)

    for node_itr, xyz, (u,v) in sorted_point_pairs[face_itr]:
        shape_function.Compute(surface_geometry.KnotsU, surface_geometry.KnotsV, u, v)

        eq_id = get_eq_id((face_itr,node_itr))

        for j, (r, s) in enumerate(shape_function.NonzeroPoleIndices):
            dof_id = get_dof_id((face_itr,r,s))

num_equations = len(eqs)
num_relevant_control_points = len(dofs)

print("> Number of equations: ", num_equations)
print("> Number of active control points: ", num_relevant_control_points)

# Fill matrix and RHS
A = np.zeros((num_equations, num_relevant_control_points))
rhs_x = np.zeros(num_equations)
rhs_y = np.zeros(num_equations)
rhs_z = np.zeros(num_equations)

for face_itr, face_i in enumerate(model.of_type('BrepFace')):

    surface_geometry = face_i.surface_geometry_3d().geometry
    shape_function = an.SurfaceShapeEvaluator(DegreeU=surface_geometry.DegreeU, DegreeV=surface_geometry.DegreeV, Order=0)

    for node_itr, xyz, (u,v) in sorted_point_pairs[face_itr]:
        shape_function.Compute(surface_geometry.KnotsU, surface_geometry.KnotsV, u, v)

        eq_id = get_eq_id((face_itr,node_itr))

        rhs_x[eq_id] = node_updates[node_itr,0]
        rhs_y[eq_id] = node_updates[node_itr,1]
        rhs_z[eq_id] = node_updates[node_itr,2]

        for j, (r, s) in enumerate(shape_function.NonzeroPoleIndices):

            dof_id = get_dof_id((face_itr,r,s))
            A[eq_id, dof_id] = shape_function(0, j)

print("> Finished assembly in " ,round( time.time()-start_time_subprocess, 3 ), " s.")

# Solution
print("> Starting solution....")

# xs, r, _, _ = la.lstsq(A, rhs_x, rcond=None)
# ys, r, _, _ = la.lstsq(A, rhs_y, rcond=None)
# zs, r, _, _ = la.lstsq(A, rhs_z, rcond=None)

AT = np.transpose(A)

rhs_modified_x = np.dot(AT,rhs_x)
rhs_modified_y = np.dot(AT,rhs_y)
rhs_modified_z = np.dot(AT,rhs_z)

B = np.dot(AT, A)

B_diag = np.diag(B)

for i in range(B_diag.shape[0]):
    entry = B[i,i]

    # regularization
    B[i,i] += 0.001


    if entry == 0:
        raise RuntimeError
        print(i)
        print("Zero on main diagonal found")


xs = la.solve(B,rhs_modified_x)
ys = la.solve(B,rhs_modified_y)
zs = la.solve(B,rhs_modified_z)

pole_update = list(zip(xs, ys, zs))

print("> Finished solution in " ,round( time.time()-start_time_subprocess, 3 ), " s.")

# Upate model
for face_itr, face_i in enumerate(model.of_type('BrepFace')):
    surface_geometry = face_i.surface_geometry_3d().geometry

    for r in range(surface_geometry.NbPolesU):
        for s in range(surface_geometry.NbPolesV):
            dof_i = (face_itr,r,s)
            if dof_i in dof_ids.keys():
                dof_id = get_dof_id(dof_i)
                pole_coords = surface_geometry.Pole(r,s)
                new_pole_coords = pole_coords + pole_update[dof_id]
                surface_geometry.SetPole(r,s,new_pole_coords)

# # testing
# print(point_pairs[0])
# print("-------------")
# for i in range(num_control_points):
#     if A[0,i] != 0:
#         print("###")
#         print(i)
#         print(get_dof(i))
#         print(A[0,i])

# surface_geometry = model.of_type('BrepFace')[8].surface_geometry_3d().geometry
# shape_function = an.SurfaceShapeEvaluator(DegreeU=surface_geometry.DegreeU, DegreeV=surface_geometry.DegreeV, Order=0)
# u = point_pairs[0][0][1][0]
# v = point_pairs[0][0][1][1]
# shape_function.Compute(surface_geometry.KnotsU, surface_geometry.KnotsV, u, v)

# nonzeros = shape_function.NonzeroPoleIndices
# print(nonzeros)


# ======================================================================================================================================
# Output
# ======================================================================================================================================
output_folder = parameters["output_parameters"]["output_folder"].GetString()
output_filename = parameters["output_parameters"]["result_geometry_filename"].GetString()
output_filename_with_path = os.path.join(output_folder,output_filename)
model.save(output_filename_with_path)

print("\n========================================================================================================")
print("> Finished reconstruction in " ,round( time.time()-start_time, 3 ), " s.")
print("========================================================================================================")