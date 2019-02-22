# ==============================================================================
#  KratosShapeOptimizationApplication
#
#  License:         BSD License
#                   license: ShapeOptimizationApplication/license.txt
#
#  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
#                   Oberbichler Thomas, https://github.com/oberbichler
#
# ==============================================================================

# importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.ShapeOptimizationApplication as KratosShape
import KratosMultiphysics.MappingApplication as KratosMapping

# Import helper classes
import cad_reconstruction_conditions as clib

# Additional imports
import ANurbs as an
import numpy as np
import numpy.linalg as la
import time, shutil

# ==============================================================================
class ConditionFactory:
    # --------------------------------------------------------------------------
    def __init__(self, fe_model_part, cad_model, parameters):
        self.fe_model_part = fe_model_part
        self.cad_model = cad_model
        self.parameters = parameters

        self.bounding_box_tolerance = self.parameters["drawing_parameters"]["patch_bounding_box_tolerance"].GetDouble()
        self.boundary_tessellation_tolerance = self.parameters["drawing_parameters"]["boundary_tessellation_tolerance"].GetDouble()

        self.boundary_polygons = None
        self.fe_point_parametrization = None

    # --------------------------------------------------------------------------
    def CreateConditions(self):
        conditions = {}
        for face_i in self.cad_model.GetByType('BrepFace'):
            conditions[face_i.Key()] = []

        if self.parameters["conditions"]["general"]["apply_integral_method"].GetBool():
            self.AddDistanceMinimizationWithIntegrationConditions(conditions)
        else:
            self.AddDistanceMinimizationConditions(conditions)
        if self.parameters["conditions"]["faces"]["curvature"]["apply_curvature_minimization"].GetBool() or \
           self.parameters["conditions"]["faces"]["mechanical"]["apply_KL_shell"].GetBool():
            self.AddFaceConditions(conditions)
        if self.parameters["conditions"]["faces"]["rigid"]["apply_rigid_conditions"].GetBool():
            self.AddRigidConditions(conditions)
        if self.parameters["conditions"]["edges"]["fe_based"]["apply_enforcement_conditions"].GetBool():
            self.AddEnforcementConditions(conditions)
        if self.parameters["conditions"]["edges"]["fe_based"]["apply_corner_enforcement_conditions"].GetBool():
            self.AddCornerEnforcementConditions(conditions)
        if self.parameters["conditions"]["edges"]["coupling"]["apply_coupling_conditions"].GetBool():
            self.AddCouplingConditions(conditions)
        if self.parameters["regularization"]["alpha"].GetDouble() != 0:
            self.AddAlphaRegularizationConditions(conditions)

        return conditions

    # --------------------------------------------------------------------------
    def AddDistanceMinimizationConditions(self, conditions):
        fe_point_parametric = self.GetFEPointParametrization()

        for entry in fe_point_parametric:

            node = entry["node"]
            node_coords = [node.X0, node.Y0, node.Z0]
            list_of_faces = entry["faces"]
            list_of_parameters = entry["parameters"]

            if len(list_of_faces) == 0:
                pass

            for face, (u,v) in zip(list_of_faces, list_of_parameters):

                surface_geometry = face.Data().Geometry()
                surface_geometry_data = surface_geometry.Data()

                shape_function = an.SurfaceShapeEvaluator(degreeU=surface_geometry_data.DegreeU(), degreeV=surface_geometry_data.DegreeV(), order=0)
                shape_function.Compute(surface_geometry_data.KnotsU(), surface_geometry_data.KnotsV(), u, v)
                nonzero_pole_indices = shape_function.NonzeroPoleIndices()

                shape_function_values = np.empty(shape_function.NbNonzeroPoles(), dtype=float)
                for i in range(shape_function.NbNonzeroPoles()):
                    shape_function_values[i] = shape_function(0,i)

                # One might introduce a possible penalty factor for nodes on boundary
                weight = 1.0
                # if is_on_boundary:
                # point_ptr = self.cad_model.Add(an.Point3D(location=node_coords))
                # point_ptr.Attributes().SetLayer('FEPointsInside')

                new_condition = clib.DistanceMinimizationCondition(node, surface_geometry, nonzero_pole_indices, shape_function_values, KratosShape.SHAPE_CHANGE, weight)
                conditions[face.Key()].append(new_condition)

        return conditions

    # --------------------------------------------------------------------------
    def AddDistanceMinimizationWithIntegrationConditions(self, conditions):
        drawing_tolerance = self.parameters["drawing_parameters"]["cad_drawing_tolerance"].GetDouble()
        total_area = 0

        for face_i in self.cad_model.GetByType('BrepFace'):

            print("> Processing face ",face_i.Key())

            surface_geometry = face_i.Data().Geometry()
            surface_geometry_data = face_i.Data().Geometry().Data()

            shape_function = an.SurfaceShapeEvaluator(degreeU=surface_geometry_data.DegreeU(), degreeV=surface_geometry_data.DegreeV(), order=0)

            list_of_points, list_of_parameters, list_of_integration_weights = self.__CreateIntegrationPointsForFace(face_i, drawing_tolerance)

            # Collect integration points in model part
            temp_model = KratosMultiphysics.Model()
            destination_mdpa = temp_model.CreateModelPart("temp_model_part")
            destination_mdpa.AddNodalSolutionStepVariable(KratosShape.NORMALIZED_SURFACE_NORMAL)
            destination_mdpa.AddNodalSolutionStepVariable(KratosShape.SHAPE_CHANGE)

            for itr, [x,y,z] in enumerate(list_of_points):
                destination_mdpa.CreateNewNode(itr, x, y, z)

            # Map information from fem to integration points using element based mapper
            mapper = KratosMapping.MapperFactory.CreateMapper( self.fe_model_part, destination_mdpa, self.parameters["conditions"]["general"]["mapping_cad_fem"].Clone() )
            mapper.Map( KratosShape.SHAPE_CHANGE, KratosShape.SHAPE_CHANGE )

            # Create conditions
            for itr, node_i in enumerate(destination_mdpa.Nodes):

                [u,v] = list_of_parameters[itr]
                weight = list_of_integration_weights[itr]

                total_area += weight

                shape_function.Compute(surface_geometry_data.KnotsU(), surface_geometry_data.KnotsV(), u, v)
                nonzero_pole_indices = shape_function.NonzeroPoleIndices()

                shape_function_values = np.empty(shape_function.NbNonzeroPoles(), dtype=float)
                for i in range(shape_function.NbNonzeroPoles()):
                    shape_function_values[i] = shape_function(0,i)

                new_condition = clib.DistanceMinimizationCondition(node_i, surface_geometry, nonzero_pole_indices, shape_function_values, KratosShape.SHAPE_CHANGE, weight)
                conditions[face_i.Key()].append(new_condition)

        print("> Total area of cad surface = ",total_area,"\n")

        return conditions

    # --------------------------------------------------------------------------
    def AddFaceConditions(self, conditions):
        drawing_tolerance = self.parameters["drawing_parameters"]["cad_drawing_tolerance"].GetDouble()
        apply_curvature_min = self.parameters["conditions"]["faces"]["curvature"]["apply_curvature_minimization"].GetBool()
        curvature_penalty_fac = self.parameters["conditions"]["faces"]["curvature"]["penalty_factor"].GetDouble()
        apply_kl_shell = self.parameters["conditions"]["faces"]["mechanical"]["apply_KL_shell"].GetBool()
        shell_penalty_fac = self.parameters["conditions"]["faces"]["mechanical"]["penalty_factor"].GetDouble()

        list_of_exclusive_faces = []
        for itr in range(self.parameters["conditions"]["faces"]["mechanical"]["exclusive_face_list"].size()):
            list_of_exclusive_faces.append(self.parameters["conditions"]["faces"]["mechanical"]["exclusive_face_list"][itr].GetString())

        for face_i in self.cad_model.GetByType('BrepFace'):

            print("> Processing face ",face_i.Key())

            # Skip faces if exclusive face list is specified (if list is empty, use all faces)
            if len(list_of_exclusive_faces) > 0 and face_i.Key() not in list_of_exclusive_faces:
                continue

            surface_geometry = face_i.Data().Geometry()
            surface_geometry_data = face_i.Data().Geometry().Data()

            shape_function = an.SurfaceShapeEvaluator(degreeU=surface_geometry_data.DegreeU(), degreeV=surface_geometry_data.DegreeV(), order=2)

            list_of_points, list_of_parameters, list_of_integration_weights = self.__CreateIntegrationPointsForFace(face_i, drawing_tolerance)

            # Create conditions
            for itr in range(len(list_of_points)):

                [u,v] = list_of_parameters[itr]
                weight = list_of_integration_weights[itr]

                shape_function.Compute(surface_geometry_data.KnotsU(), surface_geometry_data.KnotsV(), u, v)
                nonzero_pole_indices = shape_function.NonzeroPoleIndices()

                shape_function_values = np.empty(shape_function.NbNonzeroPoles(), dtype=float)
                shape_function_derivatives_u = np.empty(shape_function.NbNonzeroPoles(), dtype=float)
                shape_function_derivatives_v = np.empty(shape_function.NbNonzeroPoles(), dtype=float)
                shape_function_derivatives_uu = np.empty(shape_function.NbNonzeroPoles(), dtype=float)
                shape_function_derivatives_uv = np.empty(shape_function.NbNonzeroPoles(), dtype=float)
                shape_function_derivatives_vv = np.empty(shape_function.NbNonzeroPoles(), dtype=float)
                for i in range(shape_function.NbNonzeroPoles()):
                    shape_function_values[i] = shape_function(0,i)
                    shape_function_derivatives_u[i] = shape_function(1,i)
                    shape_function_derivatives_v[i] = shape_function(2,i)
                    shape_function_derivatives_uu[i] = shape_function(3,i)
                    shape_function_derivatives_uv[i] = shape_function(4,i)
                    shape_function_derivatives_vv[i] = shape_function(5,i)

                if apply_curvature_min:
                    weight = curvature_penalty_fac * weight
                    new_condition = clib.CurvatureMinimizationConditionWithAD(surface_geometry, nonzero_pole_indices, shape_function_derivatives_u, shape_function_derivatives_v, shape_function_derivatives_uu, shape_function_derivatives_uv, shape_function_derivatives_vv, weight)
                    conditions[face_i.Key()].append(new_condition)

                if apply_kl_shell:
                    weight = shell_penalty_fac * weight
                    new_condition = clib.KLShellConditionWithAD(surface_geometry, nonzero_pole_indices, shape_function_derivatives_u, shape_function_derivatives_v, shape_function_derivatives_uu, shape_function_derivatives_uv, shape_function_derivatives_vv, weight)
                    conditions[face_i.Key()].append(new_condition)

        return conditions

    # --------------------------------------------------------------------------
    def AddRigidConditions(self, conditions):
        list_of_exclusive_faces = []
        for itr in range(self.parameters["conditions"]["faces"]["rigid"]["exclusive_face_list"].size()):
            list_of_exclusive_faces.append(self.parameters["conditions"]["faces"]["rigid"]["exclusive_face_list"][itr].GetString())

        penalty_fac = self.parameters["conditions"]["faces"]["rigid"]["penalty_factor"].GetDouble()

        boundary_polygons = self.GetBoundaryPolygons()

        relevant_fe_points = []
        relevant_fe_points_displaced = []
        relevant_cad_uvs = []

        for face_i in self.cad_model.GetByType('BrepFace'):

            print("> Processing face ",face_i.Key())

            # Skip faces if exclusive face list is specified (if list is empty, use all faces)
            if face_i.Key() not in list_of_exclusive_faces:
                continue

            surface_geometry = face_i.Data().Geometry()
            surface_geometry_data = face_i.Data().Geometry().Data()

            shape_function = an.SurfaceShapeEvaluator(degreeU=surface_geometry_data.DegreeU(), degreeV=surface_geometry_data.DegreeV(), order=0)

            surface = an.Surface3D(face_i.Data().Geometry())
            projection = an.PointOnSurfaceProjection3D(surface)

            # Bounding box and scaling as the bounding box is based on a simple pointwise discretization of surface
            min_x, min_y, min_z, max_x, max_y, max_z = projection.BoundingBox()
            min_x -= self.bounding_box_tolerance
            min_y -= self.bounding_box_tolerance
            min_z -= self.bounding_box_tolerance
            max_x += self.bounding_box_tolerance
            max_y += self.bounding_box_tolerance
            max_z += self.bounding_box_tolerance

            for node_itr, node_i in enumerate(self.fe_model_part.Nodes):

                node_coords_i = np.array([node_i.X0, node_i.Y0, node_i.Z0])

                # Points outside bounding box are not considered
                if node_coords_i[0] < min_x or max_x < node_coords_i[0]:
                    continue
                if node_coords_i[1] < min_y or max_y < node_coords_i[1]:
                    continue
                if node_coords_i[2] < min_z or max_z < node_coords_i[2]:
                    continue

                projection.Compute(point=node_coords_i)
                projected_point_uv = np.array([projection.ParameterU(), projection.ParameterV()])

                is_inside, is_on_boundary = self.Contains(projected_point_uv, boundary_polygons[face_i.Key()], self.boundary_tessellation_tolerance*1.1)
                if is_inside:
                    relevant_fe_points.append(node_coords_i)
                    relevant_fe_points_displaced.append( node_coords_i + np.array(node_i.GetSolutionStepValue(KratosShape.SHAPE_CHANGE)) )
                    relevant_cad_uvs.append(projected_point_uv)

            # Analyze rigid body movement
            import numpy.matlib

            num_nodes = len(relevant_fe_points)

            A = np.array(relevant_fe_points)
            B = relevant_fe_points_displaced

            centroid_A = np.mean(A, axis=0)
            centroid_B = np.mean(B, axis=0)

            H = np.transpose( (A - np.matlib.repmat(centroid_A, num_nodes, 1)) ) @ (B - np.matlib.repmat(centroid_B, num_nodes, 1))

            # Note that the definition of V is different to Matlab
            # Matlab: X = U*S*V'
            # Numpy: X = U*S*V
            U, S, V = np.linalg.svd(H)

            # Rotation matrix
            R = np.transpose(V) @ np.transpose(U)

            if np.linalg.det(R)<0:
                print("> Reflection occured and was corrected!")
                R[:,2] = -R[:,2]

            # Translation vector
            t = - np.inner(R, centroid_A) + np.transpose(centroid_B)

            # Rigid motion
            rigididly_displaced_points = A @ np.transpose(R) + np.matlib.repmat(t, num_nodes, 1)

            for node_coords_disp, rigididly_displaced_point_coords, projected_point_uv in zip(relevant_fe_points_displaced, rigididly_displaced_points,relevant_cad_uvs):
                # self.cad_model.Add(an.Point3D(location=node_coords_disp))

                u = projected_point_uv[0]
                v = projected_point_uv[1]
                shape_function.Compute(surface_geometry_data.KnotsU(), surface_geometry_data.KnotsV(), u, v)
                nonzero_pole_indices = shape_function.NonzeroPoleIndices()

                shape_function_values = np.empty(shape_function.NbNonzeroPoles(), dtype=float)
                for i in range(shape_function.NbNonzeroPoles()):
                    shape_function_values[i] = shape_function(0,i)

                new_condition = clib.PositionEnforcementCondition(rigididly_displaced_point_coords, surface_geometry, nonzero_pole_indices, shape_function_values, penalty_fac)
                conditions[face_i.Key()].append(new_condition)

        return conditions

    # --------------------------------------------------------------------------
    def AddEnforcementConditions(self, conditions):
        drawing_tolerance = self.parameters["drawing_parameters"]["cad_drawing_tolerance"].GetDouble()
        min_span_length = self.parameters["drawing_parameters"]["min_span_length"].GetDouble()
        penalty_factor_tangent_enforcement = self.parameters["conditions"]["edges"]["fe_based"]["penalty_factor_tangent_enforcement"].GetDouble()
        penalty_factor_position_enforcement = self.parameters["conditions"]["edges"]["fe_based"]["penalty_factor_position_enforcement"].GetDouble()

        # Create corresponding points
        for edge_itr, edge_i in enumerate(self.cad_model.GetByType('BrepEdge')):

            print("> Processing edge ",edge_itr)

            adjacent_faces = edge_i.Data().Faces()
            num_adjacent_faces = len(adjacent_faces)

            if num_adjacent_faces == 1:
                # Skip non coupling-edges
                pass
            elif num_adjacent_faces == 2:

                face_a = adjacent_faces[0]
                face_b = adjacent_faces[1]

                surface_geometry_a = face_a.Data().Geometry()
                surface_geometry_b = face_b.Data().Geometry()
                surface_geometry_data_a = face_a.Data().Geometry().Data()
                surface_geometry_data_b = face_b.Data().Geometry().Data()

                list_of_points, list_of_parameters_a, list_of_parameters_b, _, _, list_of_integration_weights = self.__CreateIntegrationPointsForEdge(edge_i, drawing_tolerance, min_span_length)

                # Collect integration points in model part
                temp_model = KratosMultiphysics.Model()
                destination_mdpa = temp_model.CreateModelPart("temp_model_part")
                destination_mdpa.AddNodalSolutionStepVariable(KratosShape.NORMALIZED_SURFACE_NORMAL)
                destination_mdpa.AddNodalSolutionStepVariable(KratosShape.SHAPE_CHANGE)

                for itr, [x,y,z] in enumerate(list_of_points):
                    # point_ptr = self.cad_model.Add(an.Point3D(location=[x, y, z]))
                    # point_ptr.Attributes().SetLayer('CouplingPoints')
                    destination_mdpa.CreateNewNode(itr, x, y, z)

                # Map information from fem to integration points using element based mapper
                mapper = KratosMapping.MapperFactory.CreateMapper( self.fe_model_part, destination_mdpa, self.parameters["conditions"]["general"]["mapping_cad_fem"].Clone() )
                if penalty_factor_position_enforcement > 0:
                    mapper.Map( KratosShape.SHAPE_CHANGE, KratosShape.SHAPE_CHANGE )
                if penalty_factor_tangent_enforcement > 0:
                    mapper.Map( KratosShape.NORMALIZED_SURFACE_NORMAL, KratosShape.NORMALIZED_SURFACE_NORMAL )

                # Create conditions
                for itr, node in enumerate(destination_mdpa.Nodes):

                    integration_weight = list_of_integration_weights[itr]

                    shape_function_a = an.SurfaceShapeEvaluator(degreeU=surface_geometry_data_a.DegreeU(), degreeV=surface_geometry_data_a.DegreeV(), order=1)
                    shape_function_b = an.SurfaceShapeEvaluator(degreeU=surface_geometry_data_b.DegreeU(), degreeV=surface_geometry_data_b.DegreeV(), order=1)

                    # Create conditions to enforce t1 and t2 on both face a and face b
                    u_a = list_of_parameters_a[itr][0]
                    v_a = list_of_parameters_a[itr][1]
                    shape_function_a.Compute(surface_geometry_data_a.KnotsU(), surface_geometry_data_a.KnotsV(), u_a, v_a)
                    nonzero_pole_indices_a = shape_function_a.NonzeroPoleIndices()

                    shape_function_values_a = np.empty(shape_function_a.NbNonzeroPoles(), dtype=float)
                    shape_function_derivatives_u_a = np.empty(shape_function_a.NbNonzeroPoles(), dtype=float)
                    shape_function_derivatives_v_a = np.empty(shape_function_a.NbNonzeroPoles(), dtype=float)
                    for i in range(shape_function_a.NbNonzeroPoles()):
                        shape_function_values_a[i] = shape_function_a(0,i)
                        shape_function_derivatives_u_a[i] = shape_function_a(1,i)
                        shape_function_derivatives_v_a[i] = shape_function_a(2,i)

                    u_b = list_of_parameters_b[itr][0]
                    v_b = list_of_parameters_b[itr][1]
                    shape_function_b.Compute(surface_geometry_data_b.KnotsU(), surface_geometry_data_b.KnotsV(), u_b, v_b)
                    nonzero_pole_indices_b = shape_function_b.NonzeroPoleIndices()

                    shape_function_values_b = np.empty(shape_function_b.NbNonzeroPoles(), dtype=float)
                    shape_function_derivatives_u_b = np.empty(shape_function_b.NbNonzeroPoles(), dtype=float)
                    shape_function_derivatives_v_b = np.empty(shape_function_b.NbNonzeroPoles(), dtype=float)
                    for i in range(shape_function_b.NbNonzeroPoles()):
                        shape_function_values_b[i] = shape_function_b(0,i)
                        shape_function_derivatives_u_b[i] = shape_function_b(1,i)
                        shape_function_derivatives_v_b[i] = shape_function_b(2,i)

                    # Positions enforcement
                    if penalty_factor_position_enforcement > 0:
                        target_displacement = node.GetSolutionStepValue(KratosShape.SHAPE_CHANGE)
                        target_position = np.array([node.X+target_displacement[0], node.Y+target_displacement[1], node.Z+target_displacement[2]])
                        weight = penalty_factor_position_enforcement * integration_weight

                        new_condition_a = clib.PositionEnforcementCondition(target_position, surface_geometry_a, nonzero_pole_indices_a, shape_function_values_a, weight)
                        conditions[face_a.Key()].append(new_condition_a)

                        new_condition_b = clib.PositionEnforcementCondition(target_position, surface_geometry_b, nonzero_pole_indices_b, shape_function_values_b, weight)
                        conditions[face_b.Key()].append(new_condition_b)

                    # Tangents enforcement
                    if penalty_factor_tangent_enforcement > 0:
                        target_normal = node.GetSolutionStepValue(KratosShape.NORMALIZED_SURFACE_NORMAL)
                        weight = penalty_factor_tangent_enforcement * integration_weight

                        new_condition_a = clib.TangentEnforcementCondition(target_normal, surface_geometry_a, nonzero_pole_indices_a, shape_function_derivatives_u_a, shape_function_derivatives_v_a, weight)
                        conditions[face_a.Key()].append(new_condition_a)

                        new_condition_b = clib.TangentEnforcementCondition(target_normal, surface_geometry_b, nonzero_pole_indices_b, shape_function_derivatives_u_b, shape_function_derivatives_v_b, weight)
                        conditions[face_b.Key()].append(new_condition_b)
            else:
                raise RuntimeError("Max number of adjacent has to be 2!!")

    # --------------------------------------------------------------------------
    def AddCornerEnforcementConditions(self, conditions):
        # This conditions assumes an integration weight of 1 and other than that uses the penalty factors from the tangent and position enforcement
        penalty_factor = self.parameters["conditions"]["edges"]["fe_based"]["penalty_factor_corner_enforcement"].GetDouble()

        corner_points = []
        corner_point_faces = []
        corner_point_parameters = []
        considered_point_ids = []

        # Identifying corner points on each edge on each face and corresponding geometry information (avoiding duplications)
        for edge_itr, edge_i in enumerate(self.cad_model.GetByType('BrepEdge')):

            print("> Processing edge ",edge_itr)

            adjacent_faces = edge_i.Data().Faces()
            num_adjacent_faces = len(adjacent_faces)

            if num_adjacent_faces == 1:
                # Skip non coupling-edges
                pass
            elif num_adjacent_faces == 2:
                face_a = adjacent_faces[0]
                face_b = adjacent_faces[1]

                # Collect corner points
                corner_points_a, corner_points_b, parameters_a, parameters_b = self.__IdentifyStartAndEndPoint(edge_i)

                coordinate_tolerance = 10

                # Rounded values for comparison
                start_point_a = np.around(corner_points_a[0], decimals=coordinate_tolerance).tolist()
                end_point_a = np.around(corner_points_a[1], decimals=coordinate_tolerance).tolist()

                start_point_b = np.around(corner_points_b[0], decimals=coordinate_tolerance).tolist()
                end_point_b = np.around(corner_points_b[1], decimals=coordinate_tolerance).tolist()

                # Start point on face a
                considered_point_1_id = (start_point_a[0], start_point_a[1], start_point_a[2], face_a.Key())

                if considered_point_1_id not in considered_point_ids:
                    considered_point_ids.append(considered_point_1_id)
                    corner_points.append(corner_points_a[0])
                    corner_point_faces.append(face_a)
                    corner_point_parameters.append(parameters_a[0])

                # End point on face a
                considered_point_3_id = (end_point_a[0], end_point_a[1], end_point_a[2], face_a.Key())

                if considered_point_3_id not in considered_point_ids:
                    considered_point_ids.append(considered_point_3_id)
                    corner_points.append(corner_points_a[1])
                    corner_point_faces.append(face_a)
                    corner_point_parameters.append(parameters_a[1])

                # Start point on face b
                considered_point_2_id = (start_point_b[0], start_point_b[1], start_point_b[2], face_b.Key())

                if considered_point_2_id not in considered_point_ids:
                    considered_point_ids.append(considered_point_2_id)
                    corner_points.append(corner_points_b[0])
                    corner_point_faces.append(face_b)
                    corner_point_parameters.append(parameters_b[0])

                # End point on face b
                considered_point_4_id = (end_point_b[0], end_point_b[1], end_point_b[2], face_b.Key())

                if considered_point_4_id not in considered_point_ids:
                    considered_point_ids.append(considered_point_4_id)
                    corner_points.append(corner_points_b[1])
                    corner_point_faces.append(face_b)
                    corner_point_parameters.append(parameters_b[1])
            else:
                raise RuntimeError("Max number of adjacent has to be 2!!")

        # Collect corner points in model part
        temp_model = KratosMultiphysics.Model()
        destination_mdpa = temp_model.CreateModelPart("temp_model_part")
        destination_mdpa.AddNodalSolutionStepVariable(KratosShape.NORMALIZED_SURFACE_NORMAL)
        destination_mdpa.AddNodalSolutionStepVariable(KratosShape.SHAPE_CHANGE)

        for itr, [x,y,z] in enumerate(corner_points):
            point_ptr = self.cad_model.Add(an.Point3D(location=[x, y, z]))
            point_ptr.Attributes().SetLayer('CornerPoints')
            destination_mdpa.CreateNewNode(itr+1, x, y, z)

        # Map information from fem to corner points using element based mapper
        mapper = KratosMapping.MapperFactory.CreateMapper( self.fe_model_part, destination_mdpa, self.parameters["conditions"]["general"]["mapping_cad_fem"].Clone() )
        mapper.Map( KratosShape.SHAPE_CHANGE, KratosShape.SHAPE_CHANGE )
        mapper.Map( KratosShape.NORMALIZED_SURFACE_NORMAL, KratosShape.NORMALIZED_SURFACE_NORMAL )

        # Create conditions for each corner point
        for itr, node in enumerate(destination_mdpa.Nodes):
            integration_weight = 1

            face = corner_point_faces[itr]

            surface_geometry = face.Data().Geometry()
            surface_geometry_data = face.Data().Geometry().Data()

            shape_function = an.SurfaceShapeEvaluator(degreeU=surface_geometry_data.DegreeU(), degreeV=surface_geometry_data.DegreeV(), order=1)

            # Create conditions to enforce t1 and t2 on both face a and face b
            u = corner_point_parameters[itr][0]
            v = corner_point_parameters[itr][1]

            shape_function.Compute(surface_geometry_data.KnotsU(), surface_geometry_data.KnotsV(), u, v)
            nonzero_pole_indices = shape_function.NonzeroPoleIndices()

            shape_function_values = np.empty(shape_function.NbNonzeroPoles(), dtype=float)
            shape_function_derivatives_u = np.empty(shape_function.NbNonzeroPoles(), dtype=float)
            shape_function_derivatives_v = np.empty(shape_function.NbNonzeroPoles(), dtype=float)
            for i in range(shape_function.NbNonzeroPoles()):
                shape_function_values[i] = shape_function(0,i)
                shape_function_derivatives_u[i] = shape_function(1,i)
                shape_function_derivatives_v[i] = shape_function(2,i)

            # Tangents enforcement
            target_normal = node.GetSolutionStepValue(KratosShape.NORMALIZED_SURFACE_NORMAL)
            weight = penalty_factor * integration_weight

            new_condition = clib.TangentEnforcementCondition(target_normal, surface_geometry, nonzero_pole_indices, shape_function_derivatives_u, shape_function_derivatives_v, weight)
            conditions[face.Key()].append(new_condition)

            # Positions enforcement
            target_displacement = node.GetSolutionStepValue(KratosShape.SHAPE_CHANGE)
            target_position = np.array([node.X+target_displacement[0], node.Y+target_displacement[1], node.Z+target_displacement[2]])
            weight = penalty_factor * integration_weight

            new_condition = clib.PositionEnforcementCondition(target_position, surface_geometry, nonzero_pole_indices, shape_function_values, weight)
            conditions[face.Key()].append(new_condition)

    # --------------------------------------------------------------------------
    def AddCouplingConditions(self, conditions):
        drawing_tolerance = self.parameters["drawing_parameters"]["cad_drawing_tolerance"].GetDouble()
        min_span_length = self.parameters["drawing_parameters"]["min_span_length"].GetDouble()
        penalty_factor_displacement = self.parameters["conditions"]["edges"]["coupling"]["penalty_factor_displacement_coupling"].GetDouble()
        penalty_factor_rotation = self.parameters["conditions"]["edges"]["coupling"]["penalty_factor_rotation_coupling"].GetDouble()

        # Create corresponding points
        for edge_itr, edge_i in enumerate(self.cad_model.GetByType('BrepEdge')):

            adjacent_faces = edge_i.Data().Faces()
            num_adjacent_faces = len(adjacent_faces)

            if num_adjacent_faces == 1:
                # Skip non coupling-edges
                pass
            elif num_adjacent_faces == 2:

                face_a = adjacent_faces[0]
                face_b = adjacent_faces[1]

                surface_geometry_a = face_a.Data().Geometry()
                surface_geometry_b = face_b.Data().Geometry()
                surface_geometry_data_a = face_a.Data().Geometry().Data()
                surface_geometry_data_b = face_b.Data().Geometry().Data()

                trim_a, _ = edge_i.Data().Trims()
                edge_curve_a = an.CurveOnSurface3D(trim_a.Data().Geometry().Data(), trim_a.Data().Face().Data().Geometry().Data(), trim_a.Data().Geometry().Data().Domain())

                list_of_points, list_of_parameters_a, list_of_parameters_b, list_of_curve_parameters_a, _, list_of_integration_weights = self.__CreateIntegrationPointsForEdge(edge_i, drawing_tolerance, min_span_length)

                # for point in list_of_points:
                #     point_ptr = self.cad_model.Add(an.Point3D(location=point))
                #     point_ptr.Attributes().SetLayer('CouplingPoints')

                # Create conditions
                for (u_a, v_a), (u_b, v_b), t_a, integration_weight in zip(list_of_parameters_a, list_of_parameters_b, list_of_curve_parameters_a, list_of_integration_weights):

                    shape_function_a = an.SurfaceShapeEvaluator(degreeU=surface_geometry_data_a.DegreeU(), degreeV=surface_geometry_data_a.DegreeV(), order=1)
                    shape_function_b = an.SurfaceShapeEvaluator(degreeU=surface_geometry_data_b.DegreeU(), degreeV=surface_geometry_data_b.DegreeV(), order=1)

                    # Create conditions to enforce t1 and t2 on both face a and face b
                    shape_function_a.Compute(surface_geometry_data_a.KnotsU(), surface_geometry_data_a.KnotsV(), u_a, v_a)
                    nonzero_pole_indices_a = shape_function_a.NonzeroPoleIndices()

                    shape_function_values_a = np.empty(shape_function_a.NbNonzeroPoles(), dtype=float)
                    shape_function_derivatives_u_a = np.empty(shape_function_a.NbNonzeroPoles(), dtype=float)
                    shape_function_derivatives_v_a = np.empty(shape_function_a.NbNonzeroPoles(), dtype=float)
                    for i in range(shape_function_a.NbNonzeroPoles()):
                        shape_function_values_a[i] = shape_function_a(0,i)
                        shape_function_derivatives_u_a[i] = shape_function_a(1,i)
                        shape_function_derivatives_v_a[i] = shape_function_a(2,i)

                    shape_function_b.Compute(surface_geometry_data_b.KnotsU(), surface_geometry_data_b.KnotsV(), u_b, v_b)
                    nonzero_pole_indices_b = shape_function_b.NonzeroPoleIndices()

                    shape_function_values_b = np.empty(shape_function_b.NbNonzeroPoles(), dtype=float)
                    shape_function_derivatives_u_b = np.empty(shape_function_b.NbNonzeroPoles(), dtype=float)
                    shape_function_derivatives_v_b = np.empty(shape_function_b.NbNonzeroPoles(), dtype=float)
                    for i in range(shape_function_b.NbNonzeroPoles()):
                        shape_function_values_b[i] = shape_function_b(0,i)
                        shape_function_derivatives_u_b[i] = shape_function_b(1,i)
                        shape_function_derivatives_v_b[i] = shape_function_b(2,i)

                    _, T2_edge = edge_curve_a.DerivativesAt(t_a, order=1)
                    T2_edge /= np.linalg.norm(T2_edge)

                    # displacement coupling condition
                    if penalty_factor_displacement > 0:
                        weight = penalty_factor_displacement * integration_weight

                        new_condition = clib.DisplacementCouplingCondition( surface_geometry_a,
                                                                            surface_geometry_b,
                                                                            (u_a, v_a),
                                                                            (u_b, v_b),
                                                                            nonzero_pole_indices_a,
                                                                            nonzero_pole_indices_b,
                                                                            shape_function_values_a,
                                                                            shape_function_values_b,
                                                                            weight )
                        conditions[face_a.Key()].append(new_condition)

                    # rotation coupling condition
                    if penalty_factor_rotation > 0:
                        weight = penalty_factor_rotation * integration_weight

                        new_condition = clib.RotationCouplingConditionWithAD( surface_geometry_a,
                                                                              surface_geometry_b,
                                                                              (u_a, v_a),
                                                                              (u_b, v_b),
                                                                              T2_edge,
                                                                              nonzero_pole_indices_a,
                                                                              nonzero_pole_indices_b,
                                                                              shape_function_values_a,
                                                                              shape_function_values_b,
                                                                              shape_function_derivatives_u_a,
                                                                              shape_function_derivatives_u_b,
                                                                              shape_function_derivatives_v_a,
                                                                              shape_function_derivatives_v_b,
                                                                              weight )
                        conditions[face_a.Key()].append(new_condition)
            else:
                raise RuntimeError("Max number of adjacent has to be 2!!")

    # --------------------------------------------------------------------------
    def AddAlphaRegularizationConditions(self, conditions):
        alpha = self.parameters["regularization"]["alpha"].GetDouble()

        boundary_polygons = self.GetBoundaryPolygons()

        for face_i in self.cad_model.GetByType('BrepFace'):

            # Skipp embedded faces to not have two identical contributions from the same unknowns (embedded faces share the unkonws of a given geometry)
            if face_i.Attributes().HasTag('Embedded'):
                print(f'Skip {face_i.Key()}')
                continue

            surface_geometry = face_i.Data().Geometry()
            surface_geometry_data = face_i.Data().Geometry().Data()

            shape_function = an.SurfaceShapeEvaluator(degreeU=surface_geometry_data.DegreeU(), degreeV=surface_geometry_data.DegreeV(), order=0)

            deg_u = surface_geometry_data.DegreeU()
            deg_v = surface_geometry_data.DegreeV()
            knot_vec_u = surface_geometry_data.KnotsU()
            knot_vec_v = surface_geometry_data.KnotsV()

            # Compute Grevillle abscissa
            pole_indices = []
            greville_points = []
            greville_abscissa_parameters = []

            for r in range(surface_geometry_data.NbPolesU()):
                for s in range(surface_geometry_data.NbPolesV()):
                    u_value = 0.0
                    v_value = 0.0

                    for p_index in range(0,deg_u):
                        u_value += knot_vec_u[r+p_index]
                    u_value /= deg_u

                    for q_index in range(0,deg_v):
                        v_value += knot_vec_v[s+q_index]
                    v_value /= deg_v

                    is_inside, is_on_boundary = self.Contains((u_value,v_value), boundary_polygons[face_i.Key()], self.boundary_tessellation_tolerance*1.1)

                    # Only control points within the visible surface shall be considered
                    if is_inside or self.parameters["regularization"]["include_all_poles"].GetBool():
                        pole_indices.append((r,s))
                        greville_abscissa_parameters.append((u_value,v_value))
                        greville_points.append(surface_geometry_data.PointAt(u_value, v_value))
                        # point_ptr = self.cad_model.Add(an.Point3D(location=grevile_point))
                        # point_ptr.Attributes().SetLayer('GrevillePoints')

            # Create condition for each of the relevant control points
            for pole_id, grevile_point, greville_params in zip(pole_indices, greville_points, greville_abscissa_parameters):
                u = greville_params[0]
                v = greville_params[1]

                shape_function.Compute(surface_geometry_data.KnotsU(), surface_geometry_data.KnotsV(), u, v)
                nonzero_pole_indices = shape_function.NonzeroPoleIndices()

                shape_function_values = np.empty(shape_function.NbNonzeroPoles(), dtype=float)
                for i in range(shape_function.NbNonzeroPoles()):
                    shape_function_values[i] = shape_function(0,i)

                new_condition = clib.AlphaRegularizationCondition(pole_id, greville_params, surface_geometry, nonzero_pole_indices, shape_function_values, alpha)
                conditions[face_i.Key()].append(new_condition)

        return conditions

    # --------------------------------------------------------------------------
    def GetFEPointParametrization(self):
        if self.fe_point_parametrization is not None:
            return self.fe_point_parametrization
        else:
            self.fe_point_parametrization = self.__CreateFEPointParametrization()
            return self.fe_point_parametrization

    # --------------------------------------------------------------------------
    def GetBoundaryPolygons(self):
        if self.boundary_polygons is not None:
            return self.boundary_polygons
        else:
            self.boundary_polygons = self.__CreateBoundaryPolygons(self.cad_model, self.boundary_tessellation_tolerance)
            return self.boundary_polygons

    # --------------------------------------------------------------------------
    def __CreateFEPointParametrization(self):
        print("> Starting to parametrize FE points...")
        start_time = time.time()

        fe_point_parametrization = []
        for node_i in self.fe_model_part.Nodes:
            fe_point_parametrization.append({"node": node_i , "faces": [], "parameters": [], "is_on_boundary": False})

        boundary_polygons = self.GetBoundaryPolygons()

        for face_i in self.cad_model.GetByType('BrepFace'):

            # Skipp embedded faces to not have two identical contributions from the same unknowns (embedded faces share the unkonws of a given geometry)
            if face_i.Attributes().HasTag('Embedded'):
                print(f'Skip {face_i.Key()}')
                continue

            surface = an.Surface3D(face_i.Data().Geometry())
            projection = an.PointOnSurfaceProjection3D(surface)

            # Bounding box and scaling as the bounding box is based on a simple pointwise discretization of surface
            min_x, min_y, min_z, max_x, max_y, max_z = projection.BoundingBox()
            min_x -= self.bounding_box_tolerance
            min_y -= self.bounding_box_tolerance
            min_z -= self.bounding_box_tolerance
            max_x += self.bounding_box_tolerance
            max_y += self.bounding_box_tolerance
            max_z += self.bounding_box_tolerance

            for node_itr, node_i in enumerate(self.fe_model_part.Nodes):

                node_coords_i = [node_i.X0, node_i.Y0, node_i.Z0]

                # Points outside bounding box are not considered
                if node_coords_i[0] < min_x or max_x < node_coords_i[0]:
                    continue
                if node_coords_i[1] < min_y or max_y < node_coords_i[1]:
                    continue
                if node_coords_i[2] < min_z or max_z < node_coords_i[2]:
                    continue

                projection.Compute(point=node_coords_i)
                projected_point_uv = np.array([projection.ParameterU(), projection.ParameterV()])

                is_inside, is_on_boundary = self.Contains(projected_point_uv, boundary_polygons[face_i.Key()], self.boundary_tessellation_tolerance*1.1)
                if is_inside:
                    fe_point_parametrization[node_itr]["faces"].append(face_i)
                    fe_point_parametrization[node_itr]["parameters"].append(projected_point_uv)
                    fe_point_parametrization[node_itr]["is_on_boundary"] = is_on_boundary

        # Check results
        for entry in fe_point_parametrization:
            node = entry["node"]
            list_of_faces = entry["faces"]
            if len(list_of_faces) == 0:
                print("> WARNING: Missing point pair for point: ", node.Id)
                point_ptr = self.cad_model.Add(an.Point3D(location=[node.X, node.Y, node.Z]))
                point_ptr.Attributes().SetLayer('FEPointsWithNoCADPartner')

        print("> Finished parametrization of fe points in" ,round( time.time()-start_time, 3 ), " s.")

        return fe_point_parametrization

    # --------------------------------------------------------------------------
    @staticmethod
    def __CreateBoundaryPolygons(cad_model, boundary_tessellation_tolerance):
        boundary_polygons = {face.Key(): [] for face in cad_model.GetByType('BrepFace')}

        tessellation = an.CurveTessellation2D()
        for face_i in cad_model.GetByType('BrepFace'):

            for trim in face_i.Data().Trims():
                tessellation.Compute(an.Curve2D(trim.Data().Geometry()), boundary_tessellation_tolerance)
                for i in range(tessellation.NbPoints()):
                    boundary_polygons[face_i.Key()].append(tessellation.Point(i))

                    # (u,v) = tessellation.Point(i)
                    # point = face_i.Data().Geometry().Data().PointAt(u,v)
                    # point_ptr = cad_model.Add(an.Point3D(location=point))

        return boundary_polygons

    # --------------------------------------------------------------------------
    @staticmethod
    def __CreateIntegrationPointsForEdge(edge, drawing_tolerance, min_span_length):
        projection_tolerance = drawing_tolerance * 10

        trim_a, trim_b = edge.Data().Trims()

        curve_3d_a = an.CurveOnSurface3D(trim_a.Data().Geometry().Data(), trim_a.Data().Face().Data().Geometry().Data(), trim_a.Data().Geometry().Data().Domain())
        curve_3d_b = an.CurveOnSurface3D(trim_b.Data().Geometry().Data(), trim_b.Data().Face().Data().Geometry().Data(), trim_b.Data().Geometry().Data().Domain())

        projection_a = an.PointOnCurveProjection3D(curve=curve_3d_a, tolerance=projection_tolerance)
        projection_b = an.PointOnCurveProjection3D(curve=curve_3d_b, tolerance=projection_tolerance)

        spans_on_curve_b = [span_b.T0() for span_b in curve_3d_b.Spans()]

        for span_a in curve_3d_a.Spans():
            t_a = span_a.T0()
            point = curve_3d_a.PointAt(t_a)
            projection_b.Compute(point)
            t_b = projection_b.Parameter()

            spans_on_curve_b.append(t_b)

        spans_on_curve_b.append(curve_3d_b.Domain().T1())

        spans_on_curve_b.sort()

        face_a, face_b = edge.Data().Faces()

        surface_3d_a = face_a.Data().Geometry().Data()
        surface_3d_b = face_b.Data().Geometry().Data()

        degree = max(surface_3d_a.DegreeU(), surface_3d_a.DegreeV(), surface_3d_b.DegreeU(), surface_3d_b.DegreeV()) + 1

        list_of_points = []
        list_of_parameters_a = []
        list_of_parameters_b = []
        list_of_curve_parameters_a = []
        list_of_curve_parameters_b = []
        list_of_weights = []

        for t0_b, t1_b in zip(spans_on_curve_b, spans_on_curve_b[1:]):
            span_b = an.Interval(t0_b, t1_b)

            if span_b.Length() < min_span_length:
                continue

            for t_b, weight in an.IntegrationPoints.Points1D(degree, span_b):
                point, a1 = curve_3d_b.DerivativesAt(t_b, 1)
                list_of_points.append(point)

                projection_a.Compute(point)
                t_a = projection_a.Parameter()

                u_a, v_a = trim_a.Data().Geometry().Data().PointAt(t_a)
                u_b, v_b = trim_b.Data().Geometry().Data().PointAt(t_b)

                list_of_parameters_a.append((u_a, v_a))
                list_of_parameters_b.append((u_b, v_b))
                list_of_curve_parameters_a.append(t_a)
                list_of_curve_parameters_b.append(t_b)
                list_of_weights.append(weight * la.norm(a1))

        return list_of_points, list_of_parameters_a, list_of_parameters_b, list_of_curve_parameters_a, list_of_curve_parameters_b, list_of_weights

    # --------------------------------------------------------------------------
    @staticmethod
    def __CreateIntegrationPointsForFace(face, drawing_tolerance):
        surface_geometry = face.Data().Geometry().Data()

        clipper = an.TrimmedSurfaceClipping(tolerance=drawing_tolerance*10, scale=drawing_tolerance/10.0)
        integration_points = an.PolygonIntegrationPoints()

        clipper.Clear()

        for loop in face.Data().Loops():
            clipper.BeginLoop()

            for trim in loop.Data().Trims():
                clipper.AddCurve(an.Curve2D(trim.Data().Geometry()))

            clipper.EndLoop()

        clipper.Compute(surface_geometry.SpansU(), surface_geometry.SpansV())

        degree_u = surface_geometry.DegreeU()
        degree_v = surface_geometry.DegreeV()

        degree = max(degree_u, degree_v) + 1

        list_of_points = []
        list_of_parameters = []
        list_of_weights = []

        for i in range(clipper.NbSpansU()):
            for j in range(clipper.NbSpansV()):
                if clipper.SpanTrimType(i, j) == an.Empty:
                    continue

                if clipper.SpanTrimType(i, j) == an.Full:
                    for u, v, weight in an.IntegrationPoints.Points2D(degree_u+1, degree_v+1, clipper.SpanU(i), clipper.SpanV(j)):
                        [x,y,z], a1, a2  = surface_geometry.DerivativesAt(u, v, 1)

                        list_of_points.append([x,y,z])
                        list_of_parameters.append((u, v))
                        list_of_weights.append(weight * la.norm(np.cross(a1,a2)))
                else:
                    for polygon_i, polygon in enumerate(clipper.SpanPolygons(i, j)):
                        integration_points.Compute(degree, polygon)

                        for k in range(integration_points.NbIntegrationPoints()):
                            u, v, weight = integration_points.IntegrationPoint(k)

                            [x,y,z], a1, a2  = surface_geometry.DerivativesAt(u, v, 1)

                            list_of_points.append([x,y,z])
                            list_of_parameters.append((u, v))
                            list_of_weights.append(weight * la.norm(np.cross(a1,a2)))

        return list_of_points, list_of_parameters, list_of_weights

    # --------------------------------------------------------------------------
    @staticmethod
    def __IdentifyStartAndEndPoint(edge):
        trim_a, trim_b = edge.Data().Trims()

        curve_3d_a = an.CurveOnSurface3D(trim_a.Data().Geometry().Data(), trim_a.Data().Face().Data().Geometry().Data(), trim_a.Data().Geometry().Data().Domain())
        curve_3d_b = an.CurveOnSurface3D(trim_b.Data().Geometry().Data(), trim_b.Data().Face().Data().Geometry().Data(), trim_b.Data().Geometry().Data().Domain())

        t0_a = trim_a.Data().Geometry().Data().Domain().T0()
        t1_a = trim_a.Data().Geometry().Data().Domain().T1()

        t0_b = trim_b.Data().Geometry().Data().Domain().T0()
        t1_b = trim_b.Data().Geometry().Data().Domain().T1()

        start_point_a = curve_3d_a.PointAt(t0_a)
        end_point_a = curve_3d_a.PointAt(t1_a)

        start_point_b = curve_3d_b.PointAt(t0_b)
        end_point_b = curve_3d_b.PointAt(t1_b)

        list_of_parameters_a = []
        list_of_parameters_b = []

        # Parameters for start point
        u_start_a, v_start_a = trim_a.Data().Geometry().Data().PointAt(t0_a)
        u_start_b, v_start_b = trim_b.Data().Geometry().Data().PointAt(t0_b)

        list_of_parameters_a.append((u_start_a, v_start_a))
        list_of_parameters_b.append((u_start_b, v_start_b))

        # Parameters for end point
        u_end_a, v_end_a = trim_a.Data().Geometry().Data().PointAt(t1_a)
        u_end_b, v_end_b = trim_b.Data().Geometry().Data().PointAt(t1_b)

        list_of_parameters_a.append((u_end_a, v_end_a))
        list_of_parameters_b.append((u_end_b, v_end_b))

        return  [start_point_a, end_point_a], [start_point_b, end_point_b], list_of_parameters_a, list_of_parameters_b

    # --------------------------------------------------------------------------
    @classmethod
    def Contains(cls, point, polygon, tolerance):
        inside = False
        on_boundary = False

        i = 0
        j = len(polygon) - 1

        while i < len(polygon):
            U0 = polygon[i]
            U1 = polygon[j]

            if cls.__IsPointOnLine(point, U0, U1, tolerance):
                return True, True

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

        return inside, on_boundary

    # --------------------------------------------------------------------------
    @staticmethod
    def __IsPointOnLine(point, line_a, line_b, tolerance):
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

# ==============================================================================