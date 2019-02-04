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
import KratosMultiphysics.MeshingApplication as KratosMeshingApp
import KratosMultiphysics.StructuralMechanicsApplication as KratosCSM
import KratosMultiphysics.MappingApplication as KratosMapping

# Import ANurbs library
import ANurbs as an
import numpy as np
import numpy.linalg as la

# Additional imports
import time
import os

# ==============================================================================
class CADMapper:
    # --------------------------------------------------------------------------
    def __init__(self, fe_model, cad_model, parameters):
        self.fe_model = fe_model
        self.cad_model = cad_model
        self.parameters = parameters

        self.conditions = []
        self.fe_model_part = None
        self.assembler = None

    # --------------------------------------------------------------------------
    def Initialize(self):
        print("\n> Starting preprocessing...")
        start_time = time.time()

        # Create results folder
        output_dir = self.parameters["output"]["results_directory"].GetString()
        if not os.path.exists( output_dir ):
            os.makedirs( output_dir )

        # Read FE data
        fe_model_part_name = "origin_part"
        name_variable_to_map = self.parameters["input"]["variable_to_map"].GetString()
        variable_to_map = KratosMultiphysics.KratosGlobals.GetVariable(name_variable_to_map)

        if self.fe_model.HasModelPart(fe_model_part_name):
            self.fe_model_part = self.fe_model.GetModelPart(fe_model_part_name)
        else:

            self.fe_model_part = self.fe_model.CreateModelPart(fe_model_part_name)
            self.fe_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 3)
            self.fe_model_part.AddNodalSolutionStepVariable(variable_to_map)
            self.fe_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)
            self.fe_model_part.AddNodalSolutionStepVariable(KratosShape.NORMALIZED_SURFACE_NORMAL)

            fem_input_filename = self.parameters["input"]["fem_filename"].GetString()
            model_part_io = KratosMultiphysics.ModelPartIO(fem_input_filename[:-5])
            model_part_io.ReadModelPart(self.fe_model_part)

        # Refine if specified
        prop_id = 0
        prop = self.fe_model_part.Properties[prop_id]
        mat = KratosCSM.LinearElasticPlaneStress2DLaw()
        prop.SetValue(KratosMultiphysics.CONSTITUTIVE_LAW, mat.Clone())

        refinement_level = self.parameters["input"]["fe_refinement_level"].GetInt()
        for refinement_level in range(0,refinement_level):

            number_of_avg_elems = 10
            number_of_avg_nodes = 10
            nodal_neighbour_search = KratosMultiphysics.FindNodalNeighboursProcess(self.fe_model_part, number_of_avg_elems, number_of_avg_nodes)
            neighbour_calculator = KratosMultiphysics.FindElementalNeighboursProcess(self.fe_model_part,2,10)
            nodal_neighbour_search.Execute()
            neighbour_calculator.Execute()

            for elem in self.fe_model_part.Elements:
                elem.SetValue(KratosMultiphysics.SPLIT_ELEMENT,True)

            refine_on_reference = False
            interpolate_internal_variables = True
            Refine = KratosMeshingApp.LocalRefineTriangleMesh(self.fe_model_part)
            Refine.LocalRefineMesh(refine_on_reference, interpolate_internal_variables)

        # Compute average surface normals of target design
        for node in self.fe_model_part.Nodes:
            shape_update = node.GetSolutionStepValue(variable_to_map)
            node.X += shape_update[0]
            node.Y += shape_update[1]
            node.Z += shape_update[2]

        KratosShape.GeometryUtilities(self.fe_model_part).ComputeUnitSurfaceNormals(True)

        for node in self.fe_model_part.Nodes:
            shape_update = node.GetSolutionStepValue(variable_to_map)
            node.X -= shape_update[0]
            node.Y -= shape_update[1]
            node.Z -= shape_update[2]

        # Output FE-data with projected points
        output_dir = self.parameters["output"]["results_directory"].GetString()
        fem_output_filename = "fe_model_used_for_reconstruction"
        fem_output_filename_with_path = os.path.join(output_dir,fem_output_filename)
        nodal_variables = [self.parameters["input"]["variable_to_map"].GetString(), "NORMAL", "NORMALIZED_SURFACE_NORMAL"]
        self.__OutputFEData(self.fe_model_part, fem_output_filename_with_path, nodal_variables)

        # Read CAD data
        cad_filename = self.parameters["input"]["cad_filename"].GetString()
        self.cad_model = an.Model()
        self.cad_model.Load(cad_filename)

        print("> Preprocessing finished in" ,round( time.time()-start_time, 3 ), " s.")
        print("\n> Starting creation of conditions...")
        start_time = time.time()

        # Create conditions
        for face_itr, face_i in enumerate(self.cad_model.GetByType('BrepFace')):
            self.conditions.append([])

        condition_factory = ConditionsFactory(self.fe_model_part, self.cad_model, self.parameters)

        if self.parameters["conditions"]["apply_integral_method"].GetBool():
            condition_factory.CreateDistanceMinimizationWithIntegrationConditions(self.conditions)
        else:
            condition_factory.CreateDistanceMinimizationConditions(self.conditions)
        if self.parameters["conditions"]["faces"]["curvature"]["apply_curvature_minimization"].GetBool() or \
           self.parameters["conditions"]["faces"]["mechanical"]["apply_KL_shell"].GetBool():
            condition_factory.CreateFaceConditions(self.conditions)
        if self.parameters["conditions"]["faces"]["rigid"]["apply_rigid_conditions"].GetBool():
            condition_factory.CreateRigidConditions(self.conditions)
        if self.parameters["conditions"]["edges"]["fe_based"]["apply_enforcement_conditions"].GetBool():
            condition_factory.CreateEnforcementConditions(self.conditions)

        print("> Finished creation of conditions in" ,round( time.time()-start_time, 3 ), " s.")
        print("\n> Initializing assembly...")
        start_time = time.time()

        # Initialize Assembly
        self.assembler = Assembler(self.fe_model_part.Nodes, self.cad_model, self.conditions)
        self.assembler.Initialize()

        print("> Initialization of assembly finished in" ,round( time.time()-start_time, 3 ), " s.")

    # --------------------------------------------------------------------------
    def Map(self):

        # Nonlinear solution iterations
        for solution_itr in range(1,self.parameters["solution"]["iterations"].GetInt()+1):

            # Assemble
            lhs, rhs = self.assembler.AssembleSystem()

            if solution_itr == 1:
                print("> Number of equations: ", lhs.shape[0])
                print("> Number of relevant control points: ", lhs.shape[1])

            print("\n> ----------------------------------------------------")
            print("> Starting solution iteration", solution_itr,"...")
            start_time_iteration = time.time()

            print("\n> Starting system solution ....")
            start_time_solution = time.time()

            # Beta regularization
            beta = self.parameters["regularization"]["beta"].GetDouble()
            lhs_diag = np.diag(lhs)
            for i in range(lhs_diag.shape[0]):
                entry = lhs[i,i]
                # if abs(entry) < 1e-12:
                #     print("WARNING!!!!Zero on main diagonal found at position",i,". Make sure to include beta regularization.")
                lhs[i,i] += beta

            solution = la.solve(lhs,rhs)

            print("> Finished system solution in" ,round( time.time()-start_time_solution, 3 ), " s.")

            self.__UpdateCADModel(solution)

            if self.parameters["solution"]["test_solution"].GetBool():

                # Test solution quality
                test_rhs = np.zeros(rhs.shape)
                test_rhs[:] = lhs.dot(solution)

                delta = rhs-test_rhs
                error_norm = la.norm(delta)
                print("\n> Error in linear solution = ",error_norm)

                # Test residuals
                error_norm = la.norm(rhs)
                print("> RHS before current solution iteration = ",error_norm)

                # Test rhs after update
                rhs = self.assembler.AssembleRHS()

                error_norm = la.norm(rhs)
                print("\n> RHS after current solution iteration = ",error_norm)

                # Varying contribution of beta regularization is neglected as each solution iteration may be seen indendently
                # in terms of minimization of the control point displacement

                # dof_ids = self.assembler.GetDofIds()
                # dofs = self.assembler.GetDofs()

                # for face_itr, face_i in enumerate(self.cad_model.GetByType('BrepFace')):
                #     surface_geometry = face_i.Data().Geometry().Data()

                #     for r in range(surface_geometry.NbPolesU):
                #         for s in range(surface_geometry.NbPolesV):
                #             dof_i = (face_itr,r,s)
                #             if dof_i in dofs:
                #                 dof_id = dof_ids[dof_i]
                #                 rhs[dof_id,:] += beta*absolute_pole_update[dof_id,:]

            print("\n> Finished solution iteration in" ,round( time.time()-start_time_iteration, 3 ), " s.")
            print("> ----------------------------------------------------")

    # --------------------------------------------------------------------------
    def Finalize(self):
        print("\n> Finalizing mapping....")
        start_time = time.time()

        # Output cad model
        output_dir = self.parameters["output"]["results_directory"].GetString()
        output_filename = self.parameters["output"]["resulting_geometry_filename"].GetString()
        output_filename_with_path = os.path.join(output_dir,output_filename)
        self.cad_model.Save(output_filename_with_path)

        print("> Finished finalization of mapping in" ,round( time.time()-start_time, 3 ), " s.")

    # --------------------------------------------------------------------------
    def __UpdateCADModel(self, solution):
        print("\n> Updating cad database....")
        start_time = time.time()

        dof_ids = self.assembler.GetDofIds()
        dofs = self.assembler.GetDofs()

        for surface_i in self.cad_model.GetByType('SurfaceGeometry3D'):
            surface_geometry = surface_i.Data()
            surface_geometry_key = surface_i.Key()

            for r in range(surface_geometry.NbPolesU()):
                for s in range(surface_geometry.NbPolesV()):
                    dof_i_x = (surface_geometry_key,r,s,"x")
                    dof_i_y = (surface_geometry_key,r,s,"y")
                    dof_i_z = (surface_geometry_key,r,s,"z")

                    if dof_i_x in dofs:
                        dof_id_x = dof_ids[dof_i_x]
                        dof_id_y = dof_ids[dof_i_y]
                        dof_id_z = dof_ids[dof_i_z]

                        pole_coords = surface_geometry.Pole(r,s)
                        pole_update = np.array([solution[dof_id_x], solution[dof_id_y], solution[dof_id_z]])

                        new_pole_coords = pole_coords + pole_update
                        surface_geometry.SetPole(r,s,new_pole_coords)

        print("> Finished updating cad database in" ,round( time.time()-start_time, 3 ), " s.")

    # --------------------------------------------------------------------------
    @staticmethod
    def __OutputFEData(model_part, fem_output_filename, nodal_variables):
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

# ==============================================================================
class ConditionsFactory:
    # --------------------------------------------------------------------------
    def __init__(self, fe_model_part, cad_model, parameters):
        self.fe_model_part = fe_model_part
        self.cad_model = cad_model
        self.parameters = parameters

    # --------------------------------------------------------------------------
    def CreateDistanceMinimizationConditions(self, conditions):
        from cad_reconstruction_conditions import DistanceMinimizationCondition

        name_variable_to_map = self.parameters["input"]["variable_to_map"].GetString()
        variable_to_map = KratosMultiphysics.KratosGlobals.GetVariable(name_variable_to_map)

        tesselation_tolerance = self.parameters["points_projection"]["boundary_tessellation_tolerance"].GetDouble()
        bounding_box_tolerance = self.parameters["points_projection"]["patch_bounding_box_tolerance"].GetDouble()

        point_pairs = []
        for node_i in self.fe_model_part.Nodes:
            point_pairs.append([])

        tessellation = an.CurveTessellation2D()

        for face_itr, face_i in enumerate(self.cad_model.GetByType('BrepFace')):

            # Skipp embedded faces to not have two identical contributions from the same unknowns (embedded faces share the unkonws of a given geometry)
            if face_i.Attributes().HasTag('Embedded'):
                print(f'Skip {face_i.Key()}')
                continue

            print("> Processing face ",face_itr)

            surface_geometry = face_i.Data().Geometry().Data()
            shape_function = an.SurfaceShapeEvaluator(degreeU=surface_geometry.DegreeU(), degreeV=surface_geometry.DegreeV(), order=0)

            surface = an.Surface3D(face_i.Data().Geometry())
            projection = an.PointOnSurfaceProjection3D(surface)

            boundary_polygon = []
            for trim in face_i.Data().Trims():
                tessellation.Compute(an.Curve2D(trim.Data().Geometry()), tesselation_tolerance)
                for i in range(tessellation.NbPoints()):
                    boundary_polygon.append(tessellation.Point(i))

            # Bounding box and scaling as the bounding box is based on a simple pointwise discretization of surface
            min_x, min_y, min_z, max_x, max_y, max_z = projection.BoundingBox()
            min_x -= bounding_box_tolerance
            min_y -= bounding_box_tolerance
            min_z -= bounding_box_tolerance
            max_x += bounding_box_tolerance
            max_y += bounding_box_tolerance
            max_z += bounding_box_tolerance

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

                is_inside, is_on_boundary = self.__Contains(projected_point_uv, boundary_polygon, tesselation_tolerance*1.1)
                if is_inside:
                    u = projected_point_uv[0]
                    v = projected_point_uv[1]
                    shape_function.Compute(surface_geometry.KnotsU(), surface_geometry.KnotsV(), u, v)
                    nonzero_pole_indices = shape_function.NonzeroPoleIndices()

                    shape_function_values = np.empty(shape_function.NbNonzeroPoles(), dtype=float)
                    for i in range(shape_function.NbNonzeroPoles()):
                        shape_function_values[i] = shape_function(0,i)

                    # Introduce a possible penalty factor for nodes on boundary
                    weight = 1.0
                    # if is_on_boundary:
                    #     self.cad_model.Add(an.Point3D(location=node_coords_i))

                    new_condition = DistanceMinimizationCondition(node_i, surface_geometry, nonzero_pole_indices, shape_function_values, variable_to_map, weight)
                    conditions[face_itr].append(new_condition)

                    projected_point_i = projection.Point
                    point_pairs[node_itr].append(projected_point_i)

        # Check results
        for itr, entry in enumerate(point_pairs):
            if len(entry) == 0:
                print("> WARNING: Missing point pair for point: ", itr)

        # i = 0
        # for itr, entry in enumerate(point_pairs):
        #     if len(entry) >= 1:
        #         for (x, y, z), _, face_itr in entry:
        #             self.cad_model.add({
        #                 'Key': f'Point3D<{i}>',
        #                 'Type': 'Point3D',
        #                 'Location': [x, y, z],
        #                 'Layer': f'PointPairs{face_itr}',
        #             })
        #             i += 1
        #     if len(entry) == 0:
        #         print("> WARNING: Missing point pair for point: ", itr)

        # for face_itr, entry in enumerate(conditions):
        #         for node_itr, (x, y, z), _ in entry:
        #             self.cad_model.add({
        #                 'Key': f'Point3D<{i}>',
        #                 'Type': 'Point3D',
        #                 'Location': [x, y, z],
        #                 'Layer': f'PointPairs{face_itr}',
        #             })
        #             i += 1

        # Some additional output
        # # Output FE-data with projected points
        # for node_itr, node_i in enumerate(self.fe_model_part.Nodes):
        #     node_coords_i = node_coords[node_itr]
        #     projected_point_i = point_pairs[node_itr][0]
        #     distance = projected_point_i - node_coords_i
        #     node_i.SetSolutionStepValue(CONTROL_POINT_UPDATE,[distance[0],distance[1],distance[2]])

        # output_dir = self.parameters["output"]["results_directory"].GetString()
        # fem_output_filename = "fe_model_used_for_reconstruction"
        # fem_output_filename_with_path = os.path.join(output_dir,fem_output_filename)
        # nodal_variables = [self.parameters["input"]["update_variable_name"].GetString(), "CONTROL_POINT_UPDATE"]
        # OutputFEData(self.fe_model_part, fem_output_filename_with_path, nodal_variables)

        return conditions

    # --------------------------------------------------------------------------
    def CreateDistanceMinimizationWithIntegrationConditions(self, conditions):
        from cad_reconstruction_conditions import DistanceMinimizationCondition

        name_variable_to_map = self.parameters["input"]["variable_to_map"].GetString()
        variable_to_map = KratosMultiphysics.KratosGlobals.GetVariable(name_variable_to_map)

        total_area = 0

        for face_itr, face_i in enumerate(self.cad_model.GetByType('BrepFace')):

            print("> Processing face ",face_itr)

            surface_geometry = face_i.Data().Geometry().Data()
            shape_function = an.SurfaceShapeEvaluator(degreeU=surface_geometry.DegreeU(), degreeV=surface_geometry.DegreeV(), order=0)

            list_of_points, list_of_parameters, list_of_integration_weights = self.__CreateIntegrationPointsForFace(face_i)

            # Collect integration points in model part
            temp_model = KratosMultiphysics.Model()
            destination_mdpa = temp_model.CreateModelPart("temp_model_part")
            destination_mdpa.AddNodalSolutionStepVariable(KratosShape.NORMALIZED_SURFACE_NORMAL)
            destination_mdpa.AddNodalSolutionStepVariable(KratosShape.SHAPE_CHANGE)

            for itr, [x,y,z] in enumerate(list_of_points):
                # self.cad_model.add({
                #     'Key': f'FaceIntegrationPoint3D<{face_itr},{itr}>',
                #     'Type': 'Point3D',
                #     'Location': [x, y, z],
                #     'Layer': 'FaceIntegrationPoints',
                # })
                destination_mdpa.CreateNewNode(itr, x, y, z)

            # Map information from fem to integration points using element based mapper
            mapper_parameters = KratosMultiphysics.Parameters("""{
                "mapper_type" : "nearest_element",
                "search_radius" : 1.0
            }""")

            mapper = KratosMapping.MapperFactory.CreateMapper( self.fe_model_part, destination_mdpa, mapper_parameters )
            mapper.Map( variable_to_map, variable_to_map )

            # Create conditions
            for itr, node_i in enumerate(destination_mdpa.Nodes):

                [u,v] = list_of_parameters[itr]
                weight = list_of_integration_weights[itr]

                total_area += weight

                shape_function.Compute(surface_geometry.KnotsU(), surface_geometry.KnotsV(), u, v)
                nonzero_pole_indices = shape_function.NonzeroPoleIndices()

                shape_function_values = np.empty(shape_function.NbNonzeroPoles(), dtype=float)
                for i in range(shape_function.NbNonzeroPoles()):
                    shape_function_values[i] = shape_function(0,i)

                new_condition = DistanceMinimizationCondition(node_i, surface_geometry, nonzero_pole_indices, shape_function_values, variable_to_map, weight)
                conditions[face_itr].append(new_condition)

        print("> Total area of cad surface = ",total_area,"\n")

        return conditions

    # --------------------------------------------------------------------------
    def CreateFaceConditions(self, conditions):
        from cad_reconstruction_conditions import CurvatureMinimizationConditionWithAD,  KLShellConditionWithAD

        apply_curvature_min = self.parameters["conditions"]["faces"]["curvature"]["apply_curvature_minimization"].GetBool()
        curvature_penalty_fac = self.parameters["conditions"]["faces"]["curvature"]["penalty_factor"].GetDouble()

        apply_kl_shell = self.parameters["conditions"]["faces"]["mechanical"]["apply_KL_shell"].GetBool()
        shell_penalty_fac = self.parameters["conditions"]["faces"]["mechanical"]["penalty_factor"].GetDouble()

        list_of_exclusive_faces = []
        for itr in range(self.parameters["conditions"]["faces"]["mechanical"]["exclusive_face_list"].size()):
            list_of_exclusive_faces.append(self.parameters["conditions"]["faces"]["mechanical"]["exclusive_face_list"][itr].GetString())

        for face_itr, face_i in enumerate(self.cad_model.GetByType('BrepFace')):

            print("> Processing face ",face_itr)

            # Skip faces if exclusive face list is specified (if list is empty, use all faces)
            if face_i.Key() not in list_of_exclusive_faces:
                continue

            surface_geometry = face_i.Data().Geometry().Data()
            shape_function = an.SurfaceShapeEvaluator(degreeU=surface_geometry.DegreeU(), degreeV=surface_geometry.DegreeV(), order=2)

            list_of_points, list_of_parameters, list_of_integration_weights = self.__CreateIntegrationPointsForFace(face_i)

            # Create conditions
            for itr in range(len(list_of_points)):

                [u,v] = list_of_parameters[itr]
                weight = list_of_integration_weights[itr]

                shape_function.Compute(surface_geometry.KnotsU(), surface_geometry.KnotsV(), u, v)
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
                    new_condition = CurvatureMinimizationConditionWithAD(surface_geometry, nonzero_pole_indices, shape_function_derivatives_u, shape_function_derivatives_v, shape_function_derivatives_uu, shape_function_derivatives_uv, shape_function_derivatives_vv, weight)
                    conditions[face_itr].append(new_condition)

                if apply_kl_shell:
                    weight = shell_penalty_fac * weight
                    new_condition = KLShellConditionWithAD(surface_geometry, nonzero_pole_indices, shape_function_derivatives_u, shape_function_derivatives_v, shape_function_derivatives_uu, shape_function_derivatives_uv, shape_function_derivatives_vv, weight)
                    conditions[face_itr].append(new_condition)

        return conditions

 # --------------------------------------------------------------------------
    def CreateRigidConditions(self, conditions):
        from cad_reconstruction_conditions import PositionEnforcementCondition

        list_of_exclusive_faces = []
        for itr in range(self.parameters["conditions"]["faces"]["rigid"]["exclusive_face_list"].size()):
            list_of_exclusive_faces.append(self.parameters["conditions"]["faces"]["rigid"]["exclusive_face_list"][itr].GetString())

        penalty_fac = self.parameters["conditions"]["faces"]["rigid"]["penalty_factor"].GetDouble()
        tesselation_tolerance = self.parameters["points_projection"]["boundary_tessellation_tolerance"].GetDouble()
        bounding_box_tolerance = self.parameters["points_projection"]["patch_bounding_box_tolerance"].GetDouble()
        name_variable_to_map = self.parameters["input"]["variable_to_map"].GetString()
        variable_to_map = KratosMultiphysics.KratosGlobals.GetVariable(name_variable_to_map)

        relevant_fe_points = []
        relevant_fe_points_displaced = []
        relevant_cad_uvs = []

        tessellation = an.CurveTessellation2D()

        for face_itr, face_i in enumerate(self.cad_model.GetByType('BrepFace')):

            print("> Processing face ",face_itr)

            # Skip faces if exclusive face list is specified (if list is empty, use all faces)
            if face_i.Key() not in list_of_exclusive_faces:
                continue

            surface_geometry = face_i.Data().Geometry().Data()
            shape_function = an.SurfaceShapeEvaluator(degreeU=surface_geometry.DegreeU(), degreeV=surface_geometry.DegreeV(), order=0)

            surface = an.Surface3D(face_i.Data().Geometry())
            projection = an.PointOnSurfaceProjection3D(surface)

            boundary_polygon = []
            for trim in face_i.Data().Trims():
                tessellation.Compute(an.Curve2D(trim.Data().Geometry()), tesselation_tolerance)
                for i in range(tessellation.NbPoints()):
                    boundary_polygon.append(tessellation.Point(i))

            # Bounding box and scaling as the bounding box is based on a simple pointwise discretization of surface
            min_x, min_y, min_z, max_x, max_y, max_z = projection.BoundingBox()
            min_x -= bounding_box_tolerance
            min_y -= bounding_box_tolerance
            min_z -= bounding_box_tolerance
            max_x += bounding_box_tolerance
            max_y += bounding_box_tolerance
            max_z += bounding_box_tolerance

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

                is_inside, is_on_boundary = self.__Contains(projected_point_uv, boundary_polygon, tesselation_tolerance*1.1)
                if is_inside:

                    relevant_fe_points.append(node_coords_i)
                    relevant_fe_points_displaced.append( node_coords_i + np.array(node_i.GetSolutionStepValue(variable_to_map)) )
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
                self.cad_model.Add(an.Point3D(location=node_coords_disp))

                u = projected_point_uv[0]
                v = projected_point_uv[1]
                shape_function.Compute(surface_geometry.KnotsU(), surface_geometry.KnotsV(), u, v)
                nonzero_pole_indices = shape_function.NonzeroPoleIndices()

                shape_function_values = np.empty(shape_function.NbNonzeroPoles(), dtype=float)
                for i in range(shape_function.NbNonzeroPoles()):
                    shape_function_values[i] = shape_function(0,i)

                new_condition = PositionEnforcementCondition(rigididly_displaced_point_coords, surface_geometry, nonzero_pole_indices, shape_function_values, penalty_fac)
                conditions[face_itr].append(new_condition)

        return conditions

 # --------------------------------------------------------------------------
    def CreateEnforcementConditions(self, conditions):
        from cad_reconstruction_conditions import TangentEnforcementCondition, PositionEnforcementCondition

        penalty_factor_tangent_enforcement = self.parameters["conditions"]["edges"]["fe_based"]["penalty_factor_tangent_enforcement"].GetDouble()
        penalty_factor_position_enforcement = self.parameters["conditions"]["edges"]["fe_based"]["penalty_factor_position_enforcement"].GetDouble()

        face_id_to_itr = {}
        for face_itr, face_i in enumerate(self.cad_model.GetByType('BrepFace')):
            face_id_to_itr[face_i.Key()] = face_itr

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

                face_a_itr = face_id_to_itr[adjacent_faces[0].Key()]
                face_b_itr = face_id_to_itr[adjacent_faces[1].Key()]

                surface_geometry_a = face_a.Data().Geometry().Data()
                surface_geometry_b = face_b.Data().Geometry().Data()

                list_of_points, list_of_parameters_a, list_of_parameters_b, list_of_integration_weights = self.__CreateIntegrationPointsForEdge(edge_i)

                # Collect integration points in model part
                temp_model = KratosMultiphysics.Model()
                destination_mdpa = temp_model.CreateModelPart("temp_model_part")
                destination_mdpa.AddNodalSolutionStepVariable(KratosShape.NORMALIZED_SURFACE_NORMAL)
                destination_mdpa.AddNodalSolutionStepVariable(KratosShape.SHAPE_CHANGE)

                for itr, [x,y,z] in enumerate(list_of_points):
                    # self.cad_model.add({
                    #     'Key': f'CouplingPoint3D<{face_itr},{itr}>',
                    #     'Type': 'Point3D',
                    #     'Location': [x, y, z],
                    #     'Layer': 'CouplingPoints',
                    # })
                    destination_mdpa.CreateNewNode(itr, x, y, z)

                # Map information from fem to integration points using element based mapper
                mapper_parameters = KratosMultiphysics.Parameters("""{
                    "mapper_type" : "nearest_element",
                    "search_radius" : 1.0
                }""")

                mapper = KratosMapping.MapperFactory.CreateMapper( self.fe_model_part, destination_mdpa, mapper_parameters )
                mapper.Map( KratosShape.NORMALIZED_SURFACE_NORMAL, KratosShape.NORMALIZED_SURFACE_NORMAL )
                mapper.Map( KratosShape.SHAPE_CHANGE, KratosShape.SHAPE_CHANGE )

                # Create conditions
                for itr, node in enumerate(destination_mdpa.Nodes):

                    integration_weight = list_of_integration_weights[itr]

                    shape_function_a = an.SurfaceShapeEvaluator(degreeU=surface_geometry_a.DegreeU(), degreeV=surface_geometry_a.DegreeV(), order=1)
                    shape_function_b = an.SurfaceShapeEvaluator(degreeU=surface_geometry_b.DegreeU(), degreeV=surface_geometry_b.DegreeV(), order=1)

                    # Create conditions to enforce t1 and t2 on both face a and face b
                    u_a = list_of_parameters_a[itr][0]
                    v_a = list_of_parameters_a[itr][1]
                    shape_function_a.Compute(surface_geometry_a.KnotsU(), surface_geometry_a.KnotsV(), u_a, v_a)
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
                    shape_function_b.Compute(surface_geometry_b.KnotsU(), surface_geometry_b.KnotsV(), u_b, v_b)
                    nonzero_pole_indices_b = shape_function_b.NonzeroPoleIndices()

                    shape_function_values_b = np.empty(shape_function_b.NbNonzeroPoles(), dtype=float)
                    shape_function_derivatives_u_b = np.empty(shape_function_b.NbNonzeroPoles(), dtype=float)
                    shape_function_derivatives_v_b = np.empty(shape_function_b.NbNonzeroPoles(), dtype=float)
                    for i in range(shape_function_b.NbNonzeroPoles()):
                        shape_function_values_b[i] = shape_function_b(0,i)
                        shape_function_derivatives_u_b[i] = shape_function_b(1,i)
                        shape_function_derivatives_v_b[i] = shape_function_b(2,i)

                    # Tangents enforcement
                    target_normal = node.GetSolutionStepValue(KratosShape.NORMALIZED_SURFACE_NORMAL)
                    weight = penalty_factor_tangent_enforcement * integration_weight

                    new_condition_a = TangentEnforcementCondition(target_normal, surface_geometry_a, nonzero_pole_indices_a, shape_function_derivatives_u_a, shape_function_derivatives_v_a, weight)
                    conditions[face_a_itr].append(new_condition_a)

                    new_condition_b = TangentEnforcementCondition(target_normal, surface_geometry_b, nonzero_pole_indices_b, shape_function_derivatives_u_b, shape_function_derivatives_v_b, weight)
                    conditions[face_b_itr].append(new_condition_b)

                    # Positions enforcement
                    target_displacement = node.GetSolutionStepValue(KratosShape.SHAPE_CHANGE)
                    target_position = np.array([node.X+target_displacement[0], node.Y+target_displacement[1], node.Z+target_displacement[2]])
                    weight = penalty_factor_position_enforcement * integration_weight

                    new_condition_a = PositionEnforcementCondition(target_position, surface_geometry_a, nonzero_pole_indices_a, shape_function_values_a, weight)
                    conditions[face_a_itr].append(new_condition_a)

                    new_condition_b = PositionEnforcementCondition(target_position, surface_geometry_b, nonzero_pole_indices_b, shape_function_values_b, weight)
                    conditions[face_b_itr].append(new_condition_b)

                    # point_counter += 1
                    # self.cad_model.add({
                    #     'Key': f'target_position<{point_counter}>',
                    #     'Type': 'Point3D',
                    #     'Location': [target_position[0],target_position[1],target_position[2]],
                    #     'Layer': 'boundary_target_pos',
                    # })


                    # pole_coords = np.zeros((len(nonzero_pole_indices_a), 3))
                    # for i, (r,s) in enumerate(nonzero_pole_indices_a):
                    #     pole_coords[i,:] = surface_geometry_a.Pole(r,s)

                    # a1 = shape_function_derivatives_u_a @ pole_coords

                    # pole_coords = np.zeros((len(nonzero_pole_indices_b), 3))
                    # for i, (r,s) in enumerate(nonzero_pole_indices_b):
                    #     pole_coords[i,:] = surface_geometry_b.Pole(r,s)

                    # a2 = shape_function_derivatives_u_b @ pole_coords


                    # self.cad_model.add({
                    #     'Key': f'a1_a<{points_counter}>',
                    #     'Type': 'Line3D',
                    #     'Start': [x,y,z],
                    #     'End': [x + a1[0],y + a1[1], z + a1[2]],
                    #     'Layer': 'a1_a',
                    # })

            else:
                raise RuntimeError("Max number of adjacent has to be 2!!")

    # --------------------------------------------------------------------------
    @staticmethod
    def __CreateIntegrationPointsForEdge(edge):

        projection_tolerance = 0.01 # Sollte Zeichengenaigkeit sein

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
        list_of_weights = []

        for t0_b, t1_b in zip(spans_on_curve_b, spans_on_curve_b[1:]):
            span_b = an.Interval(t0_b, t1_b)

            print("Here some parameter is defined!!!!!!!!!")
            if span_b.Length() < 1e-7:
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
                list_of_weights.append(weight * la.norm(a1))

        return list_of_points, list_of_parameters_a, list_of_parameters_b, list_of_weights

    # --------------------------------------------------------------------------
    @staticmethod
    def __CreateIntegrationPointsForFace(face):
        surface_geometry = face.Data().Geometry().Data()

        print("Here there is some parameter!!!!")
        drawing_tolerance = 0.01

        clipper = an.TrimmedSurfaceClipping(tolerance=drawing_tolerance, scale=drawing_tolerance/100.0)
        integration_points = an.PolygonIntegrationPoints()
        tessellation = an.PolygonTessellation3D()

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

                        # self.cad_model.add({
                        #     'Key': f'IntegrationPoint{face_itr}{pt_i}',
                        #     'Type': 'Point3D',
                        #     'Location': [x, y, z],
                        #     'Color': '#0000ff',
                        #     'Layer': 'Full',
                        # })

                        list_of_points.append([x,y,z])
                        list_of_parameters.append((u, v))
                        list_of_weights.append(weight * la.norm(np.cross(a1,a2)))
                else:
                    for polygon_i, polygon in enumerate(clipper.SpanPolygons(i, j)):
                        integration_points.Compute(degree, polygon)

                        for k in range(integration_points.NbIntegrationPoints()):
                            u, v, weight = integration_points.IntegrationPoint(k)

                            [x,y,z], a1, a2  = surface_geometry.DerivativesAt(u, v, 1)

                            # self.cad_model.add({
                            #     'Key': f'IntegrationPoint{face_itr}{pt_i}',
                            #     'Type': 'Point3D',
                            #     'Location': [x, y, z],
                            #     'Color': '#ff0000',
                            #     'Layer': 'Quad3D',
                            # })

                            list_of_points.append([x,y,z])
                            list_of_parameters.append((u, v))
                            list_of_weights.append(weight * la.norm(np.cross(a1,a2)))

        return list_of_points, list_of_parameters, list_of_weights

    # --------------------------------------------------------------------------
    @classmethod
    def __Contains(cls, point, polygon, tolerance):
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
class Assembler():
    # --------------------------------------------------------------------------
    def __init__(self, fe_node_set, cad_model, conditions):
        self.fe_node_set = fe_node_set
        self.cad_model = cad_model
        self.conditions = conditions

        self.dof_ids ={}
        self.dofs = []
        self.eq_ids ={}
        self.eqs = []

        self.lhs = np.zeros((0, 0))
        self.rhs = np.zeros((0, 0))

    # --------------------------------------------------------------------------
    def Initialize(self):
        # Assign dof and equation ids
        for face_itr, face_i in enumerate(self.cad_model.GetByType('BrepFace')):
            surface_geometry_key = face_i.Data().Geometry().Key()
            for condition in self.conditions[face_itr]:
                for (r, s) in condition.nonzero_pole_indices:
                    _ = self.__GetDofId((surface_geometry_key,r,s,"x"))
                    _ = self.__GetDofId((surface_geometry_key,r,s,"y"))
                    _ = self.__GetDofId((surface_geometry_key,r,s,"z"))

        num_dofs = len(self.dofs)

        # Initilize equation system
        self.lhs = np.zeros((num_dofs, num_dofs))
        self.rhs = np.zeros(num_dofs)

        # # testing
        # print(point_pairs[0])
        # print("-------------")
        # for i in range(num_control_points):
        #     if lhs[0,i] != 0:
        #         print("###")
        #         print(i)
        #         print(GetDof(i))
        #         print(lhs[0,i])

        # surface_geometry = model.GetByType('BrepFace')[8].Data().Geometry().Data()
        # shape_function = an.SurfaceShapeEvaluator(DegreeU=surface_geometry.DegreeU, DegreeV=surface_geometry.DegreeV, Order=0)
        # u = point_pairs[0][0][1][0]
        # v = point_pairs[0][0][1][1]
        # shape_function.Compute(surface_geometry.KnotsU(), surface_geometry.KnotsV(), u, v)

    # --------------------------------------------------------------------------
    def AssembleSystem(self):
        print("\n> Starting assembly....")
        start_time = time.time()

        self.lhs.fill(0)
        self.rhs.fill(0)

        total_num_conditions = 0

        for face_itr, face_i in enumerate(self.cad_model.GetByType('BrepFace')):

            print("Processing face", face_itr, "with", len(self.conditions[face_itr]), "conditions.")
            surface_geometry_key = face_i.Data().Geometry().Key()

            for condition in self.conditions[face_itr]:
                total_num_conditions += 1

                local_lhs, local_rhs = condition.CalculateLocalSystem()

                nonzero_pole_indices = condition.nonzero_pole_indices

                global_dof_ids = []
                for j, (r, s) in enumerate(nonzero_pole_indices):
                    global_dof_ids.append(self.__GetDofId((surface_geometry_key,r,s,"x")))

                for j, (r, s) in enumerate(nonzero_pole_indices):
                    global_dof_ids.append(self.__GetDofId((surface_geometry_key,r,s,"y")))

                for j, (r, s) in enumerate(nonzero_pole_indices):
                    global_dof_ids.append(self.__GetDofId((surface_geometry_key,r,s,"z")))


                self.lhs[np.ix_(global_dof_ids, global_dof_ids)] += local_lhs
                self.rhs[global_dof_ids] += local_rhs

        print("Total number of conditions = ",total_num_conditions)
        print("> Finished assembly in" ,round( time.time()-start_time, 3 ), " s.")

        return self.lhs, self.rhs

    # --------------------------------------------------------------------------
    def AssembleRHS(self):
        print("\n> Starting to assemble RHS....")
        start_time = time.time()

        self.rhs.fill(0)

        for face_itr, face_i in enumerate(self.cad_model.GetByType('BrepFace')):

            print("Processing face", face_itr, "with", len(self.conditions[face_itr]), "conditions.")
            surface_geometry_key = face_i.Data().Geometry().Key()

            for condition in self.conditions[face_itr]:

                local_rhs = condition.CalculateRHS()

                nonzero_pole_indices = condition.nonzero_pole_indices

                global_dof_ids = []
                for j, (r, s) in enumerate(nonzero_pole_indices):
                    global_dof_ids.append(self.__GetDofId((surface_geometry_key,r,s,"x")))

                for j, (r, s) in enumerate(nonzero_pole_indices):
                    global_dof_ids.append(self.__GetDofId((surface_geometry_key,r,s,"y")))

                for j, (r, s) in enumerate(nonzero_pole_indices):
                    global_dof_ids.append(self.__GetDofId((surface_geometry_key,r,s,"z")))

                self.rhs[global_dof_ids] += local_rhs

        print("> Finished assembling RHS in" ,round( time.time()-start_time, 3 ), " s.")

        return self.rhs

    # --------------------------------------------------------------------------
    def GetDofs(self):
        if self.dofs == None:
            raise Exception("> Assembler:: No dofs specified yet! First initialize assembler.")
        return self.dofs

    # --------------------------------------------------------------------------
    def GetDofIds(self):
        if self.dofs == None:
            raise Exception("> Assembler:: No dof ids specified yet! First initialize assembler.")
        return self.dof_ids

    # --------------------------------------------------------------------------
    def __GetDofId(self, dof):
        if dof not in self.dof_ids:
            self.dof_ids[dof] = len(self.dof_ids)
            self.dofs.append(dof)
        return self.dof_ids[dof]

    # --------------------------------------------------------------------------
    def __GetDof(self, index):
        return self.dofs[index]

# ==============================================================================