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

# check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

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
        fem_input_filename = self.parameters["inpute"]["fem_filename"].GetString()
        name_variable_to_map = self.parameters["inpute"]["variable_to_map"].GetString()
        variable_to_map = KratosMultiphysics.KratosGlobals.GetVariable(name_variable_to_map)

        self.fe_model_part = self.fe_model.CreateModelPart("origin_part")
        self.fe_model_part.AddNodalSolutionStepVariable(variable_to_map)

        model_part_io = KratosMultiphysics.ModelPartIO(fem_input_filename[:-5])
        model_part_io.ReadModelPart(self.fe_model_part)

        # Refine if specified
        prop_id = 1
        prop = self.fe_model_part.Properties[prop_id]
        mat = KratosCSM.LinearElasticPlaneStress2DLaw()
        prop.SetValue(KratosMultiphysics.CONSTITUTIVE_LAW, mat.Clone())

        refinement_level = self.parameters["inpute"]["fe_refinement_level"].GetInt()
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

        # Output FE-data with projected points
        output_dir = self.parameters["output"]["results_directory"].GetString()
        fem_output_filename = "fe_model_used_for_reconstruction"
        fem_output_filename_with_path = os.path.join(output_dir,fem_output_filename)
        nodal_variables = [self.parameters["inpute"]["variable_to_map"].GetString()]
        self.__OutputFEData(self.fe_model_part, fem_output_filename_with_path, nodal_variables)

        # Read CAD data
        cad_filename = self.parameters["inpute"]["cad_filename"].GetString()
        self.cad_model = an.Model.open(cad_filename)

        print("> Preprocessing finished in" ,round( time.time()-start_time, 3 ), " s.")
        print("\n> Starting creation of conditions...")
        start_time = time.time()

        # Create conditions
        for face_itr, face_i in enumerate(self.cad_model.of_type('BrepFace')):
            self.conditions.append([])

        condition_factory = ConditionsFactory(self.fe_model_part.Nodes, self.cad_model, self.parameters)
        condition_factory.CreateDisplacementMappingConditions(self.conditions)

        print("> Finished creation of conditions in" ,round( time.time()-start_time, 3 ), " s.")
        print("\n> Initializing assembly...")
        start_time = time.time()

        # Initialize Assembly
        self.assembler = Assembler(self.fe_model_part.Nodes, self.cad_model, self.conditions)
        self.assembler.Initialize()

        print("> Initialization of assembly finished in" ,round( time.time()-start_time, 3 ), " s.")

    # --------------------------------------------------------------------------
    def Map(self):
        # Assemble
        lhs, rhs = self.assembler.AssembleSystem()

        print("> Number of equations: ", lhs.shape[0])
        print("> Number of relevant control points: ", lhs.shape[1])

        # Beta regularization
        lhs_diag = np.diag(lhs)
        for i in range(lhs_diag.shape[0]):
            entry = lhs[i,i]

            # regularization
            lhs[i,i] += 0.001

            if entry == 0:
                raise RuntimeError
                print(i)
                print("Zero on main diagonal found")

        # Nonlinear solution iterations
        for solution_itr in range(1,self.parameters["solution"]["iterations"].GetInt()+1):
            print("\n> ----------------------------------------------------")
            print("> Starting solution iteration", solution_itr,"...")
            start_time = time.time()

            print("\n> Starting system solution ....")
            start_time = time.time()

            solution_x = la.solve(lhs,rhs[:,0])
            solution_y = la.solve(lhs,rhs[:,1])
            solution_z = la.solve(lhs,rhs[:,2])

            print("> Finished system solution in" ,round( time.time()-start_time, 3 ), " s.")

            if self.parameters["solution"]["test_solution"].GetBool():
                # Test solution quality
                test_rhs = np.zeros(rhs.shape)
                test_rhs[:,0] = lhs.dot(solution_x)
                test_rhs[:,1] = lhs.dot(solution_y)
                test_rhs[:,2] = lhs.dot(solution_z)

                delta = rhs-test_rhs
                delta_combined = np.append(delta[:,0], [delta[:,1], delta[:,2]])
                error_norm = la.norm(delta_combined)
                print("\n> Error in linear solution = ",error_norm)

                # Test residuals
                rhs_combined = np.append(rhs[:,0], [rhs[:,1], rhs[:,2]])
                error_norm = la.norm(rhs_combined)
                print("> RHS before current solution iteration = ",error_norm)

            self.__UpdateCADModel(solution_x, solution_y, solution_z)

            if self.parameters["solution"]["test_solution"].GetBool() or self.parameters["solution"]["iterations"].GetInt()>1:
                rhs = self.assembler.AssembleRHS()

            if self.parameters["solution"]["test_solution"].GetBool():
                rhs_combined = np.append(rhs[:,0], [rhs[:,1], rhs[:,2]])
                error_norm = la.norm(rhs_combined)
                print("\n> RHS after current solution iteration = ",error_norm)

            print("\n> Finished solution iteration in" ,round( time.time()-start_time, 3 ), " s.")
            print("> ----------------------------------------------------")

    # --------------------------------------------------------------------------
    def Finalize(self):
        print("\n> Finalizing mapping....")
        start_time = time.time()

        # Output cad model
        output_dir = self.parameters["output"]["results_directory"].GetString()
        output_filename = self.parameters["output"]["resulting_geometry_filename"].GetString()
        output_filename_with_path = os.path.join(output_dir,output_filename)
        self.cad_model.save(output_filename_with_path)

        print("> Finished finalization of mapping in" ,round( time.time()-start_time, 3 ), " s.")

    # --------------------------------------------------------------------------
    def __UpdateCADModel(self, solution_x, solution_y, solution_z):
        print("\n> Updating cad database....")
        start_time = time.time()

        pole_update = list(zip(solution_x, solution_y, solution_z))
        dof_ids = self.assembler.GetDofIds()
        dofs = self.assembler.GetDofs()

        for face_itr, face_i in enumerate(self.cad_model.of_type('BrepFace')):
            surface_geometry = face_i.surface_geometry_3d().geometry

            for r in range(surface_geometry.NbPolesU):
                for s in range(surface_geometry.NbPolesV):
                    dof_i = (face_itr,r,s)
                    if dof_i in dofs:
                        dof_id = dof_ids[dof_i]
                        pole_coords = surface_geometry.Pole(r,s)
                        new_pole_coords = pole_coords + pole_update[dof_id]
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
    def __init__(self, fe_node_set, cad_model, parameters):
        self.fe_node_set = fe_node_set
        self.cad_model = cad_model

        name_variable_to_map = parameters["inpute"]["variable_to_map"].GetString()
        self.variable_to_map = KratosMultiphysics.KratosGlobals.GetVariable(name_variable_to_map)

        self.tesselation_tolerance = parameters["points_projection"]["boundary_tessellation_tolerance"].GetDouble()
        self.bounding_box_tolerance = parameters["points_projection"]["patch_bounding_box_tolerance"].GetDouble()

    # --------------------------------------------------------------------------
    def CreateDisplacementMappingConditions(self, conditions):
        from cad_reconstruction_conditions import DisplacementMappingCondition

        point_pairs = []
        for node_i in self.fe_node_set:
            point_pairs.append([])

        tessellation = an.CurveTessellation2D()

        for face_itr, face_i in enumerate(self.cad_model.of_type('BrepFace')):

            print("> Processing face ",face_itr)

            surface_geometry = face_i.surface_geometry_3d().geometry
            shape_function = an.SurfaceShapeEvaluator(DegreeU=surface_geometry.DegreeU, DegreeV=surface_geometry.DegreeV, Order=0)

            surface = face_i.surface_3d()
            projection = an.PointOnSurfaceProjection3D(surface)

            boundary_polygon = []
            for trim in face_i.trims():
                tessellation.Compute(trim.curve_2d(), self.tesselation_tolerance)
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

            for node_itr, node_i in enumerate(self.fe_node_set):

                node_coords_i = [node_i.X0, node_i.Y0, node_i.Z0]

                # Points outside bounding box are not considered
                if node_coords_i[0] < min_x or max_x < node_coords_i[0]:
                    continue
                if node_coords_i[1] < min_y or max_y < node_coords_i[1]:
                    continue
                if node_coords_i[2] < min_z or max_z < node_coords_i[2]:
                    continue

                projection.Compute(Point=node_coords_i)
                projected_point_uv = np.array([projection.ParameterU, projection.ParameterV])

                is_inside = self.__Contains(projected_point_uv, boundary_polygon, self.tesselation_tolerance*1.1)
                if is_inside:
                    u = projected_point_uv[0]
                    v = projected_point_uv[1]
                    shape_function.Compute(surface_geometry.KnotsU, surface_geometry.KnotsV, u, v)
                    nonzero_pole_indices = shape_function.NonzeroPoleIndices

                    shape_function_values = []
                    for i in range(len(nonzero_pole_indices)):
                        shape_function_values.append(shape_function(0,i))

                    new_condition = DisplacementMappingCondition(node_i, surface_geometry, nonzero_pole_indices, shape_function_values, self.variable_to_map)
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
        # for node_itr, node_i in enumerate(self.fe_node_set):
        #     node_coords_i = node_coords[node_itr]
        #     projected_point_i = point_pairs[node_itr][0]
        #     distance = projected_point_i - node_coords_i
        #     node_i.SetSolutionStepValue(CONTROL_POINT_UPDATE,[distance[0],distance[1],distance[2]])

        # output_dir = self.parameters["output"]["results_directory"].GetString()
        # fem_output_filename = "fe_model_used_for_reconstruction"
        # fem_output_filename_with_path = os.path.join(output_dir,fem_output_filename)
        # nodal_variables = [self.parameters["inpute"]["update_variable_name"].GetString(), "CONTROL_POINT_UPDATE"]
        # OutputFEData(self.fe_model_part, fem_output_filename_with_path, nodal_variables)

        return conditions

    # --------------------------------------------------------------------------
    @classmethod
    def __Contains(cls, point, polygon, tolerance):
        inside = False

        i = 0
        j = len(polygon) - 1

        while i < len(polygon):
            U0 = polygon[i]
            U1 = polygon[j]

            if cls.__IsPointOnLine(point, U0, U1, tolerance):
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
        for face_itr, face_i in enumerate(self.cad_model.of_type('BrepFace')):
            for condition in self.conditions[face_itr]:
                for (r, s) in condition.nonzero_pole_indices:
                    _ = self.__GetDofId((face_itr,r,s))

        num_dofs = len(self.dofs)

        # Initilize equation system
        self.lhs = np.zeros((num_dofs, num_dofs))
        self.rhs = np.zeros((num_dofs, 3))

        # # testing
        # print(point_pairs[0])
        # print("-------------")
        # for i in range(num_control_points):
        #     if lhs[0,i] != 0:
        #         print("###")
        #         print(i)
        #         print(GetDof(i))
        #         print(lhs[0,i])

        # surface_geometry = model.of_type('BrepFace')[8].surface_geometry_3d().geometry
        # shape_function = an.SurfaceShapeEvaluator(DegreeU=surface_geometry.DegreeU, DegreeV=surface_geometry.DegreeV, Order=0)
        # u = point_pairs[0][0][1][0]
        # v = point_pairs[0][0][1][1]
        # shape_function.Compute(surface_geometry.KnotsU, surface_geometry.KnotsV, u, v)

    # --------------------------------------------------------------------------
    def AssembleSystem(self):
        print("\n> Starting assembly....")
        start_time = time.time()

        self.lhs.fill(0)

        for face_itr, face_i in enumerate(self.cad_model.of_type('BrepFace')):

            print("Processing face", face_itr, "with", len(self.conditions[face_itr]), "conditions.")

            for condition in self.conditions[face_itr]:

                local_lhs = condition.CalculateLHS()
                local_rhs = condition.CalculateRHS()

                nonzero_pole_indices = condition.nonzero_pole_indices

                global_dof_ids = []
                for j, (r, s) in enumerate(nonzero_pole_indices):
                    global_dof_ids.append(self.__GetDofId((face_itr,r,s)))

                for i, dof_id_i in enumerate(global_dof_ids):
                    self.rhs[dof_id_i,:] += local_rhs[i,:]
                    for j, dof_id_j in enumerate(global_dof_ids):
                        self.lhs[dof_id_i, dof_id_j] += local_lhs[i,j]

        print("> Finished assembly in" ,round( time.time()-start_time, 3 ), " s.")

        return self.lhs, self.rhs

    # --------------------------------------------------------------------------
    def AssembleRHS(self):
        print("\n> Starting to assemble RHS....")
        start_time = time.time()

        self.rhs.fill(0)

        for face_itr, face_i in enumerate(self.cad_model.of_type('BrepFace')):

            print("Processing face", face_itr, "with", len(self.conditions[face_itr]), "conditions.")

            for condition in self.conditions[face_itr]:

                local_lhs = condition.CalculateLHS()
                local_rhs = condition.CalculateRHS()

                nonzero_pole_indices = condition.nonzero_pole_indices

                global_dof_ids = []
                for j, (r, s) in enumerate(nonzero_pole_indices):
                    global_dof_ids.append(self.__GetDofId((face_itr,r,s)))

                for i, dof_id_i in enumerate(global_dof_ids):
                    self.rhs[dof_id_i,:] += local_rhs[i,:]

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