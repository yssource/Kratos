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
import time, os, shutil

# ==============================================================================
class CADMapper:
    # --------------------------------------------------------------------------
    def __init__(self, fe_model, cad_model, parameters):
        default_parameters = KratosMultiphysics.Parameters("""
        {
            "input" :
            {
                "cad_filename"                  : "name.iga",
                "fem_filename"                  : "name.mdpa",
                "fe_refinement_level"           : 0
            },
            "conditions" :
            {
                "general" :
                {
                    "apply_integral_method" : false,
                    "mapping_cad_fem"       :
                    {
                        "mapper_type"   : "nearest_element",
                        "search_radius" : -1.0,
                        "echo_level"    : 0
                    }
                },
                "faces" :
                {
                    "curvature" :
                    {
                        "apply_curvature_minimization" : false,
                        "penalty_factor"               : 1e-1
                    },
                    "mechanical" :
                    {
                        "apply_KL_shell"      : false,
                        "exclusive_face_list" : [],
                        "penalty_factor"      : 1e3
                    },
                    "rigid" :
                    {
                        "apply_rigid_conditions" : false,
                        "exclusive_face_list"    : [],
                        "penalty_factor"         : 1e3
                    }
                },
                "edges" :
                {
                    "fe_based" :
                    {
                        "apply_enforcement_conditions"        : false,
                        "penalty_factor_position_enforcement" : 1e3,
                        "penalty_factor_tangent_enforcement"  : 1e3,
                        "apply_corner_enforcement_conditions" : false,
                        "penalty_factor_corner_enforcement"   : 1e4
                    },
                    "coupling" :
                    {
                        "apply_coupling_conditions"            : false,
                        "penalty_factor_displacement_coupling" : 1e3,
                        "penalty_factor_rotation_coupling"     : 1e3
                    }
                }
            },
            "drawing_parameters" :
            {
                "cad_drawing_tolerance"           : 1e-3,
                "boundary_tessellation_tolerance" : 1e-2,
                "patch_bounding_box_tolerance"    : 1.0,
                "min_span_length"                 : 1e-7
            },
            "solution" :
            {
                "iterations"    : 1,
                "test_solution" : true
            },
            "regularization" :
            {
                "alpha"             : 0.1,
                "beta"              : 0.001,
                "include_all_poles" : false
            },
            "refinement" :
            {
                "a_posteriori" :
                {
                    "apply_a_posteriori_refinement" : false,
                    "max_levels_of_refinement"      : 3,
                    "mininimum_knot_distance"       : 1.0,
                    "fe_point_distance_tolerance"   : 0.01,
                    "disp_coupling_tolerance"       : 0.01,
                    "rot_coupling_tolerance"        : 0.5
                },
                "a_priori" :
                {
                    "apply_a_priori_refinement" : false
                }
            },
            "output":
            {
                "results_directory"           : "01_Results",
                "resulting_geometry_filename" : "reconstructed_geometry.iga"
            }
        }""")
        parameters.RecursivelyValidateAndAssignDefaults(default_parameters)

        self.fe_model = fe_model
        self.cad_model = cad_model
        self.parameters = parameters

        fe_model_part_name = "origin_part"
        if self.fe_model.HasModelPart(fe_model_part_name):
            self.fe_model_part_has_to_be_read = False
            self.fe_model_part = self.fe_model.GetModelPart(fe_model_part_name)
        else:
            self.fe_model_part_has_to_be_read = True
            self.fe_model_part = self.fe_model.CreateModelPart(fe_model_part_name)
            self.fe_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 3)
            self.fe_model_part.AddNodalSolutionStepVariable(KratosShape.SHAPE_CHANGE)
            self.fe_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)
            self.fe_model_part.AddNodalSolutionStepVariable(KratosShape.NORMALIZED_SURFACE_NORMAL)
            self.fe_model_part.AddNodalSolutionStepVariable(KratosShape.FITTING_ERROR)

        self.filename_fem_for_reconstruction = "fe_model_used_for_reconstruction"
        self.filename_fem_for_quality_evaluation = "fe_model_with_reconstruction_quality"

        self.assembler = None
        self.conditions = {}
        self.fe_point_parametric = None
        self.absolute_pole_displacement = None

    # --------------------------------------------------------------------------
    def RunMappingProcess(self):
        output_dir = self.parameters["output"]["results_directory"].GetString()
        if os.path.exists(output_dir):
            shutil.rmtree(output_dir)

        self.ReadModelData()

        apply_a_posteriori_refinement = self.parameters["refinement"]["a_posteriori"]["apply_a_posteriori_refinement"].GetBool()
        max_iterations = self.parameters["refinement"]["a_posteriori"]["max_levels_of_refinement"].GetInt()
        resulting_cad_filename = self.parameters["output"]["resulting_geometry_filename"].GetString()

        if apply_a_posteriori_refinement:
            for itr in range(max_iterations):

                print("\n\n========================================================================================================")
                print("> Starting a posteriori refinement level " +  str(itr+1) + "...")
                print("========================================================================================================")

                self.Initialize()
                self.Map()
                self.Finalize()

                cad_model, nothing_to_refine = self.ResetDisplacementsAndRefineCadModel()

                # Rename result files of this iteration
                if max_iterations > 1:
                    case_directory = os.getcwd()
                    os.chdir(self.parameters["output"]["results_directory"].GetString())
                    os.rename(self.filename_fem_for_quality_evaluation+".post.bin", self.filename_fem_for_quality_evaluation+"_level_"+str(itr+1)+".post.bin")
                    os.rename(resulting_cad_filename, resulting_cad_filename.replace(".iga","_level_"+str(itr+1)+".iga"))
                    os.chdir(case_directory)

                if nothing_to_refine:
                    break
                else:
                    self.__SaveCadModel("a_posteriori_refinement_level_"+str(itr+1)+".iga")
        else:
            self.Initialize()
            self.Map()
            self.Finalize()

    # --------------------------------------------------------------------------
    def ReadModelData(self):
        print("\n> Starting to read model data...")

        # Read FEM data
        if self.fe_model_part_has_to_be_read:
            fem_input_filename = self.parameters["input"]["fem_filename"].GetString()
            model_part_io = KratosMultiphysics.ModelPartIO(fem_input_filename[:-5])
            model_part_io.ReadModelPart(self.fe_model_part)

        # Refine if specified
        prop_id = 1
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
            shape_update = node.GetSolutionStepValue(KratosShape.SHAPE_CHANGE)
            node.X += shape_update[0]
            node.Y += shape_update[1]
            node.Z += shape_update[2]

        KratosShape.GeometryUtilities(self.fe_model_part).ComputeUnitSurfaceNormals(True)

        for node in self.fe_model_part.Nodes:
            shape_update = node.GetSolutionStepValue(KratosShape.SHAPE_CHANGE)
            node.X -= shape_update[0]
            node.Y -= shape_update[1]
            node.Z -= shape_update[2]

        # Read CAD data
        cad_filename = self.parameters["input"]["cad_filename"].GetString()
        if len(self.cad_model.GetByType('BrepFace')) == 0:
            self.cad_model.Load(cad_filename)

        start_time = time.time()
        print("> Finished reading of model data in" ,round( time.time()-start_time, 3 ), " s.")

    # --------------------------------------------------------------------------
    def Initialize(self):
        # Initialize (reset) class attributes
        self.assembler = None
        self.conditions = {}
        self.fe_point_parametric = None
        self.absolute_pole_displacement = None

        # Create results folder
        output_dir = self.parameters["output"]["results_directory"].GetString()
        if not os.path.exists( output_dir ):
            os.makedirs( output_dir )

        # Output FE-data with projected points
        nodal_variables = ["SHAPE_CHANGE", "NORMAL", "NORMALIZED_SURFACE_NORMAL"]
        self.__OutputFEData(self.fe_model_part, self.filename_fem_for_reconstruction, nodal_variables)

        print("\n> Starting creation of conditions...")
        start_time = time.time()

        self.__CreateConditions()

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
            self.__ComputeAbsolutePoleDisplacements(solution)

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

            print("\n> Finished solution iteration in" ,round( time.time()-start_time_iteration, 3 ), " s.")
            print("> ----------------------------------------------------")

    # --------------------------------------------------------------------------
    def Finalize(self):
        print("\n> Finalizing mapping....")
        start_time = time.time()

        # Output cad model
        output_filename = self.parameters["output"]["resulting_geometry_filename"].GetString()
        self.__SaveCadModel(output_filename)

        print("> Finished finalization of mapping in" ,round( time.time()-start_time, 3 ), " s.")

    # --------------------------------------------------------------------------
    def ResetDisplacementsAndRefineCadModel(self):
        fe_point_distance_tolerance = self.parameters["refinement"]["a_posteriori"]["fe_point_distance_tolerance"].GetDouble() # in given length unit
        disp_coupling_tolerance = self.parameters["refinement"]["a_posteriori"]["disp_coupling_tolerance"].GetDouble() # in given length unit
        rot_coupling_tolerance = self.parameters["refinement"]["a_posteriori"]["rot_coupling_tolerance"].GetDouble() # in degree

        exist_intervals_to_refine = False
        is_fe_point_distance_satisfied = True
        is_disp_coupling_satisfied = True
        is_rot_coupling_satisfied = True

        # Identify spots to refine
        intervals_along_u_to_refine = {geometry.Key(): [] for geometry in self.cad_model.GetByType('SurfaceGeometry3D')}
        intervals_along_v_to_refine = {geometry.Key(): [] for geometry in self.cad_model.GetByType('SurfaceGeometry3D')}

        # First identify spots where distance to fe points exceeds limit)
        for entry in self.fe_point_parametric:

            node = entry["node"]
            shape_update = node.GetSolutionStepValue(KratosShape.SHAPE_CHANGE)
            node_target_position = np.array([node.X0 + shape_update[0], node.Y0 + shape_update[1], node.Z0 + shape_update[2]])

            list_of_faces = entry["faces"]
            list_of_parameters = entry["parameters"]

            for face, (u,v) in zip(list_of_faces, list_of_parameters):
                geometry = face.Data().Geometry()
                geometry_data = geometry.Data()

                distance = geometry_data.PointAt(u, v) - node_target_position

                node.SetSolutionStepValue(KratosShape.FITTING_ERROR, distance.tolist())

                if np.linalg.norm(distance) > fe_point_distance_tolerance:
                    is_fe_point_distance_satisfied = False
                    exist_intervals_to_refine = self.__AddUIntervalToList(geometry, u, v, intervals_along_u_to_refine)
                    exist_intervals_to_refine = self.__AddVIntervalToList(geometry, u, v, intervals_along_v_to_refine)

          # Output fitting error
        nodal_variables = ["SHAPE_CHANGE", "FITTING_ERROR"]
        self.__OutputFEData(self.fe_model_part, self.filename_fem_for_quality_evaluation, nodal_variables)

        # Then identify spots where coupling or enforcement conditions are not met
        from cad_reconstruction_conditions import DisplacementCouplingCondition, RotationCouplingConditionWithAD
        for conditions_face_i in self.conditions.values():
            for condition in conditions_face_i:
                if isinstance(condition, DisplacementCouplingCondition) or isinstance(condition, RotationCouplingConditionWithAD):
                    geometry_a = condition.geometry_a
                    geometry_b = condition.geometry_b
                    (u_a,v_a) = condition.parameters_a
                    (u_b,v_b) = condition.parameters_b

                    is_spot_to_be_refined = False

                    if isinstance(condition, DisplacementCouplingCondition):
                        delta_disp = np.linalg.norm(condition.CalculateQualityIndicator())
                        if delta_disp > disp_coupling_tolerance:
                            is_disp_coupling_satisfied = False
                            is_spot_to_be_refined = True
                            print("delta_disp =", delta_disp)

                    if isinstance(condition, RotationCouplingConditionWithAD):
                        delta_rot = condition.CalculateQualityIndicator() * 180 / np.pi
                        if delta_rot > rot_coupling_tolerance:
                            is_rot_coupling_satisfied = False
                            is_spot_to_be_refined = True
                            print("delta_rot =", delta_rot)

                    if is_spot_to_be_refined:
                        # Geometry a
                        exist_intervals_to_refine = self.__AddUIntervalToList(geometry_a, u_a, v_a, intervals_along_u_to_refine)
                        exist_intervals_to_refine = self.__AddVIntervalToList(geometry_a, u_a, v_a, intervals_along_v_to_refine)

                        # Geometry b
                        exist_intervals_to_refine = self.__AddUIntervalToList(geometry_b, u_b, v_b, intervals_along_u_to_refine)
                        exist_intervals_to_refine = self.__AddVIntervalToList(geometry_b, u_b, v_b, intervals_along_v_to_refine)

        self.ResetPoleDisplacements()

        # Read original cad file
        nothing_to_refine = False

        if is_fe_point_distance_satisfied and is_disp_coupling_satisfied and is_rot_coupling_satisfied:
            print("\n> Refinement reached convergence! Nothing to refine anymore.")
            nothing_to_refine = True

        elif exist_intervals_to_refine == False:
            print("\n> WARNING!!!! Refinement did not converge but finished as knot tolerance was reached for U- and V-Direction.")
            nothing_to_refine = True

        else:
            for face in self.cad_model.GetByType('BrepFace'):
                print("> Refining face ",face.Key())
                geometry = face.Data().Geometry()

                u_intervals = intervals_along_u_to_refine[geometry.Key()]
                if len(u_intervals) > 0:
                    u_interval_centers = [(interval[0]+interval[1])/2.0 for interval in u_intervals]
                    u_interval_centers = list(set(u_interval_centers)) # Remove duplicated entries
                    u_interval_centers.sort() # Sort in ascending order

                    print(u_interval_centers)

                    refined_geometry = an.KnotRefinement.InsertKnotsU(geometry.Data(), u_interval_centers)
                    self.cad_model.Replace(geometry.Key(), refined_geometry)
                else:
                    print("> Nothing to refine in U-driection!")

                v_intervals = intervals_along_v_to_refine[geometry.Key()]
                if len(v_intervals) > 0:
                    v_interval_centers = [(interval[0]+interval[1])/2.0 for interval in v_intervals]
                    v_interval_centers = list(set(v_interval_centers)) # Remove duplicated entries
                    v_interval_centers.sort() # Sort in ascending order

                    print(v_interval_centers)

                    refined_geometry = an.KnotRefinement.InsertKnotsV(geometry.Data(), v_interval_centers)
                    self.cad_model.Replace(geometry.Key(), refined_geometry)
                else:
                    print("> Nothing to refine in V-driection!")

        return self.cad_model, nothing_to_refine

    # --------------------------------------------------------------------------
    def ResetPoleDisplacements(self):
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
                        pole_update = np.array([self.absolute_pole_displacement[dof_id_x], self.absolute_pole_displacement[dof_id_y], self.absolute_pole_displacement[dof_id_z]])

                        new_pole_coords = pole_coords - pole_update
                        surface_geometry.SetPole(r,s,new_pole_coords)

    # --------------------------------------------------------------------------
    def __CreateConditions(self):
        for face_i in self.cad_model.GetByType('BrepFace'):
            self.conditions[face_i.Key()] = []

        condition_factory = ConditionsFactory(self.fe_model_part, self.cad_model, self.parameters)

        if self.parameters["conditions"]["general"]["apply_integral_method"].GetBool():
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
        if self.parameters["conditions"]["edges"]["fe_based"]["apply_corner_enforcement_conditions"].GetBool():
            condition_factory.CreateCornerEnforcementConditions(self.conditions)
        if self.parameters["conditions"]["edges"]["coupling"]["apply_coupling_conditions"].GetBool():
            condition_factory.CreateCouplingConditions(self.conditions)
        if self.parameters["regularization"]["alpha"].GetDouble() != 0:
            condition_factory.CreateAlphaRegularizationConditions(self.conditions)

        self.fe_point_parametric = condition_factory.GetFEPointParametrization()

    # --------------------------------------------------------------------------
    def __SaveCadModel(self, filename):
        output_dir = self.parameters["output"]["results_directory"].GetString()
        output_filename_with_path = os.path.join(output_dir,filename)
        self.cad_model.Save(output_filename_with_path)

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
    def __ComputeAbsolutePoleDisplacements(self, delta):
        if self.absolute_pole_displacement is None:
            self.absolute_pole_displacement = delta
        else:
            self.absolute_pole_displacement += delta

    # --------------------------------------------------------------------------
    def __AddUIntervalToList(self, geometry, u, v, interval_list):
        intervall_added = False
        min_dist = self.parameters["refinement"]["a_posteriori"]["mininimum_knot_distance"].GetDouble()

        start_span = an.Knots.UpperSpan(geometry.Data().DegreeU(), geometry.Data().KnotsU(), u)
        end_span = start_span+1
        u_start = geometry.Data().KnotsU()[start_span]
        u_end = geometry.Data().KnotsU()[end_span]

        # Avoid refining knots that are close together within a specified tolerance
        # comparisons are made in given length unit
        p_start = geometry.Data().PointAt(u_start, v)
        p_end = geometry.Data().PointAt(u_end, v)
        distance = p_end - p_start
        if np.dot(distance, distance) >= min_dist**2 :
            intervall_added = True
            u_interval = (u_start, u_end)
            interval_list[geometry.Key()].append(u_interval)

        return intervall_added

    # --------------------------------------------------------------------------
    def __AddVIntervalToList(self, geometry, u , v, interval_list):
        intervall_added = False
        min_dist = self.parameters["refinement"]["a_posteriori"]["mininimum_knot_distance"].GetDouble()

        start_span = an.Knots.UpperSpan(geometry.Data().DegreeV(), geometry.Data().KnotsV(), v)
        end_span = start_span+1
        v_start = geometry.Data().KnotsV()[start_span]
        v_end = geometry.Data().KnotsV()[end_span]

        # Avoid refining knots that are close together within a specified tolerance
        # comparisons are made in given length unit
        p_start = geometry.Data().PointAt(u, v_start)
        p_end = geometry.Data().PointAt(u, v_end)
        distance = p_end - p_start
        if np.dot(distance, distance) >= min_dist**2 :
            intervall_added = True
            v_interval = (v_start, v_end)
            interval_list[geometry.Key()].append(v_interval)

        return intervall_added

    # --------------------------------------------------------------------------
    def __OutputFEData(self, model_part, fem_output_filename, nodal_variables):
        output_dir = self.parameters["output"]["results_directory"].GetString()
        fem_output_filename_with_path = os.path.join(output_dir,fem_output_filename)

        from gid_output import GiDOutput
        nodal_results=nodal_variables
        gauss_points_results=[]
        VolumeOutput = True
        GiDPostMode = "Binary"
        GiDWriteMeshFlag = False
        GiDWriteConditionsFlag = True
        GiDMultiFileFlag = "Single"

        gig_io = GiDOutput(fem_output_filename_with_path, VolumeOutput, GiDPostMode, GiDMultiFileFlag, GiDWriteMeshFlag, GiDWriteConditionsFlag)
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

        self.bounding_box_tolerance = self.parameters["drawing_parameters"]["patch_bounding_box_tolerance"].GetDouble()
        self.boundary_tessellation_tolerance = self.parameters["drawing_parameters"]["boundary_tessellation_tolerance"].GetDouble()
        self.boundary_polygons = None
        self.fe_point_parametrization = None

    # --------------------------------------------------------------------------
    def CreateDistanceMinimizationConditions(self, conditions):
        from cad_reconstruction_conditions import DistanceMinimizationCondition

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

                new_condition = DistanceMinimizationCondition(node, surface_geometry, nonzero_pole_indices, shape_function_values, KratosShape.SHAPE_CHANGE, weight)
                conditions[face.Key()].append(new_condition)

        return conditions

    # --------------------------------------------------------------------------
    def CreateDistanceMinimizationWithIntegrationConditions(self, conditions):
        from cad_reconstruction_conditions import DistanceMinimizationCondition

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

                new_condition = DistanceMinimizationCondition(node_i, surface_geometry, nonzero_pole_indices, shape_function_values, KratosShape.SHAPE_CHANGE, weight)
                conditions[face_i.Key()].append(new_condition)

        print("> Total area of cad surface = ",total_area,"\n")

        return conditions

    # --------------------------------------------------------------------------
    def CreateFaceConditions(self, conditions):
        from cad_reconstruction_conditions import CurvatureMinimizationConditionWithAD,  KLShellConditionWithAD

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
                    new_condition = CurvatureMinimizationConditionWithAD(surface_geometry, nonzero_pole_indices, shape_function_derivatives_u, shape_function_derivatives_v, shape_function_derivatives_uu, shape_function_derivatives_uv, shape_function_derivatives_vv, weight)
                    conditions[face_i.Key()].append(new_condition)

                if apply_kl_shell:
                    weight = shell_penalty_fac * weight
                    new_condition = KLShellConditionWithAD(surface_geometry, nonzero_pole_indices, shape_function_derivatives_u, shape_function_derivatives_v, shape_function_derivatives_uu, shape_function_derivatives_uv, shape_function_derivatives_vv, weight)
                    conditions[face_i.Key()].append(new_condition)

        return conditions

    # --------------------------------------------------------------------------
    def CreateRigidConditions(self, conditions):
        from cad_reconstruction_conditions import PositionEnforcementCondition

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

                is_inside, is_on_boundary = self.__Contains(projected_point_uv, boundary_polygons[face_i.Key()], self.boundary_tessellation_tolerance*1.1)
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

                new_condition = PositionEnforcementCondition(rigididly_displaced_point_coords, surface_geometry, nonzero_pole_indices, shape_function_values, penalty_fac)
                conditions[face_i.Key()].append(new_condition)

        return conditions

    # --------------------------------------------------------------------------
    def CreateEnforcementConditions(self, conditions):
        from cad_reconstruction_conditions import TangentEnforcementCondition, PositionEnforcementCondition

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

                        new_condition_a = PositionEnforcementCondition(target_position, surface_geometry_a, nonzero_pole_indices_a, shape_function_values_a, weight)
                        conditions[face_a.Key()].append(new_condition_a)

                        new_condition_b = PositionEnforcementCondition(target_position, surface_geometry_b, nonzero_pole_indices_b, shape_function_values_b, weight)
                        conditions[face_b.Key()].append(new_condition_b)

                    # Tangents enforcement
                    if penalty_factor_tangent_enforcement > 0:
                        target_normal = node.GetSolutionStepValue(KratosShape.NORMALIZED_SURFACE_NORMAL)
                        weight = penalty_factor_tangent_enforcement * integration_weight

                        new_condition_a = TangentEnforcementCondition(target_normal, surface_geometry_a, nonzero_pole_indices_a, shape_function_derivatives_u_a, shape_function_derivatives_v_a, weight)
                        conditions[face_a.Key()].append(new_condition_a)

                        new_condition_b = TangentEnforcementCondition(target_normal, surface_geometry_b, nonzero_pole_indices_b, shape_function_derivatives_u_b, shape_function_derivatives_v_b, weight)
                        conditions[face_b.Key()].append(new_condition_b)
            else:
                raise RuntimeError("Max number of adjacent has to be 2!!")

    # --------------------------------------------------------------------------
    def CreateCornerEnforcementConditions(self, conditions):
        # This conditions assumes an integration weight of 1 and other than that uses the penalty factors from the tangent and position enforcement
        from cad_reconstruction_conditions import TangentEnforcementCondition, PositionEnforcementCondition

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

            new_condition = TangentEnforcementCondition(target_normal, surface_geometry, nonzero_pole_indices, shape_function_derivatives_u, shape_function_derivatives_v, weight)
            conditions[face.Key()].append(new_condition)

            # Positions enforcement
            target_displacement = node.GetSolutionStepValue(KratosShape.SHAPE_CHANGE)
            target_position = np.array([node.X+target_displacement[0], node.Y+target_displacement[1], node.Z+target_displacement[2]])
            weight = penalty_factor * integration_weight

            new_condition = PositionEnforcementCondition(target_position, surface_geometry, nonzero_pole_indices, shape_function_values, weight)
            conditions[face.Key()].append(new_condition)

    # --------------------------------------------------------------------------
    def CreateCouplingConditions(self, conditions):
        from cad_reconstruction_conditions import DisplacementCouplingCondition, RotationCouplingConditionWithAD

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

                        new_condition = DisplacementCouplingCondition( surface_geometry_a,
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

                        new_condition = RotationCouplingConditionWithAD( surface_geometry_a,
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
    def CreateAlphaRegularizationConditions(self, conditions):
        from cad_reconstruction_conditions import AlphaRegularizationCondition

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

                    is_inside, is_on_boundary = self.__Contains((u_value,v_value), boundary_polygons[face_i.Key()], self.boundary_tessellation_tolerance*1.1)

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

                new_condition = AlphaRegularizationCondition(pole_id, greville_params, surface_geometry, nonzero_pole_indices, shape_function_values, alpha)
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
    def __CreateFEPointParametrization(self):
        self.fe_point_parametrization = []
        for node_i in self.fe_model_part.Nodes:
            self.fe_point_parametrization.append({"node": node_i , "faces": [], "parameters": [], "is_on_boundary": False})

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

                is_inside, is_on_boundary = self.__Contains(projected_point_uv, boundary_polygons[face_i.Key()], self.boundary_tessellation_tolerance*1.1)
                if is_inside:
                    self.fe_point_parametrization[node_itr]["faces"].append(face_i)
                    self.fe_point_parametrization[node_itr]["parameters"].append(projected_point_uv)
                    self.fe_point_parametrization[node_itr]["is_on_boundary"] = is_on_boundary

        # Check results
        for entry in self.fe_point_parametrization:
            node = entry["node"]
            list_of_faces = entry["faces"]
            if len(list_of_faces) == 0:
                print("> WARNING: Missing point pair for point: ", node.Id)
                point_ptr = self.cad_model.Add(an.Point3D(location=[node.X, node.Y, node.Z]))
                point_ptr.Attributes().SetLayer('FEPointsWithNoCADPartner')

        return self.fe_point_parametrization

    # --------------------------------------------------------------------------
    def GetBoundaryPolygons(self):
        if self.boundary_polygons is not None:
            return self.boundary_polygons
        else:
            self.boundary_polygons = self.__CreateBoundaryPolygons()
            return self.boundary_polygons

    # --------------------------------------------------------------------------
    def __CreateBoundaryPolygons(self):
        self.boundary_polygons = {face.Key(): [] for face in self.cad_model.GetByType('BrepFace')}

        tessellation = an.CurveTessellation2D()
        for face_i in self.cad_model.GetByType('BrepFace'):

            for trim in face_i.Data().Trims():
                tessellation.Compute(an.Curve2D(trim.Data().Geometry()), self.boundary_tessellation_tolerance)
                for i in range(tessellation.NbPoints()):
                    self.boundary_polygons[face_i.Key()].append(tessellation.Point(i))

                    # (u,v) = tessellation.Point(i)
                    # point = face_i.Data().Geometry().Data().PointAt(u,v)
                    # point_ptr = self.cad_model.Add(an.Point3D(location=point))

        return self.boundary_polygons

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
        self.eqs = []

        self.lhs = np.zeros((0, 0))
        self.rhs = np.zeros((0, 0))

    # --------------------------------------------------------------------------
    def Initialize(self):
        # Assign dof ids
        for conditions_face_i in self.conditions.values():
            for condition in conditions_face_i:
                for dof in condition.dof_list:
                    _ = self.__GetDofId(dof)

        num_dofs = len(self.dofs)

        # Initilize equation system
        self.lhs = np.zeros((num_dofs, num_dofs))
        self.rhs = np.zeros(num_dofs)

    # --------------------------------------------------------------------------
    def AssembleSystem(self):
        print("\n> Starting assembly....")
        start_time = time.time()

        self.lhs.fill(0)
        self.rhs.fill(0)

        total_num_conditions = 0

        for face in self.cad_model.GetByType('BrepFace'):
            face_key = face.Key()

            print("Processing face", face_key, "with", len(self.conditions[face_key]), "conditions.")

            for condition in self.conditions[face_key]:
                total_num_conditions += 1

                condition_lhs, condition_rhs, condition_dof_list = condition.CalculateLocalSystem()

                global_dof_ids = []
                for dof in condition_dof_list:
                    global_dof_ids.append(self.__GetDofId(dof))

                self.lhs[np.ix_(global_dof_ids, global_dof_ids)] += condition_lhs
                self.rhs[global_dof_ids] += condition_rhs

        print("Total number of conditions = ",total_num_conditions)
        print("> Finished assembly in" ,round( time.time()-start_time, 3 ), " s.")

        return self.lhs, self.rhs

    # --------------------------------------------------------------------------
    def AssembleRHS(self):
        print("\n> Starting to assemble RHS....")
        start_time = time.time()

        self.rhs.fill(0)

        for face in self.cad_model.GetByType('BrepFace'):
            face_key = face.Key()

            print("Processing face", face_key, "with", len(self.conditions[face_key]), "conditions.")

            for condition in self.conditions[face_key]:

                local_rhs, condition_dof_list = condition.CalculateRHS()

                global_dof_ids = []
                for dof in condition_dof_list:
                    global_dof_ids.append(self.__GetDofId(dof))

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