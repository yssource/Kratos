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

# Import helper classes
import cad_reconstruction_conditions as clib
from cad_reconstruction_condition_factory import ConditionFactory

# Additional imports
import ANurbs as an
import numpy as np
import numpy.linalg as la
import time, os, shutil
import EQlib as eq

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
                "fe_refinement_level"           : 0,
                "fe_parametrization_filename"   : ""
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
                    "direct" :
                    {
                        "apply_enforcement_conditions" : false,
                        "exclusive_edge_list": [],
                        "penalty_factor_position_enforcement" : 1e3
                    },
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
                "iterations"        : 1,
                "test_solution"     : true,
                "parallel_assembly" : false
            },
            "regularization" :
            {
                "alpha"             : 0.1,
                "beta"              : 0.001
            },
            "refinement" :
            {
                "a_priori" :
                {
                    "apply_a_priori_refinement"         : false,
                    "max_levels_of_refinement"          : 3,
                    "min_knot_distance_at_max_gradient" : 1.0,
                    "exponent"                          : 2
                },
                "a_posteriori" :
                {
                    "apply_a_posteriori_refinement" : false,
                    "max_levels_of_refinement"      : 3,
                    "mininimum_knot_distance"       : 1.0,
                    "fe_point_distance_tolerance"   : 0.01,
                    "disp_coupling_tolerance"       : 0.01,
                    "rot_coupling_tolerance"        : 0.5
                }
            },
            "output":
            {
                "results_directory"                   : "01_Results",
                "resulting_geometry_filename"         : "reconstructed_geometry",
                "filename_fem_for_reconstruction"     : "fe_model_used_for_reconstruction",
                "filename_fem_for_quality_evaluation" : "fe_model_with_reconstruction_quality",
                "echo_level"                          : 0
            }
        }""")
        parameters.RecursivelyValidateAndAssignDefaults(default_parameters)

        self.fe_model = fe_model
        self.cad_model = cad_model
        self.parameters = parameters

        fe_model_part_name = "origin_part"
        if self.fe_model.HasModelPart(fe_model_part_name):
            self.fe_model_part = self.fe_model.GetModelPart(fe_model_part_name)
        else:
            self.fe_model_part = self.fe_model.CreateModelPart(fe_model_part_name)
            self.fe_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 3)
            self.fe_model_part.AddNodalSolutionStepVariable(KratosShape.SHAPE_CHANGE)
            self.fe_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)
            self.fe_model_part.AddNodalSolutionStepVariable(KratosShape.NORMALIZED_SURFACE_NORMAL)
            self.fe_model_part.AddNodalSolutionStepVariable(KratosShape.FITTING_ERROR)
            self.fe_model_part.AddNodalSolutionStepVariable(KratosShape.GRAD_SHAPE_CHANGE)
            self.fe_model_part.AddNodalSolutionStepVariable(KratosShape.SHAPE_CHANGE_ABSOLUTE)

        self.condition_factory = ConditionFactory(self.fe_model_part, self.cad_model, self.parameters)

        eq.Log.info_level = parameters["output"]["echo_level"].GetInt()

        self.conditions = []
        self.conditions_scaling_factors = []
        self.beta_scaling_factor = None
        self.fe_point_parametric = None
        self.pole_nodes = {}

    # --------------------------------------------------------------------------
    def RunMappingProcess(self):
        apply_a_priori_refinement = self.parameters["refinement"]["a_priori"]["apply_a_priori_refinement"].GetBool()
        apply_a_posteriori_refinement = self.parameters["refinement"]["a_posteriori"]["apply_a_posteriori_refinement"].GetBool()
        max_iterations = self.parameters["refinement"]["a_posteriori"]["max_levels_of_refinement"].GetInt()
        output_dir = self.parameters["output"]["results_directory"].GetString()
        fem_filename = self.parameters["output"]["filename_fem_for_quality_evaluation"].GetString()
        cad_filename = self.parameters["output"]["resulting_geometry_filename"].GetString()

        # Initialize results folder
        if os.path.exists(output_dir):
            shutil.rmtree(output_dir)

        self.ReadModelData()

        # Make sure the FE point parametrization is done on the initial mesh in case of any refinement
        if apply_a_priori_refinement or apply_a_posteriori_refinement:
            self.GetFEPointParametrization()

        if apply_a_priori_refinement:
            self.PerformAPrioriRefinement()

        if apply_a_posteriori_refinement:
            for itr in range(max_iterations):

                print("\n\n========================================================================================================")
                print("> Starting a posteriori refinement level " +  str(itr+1) + "...")
                print("========================================================================================================")

                iteration_tag = "_level_"+str(itr+1)
                self.parameters["output"]["filename_fem_for_quality_evaluation"].SetString(fem_filename + iteration_tag)
                self.parameters["output"]["resulting_geometry_filename"].SetString(cad_filename.replace(".iga", iteration_tag+".iga"))

                self.Initialize()
                self.Map()

                self.__OutputCadModel("reconstructed_geometry_refinement_level"+ iteration_tag +".iga")

                nothing_to_refine = self.ResetDisplacementsAndRefineCadModel()

                if nothing_to_refine:
                    break
                else:
                    self.__OutputCadModel("a_posteriori_refinement"+ iteration_tag +".iga")
        else:
            self.Initialize()
            self.Map()

        if self.parameters["output"]["echo_level"].GetInt() > 0:
            self.__StoreCaseFilesInResults()

    # --------------------------------------------------------------------------
    def ReadModelData(self):
        print("\n> Starting to read model data...")
        start_time = time.time()

        # Read FEM data
        if len(self.fe_model_part.Nodes) == 0:
            fem_input_filename = self.parameters["input"]["fem_filename"].GetString()
            model_part_io = KratosMultiphysics.ModelPartIO(fem_input_filename[:-5])
            model_part_io.ReadModelPart(self.fe_model_part)

        # Refine if specified
        prop_id = 1
        prop = self.fe_model_part.Properties[prop_id]
        mat = KratosCSM.LinearElasticPlaneStress2DLaw()
        prop.SetValue(KratosMultiphysics.CONSTITUTIVE_LAW, mat.Clone())

        # Refine FE if necessary
        for refinement_level in range(0,self.parameters["input"]["fe_refinement_level"].GetInt()):
            number_of_avg_elems = 10
            number_of_avg_nodes = 10

            nodal_neighbour_search = KratosMultiphysics.FindNodalNeighboursProcess(self.fe_model_part, number_of_avg_elems, number_of_avg_nodes)
            nodal_neighbour_search.Execute()

            neighbour_calculator = KratosMultiphysics.FindElementalNeighboursProcess(self.fe_model_part,2,10)
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

        print("> Finished reading of model data in" ,round( time.time()-start_time, 3 ), " s.")

    # --------------------------------------------------------------------------
    def GetFEPointParametrization(self):
        self.fe_point_parametric = self.condition_factory.GetFEPointParametrization()

    # --------------------------------------------------------------------------
    def PerformAPrioriRefinement(self):
        print("\n> Initializing prior refinement...")
        start_time = time.time()

        max_iterations = self.parameters["refinement"]["a_priori"]["max_levels_of_refinement"].GetInt()
        min_knot_distance_at_max_gradient = self.parameters["refinement"]["a_priori"]["min_knot_distance_at_max_gradient"].GetDouble()
        exponent = self.parameters["refinement"]["a_priori"]["exponent"].GetInt()

        # Compute gradient of displacment field if necessary
        for node in self.fe_model_part.Nodes:
            shape_change = np.array(node.GetSolutionStepValue(KratosShape.SHAPE_CHANGE))
            node.SetSolutionStepValue(KratosShape.SHAPE_CHANGE_ABSOLUTE, la.norm(shape_change))

        local_gradient = KratosMultiphysics.ComputeNodalGradientProcess3D(self.fe_model_part, KratosShape.SHAPE_CHANGE_ABSOLUTE, KratosShape.GRAD_SHAPE_CHANGE, KratosMultiphysics.NODAL_AREA)
        local_gradient.Execute()

        l2_norms_grads = [np.array(entry["node"].GetSolutionStepValue(KratosShape.GRAD_SHAPE_CHANGE)) for entry in self.fe_point_parametric]
        l2_norms_grads = [ la.norm(grad) for grad in l2_norms_grads]
        max_norm_grads = max(l2_norms_grads)

        # Perform iterative refinement
        for refinement_itr in range(max_iterations):

            print("> Starting a priorie refinement level " +  str(refinement_itr+1) + "...")

            # Identify spots to refine
            intervals_along_u_to_refine = {geometry.Key(): [] for geometry in self.cad_model.GetByType('SurfaceGeometry3D')}
            intervals_along_v_to_refine = {geometry.Key(): [] for geometry in self.cad_model.GetByType('SurfaceGeometry3D')}

            exist_intervals_to_refine = False

            # First identify spots where distance to fe points exceeds limit
            for entry in self.fe_point_parametric:

                node = entry["node"]

                list_of_faces = entry["faces"]
                list_of_parameters = entry["parameters"]

                for face, (u,v) in zip(list_of_faces, list_of_parameters):
                    geometry = face.Data().Geometry()
                    geometry_data = geometry.Data()

                    grad_shape_change = np.array(node.GetSolutionStepValue(KratosShape.GRAD_SHAPE_CHANGE))
                    scaled_grad_norm = la.norm(grad_shape_change) / max_norm_grads
                    scaled_grad_norm = max(scaled_grad_norm,1e-6)

                    # hyperbolic increase of min_distance: 1 / x^exp : exp = 2 means that ca 70% reduced gradient leads to ca 10x increased min distance
                    required_knot_distance = min_knot_distance_at_max_gradient * (1 / scaled_grad_norm)**exponent

                    u_added = self.__AddUIntervalToList(geometry, u, v, intervals_along_u_to_refine, required_knot_distance)
                    v_added = self.__AddVIntervalToList(geometry, u, v, intervals_along_v_to_refine, required_knot_distance)

                    exist_intervals_to_refine = exist_intervals_to_refine or u_added or v_added

            # Perform actual refinement
            if exist_intervals_to_refine == False:
                print("> Maximum refinement for specified minimum knot distances reached!")
                break
            else:
                self.__RefineAtIntervalCenters(self.cad_model, intervals_along_u_to_refine, intervals_along_v_to_refine)

        self.__OutputCadModel("cad_model_after_prior_refinement.iga")

        print("> Prior refinement finished in" ,round( time.time()-start_time, 3 ), " s.")

    # --------------------------------------------------------------------------
    def Initialize(self):
        nodal_variables = ["SHAPE_CHANGE", "GRAD_SHAPE_CHANGE", "NORMAL", "NORMALIZED_SURFACE_NORMAL"]
        filename = self.parameters["output"]["filename_fem_for_reconstruction"].GetString()
        output_dir = self.parameters["output"]["results_directory"].GetString()
        self.__OutputFEData(self.fe_model_part, filename, output_dir, nodal_variables)

        self.pole_nodes = self.__CreateNodesFromPolesForSolution(self.cad_model)

        print("\n> Starting creation of conditions...")
        start_time = time.time()

        self.conditions = self.condition_factory.CreateConditions(self.pole_nodes)

        print("> Finished creation of conditions in" ,round( time.time()-start_time, 3 ), " s.")
        print("\n> Initializing assembly...")
        start_time = time.time()

        self.system = eq.DSymmetricSystem(element_lists=self.conditions, linear_solver={'type': 'pardiso_ldlt'})

        print("> Initialization of assembly finished in" ,round( time.time()-start_time, 3 ), " s.")

    # --------------------------------------------------------------------------
    def Map(self):
        # Nonlinear solution iterations
        for solution_itr in range(1,self.parameters["solution"]["iterations"].GetInt()+1):

            print("\n> -----------------------------------------------------------------")
            print("> Starting solution iteration", solution_itr,"...")
            start_time_iteration = time.time()

            # Assemble considering scaling of constraints
            print("\n> Starting to compute RHS and LHS ....")
            self.__AssembleSystem()
            print("> Finished computing RHs and LHS in" ,round( time.time()-start_time_iteration, 3 ), " s.")

            if solution_itr == 1:
                total_num_conditions = 0
                for condition_list in self.conditions:
                    total_num_conditions += len(condition_list)
                print("\n> Number of conditions =",total_num_conditions)
                print("> Number of unkowns =", self.system.nb_dofs)

            print("\n> Starting system solution ....")
            start_time_solution = time.time()

            # Solve system and add delta to eqlib solution vector
            rhs = self.system.g
            solution = self.system.h_inv_v(rhs)
            self.system.x += solution

            # Store system data if specified
            if self.parameters["output"]["echo_level"].GetInt() > 5:
                from scipy.io import mmwrite
                mmwrite("rhs.mtx",[rhs])
                mmwrite("lhs.mtx",self.system.h)
                mmwrite("solution.mtx",[solution])

            print("> Finished system solution in" ,round( time.time()-start_time_solution, 3 ), " s.")

            print("\n> Finalizing solution step....")
            start_time = time.time()

            self.__UpdateCADModel()

            filename = self.parameters["output"]["resulting_geometry_filename"].GetString() + "_itr_"+str(solution_itr)+".iga"
            self.__OutputCadModel(filename)

            print("> Finished finalization of solution step in" ,round( time.time()-start_time, 3 ), " s.")

            if self.parameters["solution"]["test_solution"].GetBool():

                # Check fitting system
                print("\n> Max absolute control point displacement = ",la.norm(solution, np.inf))
                # print("> Computing condition number of LHS...")
                # print("> Condition number LHS = ", la.cond(self.system.h.todense()))
                # print("\n> Max absolute value in LHS = ", abs(self.system.h).max())

                # Test solution quality
                test_rhs = self.system.h_v(solution)

                delta = rhs-test_rhs
                error_norm = la.norm(delta)
                print("\n> Error in linear solution = ",error_norm)

                # Test residuals
                error_norm = la.norm(rhs)
                print("> RHS before current solution iteration = ",error_norm)

                # Test rhs after update
                self.__AssembleSystemRHS()
                rhs = self.system.g

                error_norm = la.norm(rhs)
                print("\n> RHS after current solution iteration = ",error_norm)

                # Varying contribution of beta regularization is neglected as each solution iteration may be seen indendently
                # in terms of minimization of the control point displacement

            print("\n> Finished solution iteration in" ,round( time.time()-start_time_iteration, 3 ), " s.")
            print("> -----------------------------------------------------------------")

    # --------------------------------------------------------------------------
    def ResetDisplacementsAndRefineCadModel(self):
        fe_point_distance_tolerance = self.parameters["refinement"]["a_posteriori"]["fe_point_distance_tolerance"].GetDouble() # in given length unit
        disp_coupling_tolerance = self.parameters["refinement"]["a_posteriori"]["disp_coupling_tolerance"].GetDouble() # in given length unit
        rot_coupling_tolerance = self.parameters["refinement"]["a_posteriori"]["rot_coupling_tolerance"].GetDouble() # in degree
        minimum_knot_distance = self.parameters["refinement"]["a_posteriori"]["mininimum_knot_distance"].GetDouble()

        exist_intervals_to_refine = False
        is_fe_point_distance_satisfied = True
        is_disp_coupling_satisfied = True
        is_rot_coupling_satisfied = True

        # Identify spots to refine
        intervals_along_u_to_refine = {geometry.Key(): [] for geometry in self.cad_model.GetByType('SurfaceGeometry3D')}
        intervals_along_v_to_refine = {geometry.Key(): [] for geometry in self.cad_model.GetByType('SurfaceGeometry3D')}

        # First identify spots where distance to fe points exceeds limit
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

                if la.norm(distance) > fe_point_distance_tolerance:
                    u_added = self.__AddUIntervalToList(geometry, u, v, intervals_along_u_to_refine, minimum_knot_distance)
                    v_added = self.__AddVIntervalToList(geometry, u, v, intervals_along_v_to_refine, minimum_knot_distance)

                    exist_intervals_to_refine = exist_intervals_to_refine or u_added or v_added
                    is_fe_point_distance_satisfied = False

        # Output fitting error
        nodal_variables = ["SHAPE_CHANGE", "FITTING_ERROR"]
        filename = self.parameters["output"]["filename_fem_for_quality_evaluation"].GetString()
        output_dir = self.parameters["output"]["results_directory"].GetString()
        self.__OutputFEData(self.fe_model_part, filename, output_dir, nodal_variables)

        # Then identify spots where coupling or enforcement conditions are not met
        for condition_list in self.conditions:
            for condition in condition_list:
                if isinstance(condition, clib.DisplacementCouplingCondition) or isinstance(condition, clib.DisplacementCouplingConditionWithAD) or isinstance(condition, clib.RotationCouplingConditionWithAD):
                    geometry_a = condition.geometry_a
                    geometry_b = condition.geometry_b
                    (u_a,v_a) = condition.parameters_a
                    (u_b,v_b) = condition.parameters_b

                    is_spot_to_be_refined = False

                    if isinstance(condition, clib.DisplacementCouplingCondition) or isinstance(condition, clib.DisplacementCouplingConditionWithAD):
                        delta_disp = la.norm(condition.CalculateQualityIndicator())
                        if delta_disp > disp_coupling_tolerance:
                            is_disp_coupling_satisfied = False
                            is_spot_to_be_refined = True
                            print("delta_disp =", delta_disp)

                    if isinstance(condition, clib.RotationCouplingConditionWithAD):
                        delta_rot = condition.CalculateQualityIndicator() * 180 / np.pi
                        if delta_rot > rot_coupling_tolerance:
                            is_rot_coupling_satisfied = False
                            is_spot_to_be_refined = True
                            print("delta_rot =", delta_rot)

                    if is_spot_to_be_refined:
                        # Geometry a
                        u_a_added = self.__AddUIntervalToList(geometry_a, u_a, v_a, intervals_along_u_to_refine, minimum_knot_distance)
                        v_a_added = self.__AddVIntervalToList(geometry_a, u_a, v_a, intervals_along_v_to_refine, minimum_knot_distance)

                        # Geometry b
                        u_b_added = self.__AddUIntervalToList(geometry_b, u_b, v_b, intervals_along_u_to_refine, minimum_knot_distance)
                        v_b_added = self.__AddVIntervalToList(geometry_b, u_b, v_b, intervals_along_v_to_refine, minimum_knot_distance)

                        exist_intervals_to_refine = exist_intervals_to_refine or u_a_added or v_a_added or u_b_added or v_b_added

        self.ResetPoleDisplacements()

        # Perform actual refinement
        nothing_to_refine = False

        if is_fe_point_distance_satisfied and is_disp_coupling_satisfied and is_rot_coupling_satisfied:
            print("\n> Refinement reached convergence! Nothing to refine anymore.")
            nothing_to_refine = True

        elif exist_intervals_to_refine == False:
            print("\n> WARNING!!!! Refinement did not converge but finished as knot tolerance was reached for U- and V-Direction.")
            nothing_to_refine = True

        else:
            self.__RefineAtIntervalCenters(self.cad_model, intervals_along_u_to_refine, intervals_along_v_to_refine)

        return nothing_to_refine

    # --------------------------------------------------------------------------
    def ResetPoleDisplacements(self):
        for surface_key, pole_nodes in self.pole_nodes.items():
            surface_geometry = self.cad_model.Get(surface_key)
            surface_geometry_data = surface_geometry.Data()
            for i, pole_node in enumerate(pole_nodes):
                surface_geometry_data.SetPole(i, pole_node.ref_location)

    # --------------------------------------------------------------------------
    @staticmethod
    def __CreateNodesFromPolesForSolution(cad_model):
        pole_nodes = {}
        for surface_geometry in cad_model.GetByType('SurfaceGeometry3D'):
            surface = surface_geometry.Data()
            surface_dofs = [eq.Node(x, y, z) for x, y, z in surface.Poles()]
            pole_nodes[surface_geometry.Key()] = surface_dofs
        return pole_nodes

    # --------------------------------------------------------------------------
    def __AssembleSystem(self, only_rhs=False):
        is_initial_assembly = (len(self.conditions_scaling_factors) == 0)

        # Initialize
        self.system.set_zero()

        # Assemble reference system
        self.system.compute(index=0, parallel=self.parameters["solution"]["parallel_assembly"].GetBool())
        self.system.add_working_system(factor=1)
        self.conditions_scaling_factors.append(1)

        if is_initial_assembly:
            system_trace = self.system.working_trace()
            print("> Max absolute system LHS = " + str(abs(self.system.working_h).max()))
            print("> trace_0 = "+str(system_trace))

        # Beta regularization
        beta = self.parameters["regularization"]["beta"].GetDouble()

        if is_initial_assembly:
            trace_beta = self.system.nb_dofs
            self.beta_scaling_factor = system_trace/trace_beta
            print("> beta_scaling_factor = "+str(self.beta_scaling_factor)+"\n")

        self.system.add_diagonal(self.beta_scaling_factor*beta)

        # Assembly constraints considering scaling
        for itr in range(1,len(self.conditions)):
            self.system.compute(index=itr, parallel=self.parameters["solution"]["parallel_assembly"].GetBool())

            if is_initial_assembly:
                weight_of_element_type = self.conditions[itr][0].GetWeight()
                constraint_trace = self.system.working_trace() / weight_of_element_type
                scaling_factor = system_trace/constraint_trace

                self.conditions_scaling_factors.append(scaling_factor)

                print("> trace_"+str(itr)+" = "+str(constraint_trace))
                print("> scaling_factor_"+str(itr)+" = "+str(scaling_factor))

                max_value_constraint = abs(self.system.working_h).max() / weight_of_element_type
                print("> Max absolute value constraint LHS_"+str(itr)+"= " + str(max_value_constraint)+"\n")

            self.system.add_working_system(factor=self.conditions_scaling_factors[itr])

 # --------------------------------------------------------------------------
    def __AssembleSystemRHS(self):
        if len(self.conditions_scaling_factors) == 0:
            raise RuntimeError("__AssembleSystemRHS: No scaling factors computed yet!")

        # Initialize
        self.system.set_zero()

        # Assemble reference system
        self.system.compute(index=0, order=1, parallel=self.parameters["solution"]["parallel_assembly"].GetBool())
        self.system.add_working_system(factor=1)

        # Assembly constraints considering scaling
        for itr in range(1,len(self.conditions)):
            self.system.compute(index=itr, order=1, parallel=self.parameters["solution"]["parallel_assembly"].GetBool())
            self.system.add_working_system(factor=self.conditions_scaling_factors[itr])

    # --------------------------------------------------------------------------
    def __UpdateCADModel(self):
        # Add control point displacement as field
        for face_i in self.cad_model.GetByType('BrepFace'):
            surface_key = face_i.Data().Geometry().Key()

            field = an.BrepFaceField(3)
            field.SetFace(face_i)
            for i, pole_node in enumerate(self.pole_nodes[surface_key]):
                displacement = pole_node.act_location - pole_node.ref_location
                field.SetValue(i, displacement)

            self.cad_model.Add(field)

        # Update control point position
        for surface_key, pole_nodes in self.pole_nodes.items():
            surface_geometry = self.cad_model.Get(surface_key)
            surface_geometry_data = surface_geometry.Data()
            for i, pole_node in enumerate(pole_nodes):
                # print('displacement =', pole_node.act_location - pole_node.ref_location)
                surface_geometry_data.SetPole(i, pole_node.act_location) # act_location considers the solution stored on self.system.x

    # --------------------------------------------------------------------------
    def __OutputCadModel(self, filename):
        output_dir = self.parameters["output"]["results_directory"].GetString()
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        output_filename_with_path = os.path.join(output_dir,filename)
        self.cad_model.Save(output_filename_with_path)

    # --------------------------------------------------------------------------
    def __StoreCaseFilesInResults(self):
        output_dir = self.parameters["output"]["results_directory"].GetString()
        input_cad_filename = self.parameters["input"]["cad_filename"].GetString()
        input_fem_filename = self.parameters["input"]["fem_filename"].GetString()

        import shutil, __main__
        shutil.copy(input_cad_filename, output_dir)
        shutil.copy(input_fem_filename, output_dir)
        shutil.copy(os.path.basename(__main__.__file__), output_dir)

    # --------------------------------------------------------------------------
    @staticmethod
    def __OutputFEData(model_part, fem_output_filename, output_dir, nodal_variables):
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        fem_output_filename_with_path = os.path.join(output_dir,fem_output_filename)

        from gid_output import GiDOutput
        nodal_results=nodal_variables
        gauss_points_results=[]
        VolumeOutput = True
        GiDPostMode = "Binary"
        GiDWriteMeshFlag = False
        GiDWriteConditionsFlag = True
        GiDMultiFileFlag = "Single"

        gid_io = GiDOutput(fem_output_filename_with_path, VolumeOutput, GiDPostMode, GiDMultiFileFlag, GiDWriteMeshFlag, GiDWriteConditionsFlag)
        gid_io.initialize_results(model_part)
        gid_io.write_results(1, model_part, nodal_results, gauss_points_results)
        gid_io.finalize_results()

    # --------------------------------------------------------------------------
    @staticmethod
    def __AddUIntervalToList(geometry, u, v, interval_list, mininimum_knot_distance):
        intervall_added = False

        start_span = an.Knots.UpperSpan(geometry.Data().DegreeU(), geometry.Data().KnotsU(), u)
        end_span = start_span+1
        u_start = geometry.Data().KnotsU()[start_span]
        u_end = geometry.Data().KnotsU()[end_span]

        # Avoid refining knots that are close together within a specified tolerance
        # comparisons are made in given length unit
        p_start = geometry.Data().PointAt(u_start, v)
        p_end = geometry.Data().PointAt(u_end, v)
        distance = p_end - p_start

        if np.dot(distance, distance) >= mininimum_knot_distance**2 :
            intervall_added = True
            u_interval = (u_start, u_end)
            interval_list[geometry.Key()].append(u_interval)

        return intervall_added

    # --------------------------------------------------------------------------
    @staticmethod
    def __AddVIntervalToList(geometry, u , v, interval_list, minimum_knot_distance):
        intervall_added = False

        start_span = an.Knots.UpperSpan(geometry.Data().DegreeV(), geometry.Data().KnotsV(), v)
        end_span = start_span+1
        v_start = geometry.Data().KnotsV()[start_span]
        v_end = geometry.Data().KnotsV()[end_span]

        # Avoid refining knots that are close together within a specified tolerance
        # comparisons are made in given length unit
        p_start = geometry.Data().PointAt(u, v_start)
        p_end = geometry.Data().PointAt(u, v_end)
        distance = p_end - p_start

        if np.dot(distance, distance) >= minimum_knot_distance**2 :
            intervall_added = True
            v_interval = (v_start, v_end)
            interval_list[geometry.Key()].append(v_interval)

        return intervall_added

    # --------------------------------------------------------------------------
    @staticmethod
    def __RefineAtIntervalCenters(cad_model, intervals_along_u_to_refine, intervals_along_v_to_refine):
        for face in cad_model.GetByType('BrepFace'):
                # print("> Refining face ",face.Key())
                geometry = face.Data().Geometry()

                u_intervals = intervals_along_u_to_refine[geometry.Key()]
                if len(u_intervals) > 0:
                    u_interval_centers = [(interval[0]+interval[1])/2.0 for interval in u_intervals]
                    u_interval_centers = list(set(u_interval_centers)) # Remove duplicated entries
                    u_interval_centers.sort() # Sort in ascending order

                    refined_geometry = an.KnotRefinement.InsertKnotsU(geometry.Data(), u_interval_centers)
                    cad_model.Replace(geometry.Key(), refined_geometry)
                # else:
                #     print("> Nothing to refine in U-driection!")

                v_intervals = intervals_along_v_to_refine[geometry.Key()]
                if len(v_intervals) > 0:
                    v_interval_centers = [(interval[0]+interval[1])/2.0 for interval in v_intervals]
                    v_interval_centers = list(set(v_interval_centers)) # Remove duplicated entries
                    v_interval_centers.sort() # Sort in ascending order

                    refined_geometry = an.KnotRefinement.InsertKnotsV(geometry.Data(), v_interval_centers)
                    cad_model.Replace(geometry.Key(), refined_geometry)
                # else:
                #     print("> Nothing to refine in V-driection!")

# ==============================================================================