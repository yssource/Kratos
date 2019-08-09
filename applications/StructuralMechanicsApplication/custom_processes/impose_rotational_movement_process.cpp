// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:             BSD License
//                                       license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "utilities/openmp_utils.h"
#include "custom_processes/impose_rotational_movement_process.h"
#include "utilities/transformation_utilities.h"

namespace Kratos
{
ImposeRotationalMovementProcess::ImposeRotationalMovementProcess(
    ModelPart& rThisModelPart,
    Parameters ThisParameters
    ):mrThisModelPart(rThisModelPart),
      mThisParameters(ThisParameters)
{
    Parameters default_parameters = Parameters(R"(
    {
        "model_part_name"             : "please_specify_model_part_name",
        "new_model_part_name"         : "Rigid_Movement_ModelPart",
        "transformation_settings"     : {
            "rotation_settings"         : {
                "automatic_center"         : false,
                "master_node_id"           : 0,
                "center"                   : [0,0,0],
                "axis_of_rotation"         : [0.0,0.0,0.0],
                "constrained"              : false,
                "angle_degree"             : "please give an expression in terms of the variable x, y, z, t",
                "local_axes"               : {}
            }
        }
    })" );

    mThisParameters.ValidateAndAssignDefaults(default_parameters);

    // We define the function
    mpFunction = Kratos::make_shared<PythonGenericFunctionUtility>(mThisParameters["transformation_settings"]["rotation_settings"]["angle_degree"].GetString(),  mThisParameters["transformation_settings"]["rotation_settings"]["local_axes"]);
}

/***********************************************************************************/
/***********************************************************************************/

void ImposeRotationalMovementProcess::Execute()
{
    KRATOS_TRY

    // We execute the different steps
    ExecuteInitialize();
    ExecuteInitializeSolutionStep();

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void ImposeRotationalMovementProcess::ExecuteInitialize()
{
    KRATOS_TRY

    // Getting model parts
    ModelPart& root_model_part = mrThisModelPart.GetRootModelPart();
    KRATOS_ERROR_IF_NOT(root_model_part.HasNodalSolutionStepVariable(ROTATION)) << "Please add rotation variable to your model part" << std::endl;
    ModelPart& model_part = root_model_part.GetSubModelPart(mThisParameters["model_part_name"].GetString());
    const std::string& new_model_part_name = mThisParameters["new_model_part_name"].GetString();
    ModelPart& r_rotational_model_part = new_model_part_name != model_part.Name() ? model_part.HasSubModelPart(new_model_part_name) ? model_part.GetSubModelPart(new_model_part_name) : model_part.CreateSubModelPart(new_model_part_name) : model_part;

    // Reorder nodes
    IndexType node_id = 1;
    for (auto& r_node : root_model_part.Nodes()) {
        r_node.SetId(node_id);
        ++node_id;
    }
    // Reorder constrains
    IndexType constraint_id = 1;
    for (auto& r_constrain : root_model_part.MasterSlaveConstraints()) {
        r_constrain.SetId(constraint_id);
        ++constraint_id;
    }

    // Getting index of the master node
    const int master_node_id = mThisParameters["master_node_id"].GetInt();

    // Commpute the center of the geometry
    if (mThisParameters["transformation_settings"]["rotation_settings"]["angle_degree"]["automatic_center"].GetBool()) {
        double total_area = 0.0;
        array_1d<double, 3> center = ZeroVector(3);
        if (mrThisModelPart.Conditions().size() > 0) {
            for (auto& r_cond : mrThisModelPart.Conditions()) {
                const auto& r_geometry = r_cond.GetGeometry();
                const double area = r_geometry.Area();
                total_area += area;
                noalias(center) += area * r_geometry.Center().Coordinates();
            }
            center /= total_area;
        } else { // Average nodes instead
            for (auto& r_node : mrThisModelPart.Nodes()) {
                noalias(center) += r_node.Coordinates();
            }
            center /= static_cast<double>(mrThisModelPart.Nodes().size());
        }
        mpMasterNode = mrThisModelPart.CreateNewNode(node_id, center[0], center[1], center[2]);
    } else  {
        // If we master node ID is zero then we create a new model part
        if (master_node_id == 0) {
            const auto& r_center = mThisParameters["transformation_settings"]["rotation_settings"]["center"].GetVector();
            mpMasterNode = mrThisModelPart.CreateNewNode(node_id, r_center[0], r_center[1], r_center[2]);
        } else {
            mpMasterNode = root_model_part.pGetNode(master_node_id);
        }
    }
    KRATOS_ERROR_IF_NOT(mpMasterNode->HasDofFor(ROTATION_X)) << "Please add rotation dofs to your model part" << std::endl;

    // The axis of rotation
    const auto& r_axis_of_rotation = mThisParameters["transformation_settings"]["rotation_settings"]["axis_of_rotation"].GetVector();

    // We iterate over the nodes of the rigid model part
    auto& r_nodes_array = r_rotational_model_part.Nodes();
    const int number_of_nodes = static_cast<int>(r_nodes_array.size());
    const auto it_node_begin = r_nodes_array.begin();

    // Reference constraint
    const auto& r_clone_constraint = KratosComponents<MasterSlaveConstraint>::Get("LinearMasterSlaveConstraint");

    #pragma omp parallel
    {
        ConstraintContainerType constraints_buffer;

        DofPointerVectorType master_dofs(1);
        master_dofs[0] = mpMasterNode->pGetDof(ROTATION_X);
        DofPointerVectorType slave_dofs(3);

        MatrixType relation_matrix(3, 1);
        VectorType constant_vector = ZeroVector(3);

        #pragma omp for
        for (int i = 0; i < number_of_nodes; ++i) {
            auto it_node = it_node_begin + i;
            if (it_node->Id() != mpMasterNode->Id()) {
                slave_dofs[0] = it_node->pGetDof(DISPLACEMENT_X);
                slave_dofs[1] = it_node->pGetDof(DISPLACEMENT_Y);
                slave_dofs[2] = it_node->pGetDof(DISPLACEMENT_Z);
                auto p_constraint = r_clone_constraint.Create(constraint_id + i, master_dofs, slave_dofs, relation_matrix, constant_vector);
                (constraints_buffer).insert((constraints_buffer).begin(), p_constraint);
            }
        }

        // We transfer
        #pragma omp critical
        {
            r_rotational_model_part.AddMasterSlaveConstraints(constraints_buffer.begin(),constraints_buffer.end());
            mrThisModelPart.AddMasterSlaveConstraints(constraints_buffer.begin(),constraints_buffer.end());
        }
    }

    // Fix translation movements
    mpMasterNode->Fix(DISPLACEMENT_X);
    mpMasterNode->Fix(DISPLACEMENT_Y);
    mpMasterNode->Fix(DISPLACEMENT_Z);

    KRATOS_CATCH("")
}


/***********************************************************************************/
/***********************************************************************************/

void ImposeRotationalMovementProcess::ExecuteInitializeSolutionStep()
{
    // If the movement is constrained we rotate the master node
    if (mThisParameters["transformation_settings"]["rotation_settings"]["constrained"].GetBool()) {
        mpMasterNode->Fix(ROTATION_X);
        const double time = mrThisModelPart.GetProcessInfo()[TIME];
        const double value = mpFunction->CallFunction(mpMasterNode->X(), mpMasterNode->Y(), mpMasterNode->Z(), time);
        mpMasterNode->FastGetSolutionStepValue(ROTATION_X) = value;
    }
}

// class ImposeRotationalMovementProcess
} // namespace Kratos
