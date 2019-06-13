//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics FemDem Application
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Sergio Jimenez Reyes
//

#include "custom_processes/advance_in_time_strategy_high_cycle_fatigue_process.h"

namespace Kratos
{

/***********************************************************************************/
/***********************************************************************************/

// Constructor
AdvanceInTimeStrategyHighCycleFatigueProcess::AdvanceInTimeStrategyHighCycleFatigueProcess(
    ModelPart &rModelPart, 
    Parameters ThisParameters) : mrModelPart(rModelPart)
{
    // mNNodes = mrModelPart.NumberOfNodes();
    // Parameters default_parameters = Parameters(R"(
    // {
    //     "normalized_free_energy" : false,
    //     "correct_with_displacements" : true,
    //     "correction_factor"  : 1.0 
    // })" );
    // ThisParameters.RecursivelyValidateAndAssignDefaults(default_parameters);

    // mComputeNormalizedFreeEnergy = ThisParameters["normalized_free_energy"].GetBool();
    // mCorrectWithDisplacements = ThisParameters["correct_with_displacements"].GetBool();
    // mCorrectionFactor = ThisParameters["correction_factor"].GetDouble();
}

/***********************************************************************************/
/***********************************************************************************/

void AdvanceInTimeStrategyHighCycleFatigueProcess::Execute()
{
    // auto& r_mat_properties = mrModelPart.Elements().GetProperties();
    auto& process_info = mrModelPart.GetProcessInfo();
    std::vector<int> cycle_identificator;

    mrModelPart.GetElement(1).GetValueOnIntegrationPoints(CYCLE_INDICATOR, cycle_identificator, process_info);
    KRATOS_WATCH(cycle_identificator);
    
    if (cycle_identificator[1] == 1) {
        bool advancing_strategy = false;
        this->StableConditionForAdvancingStrategy(advancing_strategy);

        if (advancing_strategy) {
        }
    }
    // NodeNormalizedFreeEnergy *NodeNormalizedFreeEnergyVector = new NodeNormalizedFreeEnergy[max_id];
    // this->NormalizedFreeEnergyExtrapolation(NodeNormalizedFreeEnergyVector);
    // delete[] NodeNormalizedFreeEnergyVector;
}

/***********************************************************************************/
/***********************************************************************************/

// void AdvanceInTimeStrategyHighCycleFatigueProcess::NormalizedFreeEnergyExtrapolation(
//     NodeNormalizedFreeEnergy* pNodeNormalizedFreeEnergyVector
//     )
// {
//     // Vector gauss_point_stress;
//     // double gauss_point_damage;
//     // // Loop over elements to extrapolate the stress to the nodes
    // // for (auto& r_elem : mrModelPart.Elements()) {
//     //     auto& r_geometry = r_elem.GetGeometry();
//     //     auto& r_mat_properties = r_elem.GetProperties();

//     //     const Vector& r_strain_vector = r_elem.GetValue(STRAIN_VECTOR);
//     //     const Vector& r_stress_vector = r_elem.GetValue(STRESS_VECTOR);
//     //     const double damage = r_elem.GetValue(DAMAGE_ELEMENT);

//     //     // Compute Normalized Free Energy on that element
        // // const double normalized_free_energy = this->StableConditionForAdvancingStrategy(r_strain_vector, r_stress_vector, damage, r_mat_properties, r_geometry);

//     //     for (unsigned int i = 0; i < r_geometry.PointsNumber(); i++) {
//     //         const std::size_t node_id = r_geometry.GetPoint(i).Id();
//     //         pNodeNormalizedFreeEnergyVector[node_id - 1].NormalizedFreeEnergy += normalized_free_energy;
//     //         pNodeNormalizedFreeEnergyVector[node_id - 1].NElems += 1;
//     //     }
//     // }

//     // // Ponderate over the elements coincident on that node
//     // for (unsigned int i = 0; i < mNNodes; i++) {
//     //     pNodeNormalizedFreeEnergyVector[i].NormalizedFreeEnergy /= pNodeNormalizedFreeEnergyVector[i].NElems;
//     // }

//     // // Loop over nodes to assign the variable
//     // for (auto& r_node : mrModelPart.Nodes()) {
//     //     const std::size_t Id = r_node.Id();
//     //     const double nodal_free_energy = pNodeNormalizedFreeEnergyVector[Id - 1].NormalizedFreeEnergy;
//     //     double& r_norm = r_node.FastGetSolutionStepValue(EQUIVALENT_NODAL_STRESS);
//     //     r_norm = nodal_free_energy;

//     //     if (mCorrectWithDisplacements) {
//     //         const double displ = norm_2(r_node.FastGetSolutionStepValue(DISPLACEMENT));
//     //         r_norm *= displ;
//     //     }
//     //     r_norm *= mCorrectionFactor;
//     // }
// }

/***********************************************************************************/
/***********************************************************************************/

void AdvanceInTimeStrategyHighCycleFatigueProcess::StableConditionForAdvancingStrategy(bool& rAdvancingStrategy) 
{
    auto& process_info = mrModelPart.GetProcessInfo();
    std::vector<double> max_stress_rel_error;
    std::vector<double> rev_factor_rel_error;
    
    double acumulated_max_stress_rel_error;
    double acumulated_rev_factor_rel_error;
    for (auto& r_elem : mrModelPart.Elements()) {
        acumulated_max_stress_rel_error = 0.0;
        acumulated_rev_factor_rel_error = 0.0;
        r_elem.GetValueOnIntegrationPoints(MAX_STRESS_RELATIVE_ERROR, max_stress_rel_error, process_info);
        r_elem.GetValueOnIntegrationPoints(REVERSION_FACTOR_RELATIVE_ERROR, rev_factor_rel_error, process_info);

        KRATOS_WATCH(max_stress_rel_error);
        KRATOS_WATCH(rev_factor_rel_error);
        for (unsigned int i=0; i < max_stress_rel_error.size(); i++) {
            acumulated_max_stress_rel_error += max_stress_rel_error[i];
            acumulated_rev_factor_rel_error += rev_factor_rel_error[i];
        }
        if (acumulated_max_stress_rel_error < 1e-4 && acumulated_rev_factor_rel_error < 1e-4) {
            rAdvancingStrategy = true;
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

// void AdvanceInTimeStrategyHighCycleFatigueProcess::AdvancingStrategyStep(double& rTimeStep) 
// {

// }


/***********************************************************************************/
/***********************************************************************************/

double AdvanceInTimeStrategyHighCycleFatigueProcess::CalculateCharacteristicLength2D(
    const Geometry<Node<3>>& rGeometry
    )
{
    // const auto& r_node_1_coordinates = rGeometry[0].Coordinates();
    // const auto& r_node_2_coordinates = rGeometry[1].Coordinates();
    // const auto& r_node_3_coordinates = rGeometry[2].Coordinates();

    // const double length1 = norm_2(r_node_1_coordinates-r_node_2_coordinates);
    // const double length2 = norm_2(r_node_2_coordinates-r_node_3_coordinates);
    // const double length3 = norm_2(r_node_3_coordinates-r_node_1_coordinates);
    // return (length1 + length2 + length3) / 3.0;
    return 1.0;
}

/***********************************************************************************/
/***********************************************************************************/

double AdvanceInTimeStrategyHighCycleFatigueProcess::CalculateCharacteristicLength3D(
    Geometry<Node<3>>& rGeometry
    )
{
    // Vector lengths = ZeroVector(6);
    // auto r_edges = rGeometry.Edges();
    // for (unsigned int edge = 0; edge < 6; edge++) { // Loop over edges
    //     const auto& r_coords_node_1 = r_edges[edge][0].Coordinates(); 
    //     const auto& r_coords_node_2 = r_edges[edge][1].Coordinates(); 
    //     lengths[edge] = std::sqrt(std::pow((r_coords_node_1[0] - r_coords_node_2[0]), 2.0) + 
    //         std::pow((r_coords_node_1[1] - r_coords_node_2[1]), 2.0) + std::pow((r_coords_node_1[2] - r_coords_node_2[2]), 2.0));
    // }
    // return (lengths[0] + lengths[1] + lengths[2] + lengths[3] + lengths[4] + lengths[5]) / 6.0;
    return 1.0;
}

/***********************************************************************************/
/***********************************************************************************/

double AdvanceInTimeStrategyHighCycleFatigueProcess::ComputeTensionFactor2D(
    const Vector& rStressVector
    )
{
    // Vector principal_stress_vector;
    // this->CalculatePrincipalStresses2D(rStressVector, principal_stress_vector);
    // double sum_a = 0.0, sum_b = 0.0, sum_c = 0.0;
    // for (unsigned int cont = 0; cont < 2; cont++) {
    //     sum_a += std::abs(principal_stress_vector[cont]);
    //     sum_b += 0.5 * (principal_stress_vector[cont]  + std::abs(principal_stress_vector[cont]));
    //     sum_c += 0.5 * (-principal_stress_vector[cont] + std::abs(principal_stress_vector[cont]));
    // }
    // if (sum_a > tolerance)
    //     return sum_b / sum_a;
    // else 
    //     return 0.0;
    return 1.0;
}

/***********************************************************************************/
/***********************************************************************************/

double AdvanceInTimeStrategyHighCycleFatigueProcess::ComputeTensionFactor3D(
    const Vector& rStressVector
    )
{
    // Vector principal_stress_vector;
    // this->CalculatePrincipalStresses3D(rStressVector, principal_stress_vector);
    // double sum_a = 0.0, sum_b = 0.0, sum_c = 0.0;
    // for (unsigned int cont = 0; cont < 3; cont++) {
    //     sum_a += std::abs(principal_stress_vector[cont]);
    //     sum_b += 0.5 * (principal_stress_vector[cont]  + std::abs(principal_stress_vector[cont]));
    //     sum_c += 0.5 * (-principal_stress_vector[cont] + std::abs(principal_stress_vector[cont]));
    // }
    // if (sum_a > tolerance)
    //     return sum_b / sum_a;
    // else 
    //     return 0.0;
    return 1.0;
}

/***********************************************************************************/
/***********************************************************************************/

void AdvanceInTimeStrategyHighCycleFatigueProcess::CalculatePrincipalStresses2D(
    const Vector& rStressVector, 
    Vector& rPrincipalStressVector
    )
{
    // rPrincipalStressVector.resize(2);
    // rPrincipalStressVector[0] = 0.5 * (rStressVector[0] + rStressVector[1]) + 
    //     std::sqrt(std::pow(0.5 * (rStressVector[0] - rStressVector[1]), 2) + std::pow(rStressVector[2], 2));
    // rPrincipalStressVector[1] = 0.5 * (rStressVector[0] + rStressVector[1]) - 
    // std::sqrt(std::pow(0.5 * (rStressVector[0] - rStressVector[1]), 2) + std::pow(rStressVector[2], 2));
}

/***********************************************************************************/
/***********************************************************************************/

void AdvanceInTimeStrategyHighCycleFatigueProcess::CalculatePrincipalStresses3D(
    const Vector& rStressVector, 
    Vector& rPrincipalStressVector
    )
{
    // rPrincipalStressVector.resize(3);
    // const double I1 = this->CalculateI1Invariant(rStressVector);
    // const double I2 = this->CalculateI2Invariant(rStressVector);
    // const double I3 = this->CalculateI3Invariant(rStressVector);
    // const double II1 = I1 * I1;

    // const double numerator = (2.0 * II1 - 9.0 * I2) * I1 + 27.0 * I3;
    // const double denominator = (II1 - 3.0 * I2);

    // if (denominator != 0.0) {
    //     double phi = numerator / (2.0 * denominator * std::sqrt(denominator));

    //     if (std::abs(phi) > 1.0) {
    //         if (phi > 0.0)
    //             phi = 1.0;
    //         else
    //             phi = -1.0;
    //     }

    //     const double acosphi = std::acos(phi);
    //     phi = acosphi / 3.0;

    //     const double aux1 = 2.0 / 3.0 * std::sqrt(II1 - 3.0 * I2);
    //     const double aux2 = I1 / 3.0;

    //     rPrincipalStressVector[0] = aux2 + aux1 * std::cos(phi);
    //     rPrincipalStressVector[1] = aux2 + aux1 * std::cos(phi - 2.09439510239);
    //     rPrincipalStressVector[2] = aux2 + aux1 * std::cos(phi - 4.18879020478);
    // } else {
    //     rPrincipalStressVector = ZeroVector(3);
    // }
}

/***********************************************************************************/
/***********************************************************************************/

void AdvanceInTimeStrategyHighCycleFatigueProcess::ObtainMaximumNodeId(
    std::size_t& rMaxId
    )
{
    // std::size_t aux = 0;
    // std::size_t id;

    // for (auto& r_node : mrModelPart.Nodes()) {
    //     id = r_node.Id();
    //     if (id > aux)
    //         aux = id;
    // }
    // rMaxId = aux;
}

/***********************************************************************************/
/***********************************************************************************/

double AdvanceInTimeStrategyHighCycleFatigueProcess::CalculateI1Invariant(
    const Vector& rStressVector
    )
{
    // return rStressVector[0] + rStressVector[1] + rStressVector[2];
    return 1.0;
}

/***********************************************************************************/
/***********************************************************************************/

double AdvanceInTimeStrategyHighCycleFatigueProcess::CalculateI2Invariant(
    const Vector& rStressVector
    )
{
    // return (rStressVector[0] + rStressVector[2]) * rStressVector[1] + rStressVector[0] * rStressVector[2] +
    //     -rStressVector[3] * rStressVector[3] - rStressVector[4] * rStressVector[4] - rStressVector[5] * rStressVector[5];
    return 1.0;
}

/***********************************************************************************/
/***********************************************************************************/

double AdvanceInTimeStrategyHighCycleFatigueProcess::CalculateI3Invariant(
    const Vector& rStressVector
    )
{
    // return (rStressVector[1] * rStressVector[2] - rStressVector[4] * rStressVector[4]) * rStressVector[0] -
    //     rStressVector[1] * rStressVector[5] * rStressVector[5] - rStressVector[2] * rStressVector[3] * rStressVector[3] +
    //     2.0 * rStressVector[3] * rStressVector[4] * rStressVector[5];
    return 1.0;
}

} // namespace Kratos
