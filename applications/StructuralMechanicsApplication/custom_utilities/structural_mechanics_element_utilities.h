//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher
//

#if !defined( KRATOS_STRUCTURAL_MECHANICS_ELEMENT_UTILITIES_H_INCLUDED )
#define  KRATOS_STRUCTURAL_MECHANICS_ELEMENT_UTILITIES_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/element.h"

namespace Kratos {
namespace StructuralMechanicsElementUtilities {

void BaseSolidElement::EquationIdVectorDisplacement(
    /*const*/ GeometryType& rGeometry,
    Element::EquationIdVectorType& rEquationIdVector,
    const ProcessInfo& rCurrentProcessInfo,
    const SizeType NumberOfNodes,
    const SizeType Dimension
    )
{
    KRATOS_TRY;

    if (rEquationIdVector.size() != Dimension*NumberOfNodes) {
        rEquationIdVector.resize(Dimension*NumberOfNodes,false);
    }

    Check in FluidElement how it is done better!

    const SizeType pos = rGeometry[0].GetDofPosition(DISPLACEMENT_X);

    if(Dimension == 2) {
        for (IndexType i = 0; i < NumberOfNodes; ++i) {
            const SizeType index = i * 2;
            rEquationIdVector[index]   = rGeometry[i].GetDof(DISPLACEMENT_X, pos).EquationId();
            rEquationIdVector[index+1] = rGeometry[i].GetDof(DISPLACEMENT_Y, pos+1).EquationId();
        }
    } else {
        for (IndexType i = 0; i < NumberOfNodes; ++i) {
            const SizeType index = i * 3;
            rEquationIdVector[index]   = rGeometry[i].GetDof(DISPLACEMENT_X, pos).EquationId();
            rEquationIdVector[index+1] = rGeometry[i].GetDof(DISPLACEMENT_Y, pos+1).EquationId();
            rEquationIdVector[index+2] = rGeometry[i].GetDof(DISPLACEMENT_Z, pos+2).EquationId();
        }
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::GetDofListDisplacement(
    Element::DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo,
    const SizeType NumberOfNodes,
    const SizeType Dimension
    )
{
    KRATOS_TRY;

    rElementalDofList.resize(0);
    rElementalDofList.reserve(Dimension*NumberOfNodes);

    if(Dimension == 2) {
        for (IndexType i = 0; i < NumberOfNodes; ++i) {
            rElementalDofList.push_back(rGeometry[i].pGetDof(DISPLACEMENT_X));
            rElementalDofList.push_back(rGeometry[i].pGetDof(DISPLACEMENT_Y));
        }
    } else {
        for (IndexType i = 0; i < NumberOfNodes; ++i) {
            rElementalDofList.push_back(rGeometry[i].pGetDof(DISPLACEMENT_X));
            rElementalDofList.push_back(rGeometry[i].pGetDof(DISPLACEMENT_Y));
            rElementalDofList.push_back(rGeometry[i].pGetDof(DISPLACEMENT_Z));
        }
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::GetValuesVectorDisplacement(
    Vector& rValues,
    int Step
    )
{
    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    const SizeType mat_size = number_of_nodes * dimension;
    if (rValues.size() != mat_size)
        rValues.resize(mat_size, false);
    for (IndexType i = 0; i < number_of_nodes; ++i)
    {
        const array_1d<double, 3 >& displacement = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT, Step);
        const SizeType index = i * dimension;
        for(unsigned int k = 0; k < dimension; ++k)
        {
            rValues[index + k] = displacement[k];
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::GetFirstDerivativesVectorDisplacement(
    Vector& rValues,
    int Step
    )
{
    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    const SizeType mat_size = number_of_nodes * dimension;
    if (rValues.size() != mat_size)
        rValues.resize(mat_size, false);
    for (IndexType i = 0; i < number_of_nodes; ++i) {
        const array_1d<double, 3 >& velocity = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, Step);
        const SizeType index = i * dimension;
        for(unsigned int k = 0; k < dimension; ++k)
            rValues[index + k] = velocity[k];
    }
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::GetSecondDerivativesVectorDisplacement(
    Vector& rValues,
    int Step
    )
{
    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    const SizeType mat_size = number_of_nodes * dimension;
    if (rValues.size() != mat_size)
        rValues.resize(mat_size, false);
    for (IndexType i = 0; i < number_of_nodes; ++i) {
        const array_1d<double, 3 >& acceleration = GetGeometry()[i].FastGetSolutionStepValue(ACCELERATION, Step);
        const SizeType index = i * dimension;
        for(unsigned int k = 0; k < dimension; ++k)
            rValues[index + k] = acceleration[k];
    }
}

bool ComputeLumpedMassMatrix(
    const Properties rProperites,
    const ProcessInfo& rCurrentProcessInfo);

void CalculateRayleighDampingMatrix(
    Element& rElement,
    Element::MatrixType& rDampingMatrix,
    /*const*/ ProcessInfo& rCurrentProcessInfo,
    const std::size_t MatrixSize);

bool HasRayleighDamping(
    const Properties rProperites,
    const ProcessInfo& rCurrentProcessInfo);

double GetRayleighAlpha(
    const Properties rProperites,
    const ProcessInfo& rCurrentProcessInfo);

double GetRayleighBeta(
    const Properties rProperites,
    const ProcessInfo& rCurrentProcessInfo);

    double CrBeamElement3D2N::CalculateShearModulus() {
  KRATOS_TRY;
  const double nu = this->GetProperties()[POISSON_RATIO];
  const double E = this->GetProperties()[YOUNG_MODULUS];
  const double G = E / (2.0 * (1.0 + nu));
  return G;
  KRATOS_CATCH("")

  void BaseSolidElement::AddExplicitContribution(
    const VectorType& rRHSVector,
    const Variable<VectorType>& rRHSVariable,
    Variable<double>& rDestinationVariable,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    auto& r_geom = this->GetGeometry();
    const SizeType dimension = r_geom.WorkingSpaceDimension();
    const SizeType number_of_nodes = r_geom.size();
    const SizeType mat_size = number_of_nodes * dimension;

    // Compiting the nodal mass
    if (rDestinationVariable == NODAL_MASS ) {
        VectorType element_mass_vector(mat_size);
        this->CalculateLumpedMassVector(element_mass_vector);

        for (IndexType i = 0; i < number_of_nodes; ++i) {
            const IndexType index = i * dimension;

            #pragma omp atomic
            r_geom[i].GetValue(NODAL_MASS) += element_mass_vector[index];
        }
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void BaseSolidElement::AddExplicitContribution(
    const VectorType& rRHSVector,
    const Variable<VectorType>& rRHSVariable,
    Variable<array_1d<double, 3>>& rDestinationVariable,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    auto& r_geom = this->GetGeometry();
    const auto& r_prop = this->GetProperties();
    const SizeType dimension = r_geom.WorkingSpaceDimension();
    const SizeType number_of_nodes = r_geom.size();
    const SizeType element_size = dimension * number_of_nodes;

    Vector damping_residual_contribution = ZeroVector(element_size);

    // Calculate damping contribution to residual -->
    if (r_prop.Has(RAYLEIGH_ALPHA) || r_prop.Has(RAYLEIGH_BETA)) {
        Vector current_nodal_velocities = ZeroVector(element_size);
        this->GetFirstDerivativesVector(current_nodal_velocities);

        Matrix damping_matrix = ZeroMatrix(element_size, element_size);
        this->CalculateDampingMatrixWithLumpedMass(damping_matrix, rCurrentProcessInfo);

        // Current residual contribution due to damping
        noalias(damping_residual_contribution) = prod(damping_matrix, current_nodal_velocities);
    }

    // Computing the force residual
    if (rRHSVariable == RESIDUAL_VECTOR && rDestinationVariable == FORCE_RESIDUAL) {
        for (IndexType i = 0; i < number_of_nodes; ++i) {
            const IndexType index = dimension * i;

            array_1d<double, 3>& r_force_residual = r_geom[i].FastGetSolutionStepValue(FORCE_RESIDUAL);

            for (IndexType j = 0; j < dimension; ++j) {
                #pragma omp atomic
                r_force_residual[j] += rRHSVector[index + j] - damping_residual_contribution[index + j];
            }
        }
    }

    KRATOS_CATCH("")
}

}

    // ?????
    if (compute_lumped_mass_matrix) {
        VectorType temp_vector(mat_size);
        CalculateLumpedMassVector(temp_vector);
        for (IndexType i = 0; i < mat_size; ++i)
            rMassMatrix(i, i) = temp_vector[i];


double CrBeamElement3D2N::CalculateReferenceLength(const GeometryType& rGeometry) {

  KRATOS_TRY;
  const double dx = this->GetGeometry()[1].X0() - this->GetGeometry()[0].X0();
  const double dy = this->GetGeometry()[1].Y0() - this->GetGeometry()[0].Y0();
  const double dz = this->GetGeometry()[1].Z0() - this->GetGeometry()[0].Z0();
  const double L = std::sqrt(dx * dx + dy * dy + dz * dz);
  return L;

    Add check if length == 0

  KRATOS_CATCH("")
}

double CrBeamElement3D2N::CalculateCurrentLength(const GeometryType& rGeometry) {

  KRATOS_TRY;
  const double du =
      this->GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT_X) -
      this->GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_X);
  const double dv =
      this->GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT_Y) -
      this->GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_Y);
  const double dw =
      this->GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT_Z) -
      this->GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT_Z);
  const double dx = this->GetGeometry()[1].X0() - this->GetGeometry()[0].X0();
  const double dy = this->GetGeometry()[1].Y0() - this->GetGeometry()[0].Y0();
  const double dz = this->GetGeometry()[1].Z0() - this->GetGeometry()[0].Z0();
  const double l = std::sqrt((du + dx) * (du + dx) + (dv + dy) * (dv + dy) +
                             (dw + dz) * (dw + dz));
  return l;
  KRATOS_CATCH("")
}



} // namespace StructuralMechanicsElementUtilities.
}  // namespace Kratos.

#endif // KRATOS_STRUCTURAL_MECHANICS_ELEMENT_UTILITIES_H_INCLUDED  defined
