// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "includes/checks.h"
#include "custom_elements/explicit_total_lagrangian_bbar.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"
#include "structural_mechanics_application_variables.h"
#include "custom_utilities/structural_mechanics_element_utilities.h"

namespace Kratos
{

template< const SizeType TDim, const SizeType TNumNodes>
Element::Pointer ExplicitTotalLagrangianBbar<TDim, TNumNodes>::Create( IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const
{
    return Kratos::make_shared<ExplicitTotalLagrangianBbar<TDim, TNumNodes>>( NewId, GetGeometry().Create( ThisNodes ), pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

template< const SizeType TDim, const SizeType TNumNodes>
Element::Pointer ExplicitTotalLagrangianBbar<TDim, TNumNodes>::Create( IndexType NewId,  GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties ) const
{
    return Kratos::make_shared<ExplicitTotalLagrangianBbar<TDim, TNumNodes>>( NewId, pGeom, pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

template< const SizeType TDim, const SizeType TNumNodes>
void ExplicitTotalLagrangianBbar<TDim, TNumNodes>::Initialize()
{
    KRATOS_TRY

    const GeometryType::IntegrationPointsArrayType& integration_points = this->GetGeometry().IntegrationPoints(this->GetIntegrationMethod());

    //Constitutive Law initialisation
    if ( mConstitutiveLawVector.size() != integration_points.size() )
        mConstitutiveLawVector.resize( integration_points.size() );

    InitializeMaterial();

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

template< const SizeType TDim, const SizeType TNumNodes>
void ExplicitTotalLagrangianBbar<TDim, TNumNodes>::InitializeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{
    KinematicVariables this_kinematic_variables = KinematicVariables();
    ConstitutiveVariables this_constitutive_variables = ConstitutiveVariables();

    // Getting geometry
    const GeometryType& r_geometry = this->GetGeometry();

    // Create constitutive law parameters:
    ConstitutiveLaw::Parameters Values(r_geometry,this->GetProperties(),rCurrentProcessInfo);

    // Set constitutive law flags:
    Flags& ConstitutiveLawOptions=Values.GetOptions();
    ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);
//
    Values.SetStrainVector(this_constitutive_variables.StrainVector);
    Values.SetStressVector(this_constitutive_variables.StressVector);
    Values.SetConstitutiveMatrix(this_constitutive_variables.D);

    // Reading integration points
    const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry.IntegrationPoints(this->GetIntegrationMethod());

    // Shape functions
    const Matrix& rNValues = r_geometry.ShapeFunctionsValues(this->GetIntegrationMethod());

    bool compute_initialization;
    for ( IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number ) {
        compute_initialization = mConstitutiveLawVector[point_number]->Has(REQUIRES_MATERIAL_INITIALIZATION) ? mConstitutiveLawVector[point_number]->GetValue(REQUIRES_MATERIAL_INITIALIZATION, compute_initialization) : true;

        if (compute_initialization) {
            // Compute element kinematics B, F, DN_DX ...
            CalculateKinematicVariables(this_kinematic_variables, point_number, this->GetIntegrationMethod(), rNValues);

            // Compute constitutive law variables
            SetConstitutiveVariables(this_kinematic_variables, this_constitutive_variables, Values, point_number, integration_points);

            // Call the constitutive law to update material variables
            mConstitutiveLawVector[point_number]->InitializeMaterialResponse(Values, ConstitutiveLaw::StressMeasure_PK2);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template< const SizeType TDim, const SizeType TNumNodes>
void ExplicitTotalLagrangianBbar<TDim, TNumNodes>::InitializeNonLinearIteration( ProcessInfo& rCurrentProcessInfo )
{
    // TODO: Add something if necessary
}

/***********************************************************************************/
/***********************************************************************************/

template< const SizeType TDim, const SizeType TNumNodes>
void ExplicitTotalLagrangianBbar<TDim, TNumNodes>::FinalizeNonLinearIteration( ProcessInfo& rCurrentProcessInfo )
{
    // TODO: Add something if necessary
}

/***********************************************************************************/
/***********************************************************************************/

template< const SizeType TDim, const SizeType TNumNodes>
void ExplicitTotalLagrangianBbar<TDim, TNumNodes>::FinalizeSolutionStep( ProcessInfo& rCurrentProcessInfo )
{
    KinematicVariables this_kinematic_variables = KinematicVariables();
    ConstitutiveVariables this_constitutive_variables = ConstitutiveVariables();

    // Create constitutive law parameters:
    ConstitutiveLaw::Parameters Values(this->GetGeometry(),this->GetProperties(),rCurrentProcessInfo);

    // Set constitutive law flags:
    Flags& ConstitutiveLawOptions=Values.GetOptions();
    ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

    Values.SetStrainVector(this_constitutive_variables.StrainVector);
    Values.SetStressVector(this_constitutive_variables.StressVector);
    Values.SetConstitutiveMatrix(this_constitutive_variables.D);

    // Reading integration points
    const GeometryType& r_geometry = this->GetGeometry();

    // Shape functions
    const Matrix& rNValues = r_geometry.ShapeFunctionsValues(this->GetIntegrationMethod());

    // Reading integration points
    const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry.IntegrationPoints(this->GetIntegrationMethod());

    bool compute_finalization;
    for ( IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number ) {
         compute_finalization = mConstitutiveLawVector[point_number]->Has(REQUIRES_MATERIAL_FINALIZATION) ? mConstitutiveLawVector[point_number]->GetValue(REQUIRES_MATERIAL_FINALIZATION, compute_finalization) : true;

        if (compute_finalization) {
            // Compute element kinematics B, F, DN_DX ...
            CalculateKinematicVariables(this_kinematic_variables, point_number, this->GetIntegrationMethod(), rNValues);

            // Compute constitutive law variables
            SetConstitutiveVariables(this_kinematic_variables, this_constitutive_variables, Values, point_number, integration_points);

            // Call the constitutive law to update material variables
            mConstitutiveLawVector[point_number]->FinalizeMaterialResponse(Values, ConstitutiveLaw::StressMeasure_PK2);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template< const SizeType TDim, const SizeType TNumNodes>
void ExplicitTotalLagrangianBbar<TDim, TNumNodes>::InitializeMaterial()
{
    KRATOS_TRY

    if ( this->GetProperties()[CONSTITUTIVE_LAW] != nullptr ) {
        const GeometryType& r_geometry = this->GetGeometry();
        const Properties& r_properties = this->GetProperties();
        const auto& N_values = r_geometry.ShapeFunctionsValues(this->GetIntegrationMethod());
        for ( IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number ) {
            mConstitutiveLawVector[point_number] = this->GetProperties()[CONSTITUTIVE_LAW]->Clone();
            mConstitutiveLawVector[point_number]->InitializeMaterial( r_properties, r_geometry, row(N_values , point_number ));
        }
    } else
        KRATOS_ERROR << "A constitutive law needs to be specified for the element with ID " << this->Id() << std::endl;

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< const SizeType TDim, const SizeType TNumNodes>
void ExplicitTotalLagrangianBbar<TDim, TNumNodes>::ResetConstitutiveLaw()
{
    KRATOS_TRY

    if ( this->GetProperties()[CONSTITUTIVE_LAW] != nullptr ) {
        const GeometryType& r_geometry = this->GetGeometry();
        const Properties& r_properties = this->GetProperties();
        const auto& N_values = r_geometry.ShapeFunctionsValues(this->GetIntegrationMethod());
        for ( IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number )
            mConstitutiveLawVector[point_number]->ResetMaterial( r_properties,  r_geometry, row( N_values, point_number ) );
    }

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

template< const SizeType TDim, const SizeType TNumNodes>
Element::Pointer ExplicitTotalLagrangianBbar<TDim, TNumNodes>::Clone (
    IndexType NewId,
    NodesArrayType const& rThisNodes
    ) const
{
    KRATOS_TRY

    ExplicitTotalLagrangianBbar<TDim, TNumNodes>::Pointer p_new_elem = Kratos::make_shared<ExplicitTotalLagrangianBbar>(NewId, this->GetGeometry().Create(rThisNodes), this->pGetProperties());
    p_new_elem->SetData(this->GetData());
    p_new_elem->Set(Flags(*this));

    // The vector containing the constitutive laws
    p_new_elem->SetConstitutiveLawVector(mConstitutiveLawVector);

    return p_new_elem;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template< const SizeType TDim, const SizeType TNumNodes>
void ExplicitTotalLagrangianBbar<TDim, TNumNodes>::GetValuesVector(
    Vector& rValues,
    int Step
    )
{
    auto& r_geom = this->GetGeometry();

    if (rValues.size() != ElementSize)
        rValues.resize(ElementSize, false);
    for (IndexType i = 0; i < TNumNodes; ++i) {
        const array_1d<double, 3 >& displacement = r_geom[i].FastGetSolutionStepValue(DISPLACEMENT, Step);
        const SizeType index = i * TDim;
        for(unsigned int k = 0; k < TDim; ++k) {
            rValues[index + k] = displacement[k];
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template< const SizeType TDim, const SizeType TNumNodes>
void ExplicitTotalLagrangianBbar<TDim, TNumNodes>::GetFirstDerivativesVector(
    Vector& rValues,
    int Step
    )
{
    auto& r_geom = this->GetGeometry();

    if (rValues.size() != ElementSize)
        rValues.resize(ElementSize, false);
    for (IndexType i = 0; i < TNumNodes; ++i) {
        const array_1d<double, 3 >& velocity = r_geom[i].FastGetSolutionStepValue(VELOCITY, Step);
        const SizeType index = i * TDim;
        for(unsigned int k = 0; k < TDim; ++k)
            rValues[index + k] = velocity[k];
    }
}

/***********************************************************************************/
/***********************************************************************************/

template< const SizeType TDim, const SizeType TNumNodes>
void ExplicitTotalLagrangianBbar<TDim, TNumNodes>::GetSecondDerivativesVector(
    Vector& rValues,
    int Step
    )
{
    auto& r_geom = this->GetGeometry();

    if (rValues.size() != ElementSize)
        rValues.resize(ElementSize, false);
    for (IndexType i = 0; i < TNumNodes; ++i) {
        const array_1d<double, 3 >& acceleration = r_geom[i].FastGetSolutionStepValue(ACCELERATION, Step);
        const SizeType index = i * TDim;
        for(unsigned int k = 0; k < TDim; ++k)
            rValues[index + k] = acceleration[k];
    }
}

/***********************************************************************************/
/***********************************************************************************/

template< const SizeType TDim, const SizeType TNumNodes>
void ExplicitTotalLagrangianBbar<TDim, TNumNodes>::AddExplicitContribution(
    const VectorType& rRHSVector,
    const Variable<VectorType>& rRHSVariable,
    Variable<double>& rDestinationVariable,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    // Auxiliar geometries and properties
    auto& r_geom = this->GetGeometry();
    const auto& r_prop = this->GetProperties();

    // Computin the vector of
    array_1d<double, ElementSize> element_mass_vector;
    this->CalculateLumpedMassVector(element_mass_vector);

    // Compiting the nodal mass
    if (rDestinationVariable == NODAL_MASS ) {
        for (IndexType i = 0; i < TNumNodes; ++i) {
            const IndexType index = i * TDim;

            #pragma omp atomic
            r_geom[i].GetValue(NODAL_MASS) += element_mass_vector[index];
        }
    }
    // Compiting the nodal damping
    if (rDestinationVariable == NODAL_DISPLACEMENT_DAMPING  || rDestinationVariable == NODAL_MASS ) {
        if (StructuralMechanicsElementUtilities::HasRayleighDamping(r_prop, rCurrentProcessInfo)) {
            const double alpha = StructuralMechanicsElementUtilities::GetRayleighAlpha(r_prop, rCurrentProcessInfo);
            for (IndexType i = 0; i < TNumNodes; ++i) {
                const IndexType index = i * TDim;

                #pragma omp atomic
                r_geom[i].GetValue(NODAL_DISPLACEMENT_DAMPING) += alpha * element_mass_vector[index];
            }
        }
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template< const SizeType TDim, const SizeType TNumNodes>
void ExplicitTotalLagrangianBbar<TDim, TNumNodes>::AddExplicitContribution(
    const VectorType& rRHSVector,
    const Variable<VectorType>& rRHSVariable,
    Variable<array_1d<double, 3>>& rDestinationVariable,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    auto& r_geom = this->GetGeometry();

    // Computing the force residual
    if (rRHSVariable == RESIDUAL_VECTOR && rDestinationVariable == FORCE_RESIDUAL) {
        for (IndexType i = 0; i < TNumNodes; ++i) {
            const IndexType index = TDim * i;

            array_1d<double, 3>& r_force_residual = r_geom[i].FastGetSolutionStepValue(FORCE_RESIDUAL);

            for (IndexType j = 0; j < TDim; ++j) {
                #pragma omp atomic
                r_force_residual[j] += rRHSVector[index + j];
            }
        }
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template< const SizeType TDim, const SizeType TNumNodes>
void ExplicitTotalLagrangianBbar<TDim, TNumNodes>::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo
    )
{
    KinematicVariables this_kinematic_variables = KinematicVariables();
    ConstitutiveVariables this_constitutive_variables = ConstitutiveVariables();

    // Resizing as needed the LHS
    const SizeType ElementSize = TNumNodes * TDim;

    // Resizing as needed the RHS if required
    if ( rRightHandSideVector.size() != ElementSize )
        rRightHandSideVector.resize( ElementSize, false );

    noalias(rRightHandSideVector) = ZeroVector( ElementSize ); //resetting RHS

    // Reading integration points
    auto& r_geometry = this->GetGeometry();
    const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry.IntegrationPoints(this->GetIntegrationMethod());

    ConstitutiveLaw::Parameters Values(r_geometry, this->GetProperties(),rCurrentProcessInfo);

    // Set constitutive law flags:
    Flags& ConstitutiveLawOptions=Values.GetOptions();
    ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

    // If strain has to be computed inside of the constitutive law with PK2
    Values.SetStrainVector(this_constitutive_variables.StrainVector); //this is the input  parameter

    // Shape functions
    const Matrix& rNValues = r_geometry.ShapeFunctionsValues(this->GetIntegrationMethod());

    // Some declarations
    const bool has_thickness = this->GetProperties().Has( THICKNESS ) ?  true : false;
    array_1d<double, 3> body_force;
    double integration_weight;

    // Computing in all integrations points
    for ( IndexType point_number = 0; point_number < integration_points.size(); ++point_number ) {
        // Contribution to external forces
        noalias(body_force) = this->GetBodyForce(integration_points, point_number);

        // Compute element kinematics B, F, DN_DX ...
        this->CalculateKinematicVariables(this_kinematic_variables, point_number, this->GetIntegrationMethod(), rNValues);

        // Compute material reponse
        this->CalculateConstitutiveVariables(this_kinematic_variables, this_constitutive_variables, Values, point_number, integration_points, ConstitutiveLaw::StressMeasure_PK2);

        // Calculating weights for integration on the reference configuration
        integration_weight = integration_points[point_number].Weight() * this_kinematic_variables.detJ0;

        if ( TDim == 2 && has_thickness)
            integration_weight *= this->GetProperties()[THICKNESS];

        this->CalculateAndAddResidualVector(rRightHandSideVector, this_kinematic_variables, rCurrentProcessInfo, body_force, this_constitutive_variables.StressVector, integration_weight);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template< const SizeType TDim, const SizeType TNumNodes>
void ExplicitTotalLagrangianBbar<TDim, TNumNodes>::CalculateOnIntegrationPoints(
    const Variable<bool>& rVariable,
    std::vector<bool>& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    const GeometryType::IntegrationPointsArrayType &integration_points = this->GetGeometry().IntegrationPoints(this->GetIntegrationMethod());

    const SizeType number_of_integration_points = integration_points.size();
    if (rOutput.size() != number_of_integration_points)
        rOutput.resize(number_of_integration_points);

    if (mConstitutiveLawVector[0]->Has( rVariable)) {
        for ( IndexType point_number = 0; point_number <number_of_integration_points; ++point_number ) {
            bool value;
            mConstitutiveLawVector[point_number]->GetValue( rVariable, value);
            rOutput[point_number] = value;
        }
    } else {
        // Create constitutive law parameters:
        ConstitutiveLaw::Parameters Values(this->GetGeometry(),this->GetProperties(),rCurrentProcessInfo);

        for ( IndexType ii = 0; ii < mConstitutiveLawVector.size(); ++ii ) {
            bool solution;
            solution = mConstitutiveLawVector[ii]->CalculateValue( Values, rVariable, solution);
            rOutput[ii] = solution;
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template< const SizeType TDim, const SizeType TNumNodes>
void ExplicitTotalLagrangianBbar<TDim, TNumNodes>::CalculateOnIntegrationPoints(
    const Variable<int>& rVariable,
    std::vector<int>& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    const GeometryType::IntegrationPointsArrayType &integration_points = this->GetGeometry().IntegrationPoints(this->GetIntegrationMethod());

    const SizeType number_of_integration_points = integration_points.size();
    if (rOutput.size() != number_of_integration_points)
        rOutput.resize(number_of_integration_points);

    if (mConstitutiveLawVector[0]->Has( rVariable)) {
        GetValueOnConstitutiveLaw(rVariable, rOutput);
    } else {
        CalculateOnConstitutiveLaw(rVariable, rOutput, rCurrentProcessInfo);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template< const SizeType TDim, const SizeType TNumNodes>
void ExplicitTotalLagrangianBbar<TDim, TNumNodes>::CalculateOnIntegrationPoints(
    const Variable<double>& rVariable,
    std::vector<double>& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    auto& r_geometry = this->GetGeometry();
    const GeometryType::IntegrationPointsArrayType &integration_points = r_geometry.IntegrationPoints(this->GetIntegrationMethod());

    if ( rOutput.size() != integration_points.size() )
        rOutput.resize( integration_points.size() );

    if (mConstitutiveLawVector[0]->Has( rVariable)) {
        GetValueOnConstitutiveLaw(rVariable, rOutput);
    } else {
        if (rVariable == INTEGRATION_WEIGHT) {
            KinematicVariables this_kinematic_variables = KinematicVariables();

            const bool has_thickness = this->GetProperties().Has(THICKNESS);
            for (IndexType point_number = 0; point_number < integration_points.size(); ++point_number) {
                this_kinematic_variables.detJ0 = CalculateDerivativesOnReferenceConfiguration(this_kinematic_variables.J0,
                                                                                    this_kinematic_variables.InvJ0,
                                                                                    this_kinematic_variables.DN_DX,
                                                                                    point_number,
                                                                                    this->GetIntegrationMethod());

                double integration_weight = integration_points[point_number].Weight() * this_kinematic_variables.detJ0;

                if (TDim == 2 && has_thickness)
                    integration_weight *= this->GetProperties()[THICKNESS];

                rOutput[point_number] = integration_weight;
            }
        } else if ( rVariable == STRAIN_ENERGY ) {
            KinematicVariables this_kinematic_variables = KinematicVariables();
            ConstitutiveVariables this_constitutive_variables = ConstitutiveVariables();

            // Create constitutive law parameters:
            ConstitutiveLaw::Parameters Values(r_geometry,this->GetProperties(),rCurrentProcessInfo);

            // Set constitutive law flags:
            Flags& ConstitutiveLawOptions=Values.GetOptions();
            ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false);
            ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, false);
            ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

            // Reading integration points
            const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry.IntegrationPoints(this->GetIntegrationMethod());

            // Shape functions
            const Matrix& rNValues = r_geometry.ShapeFunctionsValues(this->GetIntegrationMethod());

            // If strain has to be computed inside of the constitutive law with PK2
            Values.SetStrainVector(this_constitutive_variables.StrainVector); //this is the input  parameter

            for (IndexType point_number = 0; point_number < integration_points.size(); ++point_number) {
                // Compute element kinematics B, F, DN_DX ...
                CalculateKinematicVariables(this_kinematic_variables, point_number, this->GetIntegrationMethod(), rNValues);

                // Compute constitutive law variables
                SetConstitutiveVariables(this_kinematic_variables, this_constitutive_variables, Values, point_number, integration_points);

                double StrainEnergy = 0.0;

                mConstitutiveLawVector[point_number]->CalculateValue(Values, STRAIN_ENERGY, StrainEnergy);

                rOutput[point_number] = StrainEnergy;
            }
        } else if (rVariable == VON_MISES_STRESS) {
            KinematicVariables this_kinematic_variables = KinematicVariables();
            ConstitutiveVariables this_constitutive_variables = ConstitutiveVariables();

            // Shape functions
            const Matrix& rNValues = r_geometry.ShapeFunctionsValues(this->GetIntegrationMethod());

            // Create constitutive law parameters:
            ConstitutiveLaw::Parameters Values(r_geometry,this->GetProperties(),rCurrentProcessInfo);

            // Set constitutive law flags:
            Flags& ConstitutiveLawOptions=Values.GetOptions();
            ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false);
            ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
            ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

            Values.SetStrainVector(this_constitutive_variables.StrainVector);

            for (IndexType point_number = 0; point_number < integration_points.size(); ++point_number) {
                // Compute element kinematics B, F, DN_DX ...
                CalculateKinematicVariables(this_kinematic_variables, point_number, this->GetIntegrationMethod(), rNValues);

                // Compute material reponse
                CalculateConstitutiveVariables(this_kinematic_variables, this_constitutive_variables, Values, point_number, integration_points, ConstitutiveLaw::StressMeasure_PK2);

                const Matrix stress_tensor = MathUtils<double>::StressVectorToTensor( this_constitutive_variables.StressVector );

                double sigma_equivalent = 0.0;

                if (TDim == 2) {
                    sigma_equivalent = std::pow((stress_tensor(0,0) - stress_tensor(1,1)), 2.0) +
                                                3*(stress_tensor(0,1) * stress_tensor(1,0));
                } else {
                    sigma_equivalent = 0.5*(std::pow((stress_tensor(0,0) - stress_tensor(1,1)), 2.0) +
                                            std::pow((stress_tensor(1,1) - stress_tensor(2,2)), 2.0) +
                                            std::pow((stress_tensor(2,2) - stress_tensor(0,0)), 2.0) +
                                                    6*(stress_tensor(0,1) * stress_tensor(1,0) +
                                                        stress_tensor(1,2) * stress_tensor(2,1) +
                                                        stress_tensor(2,0) * stress_tensor(0,2)));
                }

                if( sigma_equivalent < 0.0 )
                    rOutput[point_number] = 0.0;
                else
                    rOutput[point_number] = std::sqrt(sigma_equivalent);
            }
        } else {
            CalculateOnConstitutiveLaw(rVariable, rOutput, rCurrentProcessInfo);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template< const SizeType TDim, const SizeType TNumNodes>
void ExplicitTotalLagrangianBbar<TDim, TNumNodes>::CalculateOnIntegrationPoints(
    const Variable<array_1d<double, 3>>& rVariable,
    std::vector<array_1d<double, 3>>& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    auto& r_geometry = this->GetGeometry();
    const GeometryType::IntegrationPointsArrayType &integration_points = r_geometry.IntegrationPoints(this->GetIntegrationMethod());

    const SizeType number_of_integration_points = integration_points.size();
    if ( rOutput.size() != number_of_integration_points )
        rOutput.resize( number_of_integration_points );

    if (mConstitutiveLawVector[0]->Has( rVariable)) {
        GetValueOnConstitutiveLaw(rVariable, rOutput);
    } else {
        if (rVariable == INTEGRATION_COORDINATES) {
            KinematicVariables this_kinematic_variables = KinematicVariables();

            for (IndexType point_number = 0; point_number < number_of_integration_points; ++point_number) {
                Point global_point;
                r_geometry.GlobalCoordinates(global_point, integration_points[point_number]);
                rOutput[point_number] = global_point.Coordinates();
            }
        } else {
            CalculateOnConstitutiveLaw(rVariable, rOutput, rCurrentProcessInfo);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template< const SizeType TDim, const SizeType TNumNodes>
void ExplicitTotalLagrangianBbar<TDim, TNumNodes>::CalculateOnIntegrationPoints(
    const Variable<array_1d<double, 6>>& rVariable,
    std::vector<array_1d<double, 6>>& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    auto& r_geometry = this->GetGeometry();
    const GeometryType::IntegrationPointsArrayType &integration_points = r_geometry.IntegrationPoints(this->GetIntegrationMethod());

    const SizeType number_of_integration_points = integration_points.size();
    if (rOutput.size() != number_of_integration_points)
        rOutput.resize(number_of_integration_points);

    if (mConstitutiveLawVector[0]->Has( rVariable)) {
        GetValueOnConstitutiveLaw(rVariable, rOutput);
    }  else {
        CalculateOnConstitutiveLaw(rVariable, rOutput, rCurrentProcessInfo);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template< const SizeType TDim, const SizeType TNumNodes>
void ExplicitTotalLagrangianBbar<TDim, TNumNodes>::CalculateOnIntegrationPoints(
    const Variable<Vector>& rVariable,
    std::vector<Vector>& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    auto& r_geometry = this->GetGeometry();
    const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry.IntegrationPoints( this->GetIntegrationMethod() );

    const SizeType number_of_integration_points = integration_points.size();
    if ( rOutput.size() != number_of_integration_points )
        rOutput.resize( number_of_integration_points );

    if (mConstitutiveLawVector[0]->Has( rVariable)) {
        GetValueOnConstitutiveLaw(rVariable, rOutput);
    } else {
        if ( rVariable == INSITU_STRESS ) {
            const SizeType StrainSize = mConstitutiveLawVector[0]->GetStrainSize();
            Vector strain_vector( StrainSize );

            for ( IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number ) {
                if ( rOutput[point_number].size() != strain_vector.size() )
                    rOutput[point_number].resize( strain_vector.size(), false );

                rOutput[point_number] = mConstitutiveLawVector[point_number]->GetValue( INSITU_STRESS, rOutput[point_number] );
            }
        } else if ( rVariable == CAUCHY_STRESS_VECTOR || rVariable == PK2_STRESS_VECTOR ) {
            // Create and initialize element variables:
            KinematicVariables this_kinematic_variables = KinematicVariables();
            ConstitutiveVariables this_constitutive_variables = ConstitutiveVariables();

            // Create constitutive law parameters:
            ConstitutiveLaw::Parameters Values(r_geometry,this->GetProperties(),rCurrentProcessInfo);

            // Set constitutive law flags:
            Flags& ConstitutiveLawOptions=Values.GetOptions();
            ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false);
            ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
            ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

            Values.SetStrainVector(this_constitutive_variables.StrainVector);

            // Shape functions
            const Matrix& rNValues = r_geometry.ShapeFunctionsValues(this->GetIntegrationMethod());

            // Reading integration points
            for ( IndexType point_number = 0; point_number < number_of_integration_points; ++point_number ) {
                // Compute element kinematics B, F, DN_DX ...
                CalculateKinematicVariables(this_kinematic_variables, point_number, this->GetIntegrationMethod(), rNValues);

                //call the constitutive law to update material variables
                if( rVariable == CAUCHY_STRESS_VECTOR) {
                    // Compute material reponse
                    CalculateConstitutiveVariables(this_kinematic_variables, this_constitutive_variables, Values, point_number, integration_points, ConstitutiveLaw::StressMeasure_Cauchy);
                } else {
                    // Compute material reponse
                    CalculateConstitutiveVariables(this_kinematic_variables, this_constitutive_variables, Values, point_number, integration_points,ConstitutiveLaw::StressMeasure_PK2);
                }

                if ( rOutput[point_number].size() != StrainSize )
                    rOutput[point_number].resize( StrainSize, false );

                rOutput[point_number] = this_constitutive_variables.StressVector;
            }
        } else if( rVariable == GREEN_LAGRANGE_STRAIN_VECTOR  || rVariable == ALMANSI_STRAIN_VECTOR ) {
            // Create and initialize element variables:
            KinematicVariables this_kinematic_variables = KinematicVariables();
            ConstitutiveVariables this_constitutive_variables = ConstitutiveVariables();

            // Create constitutive law parameters:
            ConstitutiveLaw::Parameters Values(r_geometry,this->GetProperties(),rCurrentProcessInfo);

            // Set constitutive law flags:
            Flags &ConstitutiveLawOptions=Values.GetOptions();
            ConstitutiveLawOptions.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false);
            ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_STRESS, false);
            ConstitutiveLawOptions.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

            Values.SetStrainVector(this_constitutive_variables.StrainVector);

            const ConstitutiveLaw::StressMeasure this_stress_measure = rVariable == GREEN_LAGRANGE_STRAIN_VECTOR ? ConstitutiveLaw::StressMeasure_PK2 : ConstitutiveLaw::StressMeasure_Kirchhoff;

            // Shape functions
            const Matrix& rNValues = r_geometry.ShapeFunctionsValues(this->GetIntegrationMethod());

            //reading integration points
            for ( IndexType point_number = 0; point_number < number_of_integration_points; ++point_number ) {
                // Compute element kinematics B, F, DN_DX ...
                CalculateKinematicVariables(this_kinematic_variables, point_number, this->GetIntegrationMethod(),rNValues);

                // Compute material reponse
                CalculateConstitutiveVariables(this_kinematic_variables, this_constitutive_variables, Values, point_number, integration_points, this_stress_measure);

                if ( rOutput[point_number].size() != StrainSize)
                    rOutput[point_number].resize( StrainSize, false );

                rOutput[point_number] = this_constitutive_variables.StrainVector;
            }
        } else {
            CalculateOnConstitutiveLaw(rVariable, rOutput, rCurrentProcessInfo);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template< const SizeType TDim, const SizeType TNumNodes>
void ExplicitTotalLagrangianBbar<TDim, TNumNodes>::CalculateOnIntegrationPoints(
    const Variable<Matrix>& rVariable,
    std::vector<Matrix>& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    auto& r_geometry = this->GetGeometry();
    const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry.IntegrationPoints( this->GetIntegrationMethod() );

    if ( rOutput.size() != integration_points.size() )
        rOutput.resize( integration_points.size() );

    if (mConstitutiveLawVector[0]->Has( rVariable)) {
        GetValueOnConstitutiveLaw(rVariable, rOutput);
    } else {
        if ( rVariable == CAUCHY_STRESS_TENSOR || rVariable == PK2_STRESS_TENSOR ) {
            std::vector<Vector> stress_vector;

            if( rVariable == CAUCHY_STRESS_TENSOR )
                this->CalculateOnIntegrationPoints( CAUCHY_STRESS_VECTOR, stress_vector, rCurrentProcessInfo );
            else
                this->CalculateOnIntegrationPoints( PK2_STRESS_VECTOR, stress_vector, rCurrentProcessInfo );

            // Loop integration points
            for ( IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number ) {
                if ( rOutput[point_number].size2() != TDim )
                    rOutput[point_number].resize( TDim, TDim, false );

                rOutput[point_number] = MathUtils<double>::StressVectorToTensor(stress_vector[point_number]);
            }
        }
        else if ( rVariable == GREEN_LAGRANGE_STRAIN_TENSOR  || rVariable == ALMANSI_STRAIN_TENSOR) {
            std::vector<Vector> strain_vector;
            if( rVariable == GREEN_LAGRANGE_STRAIN_TENSOR )
                CalculateOnIntegrationPoints( GREEN_LAGRANGE_STRAIN_VECTOR, strain_vector, rCurrentProcessInfo );
            else
                CalculateOnIntegrationPoints( ALMANSI_STRAIN_VECTOR, strain_vector, rCurrentProcessInfo );

            // Loop integration points
            for ( IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number ) {
                if ( rOutput[point_number].size2() != TDim )
                    rOutput[point_number].resize( TDim, TDim, false );

                rOutput[point_number] = MathUtils<double>::StrainVectorToTensor(strain_vector[point_number]);
            }
        } else if ( rVariable == DEFORMATION_GRADIENT ) { // VARIABLE SET FOR TRANSFER PURPOUSES

            // Shape functions
            const Matrix& rNValues = r_geometry.ShapeFunctionsValues(this->GetIntegrationMethod());

            // Create and initialize element variables:
            KinematicVariables this_kinematic_variables = KinematicVariables();
            ConstitutiveVariables this_constitutive_variables = ConstitutiveVariables();

            // Create constitutive law parameters:
            ConstitutiveLaw::Parameters Values(r_geometry,this->GetProperties(),rCurrentProcessInfo);

            // Reading integration points
            for ( IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number ) {
                // Compute element kinematics B, F, DN_DX ...
                CalculateKinematicVariables(this_kinematic_variables, point_number, this->GetIntegrationMethod(), rNValues);

                if( rOutput[point_number].size2() != this_kinematic_variables.F.size2() )
                    rOutput[point_number].resize( this_kinematic_variables.F.size1() , this_kinematic_variables.F.size2() , false );

                rOutput[point_number] = this_kinematic_variables.F;
            }
        }  else {
            CalculateOnConstitutiveLaw(rVariable, rOutput, rCurrentProcessInfo);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template< const SizeType TDim, const SizeType TNumNodes>
void ExplicitTotalLagrangianBbar<TDim, TNumNodes>::SetValueOnIntegrationPoints(
    const Variable<bool>& rVariable,
    std::vector<bool>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    if (mConstitutiveLawVector[0]->Has( rVariable)) {
        for ( IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number ) {
            mConstitutiveLawVector[point_number]->SetValue( rVariable,rValues[point_number], rCurrentProcessInfo);
        }
    } else {
        KRATOS_WARNING("ExplicitTotalLagrangianBbar") << "The variable " << rVariable << " is not implemented in the current ConstitutiveLaw" << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template< const SizeType TDim, const SizeType TNumNodes>
void ExplicitTotalLagrangianBbar<TDim, TNumNodes>::SetValueOnIntegrationPoints(
    const Variable<int>& rVariable,
    std::vector<int>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    if (mConstitutiveLawVector[0]->Has( rVariable)) {
        for ( IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number ) {
            mConstitutiveLawVector[point_number]->SetValue( rVariable,rValues[point_number], rCurrentProcessInfo);
        }
    } else {
        KRATOS_WARNING("ExplicitTotalLagrangianBbar") << "The variable " << rVariable << " is not implemented in the current ConstitutiveLaw" << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template< const SizeType TDim, const SizeType TNumNodes>
void ExplicitTotalLagrangianBbar<TDim, TNumNodes>::SetValueOnIntegrationPoints(
    const Variable<double>& rVariable,
    std::vector<double>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    if (mConstitutiveLawVector[0]->Has( rVariable)) {
        for ( IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number ) {
            mConstitutiveLawVector[point_number]->SetValue( rVariable,rValues[point_number], rCurrentProcessInfo);
        }
    } else {
        KRATOS_WARNING("ExplicitTotalLagrangianBbar") << "The variable " << rVariable << " is not implemented in the current ConstitutiveLaw" << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template< const SizeType TDim, const SizeType TNumNodes>
void ExplicitTotalLagrangianBbar<TDim, TNumNodes>::SetValueOnIntegrationPoints(
    const Variable<Vector>& rVariable,
    std::vector<Vector>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    if (mConstitutiveLawVector[0]->Has( rVariable)) {
        for ( IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number ) {
            mConstitutiveLawVector[point_number]->SetValue( rVariable,rValues[point_number], rCurrentProcessInfo);
        }
    } else {
        KRATOS_WARNING("ExplicitTotalLagrangianBbar") << "The variable " << rVariable << " is not implemented in the current ConstitutiveLaw" << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template< const SizeType TDim, const SizeType TNumNodes>
void ExplicitTotalLagrangianBbar<TDim, TNumNodes>::SetValueOnIntegrationPoints(
    const Variable<ConstitutiveLaw::Pointer>& rVariable,
    std::vector<ConstitutiveLaw::Pointer>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    if (rVariable == CONSTITUTIVE_LAW) {
        const SizeType integration_points_number = mConstitutiveLawVector.size();
        for ( IndexType point_number = 0; point_number < integration_points_number; ++point_number ) {
            mConstitutiveLawVector[point_number] = rValues[point_number];
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template< const SizeType TDim, const SizeType TNumNodes>
void ExplicitTotalLagrangianBbar<TDim, TNumNodes>::SetValueOnIntegrationPoints(
    const Variable<array_1d<double, 3 > >& rVariable,
    std::vector<array_1d<double, 3 > > rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    if (mConstitutiveLawVector[0]->Has( rVariable)) {
        for ( IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number ) {
            mConstitutiveLawVector[point_number]->SetValue( rVariable,rValues[point_number], rCurrentProcessInfo);
        }
    } else {
        KRATOS_WARNING("ExplicitTotalLagrangianBbar") << "The variable " << rVariable << " is not implemented in the current ConstitutiveLaw" << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template< const SizeType TDim, const SizeType TNumNodes>
void ExplicitTotalLagrangianBbar<TDim, TNumNodes>::SetValueOnIntegrationPoints(
    const Variable<array_1d<double, 6 > >& rVariable,
    std::vector<array_1d<double, 6 > > rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    if (mConstitutiveLawVector[0]->Has( rVariable)) {
        for ( IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number ) {
            mConstitutiveLawVector[point_number]->SetValue( rVariable,rValues[point_number], rCurrentProcessInfo);
        }
    } else {
        KRATOS_WARNING("ExplicitTotalLagrangianBbar") << "The variable " << rVariable << " is not implemented in the current ConstitutiveLaw" << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template< const SizeType TDim, const SizeType TNumNodes>
void ExplicitTotalLagrangianBbar<TDim, TNumNodes>::SetValueOnIntegrationPoints(
    const Variable<Matrix>& rVariable,
    std::vector<Matrix>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    if (mConstitutiveLawVector[0]->Has( rVariable)) {
        for ( IndexType point_number = 0; point_number < mConstitutiveLawVector.size(); ++point_number ) {
            mConstitutiveLawVector[point_number]->SetValue( rVariable,rValues[point_number], rCurrentProcessInfo);
        }
    } else {
        KRATOS_WARNING("ExplicitTotalLagrangianBbar") << "The variable " << rVariable << " is not implemented in the current ConstitutiveLaw" << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template< const SizeType TDim, const SizeType TNumNodes>
void ExplicitTotalLagrangianBbar<TDim, TNumNodes>::GetValueOnIntegrationPoints(
    const Variable<bool>& rVariable,
    std::vector<bool>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    this->CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
}

/***********************************************************************************/
/***********************************************************************************/

template< const SizeType TDim, const SizeType TNumNodes>
void ExplicitTotalLagrangianBbar<TDim, TNumNodes>::GetValueOnIntegrationPoints(
    const Variable<int>& rVariable,
    std::vector<int>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    this->CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
}

/***********************************************************************************/
/***********************************************************************************/

template< const SizeType TDim, const SizeType TNumNodes>
void ExplicitTotalLagrangianBbar<TDim, TNumNodes>::GetValueOnIntegrationPoints(
    const Variable<double>& rVariable,
    std::vector<double>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    this->CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
}

/***********************************************************************************/
/***********************************************************************************/

template< const SizeType TDim, const SizeType TNumNodes>
void ExplicitTotalLagrangianBbar<TDim, TNumNodes>::GetValueOnIntegrationPoints(
    const Variable<array_1d<double, 3>>& rVariable,
    std::vector<array_1d<double, 3>>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    this->CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
}

/***********************************************************************************/
/***********************************************************************************/

template< const SizeType TDim, const SizeType TNumNodes>
void ExplicitTotalLagrangianBbar<TDim, TNumNodes>::GetValueOnIntegrationPoints(
    const Variable<array_1d<double, 6>>& rVariable,
    std::vector<array_1d<double, 6>>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    this->CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
}

/***********************************************************************************/
/***********************************************************************************/

template< const SizeType TDim, const SizeType TNumNodes>
void ExplicitTotalLagrangianBbar<TDim, TNumNodes>::GetValueOnIntegrationPoints(
    const Variable<Vector>& rVariable,
    std::vector<Vector>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    this->CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
}

/***********************************************************************************/
/***********************************************************************************/

template< const SizeType TDim, const SizeType TNumNodes>
void ExplicitTotalLagrangianBbar<TDim, TNumNodes>::GetValueOnIntegrationPoints(
    const Variable<Matrix>& rVariable,
    std::vector<Matrix>& rValues,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    this->CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
}

/***********************************************************************************/
/***********************************************************************************/

template< const SizeType TDim, const SizeType TNumNodes>
void ExplicitTotalLagrangianBbar<TDim, TNumNodes>::GetValueOnIntegrationPoints(
        const Variable<ConstitutiveLaw::Pointer>& rVariable,
        std::vector<ConstitutiveLaw::Pointer>& rValues,
        const ProcessInfo& rCurrentProcessInfo
        )
{
    if (rVariable == CONSTITUTIVE_LAW) {
        const SizeType integration_points_number = mConstitutiveLawVector.size();
        if (rValues.size() != integration_points_number) {
            rValues.resize(integration_points_number);
        }
        for (IndexType point_number = 0; point_number < integration_points_number; ++point_number) {
            rValues[point_number] = mConstitutiveLawVector[point_number];
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template< const SizeType TDim, const SizeType TNumNodes>
int  ExplicitTotalLagrangianBbar<TDim, TNumNodes>::Check( const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY;

    Element::Check(rCurrentProcessInfo);

    // Verify that the variables are correctly initialized
    KRATOS_CHECK_VARIABLE_KEY(DISPLACEMENT)
    KRATOS_CHECK_VARIABLE_KEY(VELOCITY)
    KRATOS_CHECK_VARIABLE_KEY(ACCELERATION)
    KRATOS_CHECK_VARIABLE_KEY(DENSITY)
    KRATOS_CHECK_VARIABLE_KEY(VOLUME_ACCELERATION)
    KRATOS_CHECK_VARIABLE_KEY(THICKNESS)

    // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
    for ( IndexType i = 0; i < TNumNodes; i++ ) {
        const NodeType &rnode = this->GetGeometry()[i];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISPLACEMENT,rnode)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY,rnode)
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ACCELERATION,rnode)

        KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_X, rnode)
        KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Y, rnode)
        KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Z, rnode)
    }

    // Verify that the constitutive law exists
    KRATOS_ERROR_IF_NOT(this->GetProperties().Has( CONSTITUTIVE_LAW )) << "Constitutive law not provided for property " << this->GetProperties().Id() << std::endl;

    // Verify that the constitutive law has the correct TDim
    if ( TDim == 2 ) {
        KRATOS_ERROR_IF( StrainSize < 3 || StrainSize > 4) << "Wrong constitutive law used. This is a 2D element! expected strain size is 3 or 4 (el id = ) " << this->Id() << std::endl;
    } else {
        KRATOS_ERROR_IF_NOT(StrainSize == 6) << "Wrong constitutive law used. This is a 3D element! expected strain size is 6 (el id = ) "<<  this->Id() << std::endl;
    }

    // Check constitutive law
    if ( mConstitutiveLawVector.size() > 0 ) {
        return mConstitutiveLawVector[0]->Check( this->GetProperties(), this->GetGeometry(), rCurrentProcessInfo );
    }

    return 0;

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template< const SizeType TDim, const SizeType TNumNodes>
void ExplicitTotalLagrangianBbar<TDim, TNumNodes>::CalculateKinematicVariables(
    KinematicVariables& rThisKinematicVariables,
    const IndexType PointNumber,
    const GeometryType::IntegrationMethod& rIntegrationMethod,
    const Matrix& rNValues
    )
{
    // Shape functions
    auto& r_geometry = this->GetGeometry();
    rThisKinematicVariables.N = row(rNValues, PointNumber);

    rThisKinematicVariables.detJ0 = this->CalculateDerivativesOnReferenceConfiguration(rThisKinematicVariables.J0, rThisKinematicVariables.InvJ0, rThisKinematicVariables.DN_DX, PointNumber, rIntegrationMethod);
    KRATOS_ERROR_IF(rThisKinematicVariables.detJ0 < 0.0) << "WARNING:: ELEMENT ID: " << this->Id() << " INVERTED. DETJ0: " << rThisKinematicVariables.detJ0 << std::endl;

    noalias(rThisKinematicVariables.F) = ZeroMatrix(TDim, TDim);
    for (IndexType i_node = 0; i_node < TNumNodes; ++i_node) {
        const array_1d<double, 3>& r_coordinates = r_geometry[i_node].Coordinates();
        for (IndexType i_dim = 0; i_dim < TDim; ++i_dim) {
            for (IndexType j_dim = 0; j_dim < TDim; ++j_dim) {
                rThisKinematicVariables.F(i_dim, j_dim) += rThisKinematicVariables.DN_DX(i_node, j_dim) * r_coordinates[i_dim];
            }
        }
    }

    CalculateB(rThisKinematicVariables.B, rThisKinematicVariables.F, rThisKinematicVariables.DN_DX);

    rThisKinematicVariables.detF = MathUtils<double>::Det(rThisKinematicVariables.F);
}

/***********************************************************************************/
/***********************************************************************************/

template< const SizeType TDim, const SizeType TNumNodes>
void ExplicitTotalLagrangianBbar<TDim, TNumNodes>::CalculateB(
    BoundedMatrix<double, StrainSize, TDim * TNumNodes>& rB,
    const Matrix& rF,
    const BoundedMatrix<double, TNumNodes, TDim>& rDN_DX
    )
{
    KRATOS_TRY
    if (TDim == 2) {
        for (IndexType i = 0; i < TNumNodes; ++i) {
            const auto index = TDim * i;
            rB(0, index + 0) = rF(0, 0) * rDN_DX(i, 0);
            rB(0, index + 1) = rF(1, 0) * rDN_DX(i, 0);
            rB(1, index + 0) = rF(0, 1) * rDN_DX(i, 1);
            rB(1, index + 1) = rF(1, 1) * rDN_DX(i, 1);
            rB(2, index + 0) = rF(0, 0) * rDN_DX(i, 1) + rF(0, 1) * rDN_DX(i, 0);
            rB(2, index + 1) = rF(1, 0) * rDN_DX(i, 1) + rF(1, 1) * rDN_DX(i, 0);
        }
    } else {
        for (IndexType i = 0; i < TNumNodes; ++i) {
            const IndexType index = TDim * i;
            rB(0, index + 0) = rF(0, 0) * rDN_DX(i, 0);
            rB(0, index + 1) = rF(1, 0) * rDN_DX(i, 0);
            rB(0, index + 2) = rF(2, 0) * rDN_DX(i, 0);
            rB(1, index + 0) = rF(0, 1) * rDN_DX(i, 1);
            rB(1, index + 1) = rF(1, 1) * rDN_DX(i, 1);
            rB(1, index + 2) = rF(2, 1) * rDN_DX(i, 1);
            rB(2, index + 0) = rF(0, 2) * rDN_DX(i, 2);
            rB(2, index + 1) = rF(1, 2) * rDN_DX(i, 2);
            rB(2, index + 2) = rF(2, 2) * rDN_DX(i, 2);
            rB(3, index + 0) = rF(0, 0) * rDN_DX(i, 1) + rF(0, 1) * rDN_DX(i, 0);
            rB(3, index + 1) = rF(1, 0) * rDN_DX(i, 1) + rF(1, 1) * rDN_DX(i, 0);
            rB(3, index + 2) = rF(2, 0) * rDN_DX(i, 1) + rF(2, 1) * rDN_DX(i, 0);
            rB(4, index + 0) = rF(0, 1) * rDN_DX(i, 2) + rF(0, 2) * rDN_DX(i, 1);
            rB(4, index + 1) = rF(1, 1) * rDN_DX(i, 2) + rF(1, 2) * rDN_DX(i, 1);
            rB(4, index + 2) = rF(2, 1) * rDN_DX(i, 2) + rF(2, 2) * rDN_DX(i, 1);
            rB(5, index + 0) = rF(0, 2) * rDN_DX(i, 0) + rF(0, 0) * rDN_DX(i, 2);
            rB(5, index + 1) = rF(1, 2) * rDN_DX(i, 0) + rF(1, 0) * rDN_DX(i, 2);
            rB(5, index + 2) = rF(2, 2) * rDN_DX(i, 0) + rF(2, 0) * rDN_DX(i, 2);
        }
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

template< const SizeType TDim, const SizeType TNumNodes>
void ExplicitTotalLagrangianBbar<TDim, TNumNodes>::CalculateConstitutiveVariables(
    KinematicVariables& rThisKinematicVariables,
    ConstitutiveVariables& rThisConstitutiveVariables,
    ConstitutiveLaw::Parameters& rValues,
    const IndexType PointNumber,
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints,
    const ConstitutiveLaw::StressMeasure ThisStressMeasure
    )
{
    // Setting the variables for the CL
    SetConstitutiveVariables(rThisKinematicVariables, rThisConstitutiveVariables, rValues, PointNumber, IntegrationPoints);

    // Actually do the computations in the ConstitutiveLaw
    mConstitutiveLawVector[PointNumber]->CalculateMaterialResponse(rValues, ThisStressMeasure); //here the calculations are actually done
}

/***********************************************************************************/
/***********************************************************************************/

template< const SizeType TDim, const SizeType TNumNodes>
void ExplicitTotalLagrangianBbar<TDim, TNumNodes>::SetConstitutiveVariables(
    KinematicVariables& rThisKinematicVariables,
    ConstitutiveVariables& rThisConstitutiveVariables,
    ConstitutiveLaw::Parameters& rValues,
    const IndexType PointNumber,
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints
    )
{
    // Here we essentially set the input parameters
    rValues.SetShapeFunctionsValues(rThisKinematicVariables.N); // shape functions
    rValues.SetDeterminantF(rThisKinematicVariables.detF); // Assuming the determinant is computed somewhere else
    rValues.SetDeformationGradientF(rThisKinematicVariables.F); //F computed somewhere else

    // Here we set the space on which the results shall be written
    rValues.SetConstitutiveMatrix(rThisConstitutiveVariables.D); // Assuming the determinant is computed somewhere else
    rValues.SetStressVector(rThisConstitutiveVariables.StressVector); //F computed somewhere else
}

/***********************************************************************************/
/***********************************************************************************/

template< const SizeType TDim, const SizeType TNumNodes>
double ExplicitTotalLagrangianBbar<TDim, TNumNodes>::CalculateDerivativesOnReferenceConfiguration(
    BoundedMatrix<double, TDim, TDim>& rJ0,
    BoundedMatrix<double, TDim, TDim>& rInvJ0,
    BoundedMatrix<double, TNumNodes, TDim>& rDN_DX,
    const IndexType PointNumber,
    IntegrationMethod ThisIntegrationMethod
    ) const
{
    const GeometryType& r_geom = this->GetGeometry();
    GeometryUtils::DirectJacobianOnInitialConfiguration(r_geom, rJ0, PointNumber, ThisIntegrationMethod);
    double detJ0;
    MathUtils<double>::InvertMatrix(rJ0, rInvJ0, detJ0);
    const Matrix& rDN_De = this->GetGeometry().ShapeFunctionsLocalGradients(ThisIntegrationMethod)[PointNumber];
    GeometryUtils::ShapeFunctionsGradients(rDN_De, rInvJ0, rDN_DX);
    return detJ0;
}

/***********************************************************************************/
/***********************************************************************************/

template< const SizeType TDim, const SizeType TNumNodes>
double ExplicitTotalLagrangianBbar<TDim, TNumNodes>::CalculateDerivativesOnCurrentConfiguration(
    BoundedMatrix<double, TDim, TDim>& rJ,
    BoundedMatrix<double, TDim, TDim>& rInvJ,
    BoundedMatrix<double, TNumNodes, TDim>& rDN_DX,
    const IndexType PointNumber,
    IntegrationMethod ThisIntegrationMethod
    ) const
{
    double detJ;
    const GeometryType& r_geom = this->GetGeometry();
    GeometryUtils::DirectJacobianOnCurrentConfiguration(r_geom, r_geom.IntegrationPoints(ThisIntegrationMethod)[PointNumber], rJ);
    const Matrix& DN_De = this->GetGeometry().ShapeFunctionsLocalGradients(ThisIntegrationMethod)[PointNumber];
    MathUtils<double>::InvertMatrix( rJ, rInvJ, detJ );
    GeometryUtils::ShapeFunctionsGradients(DN_De, rInvJ, rDN_DX);
    return detJ;
}

/***********************************************************************************/
/***********************************************************************************/

template< const SizeType TDim, const SizeType TNumNodes>
array_1d<double, 3> ExplicitTotalLagrangianBbar<TDim, TNumNodes>::GetBodyForce(
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints,
    const IndexType PointNumber
    ) const
{
    array_1d<double, 3> body_force;
    for (IndexType i = 0; i < 3; ++i)
        body_force[i] = 0.0;

    const auto& r_geometry = this->GetGeometry();
    const auto& r_properties = this->GetProperties();
    const double density = r_properties[DENSITY];

    if (r_properties.Has( VOLUME_ACCELERATION )) {
        noalias(body_force) += density * r_properties[VOLUME_ACCELERATION];
    } else if( r_geometry[0].SolutionStepsDataHas(VOLUME_ACCELERATION) ) {
        Vector N;
        N = r_geometry.ShapeFunctionsValues(N, IntegrationPoints[PointNumber].Coordinates());
        for (IndexType i_node = 0; i_node < TNumNodes; ++i_node)
            noalias(body_force) += N[i_node] * density * r_geometry[i_node].FastGetSolutionStepValue(VOLUME_ACCELERATION);
    }

    return body_force;
}

/***********************************************************************************/
/***********************************************************************************/

template< const SizeType TDim, const SizeType TNumNodes>
void ExplicitTotalLagrangianBbar<TDim, TNumNodes>::CalculateAndAddResidualVector(
    VectorType& rRightHandSideVector,
    const KinematicVariables& rThisKinematicVariables,
    const ProcessInfo& rCurrentProcessInfo,
    const array_1d<double, 3>& rBodyForce,
    const Vector& rStressVector,
    const double IntegrationWeight
    ) const
{
    KRATOS_TRY

    // Operation performed: rRightHandSideVector += ExtForce * IntegrationWeight
    for ( IndexType i = 0; i < TNumNodes; ++i ) {
        const SizeType index = TDim * i;

        for ( IndexType j = 0; j < TDim; ++j )
            rRightHandSideVector[index + j] += IntegrationWeight * rThisKinematicVariables.N[i] * rBodyForce[j];
    }

    // Operation performed: rRightHandSideVector -= IntForce * IntegrationWeight
    noalias( rRightHandSideVector ) -= IntegrationWeight * prod( trans( rThisKinematicVariables.B ), rStressVector );

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

template< const SizeType TDim, const SizeType TNumNodes>
void ExplicitTotalLagrangianBbar<TDim, TNumNodes>::CalculateLumpedMassVector(array_1d<double, ElementSize>& rMassVector) const
{
    KRATOS_TRY;

    const auto& r_prop = this->GetProperties();

    const double density = StructuralMechanicsElementUtilities::GetDensityForMassMatrixComputation(*this);
    const double thickness = (TDim == 2 && r_prop.Has(THICKNESS)) ? r_prop[THICKNESS] : 1.0;

    // LUMPED MASS MATRIX
    const double total_mass = this->GetGeometry().DomainSize() * density * thickness;

    Vector lumping_factors;
    lumping_factors = this->GetGeometry().LumpingFactors( lumping_factors );

    for ( IndexType i = 0; i < TNumNodes; ++i ) {
        const double temp = lumping_factors[i] * total_mass;
        for ( IndexType j = 0; j < TDim; ++j ) {
            IndexType index = i * TDim + j;
            rMassVector[index] = temp;
        }
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

template< const SizeType TDim, const SizeType TNumNodes>
void ExplicitTotalLagrangianBbar<TDim, TNumNodes>::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Element );
    rSerializer.save("ConstitutiveLawVector", mConstitutiveLawVector);
}

/***********************************************************************************/
/***********************************************************************************/

template< const SizeType TDim, const SizeType TNumNodes>
void ExplicitTotalLagrangianBbar<TDim, TNumNodes>::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Element );
    rSerializer.load("ConstitutiveLawVector", mConstitutiveLawVector);
}

/***********************************************************************************/
/***********************************************************************************/

template class ExplicitTotalLagrangianBbar<3, 8>;

} // Namespace Kratos
