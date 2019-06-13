// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Riccardo Rossi
//                   Antonia Larese
//
// System includes
#include <iostream>

// External includes

// Project includes
#include "custom_constitutive/linear_plane_stress_optimization.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{

/******************************CONSTRUCTOR******************************************/
/***********************************************************************************/

LinearPlaneStressOptimization::LinearPlaneStressOptimization()
    :LinearPlaneStress()
{
}

/******************************COPY CONSTRUCTOR*************************************/
/***********************************************************************************/

LinearPlaneStressOptimization::LinearPlaneStressOptimization(const LinearPlaneStressOptimization& rOther)
    :LinearPlaneStress(rOther)
{
}

/********************************CLONE**********************************************/
/***********************************************************************************/

ConstitutiveLaw::Pointer LinearPlaneStressOptimization::Clone() const
{
    LinearPlaneStressOptimization::Pointer p_clone(new LinearPlaneStressOptimization(*this));
    return p_clone;
}

/*******************************DESTRUCTOR******************************************/
/***********************************************************************************/

LinearPlaneStressOptimization::~LinearPlaneStressOptimization()
{
}

/***********************************************************************************/
/***********************************************************************************/

bool LinearPlaneStressOptimization::Has(const Variable<double>& rThisVariable)
{
    if (rThisVariable == SCALE_FACTOR) {
        return true;
    } else {
        return false;
    }
}

/***********************************************************************************/
/***********************************************************************************/

double& LinearPlaneStressOptimization::GetValue(
    const Variable<double>& rThisVariable,
    double& rValue
    )
{
    if(rThisVariable == SCALE_FACTOR) {
        rValue = mOptimizationCoefficient;
    }

    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

void LinearPlaneStressOptimization::SetValue(
    const Variable<double>& rVariable,
    const double& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    if(rVariable == SCALE_FACTOR) {
        mOptimizationCoefficient = rValue;
    }
}

//*************************CONSTITUTIVE LAW GENERAL FEATURES *************************

void LinearPlaneStressOptimization::GetLawFeatures(Features& rFeatures)
{
    //Set the type of law
    rFeatures.mOptions.Set( PLANE_STRESS_LAW );
    rFeatures.mOptions.Set( INFINITESIMAL_STRAINS );
    rFeatures.mOptions.Set( ISOTROPIC );

    //Set strain measure required by the consitutive law
    rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);
    rFeatures.mStrainMeasures.push_back(StrainMeasure_Deformation_Gradient);

    //Set the strain size
    rFeatures.mStrainSize = 3;

    //Set the spacedimension
    rFeatures.mSpaceDimension = 2;
}

/***********************************************************************************/
/***********************************************************************************/

void LinearPlaneStressOptimization::CalculateElasticMatrix(Matrix& C, ConstitutiveLaw::Parameters& rValues)
{
    BaseType::CalculateElasticMatrix(C, rValues);

    C *= mOptimizationCoefficient;
}

/***********************************************************************************/
/***********************************************************************************/

void LinearPlaneStressOptimization::CalculatePK2Stress(
    const Vector& rStrainVector,
    Vector& rStressVector,
    ConstitutiveLaw::Parameters& rValues
    )
{
    LinearPlaneStress::CalculatePK2Stress(rStrainVector, rStressVector,rValues);

    rStressVector *= mOptimizationCoefficient;
}

/***********************************************************************************/
/***********************************************************************************/

void LinearPlaneStressOptimization::CalculateCauchyGreenStrain(Parameters& rValues, Vector& rStrainVector)
{
    //1.-Compute total deformation gradient
    const Matrix& F = rValues.GetDeformationGradientF();

    // for shells/membranes in case the DeformationGradient is of size 3x3
    BoundedMatrix<double, 2, 2> F2x2;
    for (unsigned int i = 0; i<2; ++i)
        for (unsigned int j = 0; j<2; ++j)
            F2x2(i, j) = F(i, j);

    Matrix E_tensor = prod(trans(F2x2), F2x2);

    for (unsigned int i = 0; i<2; ++i)
        E_tensor(i, i) -= 1.0;

    E_tensor *= 0.5;

    noalias(rStrainVector) = MathUtils<double>::StrainTensorToVector(E_tensor);
}

} // Namespace Kratos
