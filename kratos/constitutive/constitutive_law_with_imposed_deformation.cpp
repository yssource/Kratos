//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "includes/imposed_deformation.h"
#include "constitutive/constitutive_law_with_imposed_deformation.h"

namespace Kratos
{

/******************************CONSTRUCTOR******************************************/
/***********************************************************************************/

ConstitutiveLawWithImposedDeformation::ConstitutiveLawWithImposedDeformation()
    : ConstitutiveLaw()
{
}

/******************************COPY CONSTRUCTOR*************************************/
/***********************************************************************************/

ConstitutiveLawWithImposedDeformation::ConstitutiveLawWithImposedDeformation(const ConstitutiveLawWithImposedDeformation& rOther)
    : ConstitutiveLaw(rOther),
      mpImposedDeformation(rOther.mpImposedDeformation)
{
}

/********************************CLONE**********************************************/
/***********************************************************************************/

ConstitutiveLaw::Pointer ConstitutiveLawWithImposedDeformation::Clone() const
{
    ConstitutiveLaw::Pointer p_clone(new ConstitutiveLawWithImposedDeformation(*this));
    return p_clone;
}

/*******************************DESTRUCTOR******************************************/
/***********************************************************************************/

ConstitutiveLawWithImposedDeformation::~ConstitutiveLawWithImposedDeformation()
{
    // TODO: Add if necessary
}

/***********************************************************************************/
/***********************************************************************************/

void ConstitutiveLawWithImposedDeformation::InitializeMaterial(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues
    )
{
    // Assign imposed deformation
    if (rMaterialProperties.Has(IMPOSED_DEFORMATION)) {
        mpImposedDeformation = (&(*(rMaterialProperties.GetValue(IMPOSED_DEFORMATION)->Clone())));
    }
}

/***********************************************************************************/
/***********************************************************************************/

ImposedDeformation* ConstitutiveLawWithImposedDeformation::GetImposedDeformation (ConstitutiveLaw::Parameters& rParameterValues)
{
    return mpImposedDeformation;
}

/***********************************************************************************/
/***********************************************************************************/

int ConstitutiveLawWithImposedDeformation::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_ERROR_IF(mpImposedDeformation == NULL) << "Imposed deformation not initialized" << std::endl;

    return 0;
}

/***********************************************************************************/
/***********************************************************************************/

void ConstitutiveLawWithImposedDeformation::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, ConstitutiveLaw)
    rSerializer.save("ImposedDeformation", mpImposedDeformation);
}

/***********************************************************************************/
/***********************************************************************************/

void ConstitutiveLawWithImposedDeformation::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, ConstitutiveLaw)
    rSerializer.load("ImposedDeformation", mpImposedDeformation);
}

} // Namespace Kratos
