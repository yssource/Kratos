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

#if !defined (KRATOS_CONSTITUTIVE_LAW_WITH_IMPOSED_DEFORMATION_H_INCLUDED)
#define  KRATOS_CONSTITUTIVE_LAW_WITH_IMPOSED_DEFORMATION_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/constitutive_law.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/**
 * @class ConstitutiveLawWithImposedDeformation
 * @ingroup KratosCore
 * @brief This is a constitutive law with imposed deformation
 * @author Vicente Mataix Ferrandiz
 * @see ConstitutiveLaw
 */
class KRATOS_API(KRATOS_CORE) ConstitutiveLawWithImposedDeformation
    : public ConstitutiveLaw
{
public:
    ///@name Type Definitions
    ///@{

    /**
     * Counted pointer of ConstitutiveLawWithImposedDeformation
     */
    KRATOS_CLASS_POINTER_DEFINITION( ConstitutiveLawWithImposedDeformation );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ConstitutiveLawWithImposedDeformation();

    /// Copy constructor.
    ConstitutiveLawWithImposedDeformation (const ConstitutiveLawWithImposedDeformation& rOther);

    /// Destructor.
    ~ConstitutiveLawWithImposedDeformation() override;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Clone function (has to be implemented by any derived class)
     * @return A pointer to a new instance of this constitutive law
     * @note implementation scheme:
     *      ImposedDeformation::Pointer p_clone(ImposedDeformation);
     *      return p_clone;
     */
    ConstitutiveLaw::Pointer Clone() const override;

    /**
     * @brief This is to be called at the very beginning of the calculation (e.g. from InitializeElement) in order to initialize all relevant attributes of the constitutive law
     * @param rMaterialProperties the Properties instance of the current element
     * @param rElementGeometry the geometry of the current element
     * @param rShapeFunctionsValues the shape functions values in the current integration point
     */
    void InitializeMaterial(
        const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry,
        const Vector& rShapeFunctionsValues
        ) override;

    /**
     * @brief Gets the current imposed deformation instance
     * @param rParameterValues The needed parameters for the imposed deformation (coming from constitutive law)
     * @return The current imposed deformation instance
     */
    ImposedDeformation* GetImposedDeformation (ConstitutiveLaw::Parameters& rParameterValues) override;

    /**
     * @brief This function provides the place to perform checks on the completeness of the input.
     * @details It is designed to be called only once (or anyway, not often) typically at the beginning of the calculations, so to verify that nothing is missing from the input or that no common error is found.
     * @param rMaterialProperties: The properties of the material
     * @param rElementGeometry: The geometry of the element
     * @param rCurrentProcessInfo: The current process info instance
     */
    int Check(
        const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    ///@}

protected:

    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    ///@}
private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    ImposedDeformation* mpImposedDeformation = NULL; /// Pointer to the imposed deformation

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{
    ///@}

    ///@}
    ///@name Private  Access
    ///@{
    ///@}

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;


}; // Class ConstitutiveLawWithImposedDeformation
}  // namespace Kratos.
#endif // KRATOS_CONSTITUTIVE_LAW_WITH_IMPOSED_DEFORMATION_H_INCLUDED  defined
