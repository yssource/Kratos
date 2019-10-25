// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Alejandro Cornejo & Agustina Giuliodori
//  Collaborator:
//

#if !defined(KRATOS_NITINOL_PSEUDO_ELASTICITY_3D_H_INCLUDED)
#define KRATOS_NITINOL_PSEUDO_ELASTICITY_3D_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/constitutive_law.h"
#include "custom_constitutive/elastic_isotropic_3d.h"
#include "custom_utilities/constitutive_law_utilities.h"

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
 * @class NitinolPseudoElasticity3D
 * @ingroup StructuralMechanicsApplication
 * @brief This is a viscous law using Maxwell formulation
 * @details The definition of a maxwell material can be found in https://en.wikipedia.org/wiki/Maxwell_material
 * This definition consists in a spring and a damper in serial
 *
 *           -----^^^^^^-----------[------
 *                Spring (K)   Damper (C)
 * The Maxwell law requires the definition of the following properties:
 * - VISCOUS_PARAMETER: It is the material coefficient of viscosity. This model describes the damper as a Newtonian fluid and models the spring with Hooke's law. 
 * @param TElasticBehaviourLaw Defines the elastic behaviour of the constitutive law (can be hyperelastic or just linear elastic, or any desired elastic behaviour)
 * @author Alejandro Cornejo & Lucia Barbu
 */
template<class TElasticBehaviourLaw>
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) NitinolPseudoElasticity3D
    : public TElasticBehaviourLaw
{
  public:
    ///@name Type Definitions
    ///@{

    /// Definition of the base CL class
    typedef ConstitutiveLaw CLBaseType;
    
    /// Definition of the base class
    typedef TElasticBehaviourLaw BaseType;

    /// The index definition
    typedef std::size_t IndexType;

    /// The size definition
    typedef std::size_t SizeType;

    /// Static definition of the dimension
    static constexpr SizeType Dimension = TElasticBehaviourLaw::Dimension;
    
    /// Static definition of the VoigtSize
    static constexpr SizeType VoigtSize = TElasticBehaviourLaw::VoigtSize;
    
    /// Counted pointer of GenericYieldSurface
    KRATOS_CLASS_POINTER_DEFINITION(NitinolPseudoElasticity3D);
    
    /// The node definition
    typedef Node<3> NodeType;
    
    /// The geometry definition
    typedef Geometry<NodeType> GeometryType;
    
    /// Definition of the machine precision tolerance
    static constexpr double tolerance = std::numeric_limits<double>::epsilon();

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor.
     */
    NitinolPseudoElasticity3D();

    /**
    * @brief Clone.
    */
    ConstitutiveLaw::Pointer Clone() const override;

    /**
    * @brief Copy constructor.
    */
    NitinolPseudoElasticity3D(const NitinolPseudoElasticity3D& rOther);
    
    /**
    * @brief Destructor.
    */
    ~NitinolPseudoElasticity3D() override;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * Computes the material response in terms of 1st Piola-Kirchhoff stresses and constitutive tensor
     * @see Parameters
     */
    void CalculateMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues) override;

    /**
     * Computes the material response in terms of 2nd Piola-Kirchhoff stresses and constitutive tensor
     * @see Parameters
     */
    void CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues) override;

    /**
     * Computes the material response in terms of Kirchhoff stresses and constitutive tensor
     * @see Parameters
     */
    void CalculateMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues) override;

    /**
     * Computes the material response in terms of Cauchy stresses and constitutive tensor
     * @see Parameters
     */
    void CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues) override;

    /**
     * Finalize the material response in terms of 1st Piola-Kirchhoff stresses
     * @see Parameters
     */
    void FinalizeMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues) override;

    /**
     * Finalize the material response in terms of 2nd Piola-Kirchhoff stresses
     * @see Parameters
     */
    void FinalizeMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues) override;

    /**
     * Finalize the material response in terms of Kirchhoff stresses
     * @see Parameters
     */
    void FinalizeMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues) override;

    /**
     * Finalize the material response in terms of Cauchy stresses
     * @see Parameters
     */
    void FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues) override;

    /**
     * @brief Calculates the value of a specified variable (Vector)
     * @param rParameterValues the needed parameters for the CL calculation
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @param rValue output: the value of the specified variable
     */
    Vector& CalculateValue(
        ConstitutiveLaw::Parameters& rParameterValues,
        const Variable<Vector>& rThisVariable,
        Vector& rValue
        ) override;

    /**
     * @brief Calculates the value of a specified variable (Matrix)
     * @param rParameterValues the needed parameters for the CL calculation
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @param rValue output: the value of the specified variable
     */
    Matrix& CalculateValue(
        ConstitutiveLaw::Parameters& rParameterValues,
        const Variable<Matrix>& rThisVariable,
        Matrix& rValue
        ) override;

    /**
     * @brief This function provides the place to perform checks on the completeness of the input.
     * @details It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input or that no common error is found.
     * @param rMaterialProperties The properties of the material
     * @param rElementGeometry The geometry of the element
     * @param rCurrentProcessInfo The current process info instance
     * @return 0 if OK, 1 otherwise
     */
    int Check(
        const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    void CalculatePseudoElasticMatrix(
      Matrix &rConstitutiveMatrix, 
      ConstitutiveLaw::Parameters& rParameterValues,
      const double MartensitePercentage);

    void CheckIfLoading(
        const Matrix &rPseudoElasticMatrix,
        const Vector &rStrainVector,
        bool &rIsLoading);

    double CalculatePseudoDruckerPragerUniaxialStress(
        const array_1d<double, VoigtSize> &rStressVector,
        ConstitutiveLaw::Parameters &rValues,
        array_1d<double, VoigtSize>& rDeviator);

    double CalculateThreshold(
        ConstitutiveLaw::Parameters &rValues,
        const bool IsLoading);

    void
    IntegrateStressVector(
        array_1d<double, VoigtSize> &rStressVector,
        const double YieldCondition,
        ConstitutiveLaw::Parameters &rValues,
        const bool SaveInternalVars,
        const bool IsLoading);

    void ForwardTransformation(
        const double YieldCondition,
        ConstitutiveLaw::Parameters &rValues,
        array_1d<double, VoigtSize> &rStressVector,
        const array_1d<double, VoigtSize> &rDeviator,
        const bool SaveInternalVars);

    void BackwardTransformation(
        const double YieldCondition,
        ConstitutiveLaw::Parameters &rValues,
        array_1d<double, VoigtSize> &rStressVector,
        const array_1d<double, VoigtSize> &rDeviator,
        const bool SaveInternalVars);

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    ///@}
    ///@name Friends
    ///@{

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
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

    ///@}
  private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    double mMartensitePercentage = 0.0;
    Vector mTransformationStrain = ZeroVector(VoigtSize);
    Vector mPreviousStrain       = ZeroVector(VoigtSize);

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    // Serialization

    friend class Serializer;

    // void save(Serializer &rSerializer) const override
    // {
    //     KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw)
    //     // rSerializer.save("PrevStressVector", mPrevStressVector);
    //     // rSerializer.save("PrevStrainVector", mPrevStrainVector);
    // }

    // void load(Serializer &rSerializer) override
    // {
    //     KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw)
    //     // rSerializer.load("PrevStressVector", mPrevStressVector);
    //     // rSerializer.load("PrevStrainVector", mPrevStrainVector);
    // }

    ///@}

}; // Class GenericYieldSurface

} // namespace Kratos
#endif
