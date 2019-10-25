//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#ifndef KRATOS_IMPOSED_DEFORMATION_H_INCLUDED
#define KRATOS_IMPOSED_DEFORMATION_H_INCLUDED

// System includes


// External includes

// Project includes
#include "includes/define.h"
#include "includes/constitutive_law.h"
#include "includes/kratos_parameters.h"
#include "includes/serializer.h"

namespace Kratos 
{
///@addtogroup KratosCore
///@{

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
 * @class ImposedDeformation
 * @ingroup KratosCore
 * @brief This is the baase class defined in order to impose an initial deformation in a imposed deformation. Must be specialized in order to be able to impose different behaviours
 * @details Must be called in the imposed deformation implementation
 * @author Vicente Mataix Ferrandiz
 * @see ConstitutiveLaw
 */
class KRATOS_API(KRATOS_CORE) ImposedDeformation
{
public:
    
    ///@name Type Definitions
    ///@{
    
    /// Pointer definition of ImposedDeformation
    KRATOS_CLASS_POINTER_DEFINITION(ImposedDeformation);
    
    ///@}
    ///@name Life Cycle
    ///@{
    
    /// Default constructor.
    ImposedDeformation() 
    {

    }
    
    /// Destructor.
    virtual ~ImposedDeformation() = default;

    ///@}
    ///@name Operators
    ///@{
  
    ///@}
    ///@name Operations
    ///@{
  
  	/**
     * @brief Clone function (has to be implemented by any derived class)
     * @return A pointer to a new instance of this imposed deformation
     * @note implementation scheme:
     *      ImposedDeformation::Pointer p_clone(ImposedDeformation);
     *      return p_clone;
     */
    virtual Pointer Clone() const;

    /**
     * @brief Creates a new imposed deformation pointer
     * @param NewParameters The configuration parameters of the new imposed deformation
     * @return a Pointer to the new imposed deformation
     */
    virtual Pointer Create(Kratos::Parameters NewParameters) const;
  
    /**
     * @brief Returns whether this imposed deformation has specified variable (boolean)
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the imposed deformation
     */
    virtual bool Has(const Variable<bool>& rThisVariable);

    /**
     * @brief Returns whether this imposed deformation has specified variable (integer)
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the imposed deformation
     */
    virtual bool Has(const Variable<int>& rThisVariable);

    /**
     * @brief Returns whether this imposed deformation has specified variable (double)
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the imposed deformation
     */
    virtual bool Has(const Variable<double>& rThisVariable);

    /**
     * @brief Returns whether this imposed deformation has specified variable (Vector)
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the imposed deformation
     */
    virtual bool Has(const Variable<Vector>& rThisVariable);

    /**
     * @brief Returns whether this imposed deformation has specified variable (Matrix)
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the imposed deformation
     */
    virtual bool Has(const Variable<Matrix>& rThisVariable);

    /**
     * @brief Returns whether this imposed deformation has specified variable (array of 3 components)
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the imposed deformation
     * @note Fixed size array of 3 doubles (e.g. for 2D stresses, plastic strains, ...)
     */
    virtual bool Has(const Variable<array_1d<double, 3> >& rThisVariable);

    /**
     * @brief Returns whether this imposed deformation has specified variable (array of 6 components)
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the imposed deformation
     * @note Fixed size array of 6 doubles (e.g. for stresses, plastic strains, ...)
     */
    virtual bool Has(const Variable<array_1d<double, 6> >& rThisVariable);

    /**
     * @brief Returns the value of a specified variable (boolean)
	 * @param pConstitutiveLaw The pointer to the imposed deformation which calls this imposed deformation
     * @param rThisVariable The variable to be returned
     * @param rValue A reference to the returned value
     * @return rValue[output] The value of the specified variable
     */
    virtual bool& GetValue(
	    const ConstitutiveLaw* pConstitutiveLaw,
		const Variable<bool>& rThisVariable, 
		bool& rValue
		);

    /**
     * Returns the value of a specified variable (integer)
	 * @param pConstitutiveLaw The pointer to the imposed deformation which calls this imposed deformation
     * @param rThisVariable The variable to be returned
     * @param rValue A reference to the returned value
     * @return rValue[output] The value of the specified variable
     */
    virtual int& GetValue(
	    const ConstitutiveLaw* pConstitutiveLaw,
		const Variable<int>& rThisVariable, 
		int& rValue
		);

    /**
     * @brief Returns the value of a specified variable (double)
	 * @param pConstitutiveLaw The pointer to the imposed deformation which calls this imposed deformation
     * @param rThisVariable The variable to be returned
     * @param rValue A reference to the returned value
     * @return rValue[output] The value of the specified variable
     */
    virtual double& GetValue(
	    const ConstitutiveLaw* pConstitutiveLaw,
		const Variable<double>& rThisVariable, 
		double& rValue
		);

    /**
     * @brief Returns the value of a specified variable (Vector)
	 * @param pConstitutiveLaw The pointer to the imposed deformation which calls this imposed deformation
     * @param rThisVariable The variable to be returned
     * @param rValue A reference to the returned value
     * @return rValue[output] The value of the specified variable
     */
    virtual Vector& GetValue(
	    const ConstitutiveLaw* pConstitutiveLaw,
		const Variable<Vector>& rThisVariable, 
		Vector& rValue
		);

    /**
     * @brief Returns the value of a specified variable (Matrix)
	 * @param pConstitutiveLaw The pointer to the imposed deformation which calls this imposed deformation
     * @param rThisVariable The variable to be returned
     * @return rValue[output] The value of the specified variable
     */
    virtual Matrix& GetValue(
	    const ConstitutiveLaw* pConstitutiveLaw,
		const Variable<Matrix>& rThisVariable, 
		Matrix& rValue
		);

    /**
     * @brief Returns the value of a specified variable (array of 3 components)
	 * @param pConstitutiveLaw The pointer to the imposed deformation which calls this imposed deformation
     * @param rThisVariable The variable to be returned
     * @param rValue A reference to the returned value
     * @return rValue[output] The value of the specified variable
     */
    virtual array_1d<double, 3 > & GetValue(
	    const ConstitutiveLaw* pConstitutiveLaw,
		const Variable<array_1d<double, 3>>& rThisVariable,
        array_1d<double, 3>& rValue
		);

    /**
     * @brief Returns the value of a specified variable (array of 6 components)
	 * @param pConstitutiveLaw The pointer to the imposed deformation which calls this imposed deformation
     * @param rThisVariable The variable to be returned
     * @param rValue A reference to the returned value
     * @return the value of the specified variable
     */
    virtual array_1d<double, 6 >& GetValue(
	    const ConstitutiveLaw* pConstitutiveLaw,
		const Variable<array_1d<double, 6>>& rThisVariable,
        array_1d<double, 6>& rValue
		);

    /**
     * @brief Sets the value of a specified variable (boolean)
	 * @param pConstitutiveLaw The pointer to the imposed deformation which calls this imposed deformation
     * @param rVariable The variable to be returned
     * @param Value new value of the specified variable
     * @param rCurrentProcessInfo the process info
     */
    virtual void SetValue(
	    const ConstitutiveLaw* pConstitutiveLaw,
		const Variable<bool>& rVariable,
        const bool& Value,
        const ProcessInfo& rCurrentProcessInfo
		);

    /**
     * @brief Sets the value of a specified variable (integer)
	 * @param pConstitutiveLaw The pointer to the imposed deformation which calls this imposed deformation
     * @param rVariable The variable to be returned
     * @param Value new value of the specified variable
     * @param rCurrentProcessInfo the process info
     */
    virtual void SetValue(
	    const ConstitutiveLaw* pConstitutiveLaw,
		const Variable<int>& rVariable,
        const int& Value,
        const ProcessInfo& rCurrentProcessInfo
		);

    /**
     * @brief Sets the value of a specified variable (double)
	 * @param pConstitutiveLaw The pointer to the imposed deformation which calls this imposed deformation
     * @param rVariable The variable to be returned
     * @param rValue new value of the specified variable
     * @param rCurrentProcessInfo the process info
     */
    virtual void SetValue(
	    const ConstitutiveLaw* pConstitutiveLaw,
		const Variable<double>& rVariable,
        const double& rValue,
        const ProcessInfo& rCurrentProcessInfo
		);

    /**
     * @brief Sets the value of a specified variable (Vector)
	 * @param pConstitutiveLaw The pointer to the imposed deformation which calls this imposed deformation
     * @param rVariable The variable to be returned
     * @param rValue new value of the specified variable
     * @param rCurrentProcessInfo the process info
     */
    virtual void SetValue(
	    const ConstitutiveLaw* pConstitutiveLaw,
		const Variable<Vector>& rVariable,
        const Vector& rValue,
		const ProcessInfo& rCurrentProcessInfo
		);

    /**
     * @brief Sets the value of a specified variable (Matrix)
	 * @param pConstitutiveLaw The pointer to the imposed deformation which calls this imposed deformation
     * @param rVariable The variable to be returned
     * @param rValue new value of the specified variable
     * @param rCurrentProcessInfo the process info
     */
    virtual void SetValue(
	    const ConstitutiveLaw* pConstitutiveLaw,
		const Variable<Matrix>& rVariable,
        const Matrix& rValue,
		const ProcessInfo& rCurrentProcessInfo
		);

    /**
     * @brief Sets the value of a specified variable (array of 3 components)
	 * @param pConstitutiveLaw The pointer to the imposed deformation which calls this imposed deformation
     * @param rVariable The variable to be returned
     * @param rValue new value of the specified variable
     * @param rCurrentProcessInfo the process info
     */
    virtual void SetValue(
	    const ConstitutiveLaw* pConstitutiveLaw,
		const Variable<array_1d<double, 3>>& rVariable,
		const array_1d<double, 3>& rValue,
        const ProcessInfo& rCurrentProcessInfo
		);

    /**
     * @brief Sets the value of a specified variable (array of 6 components)
	 * @param pConstitutiveLaw The pointer to the imposed deformation which calls this imposed deformation
     * @param rVariable The variable to be returned
     * @param rValue new value of the specified variable
     * @param rCurrentProcessInfo the process info
     */
    virtual void SetValue(
	    const ConstitutiveLaw* pConstitutiveLaw,
		const Variable<array_1d<double, 6>>& rVariable,
	    const array_1d<double, 6>& rValue,
	    const ProcessInfo& rCurrentProcessInfo
		);

    /**
     * @brief Calculates the value of a specified variable (bool)
	 * @param pConstitutiveLaw The pointer to the imposed deformation which calls this imposed deformation
     * @param rParameterValues The needed parameters for the imposed deformation (coming from constitutive law)
     * @param rThisVariable The variable to be returned
     * @param rValue A reference to the returned value
     * @param rValue[output] The value of the specified variable
     */
    virtual bool& CalculateValue(
	    const ConstitutiveLaw* pConstitutiveLaw,
		Parameters& rParameterValues, 
		const Variable<bool>& rThisVariable, 
		bool& rValue
		);

    /**
     * @brief Calculates the value of a specified variable (int)
	 * @param pConstitutiveLaw The pointer to the imposed deformation which calls this imposed deformation
     * @param rParameterValues The needed parameters for the imposed deformation (coming from constitutive law)
     * @param rThisVariable The variable to be returned
     * @param rValue A reference to the returned value
     * @param rValue[output] The value of the specified variable
     */
    virtual int& CalculateValue(
	    const ConstitutiveLaw* pConstitutiveLaw,
		Parameters& rParameterValues, 
		const Variable<int>& rThisVariable, 
		int& rValue
		);

    /**
     * @brief Calculates the value of a specified variable (double)
	 * @param pConstitutiveLaw The pointer to the imposed deformation which calls this imposed deformation
     * @param rParameterValues The needed parameters for the imposed deformation (coming from constitutive law)
     * @param rThisVariable The variable to be returned
     * @param rValue A reference to the returned value
     * @param rValue[output] The value of the specified variable
     */
    virtual double& CalculateValue(
	    const ConstitutiveLaw* pConstitutiveLaw,
		ConstitutiveLaw::Parameters& rParameterValues, 
		const Variable<double>& rThisVariable, 
		double& rValue
		);

    /**
     * @brief Calculates the value of a specified variable (Vector)
	 * @param pConstitutiveLaw The pointer to the imposed deformation which calls this imposed deformation
     * @param rParameterValues The needed parameters for the imposed deformation (coming from constitutive law)
     * @param rThisVariable The variable to be returned
     * @param rValue A reference to the returned value
     * @param rValue[output] The value of the specified variable
     */
    virtual Vector& CalculateValue(
	    const ConstitutiveLaw* pConstitutiveLaw,
		ConstitutiveLaw::Parameters& rParameterValues, 
		const Variable<Vector>& rThisVariable, 
		Vector& rValue
		);

    /**
     * @brief Calculates the value of a specified variable (Matrix)
	 * @param pConstitutiveLaw The pointer to the imposed deformation which calls this imposed deformation
     * @param rParameterValues The needed parameters for the imposed deformation (coming from constitutive law)
     * @param rThisVariable The variable to be returned
     * @param rValue A reference to the returned value
     * @param rValue[output] The value of the specified variable
     */
    virtual Matrix& CalculateValue(
	    const ConstitutiveLaw* pConstitutiveLaw,
		ConstitutiveLaw::Parameters& rParameterValues, 
		const Variable<Matrix>& rThisVariable, 
		Matrix& rValue
		);

    /**
     * @brief Calculates the value of a specified variable (array of 3 components)
	 * @param pConstitutiveLaw The pointer to the imposed deformation which calls this imposed deformation
     * @param rParameterValues The needed parameters for the imposed deformation (coming from constitutive law)
     * @param rThisVariable The variable to be returned
     * @param rValue A reference to the returned value
     * @param rValue[output] The value of the specified variable
     */
    virtual array_1d<double, 3>& CalculateValue(
	    const ConstitutiveLaw* pConstitutiveLaw,
		Parameters& rParameterValues, 
		const Variable<array_1d<double, 3 > >& rVariable,
	    array_1d<double, 3 > & rValue
		);

    /**
     * returns the value of a specified variable (array of 6 components)
	 * @param pConstitutiveLaw The pointer to the imposed deformation which calls this imposed deformation
     * @param rThisVariable The variable to be returned
     * @param rValue A reference to the returned value
     * @return the value of the specified variable
     */
    virtual array_1d<double, 6>& CalculateValue(
	    const ConstitutiveLaw* pConstitutiveLaw,
		ConstitutiveLaw::Parameters& rParameterValues, 
		const Variable<array_1d<double, 6 > >& rVariable,
	    array_1d<double, 6 > & rValue
		);
  
    /**
     * @brief This is to be called at the very beginning of the calculation by the ConstitutiveLaw
	 * @details Should be called in InitializeMaterial
	 * @param pConstitutiveLaw The pointer to the imposed deformation which calls this imposed deformation
     * @param rMaterialProperties the Properties instance of the current element
     * @param rElementGeometry the geometry of the current element
     * @param rShapeFunctionsValues the shape functions values in the current integration point
     */
    virtual void Initialize(
		const ConstitutiveLaw* pConstitutiveLaw,
		const Properties& rMaterialProperties,
		const GeometryType& rElementGeometry,
		const Vector& rShapeFunctionsValues
		);
  
    /**
     * @brief Computes the material response in terms of stresses and constitutive tensor
	 * @param pConstitutiveLaw The pointer to the imposed deformation which calls this imposed deformation
	 * @param rParameterValues The needed parameters for the imposed deformation (coming from constitutive law)
	 * @param rStressMeasure The stress measure considered
     * @see ConstitutiveLaw::Parameters
     * @see ConstitutiveLaw::StressMeasure
     */
    void CalculateMaterialResponse (
		const ConstitutiveLaw* pConstitutiveLaw,
		ConstitutiveLaw::Parameters& rParameterValues,
		const ConstitutiveLaw::StressMeasure& rStressMeasure
		);

    /**
     * @brief Computes the material response in terms of 1st Piola-Kirchhoff stresses and constitutive tensor
	 * @param pConstitutiveLaw The pointer to the imposed deformation which calls this imposed deformation
	 * @param rParameterValues The needed parameters for the imposed deformation (coming from constitutive law)
     * @see ConstitutiveLaw::Parameters
     */
    virtual void CalculateMaterialResponsePK1 (
		const ConstitutiveLaw* pConstitutiveLaw,
		ConstitutiveLaw::Parameters& rParameterValues	
		);

    /**
     * @brief Computes the material response in terms of 2nd Piola-Kirchhoff stresses and constitutive tensor
	 * @param pConstitutiveLaw The pointer to the imposed deformation which calls this imposed deformation
	 * @param rParameterValues The needed parameters for the imposed deformation (coming from constitutive law)
     * @see ConstitutiveLaw::Parameters
     */
    virtual void CalculateMaterialResponsePK2 (
		const ConstitutiveLaw* pConstitutiveLaw,
		ConstitutiveLaw::Parameters& rParameterValues
		);

    /**
     * @brief Computes the material response in terms of Kirchhoff stresses and constitutive tensor
	 * @param pConstitutiveLaw The pointer to the imposed deformation which calls this imposed deformation
	 * @param rParameterValues The needed parameters for the imposed deformation (coming from constitutive law)
     * @see ConstititutiveLaw::Parameters
     */
    virtual void CalculateMaterialResponseKirchhoff (
		const ConstitutiveLaw* pConstitutiveLaw,
		ConstitutiveLaw::Parameters& rParameterValues
		);

    /**
     * @brief Computes the material response in terms of Cauchy stresses and constitutive tensor
	 * @param pConstitutiveLaw The pointer to the imposed deformation which calls this imposed deformation
	 * @param rParameterValues The needed parameters for the imposed deformation (coming from constitutive law)
     * @see ConstitutiveLaw::Parameters
     */
    virtual void CalculateMaterialResponseCauchy(
		const ConstitutiveLaw* pConstitutiveLaw,
		ConstitutiveLaw::Parameters& rParameterValues
		);

    /**
     * @brief Initialize the material response,  called by the element in FinalizeSolutionStep.
	 * @param pConstitutiveLaw The pointer to the imposed deformation which calls this imposed deformation
	 * @param rParameterValues The needed parameters for the imposed deformation (coming from constitutive law)
	 * @param rStressMeasure The stress measure considered
     * @see ConstitutiveLaw::Parameters
     * @see ConstititutiveLaw::StressMeasures
     */
    void InitializeMaterialResponse (
		const ConstitutiveLaw* pConstitutiveLaw,
		ConstitutiveLaw::Parameters& rParameterValues,
		const ConstitutiveLaw::StressMeasure& rStressMeasure
		);

    /**
     * @brief Initialize the material response in terms of 1st Piola-Kirchhoff stresses
	 * @param pConstitutiveLaw The pointer to the imposed deformation which calls this imposed deformation
	 * @param rParameterValues The needed parameters for the imposed deformation (coming from constitutive law)
     * @see ConstitutiveLaw::Parameters
     */
    virtual void InitializeMaterialResponsePK1 (
		const ConstitutiveLaw* pConstitutiveLaw,
		ConstitutiveLaw::Parameters& rParameterValues
		);

    /**
     * @brief Initialize the material response in terms of 2nd Piola-Kirchhoff stresses
	 * @param pConstitutiveLaw The pointer to the imposed deformation which calls this imposed deformation
	 * @param rParameterValues The needed parameters for the imposed deformation (coming from constitutive law)
     * @see ConstitutiveLaw::Parameters
     */
    virtual void InitializeMaterialResponsePK2 (
		const ConstitutiveLaw* pConstitutiveLaw,
		ConstitutiveLaw::Parameters& rParameterValues
		);

    /**
     * @brief Initialize the material response in terms of Kirchhoff stresses
	 * @param pConstitutiveLaw The pointer to the imposed deformation which calls this imposed deformation
	 * @param rParameterValues The needed parameters for the imposed deformation (coming from constitutive law)
     * @see ConstitutiveLaw::Parameters
     */
    virtual void InitializeMaterialResponseKirchhoff (
		const ConstitutiveLaw* pConstitutiveLaw,
		ConstitutiveLaw::Parameters& rParameterValues
		);

    /**
     * @brief Initialize the material response in terms of Cauchy stresses
	 * @param pConstitutiveLaw The pointer to the imposed deformation which calls this imposed deformation
	 * @param rParameterValues The needed parameters for the imposed deformation (coming from constitutive law)
     * @see ConstitutiveLaw::Parameters
     */
    virtual void InitializeMaterialResponseCauchy (
		const ConstitutiveLaw* pConstitutiveLaw,
		ConstitutiveLaw::Parameters& rParameterValues
		);

    /**
     * @brief Finalize the material response,  called by the element in FinalizeSolutionStep.
	 * @param pConstitutiveLaw The pointer to the imposed deformation which calls this imposed deformation
	 * @param rParameterValues The needed parameters for the imposed deformation (coming from constitutive law)
	 * @param rStressMeasure The stress measure considered
     * @see ConstitutiveLaw::Parameters
     * @see ConstitutiveLaw::StressMeasures
     */
    void FinalizeMaterialResponse (
		const ConstitutiveLaw* pConstitutiveLaw,
		ConstitutiveLaw::Parameters& rParameterValues,
		const ConstitutiveLaw::StressMeasure& rStressMeasure
		);

    /**
     * @brief Finalize the material response in terms of 1st Piola-Kirchhoff stresses
	 * @param pConstitutiveLaw The pointer to the imposed deformation which calls this imposed deformation
	 * @param rParameterValues The needed parameters for the imposed deformation (coming from constitutive law)
     * @see ConstitutiveLaw::Parameters
     */
    virtual void FinalizeMaterialResponsePK1 (
		const ConstitutiveLaw* pConstitutiveLaw,
		ConstitutiveLaw::Parameters& rParameterValues
		);

    /**
     * @brief Finalize the material response in terms of 2nd Piola-Kirchhoff stresses
	 * @param pConstitutiveLaw The pointer to the imposed deformation which calls this imposed deformation
	 * @param rParameterValues The needed parameters for the imposed deformation (coming from constitutive law)
     * @see ConstitutiveLaw::Parameters
     */
    virtual void FinalizeMaterialResponsePK2(
		const ConstitutiveLaw* pConstitutiveLaw,
		ConstitutiveLaw::Parameters& rParameterValues
		);

    /**
     * @brief Finalize the material response in terms of Kirchhoff stresses
	 * @param pConstitutiveLaw The pointer to the imposed deformation which calls this imposed deformation
	 * @param rParameterValues The needed parameters for the imposed deformation (coming from constitutive law)
     * @see ConstitutiveLaw::Parameters
     */
    virtual void FinalizeMaterialResponseKirchhoff (
		const ConstitutiveLaw* pConstitutiveLaw,
		ConstitutiveLaw::Parameters& rParameterValues
		);

    /**
     * @brief Finalize the material response in terms of Cauchy stresses
	 * @param pConstitutiveLaw The pointer to the imposed deformation which calls this imposed deformation
	 * @param rParameterValues The needed parameters for the imposed deformation (coming from constitutive law)
     * @see ConstitutiveLaw::Parameters
     */
    virtual void FinalizeMaterialResponseCauchy (
		const ConstitutiveLaw* pConstitutiveLaw,
		ConstitutiveLaw::Parameters& rParameterValues
		);

    /**
     * @brief This can be used in order to reset all internal variables of the imposed deformation
	 * @param pConstitutiveLaw The pointer to the imposed deformation which calls this imposed deformation
     * @param rMaterialProperties the Properties instance of the current element
     * @param rElementGeometry the geometry of the current element
     * @param rShapeFunctionsValues the shape functions values in the current integration point
     */
    virtual void Reset(
		const ConstitutiveLaw* pConstitutiveLaw,
		const Properties& rMaterialProperties,
	    const GeometryType& rElementGeometry,
	    const Vector& rShapeFunctionsValues
		);
  
    /**
     * @brief This function is designed to be called once to perform all the checks needed on the input provided. 
	 * @details Checks can be "expensive" as the function is designed to catch user's errors.
	 * @param pConstitutiveLaw The pointer to the imposed deformation which calls this imposed deformation
     * @param rMaterialProperties The Properties instance of the current element
     * @param rElementGeometry The geometry of the current element
     * @param rShapeFunctionsValues The shape functions values in the current integration point
     * @return 0 If all OK, 1 otherwise
     */
    virtual int Check(
		const ConstitutiveLaw* pConstitutiveLaw,
		const Properties& rMaterialProperties,
	    const GeometryType& rElementGeometry,
	    const ProcessInfo& rCurrentProcessInfo
		);
  
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
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const
    {
    }

    void load(Serializer& rSerializer)
    {
    }

    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Unaccessible methods
    ///@{
}; // Class ImposedDeformation

} // namespace Kratos.
#endif
