// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Riccardo Rossi
//                   Ruben Zorrilla
//

// System includes
#if !defined(KRATOS_FSI_LINE_LOAD_CONDITION_2D_H_INCLUDED )
#define  KRATOS_FSI_LINE_LOAD_CONDITION_2D_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/variables.h"
#include "custom_conditions/line_load_condition_2d.h"

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

/// FSI line load intended to equilibrate the incompressible fluid added mass.
/** The aim of this load is to add the FSI artificial mass appeared in the fluid
  * field due to the incompressible added mass effect. Besides, it also adds the
  * standard line load contrubition (in case it applies).
  * To do that, the next terms are added:
  *     LHS += rho*h*M
  *     RHS -= rho*h*M*(acc_i - acc_old)
  * being:
  *     rho: the fluid density
  *     h : the hypothetical fluid thickness movilized by the structure
  *     M : the consistent mass matrix
  *     acc_i: structure previous non-linear iteration acceleration
  *     acc_old: structure previous FSI non-linear iteration acceleration (it is
  *              supposed to be stored in the non-historical database)
  */

class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION)  FSILineLoadCondition2D
    : public LineLoadCondition2D
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of FSILineLoadCondition2D
    KRATOS_CLASS_POINTER_DEFINITION( FSILineLoadCondition2D );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    FSILineLoadCondition2D( IndexType NewId, GeometryType::Pointer pGeometry );
    FSILineLoadCondition2D( IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties );

    /// Destructor.
    ~FSILineLoadCondition2D() override;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    Condition::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties
        ) const override;

    Condition::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties
        ) const override;

    void CalculateMassMatrix(
        MatrixType& rMassMatrix,
        ProcessInfo& rCurrentProcessInfo) override;

    void CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
		VectorType& rRightHandSideVector,
		ProcessInfo& rCurrentProcessInfo) override;

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
//      virtual String Info() const;

    /// Print information about this object.
//      virtual void PrintInfo(std::ostream& rOStream) const;

    /// Print object's data.
//      virtual void PrintData(std::ostream& rOStream) const;


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

    // A protected default constructor necessary for serialization
    FSILineLoadCondition2D() {};

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
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save( Serializer& rSerializer ) const override {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, LineLoadCondition2D );
    }

    void load( Serializer& rSerializer ) override {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, LineLoadCondition2D );
    }

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    //FSILineLoadCondition2D& operator=(const FSILineLoadCondition2D& rOther);

    /// Copy constructor.
    //FSILineLoadCondition2D(const FSILineLoadCondition2D& rOther);

    ///@}

}; // Class FSILineLoadCondition2D

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
/*  inline std::istream& operator >> (std::istream& rIStream,
        FSILineLoadCondition2D& rThis);
*/
/// output stream function
/*  inline std::ostream& operator << (std::ostream& rOStream,
        const FSILineLoadCondition2D& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }*/
///@}

}  // namespace Kratos.

#endif // KRATOS_FSI_LINE_LOAD_CONDITION_2D_H_INCLUDED  defined
