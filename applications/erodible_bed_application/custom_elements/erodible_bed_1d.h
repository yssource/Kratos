//   
//   Project Name:        Kratos       
//   Last modified by:    $Author: it's me! $
//   Date:                $Date: 2008-08-08 23:58:38 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_ERODIBLE_BED_1D_ELEM_H_INCLUDED)
#define  KRATOS_ERODIBLE_BED_1D_ELEM_H_INCLUDED 

// System includes 


// External includes 
#include "boost/smart_ptr.hpp"


// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h" 

namespace Kratos
{

  class ErodibleBed1D
	  : public Element
   {
   public:
     
     /// Counted pointer of ErodibleBed1D
     KRATOS_CLASS_POINTER_DEFINITION(ErodibleBed1D);


    /// Default constructor.
     ErodibleBed1D(IndexType NewId, GeometryType::Pointer pGeometry);
     ErodibleBed1D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties);

     /// Destructor.
     virtual ~ ErodibleBed1D();


     Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const;

     void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);
     
     void CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo);
     
     void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo);

     void GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo);

     void InitializeSolutionStep(ProcessInfo& CurrentProcessInfo);


      
   protected:
   
   private:
	friend class Serializer;

       ErodibleBed1D() : Element()
       {
       }
       
       
   }; // Class ErodibleBed1D
}  // namespace Kratos.

#endif // KRATOS_FIC_1D_ELEM_H_INCLUDED  defined
