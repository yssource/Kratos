/*
//  KRATOS  _____________
//         /  _/ ____/   |
//         / // / __/ /| |
//       _/ // /_/ / ___ |
//      /___/\____/_/  |_| Application
//
//  Main authors:   Thomas Oberbichler
//                  Lukas Rauch
*/

#if !defined(KRATOS_IGA_BEAM_LOAD_CONDITION_H_INCLUDED)
#define KRATOS_IGA_BEAM_LOAD_CONDITION_H_INCLUDED

// System includes
#include "includes/define.h"
#include "includes/condition.h"

// External includes
#include <Eigen/Core>
#include "C:\_Masterarbeit\bibliotheken\HyperJet\HyperJet.h"

// Project includes
#include "iga_base_condition.h"


namespace Kratos
{
 class IgaBeamLoadCondition
  : public IgaBaseCondition<4>
    {
    public:
      ///@name Type Definitions
      ///@{
      
       /// Counted pointer of PointForce2D
       KRATOS_CLASS_POINTER_DEFINITION(IgaBeamLoadCondition);
      
      /// Default constructor. 
      IgaBeamLoadCondition(IndexType NewId, GeometryType::Pointer pGeometry);

      IgaBeamLoadCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

      /// Destructor.
      virtual ~IgaBeamLoadCondition();
      

      Condition::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const;
      
      void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo);

      void GetDofList(DofsVectorType& ConditionalDofList,ProcessInfo& CurrentProcessInfo);

      void CalculateAll(
          MatrixType& rLeftHandSideMatrix,
          VectorType& rRightHandSideVector,
          ProcessInfo& rCurrentProcessInfo,
          const bool ComputeLeftHandSide,
          const bool ComputeRightHandSide);

   protected:
 
 
   private:

	friend class Serializer;

	// A private default constructor necessary for serialization  
	IgaBeamLoadCondition() : IgaBaseCondition<4>()
       {
       }
        

  }; // Class PointForce2D 

} //namespace kratos 
#endif