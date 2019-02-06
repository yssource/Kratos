//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main author:     Jordi Cotela
//

#if !defined(KRATOS_TURBULENCE_K_EPSILON_PROCESS_H_INCLUDED )
#define  KRATOS_TURBULENCE_K_EPSILON_PROCESS_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"
#include "solving_strategies/convergencecriterias/displacement_criteria.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "processes/process.h"
#include "processes/variational_distance_calculation_process.h"

namespace Kratos
{
  ///@addtogroup FluidDynamicsApplication
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

  /// Auxiliary process to set Boussinesq buoyancy forces in variable temperature flows.
  /** This process modifies the BODY_FORCE variable according to the Boussinesq hypothesis
      so that the fluid element can take natural convection into account.

      This process makes use of the following data:
      - TEMPERATURE from the nodal solution step data: current temperature for the node (mandatory).
      - AMBIENT_TEMPERATURE from ProcessInfo: The reference temperature for the simulation (mandatory).
      - gravity from the Parameters passed in the constructor: an array that defines the gravity vector (mandatory).
      - thermal_expansion_coefficient from the Parameters: a double defining the thermal expansion coefficient for the fluid (optional).

      With this, the process calculates the Boussinesq force and assings it to the BODY_FORCE solution step variable of each node.
      The force is set to (1 + thermal_expansion_coefficient*(temperature - ambient_temperature) ) * g

      If the thermal expansion coefficient is not provided, it is assumed to be (1/ambient_temperature).
      This is the usual value for perfect gases (if the temperature is given in Kelvin).
   */
  template< class TSparseSpace, class TDenseSpace, class TLinearSolver >
  class TurbulenceEvmKEpsilonProcess: public Process
  {
  public:
      ///@name Type Definitions
      ///@{

      typedef TurbulenceEvmKEpsilonProcess<TSparseSpace,TDenseSpace,TLinearSolver> ModelType;

      typedef SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver> StrategyType;

      typedef Scheme< TSparseSpace, TDenseSpace> SchemeType;
      typedef typename SchemeType::Pointer SchemePointerType;
      typedef BuilderAndSolver< TSparseSpace, TDenseSpace, TLinearSolver > BuilderAndSolverType;
      typedef typename BuilderAndSolverType::Pointer BuilderAndSolverPointerType;
      typedef ConvergenceCriteria< TSparseSpace, TDenseSpace > ConvergenceCriteriaType;
      typedef typename ConvergenceCriteriaType::Pointer ConvergenceCriteriaPointerType;
      typedef ModelPart::NodesContainerType NodesArrayType;

      /// Pointer definition of TurbulenceEvmKEpsilonProcess
      KRATOS_CLASS_POINTER_DEFINITION(TurbulenceEvmKEpsilonProcess);

      ///@}
      ///@name Life Cycle
      ///@{

      /// Constructor
      TurbulenceEvmKEpsilonProcess(ModelPart& rModelPart, Parameters& rParameters, TLinearSolver& rLinearSolver);

      /// Destructor.
      ~TurbulenceEvmKEpsilonProcess() override;


      ///@}
      ///@name Operators
      ///@{


      ///@}
      ///@name Operations
      ///@{

      void Execute() override;

      void ExecuteInitialize() override;

      void ExecuteInitializeSolutionStep() override;

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
      std::string Info() const override;

      /// Print information about this object.
      void PrintInfo(std::ostream& rOStream) const override;

      /// Print object's data.
      void PrintData(std::ostream& rOStream) const override;


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

      ModelPart& mrModelPart;
      Parameters& mrParameters;
      TLinearSolver& mrLinearSolver;
      VariationalDistanceCalculationProcess<TDim, TSparseSpace,TDenseSpace, TLinearSolver> mDistanceCalculator;

      ///@}
      ///@name Private Operators
      ///@{


      ///@}
      ///@name Private Operations
      ///@{

      void CalculateWallDistances();

      // void AssignBoundaryConditions();

      ///@}
      ///@name Private  Access
      ///@{


      ///@}
      ///@name Private Inquiry
      ///@{


      ///@}
      ///@name Un accessible methods
      ///@{

      /// Assignment operator.
      TurbulenceEvmKEpsilonProcess& operator=(TurbulenceEvmKEpsilonProcess<TSparseSpace, TDenseSpace, TLinearSolver> const& rOther);

      /// Copy constructor.
      TurbulenceEvmKEpsilonProcess(TurbulenceEvmKEpsilonProcess<TSparseSpace, TDenseSpace, TLinearSolver> const& rOther);


      ///@}

    }; // Class TurbulenceEvmKEpsilonProcess

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  /// output stream function
  template< class TSparseSpace, class TDenseSpace, class TLinearSolver >
  inline std::ostream& operator << (
      std::ostream& rOStream,
      const TurbulenceEvmKEpsilonProcess<TSparseSpace, TDenseSpace, TLinearSolver>& rThis);

  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_TURBULENCE_K_EPSILON_PROCESS_H_INCLUDED  defined
