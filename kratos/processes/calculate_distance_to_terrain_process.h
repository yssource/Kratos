//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//


#if !defined(KRATOS_CALCULATE_DISTANCE_TO_TERRAIN_PROCESS_H_INCLUDED )
#define  KRATOS_CALCULATE_DISTANCE_TO_TERRAIN_PROCESS_H_INCLUDED



// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "processes/find_intersected_geometrical_objects_process.h"
#include "processes/calculate_distance_to_skin_process.h"
#include "includes/kratos_parameters.h"

namespace Kratos
{
  ///@addtogroup Kratos Core
  ///@{

  ///@name Kratos Classes
  ///@{

  /// Calculates the nodal distances using elemental discontinuous distances.
  /** This class calculates the nodal distances as a minimum elemental distances connected to it.
  */
  class KRATOS_API(KRATOS_CORE) CalculateDistanceToTerrainProcess : public CalculateDistanceToSkinProcess
    {
    public:
      ///@name Type Definitions
      ///@{

      /// Pointer definition of CalculateDistanceToTerrainProcess
      KRATOS_CLASS_POINTER_DEFINITION(CalculateDistanceToTerrainProcess);

      //TODO: These using statements have been included to make the old functions able to compile. It is still pending to update them.
      using ConfigurationType = Internals::DistanceSpatialContainersConfigure;
	  using CellType = OctreeBinaryCell<ConfigurationType>;
	  using OctreeType = OctreeBinary<CellType>;
	  using CellNodeDataType = ConfigurationType::cell_node_data_type;

      ///@}
      ///@name Life Cycle
      ///@{

	  /// Constructor to be used.
	  CalculateDistanceToTerrainProcess(ModelPart& rVolumePart, ModelPart& rTerrainPart, Parameters rParameters);

	  /// Destructor.
      ~CalculateDistanceToTerrainProcess() override;

	  ///@}
	  ///@name Deleted
	  ///@{

      /// Default constructor.
      CalculateDistanceToTerrainProcess() = delete;;

	  /// Copy constructor.
	  CalculateDistanceToTerrainProcess(CalculateDistanceToTerrainProcess const& rOther) = delete;

	  /// Assignment operator.
	  CalculateDistanceToTerrainProcess& operator=(CalculateDistanceToTerrainProcess const& rOther) = delete;

	  ///@}
      ///@name Operations
      ///@{
      virtual double DistancePositionInSpace(double* pCoords); //TODO: This method has been adapted from the previous implementation. It is still pending to update it.

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
      ///@name Un accessible methods
      ///@{



      ///@}

    }; // Class CalculateDistanceToTerrainProcess

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream,
				    CalculateDistanceToTerrainProcess& rThis);

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream,
				    const CalculateDistanceToTerrainProcess& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_CALCULATE_DISTANCE_TO_TERRAIN_PROCESS_H_INCLUDED  defined
