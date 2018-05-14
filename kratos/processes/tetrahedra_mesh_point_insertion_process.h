//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Carlos A. Roig
//

#if !defined(KRATOS_TETRAHEDRA_MESH_QUALITY_POINT_INSERTION_PROCESS_H_INCLUDED )
#define  KRATOS_TETRAHEDRA_MESH_QUALITY_POINT_INSERTION_PROCESS_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"

namespace Kratos {
///@addtogroup Kratos Core
///@{

///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
*/
class KRATOS_API(KRATOS_CORE) TetrahedraMeshQualityPointInsertionProcess {
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of TetrahedraMeshQualityPointInsertionProcess
    KRATOS_CLASS_POINTER_DEFINITION(TetrahedraMeshQualityPointInsertionProcess);

    typedef Node<3> NodeType;
    typedef WeakPointerVector< Node<3> > NeighboursVectorType;
    typedef std::vector<Point > PointsVectorType;

    ///@}
    ///@name Flags
    ///@{

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor takes the modelpart to apply smoothing to its mesh 0.
    TetrahedraMeshQualityPointInsertionProcess(ModelPart& rModelPart);

    /// Destructor.
    ~TetrahedraMeshQualityPointInsertionProcess() override;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    void Execute() override;

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

    void CalculateInsertPoint();

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
    TetrahedraMeshQualityPointInsertionProcess& operator=(TetrahedraMeshQualityWeightedSmoothingProcess const& rOther);

    /// Copy constructor.
    TetrahedraMeshQualityPointInsertionProcess(TetrahedraMeshQualityWeightedSmoothingProcess const& rOther);


    ///@}

}; // Class TetrahedraMeshQualityPointInsertionProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
          TetrahedraMeshQualityPointInsertionProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
          const TetrahedraMeshQualityPointInsertionProcess& rThis)
  {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
  }
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_TETRAHEDRA_MESH_QUALITY_POINT_INSERTION_PROCESS_H_INCLUDED  defined
