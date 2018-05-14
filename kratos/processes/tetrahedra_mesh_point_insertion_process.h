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

#if !defined(KRATOS_TETRAHEDRA_MESH_POINT_INSERTION_PROCESS_H_INCLUDED )
#define  KRATOS_TETRAHEDRA_MESH_POINT_INSERTION_PROCESS_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "processes/process.h"

namespace Kratos {
///@addtogroup Kratos Core
///@{

///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
*/
class KRATOS_API(KRATOS_CORE) TetrahedraMeshPointInsertionProcess: public Process {
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of TetrahedraMeshPointInsertionProcess
    KRATOS_CLASS_POINTER_DEFINITION(TetrahedraMeshPointInsertionProcess);

    typedef Node<3>                     NodeType;
    typedef WeakPointerVector<Node<3>>  NeighboursVectorType;
    typedef std::vector<Point>          PointsVectorType;

    struct Graph {
		std::unordered_map<std::size_t, std::size_t> nodes;
		std::unordered_map<std::size_t, std::pair<std::size_t, std::size_t>> edges;
		std::unordered_map<std::size_t, std::list<std::size_t>> connect;
	};

    ///@}
    ///@name Flags
    ///@{

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor takes the modelpart to apply smoothing to its mesh 0.
    TetrahedraMeshPointInsertionProcess(ModelPart& rModelPart);

    /// Destructor.
    ~TetrahedraMeshPointInsertionProcess() override;

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

    ModelPart & mrModelPart;


    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    void CalculateInsertPoint(Element::Pointer & pElement, Point & rInsertPoint);
    void CalculateCut(ModelPart & rModelPart, Element::Pointer & pElement, Point & rInsertPoint, Graph & meshCut, std::size_t depth);

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
    TetrahedraMeshPointInsertionProcess& operator=(TetrahedraMeshPointInsertionProcess const& rOther);

    /// Copy constructor.
    TetrahedraMeshPointInsertionProcess(TetrahedraMeshPointInsertionProcess const& rOther);


    ///@}

}; // Class TetrahedraMeshPointInsertionProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
          TetrahedraMeshPointInsertionProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
          const TetrahedraMeshPointInsertionProcess& rThis)
  {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
  }
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_TETRAHEDRA_MESH_POINT_INSERTION_PROCESS_H_INCLUDED  defined
