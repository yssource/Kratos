
//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher
//

// System includes
#include <vector>

// Project includes
#include "testing/define.h"
#include "testing/model_part.h"
#include "testing/testing.h"

namespace Kratos {
namespace Testing {
namespace MapperInterfaceInfoTeestingUtils {

typedef std::vector<std::array<double,3>> NodalCoordsContainerType;
typedef std::vector<int> EqIdVectorType;

void CreateNodesModelPart(
    const std::size_t NumNodes,
    ModelPart& rModelPart,
    NodalCoordsContainerType& rNodalCoords,
    EqIdVectorType& rEquationIds);

void ExecuteMapperInterfaceInfoTest();

}  // namespace MapperInterfaceInfoTeestingUtils
}  // namespace Testing
}  // namespace Kratos