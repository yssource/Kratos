//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Carlos Roig
//
//

// System includes


// External includes


// Project includes
#include "includes/node.h"
#include "testing/testing.h"
#include "spatial_containers/bins_static_objects.h"
#include "spatial_containers/configures/node_configure.h"

namespace Kratos {
namespace Testing {

typedef BinsObjectStatic<NodeConfigure>            ObjectStaticBins;

// typedef ObjectStaticBins::PointType                PointType;
typedef NodeConfigure::ObjectType                   ObjectType;
typedef ObjectType::Pointer                         ObjectTypePointer;
typedef ObjectStaticBins::ContainerType            ObjectContainerType;
// typedef ObjectStaticBins::PointVector              PointVector;
typedef ObjectStaticBins::IteratorType             IteratorType;
typedef std::vector<double>                         DistanceVector;
typedef std::vector<double>::iterator               DistanceIterator;

typedef ObjectStaticBins::CoordinateType           CoordinateType;
typedef ObjectStaticBins::PointerType              PointerType;
typedef ObjectStaticBins::SearchStructureType      SearchStructureType;

/**
 * @brief Test that the bins is constructed correctly
 * 
 */
KRATOS_TEST_CASE_IN_SUITE(ObjectStaticBinsDefaultConstructorBoundigBox, KratosCoreFastSuite)
{
    ObjectContainerType objects;

    for(std::size_t i = 0; i < 10; i++) {
        objects.push_back(ObjectTypePointer(new ObjectType(i, i, i, i)));
    }

    ObjectStaticBins testBins(objects.begin(), objects.end());

    // This bins uses "epsilon = (maxpoint - minpoint) * 0.01" as a safe marging in the bound box calculation
    // ( not for the static bins but we should use it )
    double epsilon = 0.0;

    KRATOS_CHECK_EQUAL(testBins.GetMinPoint()[0], 0.0 - epsilon);
    KRATOS_CHECK_EQUAL(testBins.GetMinPoint()[1], 0.0 - epsilon);
    KRATOS_CHECK_EQUAL(testBins.GetMinPoint()[2], 0.0 - epsilon);

    KRATOS_CHECK_EQUAL(testBins.GetMaxPoint()[0], 9.0 + epsilon);
    KRATOS_CHECK_EQUAL(testBins.GetMaxPoint()[1], 9.0 + epsilon);
    KRATOS_CHECK_EQUAL(testBins.GetMaxPoint()[2], 9.0 + epsilon);
}

/**
 * @brief Test that the number of cells is calculated correctly
 * 
 */
KRATOS_TEST_CASE_IN_SUITE(ObjectStaticBinsDefaultConstructorCellNumber, KratosCoreFastSuite)
{
    ObjectContainerType objects;

    for(std::size_t i = 0; i < 10; i++) {
        objects.push_back(ObjectTypePointer(new ObjectType(i, i, i, i)));
    }

    ObjectStaticBins testBins(objects.begin(), objects.end());

    KRATOS_CHECK_EQUAL(testBins.GetDivisions()[0], 3);
    KRATOS_CHECK_EQUAL(testBins.GetDivisions()[1], 3);
    KRATOS_CHECK_EQUAL(testBins.GetDivisions()[2], 3);
}

/**
 * @brief Test that the size of the cells is calculated correctly
 * 
 */
KRATOS_TEST_CASE_IN_SUITE(ObjectStaticBinsDefaultConstructorCellSize, KratosCoreFastSuite)
{
    ObjectContainerType objects;

    for(std::size_t i = 0; i < 10; i++) {
        objects.push_back(ObjectTypePointer(new ObjectType(i, i, i, i)));
    }

    ObjectStaticBins testBins(objects.begin(), objects.end());

    // This bins uses "epsilon = (maxpoint - minpoint) * 0.01" as a safe marging in the bound box calculation
    // ( not for the static bins but we should use it )
    double epsilon = 0.0;

    KRATOS_CHECK_EQUAL(testBins.GetCellSize()[0], 3.0 + epsilon * 2 / 3);
    KRATOS_CHECK_EQUAL(testBins.GetCellSize()[1], 3.0 + epsilon * 2 / 3);
    KRATOS_CHECK_EQUAL(testBins.GetCellSize()[2], 3.0 + epsilon * 2 / 3);
}

/**
 * @brief Test that the bins is constructed correctly with the cell size constructor
 * 
 */
// KRATOS_TEST_CASE_IN_SUITE(ObjectStaticBinsCellSizeConstructorBoundingBox, KratosCoreFastSuite)
// {
//     ObjectContainerType objects;

//     for(std::size_t i = 0; i < 10; i++) {
//         objects.push_back(ObjectTypePointer(new ObjectType(i, i, i, i)));
//     }

//     CoordinateType testCellSize = 3.2244;

//     ObjectStaticBins testBins(objects.begin(), objects.end(), testCellSize);

//     // This bins uses "epsilon = (maxpoint - minpoint) * 0.01" as a safe marging in the bound box calculation
//     double epsilon = 0.01 * (9 - 0);

//     KRATOS_CHECK_EQUAL(testBins.GetMinPoint()[0], 0.0 - epsilon);
//     KRATOS_CHECK_EQUAL(testBins.GetMinPoint()[1], 0.0 - epsilon);
//     KRATOS_CHECK_EQUAL(testBins.GetMinPoint()[2], 0.0 - epsilon);

//     KRATOS_CHECK_EQUAL(testBins.GetMaxPoint()[0], 9.0 + epsilon);
//     KRATOS_CHECK_EQUAL(testBins.GetMaxPoint()[1], 9.0 + epsilon);
//     KRATOS_CHECK_EQUAL(testBins.GetMaxPoint()[2], 9.0 + epsilon);
// }

/**
 * @brief Test that the number of cells is calculated correctly with the cell size constructor
 * 
 */
// KRATOS_TEST_CASE_IN_SUITE(ObjectStaticBinsCellSizeConstructorCellNumber, KratosCoreFastSuite)
// {
//     ObjectContainerType objects;

//     for(std::size_t i = 0; i < 10; i++) {
//         objects.push_back(ObjectTypePointer(new ObjectType(i, i, i, i)));
//     }

//     CoordinateType testCellSize = 3.2244;

//     ObjectStaticBins testBins(objects.begin(), objects.end(), testCellSize);

//     KRATOS_CHECK_EQUAL(testBins.GetDivisions()[0], 3);
//     KRATOS_CHECK_EQUAL(testBins.GetDivisions()[1], 3);
//     KRATOS_CHECK_EQUAL(testBins.GetDivisions()[2], 3);
// }

/**
 * @brief Test that the size of the cells is calculated correctly with the cell size constructor
 * 
 */
// KRATOS_TEST_CASE_IN_SUITE(ObjectStaticBinsCellSizeConstructorCellSize, KratosCoreFastSuite)
// {
//     ObjectContainerType objects;

//     for(std::size_t i = 0; i < 10; i++) {
//         objects.push_back(ObjectTypePointer(new ObjectType(i, i, i, i)));
//     }

//     CoordinateType testCellSize = 3.2244;

//     ObjectStaticBins testBins(objects.begin(), objects.end(), testCellSize);

//     KRATOS_CHECK_EQUAL(testBins.GetCellSize()[0], testCellSize);
//     KRATOS_CHECK_EQUAL(testBins.GetCellSize()[1], testCellSize);
//     KRATOS_CHECK_EQUAL(testBins.GetCellSize()[2], testCellSize);
// }

/**
 * @brief Searches the nearest point ( not implemented for bins_dynamic_objects )
 * 
 */
// KRATOS_TEST_CASE_IN_SUITE(ObjectStaticBinsExistPoint, KratosCoreFastSuite)
// {
//     ObjectContainerType objects;

//     for(std::size_t i = 0; i < 10; i++) {
//         objects.push_back(ObjectTypePointer(new ObjectType(i, i, i, i)));
//     }

//     ObjectStaticBins testBins(objects.begin(), objects.end());

//     PointerType nearestPoint = testBins.ExistPoint(PointerType(new ObjectType(0, 4.1, 4.1, 4.1)));

//     KRATOS_CHECK_EQUAL(nearestPoint->Id(), 4);
// }

/**
 * @brief Searches the nearest point (excluding the input)
 * 
 */
// Missing function / does not exists yet.

/**
 * @brief Searches the nearest point with (including the input) ( not implemented for bins_dynamic_objects )
 * 
 */
// KRATOS_TEST_CASE_IN_SUITE(ObjectStaticBinsNearestPoint, KratosCoreFastSuite)
// {
//     ObjectContainerType objects;

//     for(std::size_t i = 0; i < 10; i++) {
//         objects.push_back(ObjectTypePointer(new ObjectType(i, i, i, i)));
//     }

//     PointerType pointToSearch = PointerType(new ObjectType(10, 4.25, 4.25, 4.25));
//     objects.push_back(pointToSearch);
    
//     ObjectStaticBins testBins(objects.begin(), objects.end());

//     PointerType nearestPoint = testBins.SearchNearestPoint(*pointToSearch);
 
//     KRATOS_CHECK_EQUAL(nearestPoint->Id(), 10);
// }

/**
 * @brief Searches the nearest point (including the input) with distance ( not implemented for bins_dynamic_objects )
 * 
 */
// KRATOS_TEST_CASE_IN_SUITE(ObjectStaticBinsNearestPointWithDistance, KratosCoreFastSuite)
// {
//     ObjectContainerType objects;

//     for(std::size_t i = 0; i < 10; i++) {
//         objects.push_back(ObjectTypePointer(new ObjectType(i, i, i, i)));
//     }

//     PointerType pointToSearch = PointerType(new ObjectType(10, 4.25, 4.25, 4.25));
    
//     ObjectStaticBins testBins(objects.begin(), objects.end());

//     double squaredDistance = 0.0;
//     PointerType nearestPoint = testBins.SearchNearestPoint(*pointToSearch, squaredDistance);
 
//     KRATOS_CHECK_EQUAL(nearestPoint->Id(), 4);
//     KRATOS_CHECK_EQUAL(squaredDistance, 0.1875);
// }

/**
 * @brief Searches the nearest point (including the input) with distance (threadsafe) ( not implemented for bins_dynamic_objects )
 * 
 */
// KRATOS_TEST_CASE_IN_SUITE(ObjectStaticBinsNearestPointWithDistanceThreadsafe, KratosCoreFastSuite)
// {
//     ObjectContainerType objects;

//     for(std::size_t i = 0; i < 10; i++) {
//         objects.push_back(ObjectTypePointer(new ObjectType(i, i, i, i)));
//     }

//     PointerType pointToSearch = PointerType(new ObjectType(10, 4.25, 4.25, 4.25));
    
//     ObjectStaticBins testBins(objects.begin(), objects.end());
//     SearchStructureType searchBox;

//     double squaredDistance = 0.0;
//     PointerType nearestPoint = testBins.SearchNearestPoint(*pointToSearch, squaredDistance, searchBox);
 
//     KRATOS_CHECK_EQUAL(nearestPoint->Id(), 4);
//     KRATOS_CHECK_EQUAL(squaredDistance, 0.1875);
// }

    
} // namespace Testesing
} // namespace Kratos