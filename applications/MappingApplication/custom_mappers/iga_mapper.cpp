//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher, Jordi Cotela
//
// See Master-Thesis P.Bucher
// "Development and Implementation of a Parallel
//  Framework for Non-Matching Grid Mapping"

// System includes

// External includes

// Project includes
#include "iga_mapper.h"
#include "custom_utilities/mapper_typedefs.h"
#include "custom_utilities/mapping_matrix_utilities.h"
#include "custom_utilities/mapper_utilities.h"
#include "mapping_application_variables.h"
#ifdef KRATOS_USING_MPI // mpi-parallel compilation
#include "custom_searching/interface_communicator_mpi.h"
#endif

namespace Kratos {

typedef std::size_t IndexType;
typedef std::size_t SizeType;


void IgaLocalSystem::CalculateAll(MatrixType& rLocalMappingMatrix,
                    EquationIdVectorType& rOriginIds,
                    EquationIdVectorType& rDestinationIds,
                    MapperLocalSystem::PairingStatus& rPairingStatus) const
{
    ProcessInfo tmp_process_info;

    std::vector<Vector> iga_information_vector;

    // mpIgaMapperCondition->SetValue(NODE_POINTER, mpNode);

    mpIgaMapperCondition->CalculateOnIntegrationPoints(INITIAL_STRAIN, iga_information_vector, tmp_process_info);

    const IndexType num_sf_values = iga_information_vector[0].size();

    rLocalMappingMatrix.resize(1, num_sf_values, false);
    for (IndexType i=0; i<num_sf_values; ++i) {
        rLocalMappingMatrix(0,i) = iga_information_vector[0][i];
    }
    rOriginIds.resize(num_sf_values);
    for (IndexType i=0; i<num_sf_values; ++i) {
        rOriginIds[i] = iga_information_vector[1][i];
    }

    if (rDestinationIds.size() != 1) rDestinationIds.resize(1);
    rDestinationIds[0] = mpNode->GetValue(INTERFACE_EQUATION_ID);

    // if (mInterfaceInfos.size() > 0) {
    //     rPairingStatus = MapperLocalSystem::PairingStatus::InterfaceInfoFound;

    //     if (rLocalMappingMatrix.size1() != 1 || rLocalMappingMatrix.size2() != 1) {
    //         rLocalMappingMatrix.resize(1, 1, false);
    //     }
    //     if (rOriginIds.size()      != 1) rOriginIds.resize(1);
    //     if (rDestinationIds.size() != 1) rDestinationIds.resize(1);

    //     int nearest_neighbor_id;
    //     double nearest_neighbor_distance;
    //     mInterfaceInfos[0]->GetValue(nearest_neighbor_id, MapperInterfaceInfo::InfoType::Dummy);
    //     mInterfaceInfos[0]->GetValue(nearest_neighbor_distance, MapperInterfaceInfo::InfoType::Dummy);

    //     for (std::size_t i=1; i<mInterfaceInfos.size(); ++i) {
    //         // no check if this InterfaceInfo is an approximation is necessary
    //         // bcs this does not exist for NearestNeighbor
    //         double distance;
    //         mInterfaceInfos[i]->GetValue(distance, MapperInterfaceInfo::InfoType::Dummy);

    //         if (distance < nearest_neighbor_distance) {
    //             nearest_neighbor_distance = distance;
    //             mInterfaceInfos[i]->GetValue(nearest_neighbor_id, MapperInterfaceInfo::InfoType::Dummy);
    //         }
    //     }

    //     rLocalMappingMatrix(0,0) = 1.0;
    //     rOriginIds[0] = nearest_neighbor_id;
    //     KRATOS_DEBUG_ERROR_IF_NOT(mpNode) << "Members are not intitialized!" << std::endl;
    //     rDestinationIds[0] = mpNode->GetValue(INTERFACE_EQUATION_ID);
    // }
    // else ResizeToZero(rLocalMappingMatrix, rOriginIds, rDestinationIds, rPairingStatus);
}

std::string IgaLocalSystem::PairingInfo(const int EchoLevel, const int CommRank) const
{
    // KRATOS_DEBUG_ERROR_IF_NOT(mpNode) << "Members are not intitialized!" << std::endl;

    // std::stringstream buffer;
    // buffer << "IgaMapper based on " << mpNode->Info();
    // if (EchoLevel > 1) // TODO leave here?
    //     buffer << " at Coodinates " << Coordinates()[0] << " | " << Coordinates()[1] << " | " << Coordinates()[2];
    // buffer << " in rank " << CommRank;
    // return buffer.str();
}


template<class TSparseSpace, class TDenseSpace>
void IgaMapper<TSparseSpace, TDenseSpace>::InitializeInterface(Kratos::Flags MappingOptions)
{
    // CreateMapperLocalSystems(mrModelPartFem.GetCommunicator(),
    //                          mMapperLocalSystems);

    // Assumption: All Fem-stuff is local!!! => requires some mpi-exchange
    // => for mortar the geometries should be local (and then some searching is necessary)
    // => for IGA I either have to bring nodal coords + interface-eqn-ids local or
    //    geometries (nodes-coords with interface-eqn-ids) and geometry-type to reconstruct it


    // for i in range(fem_nodes): // only works for projection-based IGA-Mapper, otherwise I have to loop elements (or conditions...?)
    //     create IGA local System (iga_mapper_condition); => constructor takes node, thus I can skip the IGAInterfaceInfo

    mMapperLocalSystems.reserve(mrModelPartFem.NumberOfNodes());

    Condition* p_iga_mapper_cond = &(*mrModelPartIga.ConditionsBegin());

    for (auto& r_node : mrModelPartFem.Nodes()) {
        mMapperLocalSystems.push_back(Kratos::make_unique<IgaLocalSystem>(p_iga_mapper_cond, &r_node));
    }



    BuildMappingMatrix(MappingOptions);
}

/* Performs operations that are needed for Initialization and when the interface is updated (All cases)
I.e. Operations that can be performed several times in the livetime of the mapper
*/
template<class TSparseSpace, class TDenseSpace>
void IgaMapper<TSparseSpace, TDenseSpace>::BuildMappingMatrix(Kratos::Flags MappingOptions)
{
    AssignInterfaceEquationIds(); // Has to be done ever time in case of overlapping interfaces!

    KRATOS_ERROR_IF_NOT(mpIntefaceCommunicator) << "mpIntefaceCommunicator is a nullptr!" << std::endl;

    // const MapperInterfaceInfoUniquePointerType p_ref_interface_info = GetMapperInterfaceInfo();

    // mpIntefaceCommunicator->ExchangeInterfaceData(mrModelPartFem.GetCommunicator(),
    //                                               MappingOptions,
    //                                               p_ref_interface_info);

    const int echo_level = mMapperSettings["echo_level"].GetInt();

    MappingMatrixUtilities::BuildMappingMatrix<TSparseSpace, TDenseSpace>(
        mpMappingMatrix,
        mpInterfaceVectorContainerOrigin->pGetVector(),
        mpInterfaceVectorContainerDestination->pGetVector(),
        mpInterfaceVectorContainerOrigin->GetModelPart(),
        mpInterfaceVectorContainerDestination->GetModelPart(),
        mMapperLocalSystems,
        echo_level);

    if (echo_level > 0) {
        PrintPairingInfo(echo_level);
    }
}

template<>
void IgaMapper<MapperDefinitions::SparseSpaceType,
    MapperDefinitions::DenseSpaceType>::InitializeInterfaceCommunicator()
{
    mpIntefaceCommunicator = Kratos::make_unique<InterfaceCommunicator>(mrModelPartIga,
                                                                        mMapperLocalSystems,
                                                                        mMapperSettings);
}

#ifdef KRATOS_USING_MPI // mpi-parallel compilation
template<>
void IgaMapper<MapperDefinitions::MPISparseSpaceType,
    MapperDefinitions::DenseSpaceType>::InitializeInterfaceCommunicator()
{
    mpIntefaceCommunicator = Kratos::make_unique<InterfaceCommunicatorMPI>(mrModelPartIga,
                                                                           mMapperLocalSystems,
                                                                           mMapperSettings);
}
#endif

template<class TSparseSpace, class TDenseSpace>
void IgaMapper<TSparseSpace, TDenseSpace>::MapInternal(
    const Variable<double>& rOriginVariable,
    const Variable<double>& rDestinationVariable,
    Kratos::Flags MappingOptions)
{
    mpInterfaceVectorContainerOrigin->UpdateSystemVectorFromModelPart(rOriginVariable, MappingOptions);

    TSparseSpace::Mult(
        *mpMappingMatrix,
        mpInterfaceVectorContainerOrigin->GetVector(),
        mpInterfaceVectorContainerDestination->GetVector()); // rQd = rMdo * rQo

    mpInterfaceVectorContainerDestination->UpdateModelPartFromSystemVector(rDestinationVariable, MappingOptions);
}

template<class TSparseSpace, class TDenseSpace>
void IgaMapper<TSparseSpace, TDenseSpace>::MapInternalTranspose(
    const Variable<double>& rOriginVariable,
    const Variable<double>& rDestinationVariable,
    Kratos::Flags MappingOptions)
{
    mpInterfaceVectorContainerDestination->UpdateSystemVectorFromModelPart(rDestinationVariable, MappingOptions);

    TSparseSpace::TransposeMult(
        *mpMappingMatrix,
        mpInterfaceVectorContainerDestination->GetVector(),
        mpInterfaceVectorContainerOrigin->GetVector()); // rQo = rMdo^T * rQd

    mpInterfaceVectorContainerOrigin->UpdateModelPartFromSystemVector(rOriginVariable, MappingOptions);
}

template<class TSparseSpace, class TDenseSpace>
void IgaMapper<TSparseSpace, TDenseSpace>::MapInternal(
    const Variable<array_1d<double, 3>>& rOriginVariable,
    const Variable<array_1d<double, 3>>& rDestinationVariable,
    Kratos::Flags MappingOptions)
{
    const auto& var_x_origin = KratosComponents<ComponentVariableType>::Get(rOriginVariable.Name() + "_X");
    const auto& var_y_origin = KratosComponents<ComponentVariableType>::Get(rOriginVariable.Name() + "_Y");
    const auto& var_z_origin = KratosComponents<ComponentVariableType>::Get(rOriginVariable.Name() + "_Z");

    const auto& var_x_destination = KratosComponents<ComponentVariableType>::Get(rDestinationVariable.Name() + "_X");
    const auto& var_y_destination = KratosComponents<ComponentVariableType>::Get(rDestinationVariable.Name() + "_Y");
    const auto& var_z_destination = KratosComponents<ComponentVariableType>::Get(rDestinationVariable.Name() + "_Z");

    // X-Component
    mpInterfaceVectorContainerOrigin->UpdateSystemVectorFromModelPart(var_x_origin, MappingOptions);

    TSparseSpace::Mult(
        *mpMappingMatrix,
        mpInterfaceVectorContainerOrigin->GetVector(),
        mpInterfaceVectorContainerDestination->GetVector()); // rQd = rMdo * rQo

    mpInterfaceVectorContainerDestination->UpdateModelPartFromSystemVector(var_x_destination, MappingOptions);

    // Y-Component
    mpInterfaceVectorContainerOrigin->UpdateSystemVectorFromModelPart(var_y_origin, MappingOptions);

    TSparseSpace::Mult(
        *mpMappingMatrix,
        mpInterfaceVectorContainerOrigin->GetVector(),
        mpInterfaceVectorContainerDestination->GetVector()); // rQd = rMdo * rQo

    mpInterfaceVectorContainerDestination->UpdateModelPartFromSystemVector(var_y_destination, MappingOptions);

    // Z-Component
    mpInterfaceVectorContainerOrigin->UpdateSystemVectorFromModelPart(var_z_origin, MappingOptions);

    TSparseSpace::Mult(
        *mpMappingMatrix,
        mpInterfaceVectorContainerOrigin->GetVector(),
        mpInterfaceVectorContainerDestination->GetVector()); // rQd = rMdo * rQo

    mpInterfaceVectorContainerDestination->UpdateModelPartFromSystemVector(var_z_destination, MappingOptions);
}

template<class TSparseSpace, class TDenseSpace>
void IgaMapper<TSparseSpace, TDenseSpace>::MapInternalTranspose(
    const Variable<array_1d<double, 3>>& rOriginVariable,
    const Variable<array_1d<double, 3>>& rDestinationVariable,
    Kratos::Flags MappingOptions)
{
    const auto& var_x_origin = KratosComponents<ComponentVariableType>::Get(rOriginVariable.Name() + "_X");
    const auto& var_y_origin = KratosComponents<ComponentVariableType>::Get(rOriginVariable.Name() + "_Y");
    const auto& var_z_origin = KratosComponents<ComponentVariableType>::Get(rOriginVariable.Name() + "_Z");

    const auto& var_x_destination = KratosComponents<ComponentVariableType>::Get(rDestinationVariable.Name() + "_X");
    const auto& var_y_destination = KratosComponents<ComponentVariableType>::Get(rDestinationVariable.Name() + "_Y");
    const auto& var_z_destination = KratosComponents<ComponentVariableType>::Get(rDestinationVariable.Name() + "_Z");

    // X-Component
    mpInterfaceVectorContainerDestination->UpdateSystemVectorFromModelPart(var_x_destination, MappingOptions);

    TSparseSpace::TransposeMult(
        *mpMappingMatrix,
        mpInterfaceVectorContainerOrigin->GetVector(),
        mpInterfaceVectorContainerDestination->GetVector()); // rQo = rMdo^T * rQd

    mpInterfaceVectorContainerOrigin->UpdateModelPartFromSystemVector(var_x_origin, MappingOptions);

    // Y-Component
    mpInterfaceVectorContainerDestination->UpdateSystemVectorFromModelPart(var_y_destination, MappingOptions);

    TSparseSpace::TransposeMult(
        *mpMappingMatrix,
        mpInterfaceVectorContainerOrigin->GetVector(),
        mpInterfaceVectorContainerDestination->GetVector()); // rQo = rMdo^T * rQd

    mpInterfaceVectorContainerOrigin->UpdateModelPartFromSystemVector(var_y_origin, MappingOptions);

    // Z-Component
    mpInterfaceVectorContainerDestination->UpdateSystemVectorFromModelPart(var_z_destination, MappingOptions);

    TSparseSpace::TransposeMult(
        *mpMappingMatrix,
        mpInterfaceVectorContainerOrigin->GetVector(),
        mpInterfaceVectorContainerDestination->GetVector()); // rQo = rMdo^T * rQd

    mpInterfaceVectorContainerOrigin->UpdateModelPartFromSystemVector(var_z_origin, MappingOptions);
}

template<class TSparseSpace, class TDenseSpace>
void IgaMapper<TSparseSpace, TDenseSpace>::ValidateInput(Parameters MapperSettings)
{
    MapperUtilities::CheckInterfaceModelParts(0);
    ValidateParameters(MapperSettings);

    if (mMapperSettings["search_radius"].GetDouble() < 0.0) {
        const double search_radius = MapperUtilities::ComputeSearchRadius(
                                        mrModelPartIga,
                                        mrModelPartFem,
                                        MapperSettings["echo_level"].GetInt());
        mMapperSettings["search_radius"].SetDouble(search_radius);
    }
}

template<class TSparseSpace, class TDenseSpace>
void IgaMapper<TSparseSpace, TDenseSpace>::PrintPairingInfo(const int EchoLevel)
{
    const int comm_rank = mrModelPartFem.GetCommunicator().MyPID();
    std::stringstream warning_msg;

    for (const auto& rp_local_sys : mMapperLocalSystems)
    {
        const auto pairing_status = rp_local_sys->GetPairingStatus();

        if (pairing_status != MapperLocalSystem::PairingStatus::InterfaceInfoFound)
        {
            warning_msg << rp_local_sys->PairingInfo(EchoLevel, comm_rank);

            if (pairing_status == MapperLocalSystem::PairingStatus::Approximation)
                warning_msg << " is using an approximation";
            else if (pairing_status == MapperLocalSystem::PairingStatus::NoInterfaceInfo)
                warning_msg << " has not found a neighbor";

            KRATOS_WARNING("Mapper") << warning_msg.str() << std::endl;

            // reset the stringstream
            warning_msg.str( std::string() );
            warning_msg.clear();
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation
template class IgaMapper< MapperDefinitions::SparseSpaceType, MapperDefinitions::DenseSpaceType >;
#ifdef KRATOS_USING_MPI // mpi-parallel compilation
template class IgaMapper< MapperDefinitions::MPISparseSpaceType, MapperDefinitions::DenseSpaceType >;
#endif

}  // namespace Kratos.
