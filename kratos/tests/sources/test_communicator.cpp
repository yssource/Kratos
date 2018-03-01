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
//

// System includes

// External includes

// Project includes
#include "testing/testing.h"
#include "includes/mesh.h"
#include "includes/communicator.h"
#include "geometries/triangle_2d_3.h"

#ifdef KRATOS_USING_MPI
    #include "includes/mpi_communicator.h"
#endif

namespace Kratos {
namespace Testing {

#ifdef KRATOS_USING_MPI

	// KRATOS_TEST_CASE_IN_SUITE(TransferElement, KratosCoreFastSuite)
	// {
    //     typedef Triangle2D3<Node<3>>                        GeometryType;
    //     typedef PointerVectorSet<Element, IndexedObject>    ElementContainerType;

    //     int mpi_rank, mpi_size;

    //     MPI_Init(nullptr, nullptr);	
	//     MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    //     MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    
    //     Element::Pointer localElementLeft = Element::Pointer(new Element(
    //         0, 
    //         GeometryType::Pointer( 
    //             new GeometryType(
    //                 Node<3>::Pointer(new Node<3>(500, 0.1, 0.1, 0.1)),
    //                 Node<3>::Pointer(new Node<3>(501, 0.1, 1.1, 0.1)),
    //                 Node<3>::Pointer(new Node<3>(502, 1.1, 0.6, 0.1))
    //             )
    //         )
    //     ));

    //     Element::Pointer localElementRight = Element::Pointer(new Element(
    //         1, 
    //         GeometryType::Pointer( 
    //             new GeometryType(
    //                 Node<3>::Pointer(new Node<3>(502, 1.1, 0.6, 0.1)),
    //                 Node<3>::Pointer(new Node<3>(504, 2.1, 1.1, 0.1)),
    //                 Node<3>::Pointer(new Node<3>(503, 2.1, 0.1, 0.1))
    //             )
    //         )
    //     ));

    //     std::vector<ElementContainerType> sendObjects(mpi_size, ElementContainerType());
    //     std::vector<ElementContainerType> recvObjects(mpi_size, ElementContainerType());

    //     // Transfer the object of process 1 to 0 and 0 to 1
    //     sendObjects[0].push_back(localElementRight); 
    //     sendObjects[1].push_back(localElementLeft);

    //     MPICommunicator communicator(new VariablesList());

    //     communicator.TransferObjects(sendObjects, recvObjects);

    //     if(mpi_rank == 0) {
    //         auto recvElementGeom = recvObjects[1][1].GetGeometry(); // The received element geometry
            
    //         // First Point
    //         KRATOS_CHECK_EQUAL(recvElementGeom(0)->Id(), 502);
    //         KRATOS_CHECK_EQUAL(recvElementGeom(0)->X(), 1.1);
    //         KRATOS_CHECK_EQUAL(recvElementGeom(0)->Y(), 0.6);
    //         KRATOS_CHECK_EQUAL(recvElementGeom(0)->Z(), 0.1);

    //         // Second Point
    //         KRATOS_CHECK_EQUAL(recvElementGeom(1)->Id(), 504);
    //         KRATOS_CHECK_EQUAL(recvElementGeom(1)->X(), 2.1);
    //         KRATOS_CHECK_EQUAL(recvElementGeom(1)->Y(), 1.1);
    //         KRATOS_CHECK_EQUAL(recvElementGeom(1)->Z(), 0.1);

    //         // Third Point
    //         KRATOS_CHECK_EQUAL(recvElementGeom(2)->Id(), 503);
    //         KRATOS_CHECK_EQUAL(recvElementGeom(2)->X(), 2.1);
    //         KRATOS_CHECK_EQUAL(recvElementGeom(2)->Y(), 0.1);
    //         KRATOS_CHECK_EQUAL(recvElementGeom(2)->Z(), 0.1);

    //         // Check that check_equal can check points ( this must always hold )
    //         KRATOS_CHECK_EQUAL(localElementLeft->GetGeometry()(0), localElementLeft->GetGeometry()(0));

    //         // Check the common point is the same
    //         KRATOS_CHECK_EQUAL(recvElementGeom(0), localElementLeft->GetGeometry()(2));
    //     }

    //     if(mpi_rank == 1) {
    //         auto recvElementGeom = recvObjects[0][0].GetGeometry(); // The received element geometry
            
    //         // First Point
    //         KRATOS_CHECK_EQUAL(recvElementGeom(0)->Id(), 500);
    //         KRATOS_CHECK_EQUAL(recvElementGeom(0)->X(), 0.1);
    //         KRATOS_CHECK_EQUAL(recvElementGeom(0)->Y(), 0.1);
    //         KRATOS_CHECK_EQUAL(recvElementGeom(0)->Z(), 0.1);

    //         // Second Point
    //         KRATOS_CHECK_EQUAL(recvElementGeom(1)->Id(), 501);
    //         KRATOS_CHECK_EQUAL(recvElementGeom(1)->X(), 0.1);
    //         KRATOS_CHECK_EQUAL(recvElementGeom(1)->Y(), 1.1);
    //         KRATOS_CHECK_EQUAL(recvElementGeom(1)->Z(), 0.1);

    //         // Third Point
    //         KRATOS_CHECK_EQUAL(recvElementGeom(2)->Id(), 502);
    //         KRATOS_CHECK_EQUAL(recvElementGeom(2)->X(), 1.1);
    //         KRATOS_CHECK_EQUAL(recvElementGeom(2)->Y(), 0.6);
    //         KRATOS_CHECK_EQUAL(recvElementGeom(2)->Z(), 0.1);

    //         // Check that check_equal can check points ( this must always hold )
    //         KRATOS_CHECK_EQUAL(localElementRight->GetGeometry()(0), localElementRight->GetGeometry()(0));

    //         // Check the common point is the same
    //         KRATOS_CHECK_EQUAL(recvElementGeom(2), localElementRight->GetGeometry()(0));
    //     }

    //     MPI_Finalize();
    // }

    KRATOS_TEST_CASE_IN_SUITE(TransferMesh, KratosCoreFastSuite)
	{
        int mpi_rank, mpi_size;

        MPI_Init(nullptr, nullptr);	
	    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
        MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

        ModelPart model_part("Test");

        // Set MPICommunicator as modelpart's communicator
        model_part.AddNodalSolutionStepVariable(VELOCITY);
        VariablesList * variables_list = &model_part.GetNodalSolutionStepVariablesList();

        std::cout << mpi_rank << " - Test point B" << std::endl;
        MPI_Barrier(MPI_COMM_WORLD);

        if(mpi_rank == 0) {
            model_part.CreateNewNode(500, 0.1, 0.1, 0.1);
            model_part.CreateNewNode(501, 0.1, 1.1, 0.1);
            model_part.CreateNewNode(502, 1.1, 0.6, 0.1);
        }

        std::cout << mpi_rank << " - Test point C" << std::endl;
        MPI_Barrier(MPI_COMM_WORLD);

        if(mpi_rank == 1) {
            model_part.CreateNewNode(502, 1.1, 0.6, 0.1);
            model_part.CreateNewNode(504, 2.1, 1.1, 0.1);
            model_part.CreateNewNode(503, 2.1, 0.1, 0.1);
        }

        MPICommunicator * mpi_communicator = new MPICommunicator(variables_list);
        model_part.SetCommunicator(Communicator::Pointer(mpi_communicator));
        model_part.pGetCommunicator()->LocalMesh().Nodes() = model_part.pGetMesh(0)->Nodes();

        std::cout << mpi_rank << " - Test point D" << std::endl;
        MPI_Barrier(MPI_COMM_WORLD);

        std::vector<ModelPart::MeshType::Pointer> send_meshes(mpi_size);
        std::vector<ModelPart::MeshType::Pointer> recv_meshes(mpi_size);

        if(mpi_rank == 0) { send_meshes[1] = mpi_communicator->pLocalMesh(); }
        if(mpi_rank == 1) { send_meshes[0] = mpi_communicator->pLocalMesh(); }

        if(mpi_rank == 0) { recv_meshes[1] = ModelPart::MeshType::Pointer(new ModelPart::MeshType()); }
        if(mpi_rank == 1) { recv_meshes[0] = ModelPart::MeshType::Pointer(new ModelPart::MeshType()); }

        std::cout << mpi_rank << " - Test point E" << std::endl;
        MPI_Barrier(MPI_COMM_WORLD);

        mpi_communicator->TransferMesh(send_meshes, recv_meshes);

        std::cout << mpi_rank << " - Test point F" << std::endl;
        MPI_Barrier(MPI_COMM_WORLD);

        if(mpi_rank == 0) {
            // Check the common point is the same
            KRATOS_CHECK_EQUAL(
                &(send_meshes[1]->GetNode(502)), 
                &(recv_meshes[1]->GetNode(502))
            );
        }

        // if(mpi_rank == 1) {
        //     // Check the common point is the same
        //     KRATOS_CHECK_EQUAL(
        //         &(send_meshes[0]->GetNode(502)), 
        //         &(recv_meshes[0]->GetNode(502))
        //     );
        // }

        MPI_Finalize();
    }

#endif

} // namespace Testing
} // namespace Kratos

