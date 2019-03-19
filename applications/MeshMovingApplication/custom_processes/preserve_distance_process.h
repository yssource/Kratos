//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Andreas Winterstein
//
//


#ifndef PRESERVE_DISTANCE_PROCESS_H
#define PRESERVE_DISTANCE_PROCESS_H

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
//#include "includes/node.h"
//#include "includes/define.h"
#include "processes/process.h"
#include "utilities/math_utils.h"
//#include "includes/kratos_parameters.h"
#include "spatial_containers/spatial_containers.h"
#include "includes/model_part.h"
#include "includes/mesh_moving_variables.h"

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

/**
 * @class PreserveDistanceProcess
 *
 * @ingroup MeshMovingApplication
 *
 * @brief This method preserves the distance between two meshes
 *
 * @author Andreas Winterstein
*/

namespace Kratos
{
class KRATOS_API(MESH_MOVING_APPLICATION) PreserveDistanceProcess
    : public Process
{
  public:

    typedef Node < 3 > NodeType;
    typedef Node < 3 > ::Pointer NodeTypePointer;
    typedef std::vector<NodeTypePointer> NodeVector;
    typedef ModelPart::NodesContainerType NodesArrayType;
    typedef std::vector<NodeTypePointer>::iterator NodeIterator;
    typedef std::vector<double> DoubleVector;
    typedef DoubleVector::iterator DoubleVectorIterator;
    typedef std::size_t SizeType;
    typedef ModelPart::VariableComponentType VariableComponentType;
    typedef ModelPart::MasterSlaveConstraintType::Pointer ConstraintPointer;


    // Type definitions for tree-search
    typedef Bucket< 3, NodeType, NodeVector, NodeTypePointer, NodeIterator, DoubleVectorIterator > BucketType;
    typedef Tree< KDTreePartition<BucketType> > KDTree;


    /// Pointer definition of ApplyMultipointConstraintsProcess
    KRATOS_CLASS_POINTER_DEFINITION(PreserveDistanceProcess);

    /// Constructor.
    PreserveDistanceProcess(ModelPart &rModelPart, ModelPart &rMainModelPart,
     Parameters InputParameters):mrModelPart(rModelPart),mrMainModelPart(rMainModelPart),mParameters(InputParameters)
    {
        KRATOS_TRY;
        Parameters default_parameters = Parameters(R"(
        {
            "constraint_set_name"           : "LinearMasterSlaveConstraint",
            "master_sub_model_part_name"    : "master_connect",
            "slave_sub_model_part_name"     : "slave_connect",
            "variable_names"                : ["MESH_DISPLACEMENT_X","MESH_DISPLACEMENT_Y", "MESH_DISPLACEMENT_Z"],
            "debug_info"                    : false,
            "must_find_neighbor"            : true,
            "neighbor_search_radius"        : 0.40,
            "bucket_size"                   : 10
        })" );
        default_parameters.ValidateAndAssignDefaults(InputParameters);
        KRATOS_CATCH("")
    }

    void ExecuteInitialize() override
    {
      KRATOS_TRY;
      this->CoupleModelParts();
      std::cout<<"|||||||||||||||||||||Calling preserve distance process"<<std::endl;
      KRATOS_CATCH("");
    }



    void CoupleModelParts(){
      KRATOS_TRY;
      DoubleVector slave_nodes = {302,306,310,315,320,297};
      DoubleVector master_nodes = {143,147,152,158,162,136};

      NodesArrayType &r_nodes  = mrMainModelPart.Nodes();
      //NodesArrayType &r_nodes_slave  = mrMainModelPart.Nodes();

      int current_id = 1;
      for(unsigned int i=0;i<master_nodes.size();++i)
        {
          VariableComponentType current_dof =
          KratosComponents<VariableComponentType>::Get("MESH_DISPLACEMENT_X");
          const ConstraintPointer p_current_constraint = mrMainModelPart.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint",current_id,
                r_nodes[master_nodes[i]],current_dof,r_nodes[slave_nodes[i]],
                current_dof,1.0,0.0);
          current_id+=1;
        }

      for(unsigned int i=0;i<master_nodes.size();++i)
        {
          VariableComponentType current_dof =
          KratosComponents<VariableComponentType>::Get("MESH_DISPLACEMENT_Y");
          const ConstraintPointer p_current_constraint = mrMainModelPart.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint",current_id,
                r_nodes[master_nodes[i]],current_dof,r_nodes[slave_nodes[i]],
                current_dof,1.0,0.0);
          current_id+=1;
        }


      FindNeighbourNodes();
      KRATOS_CATCH("");
    }


     

    void CreateSubModelPartNodesVector(NodeVector& rMasterNodeVector, std::string SubModelPartName)
    {
        KRATOS_TRY;
        ModelPart &model_part = mrModelPart.GetSubModelPart(mParameters[SubModelPartName].GetString());
        NodesArrayType &r_nodes = model_part.Nodes();

        rMasterNodeVector.resize(r_nodes.size());
        auto i_begin = model_part.NodesBegin();

        for (SizeType i(0);i<r_nodes.size();++i)
        {
            NodeTypePointer pnode = *((i_begin+i).base());
            rMasterNodeVector[i] = pnode;
        }
        KRATOS_CATCH("");
    }



/*
    void CoupleModelParts()
    {
      KRATOS_TRY;
      ModelPart &master_model_part    = mrModelPart.GetSubModelPart(mParameters["master_sub_model_part_name"].GetString());
      ModelPart &slave_model_part     = mrModelPart.GetSubModelPart(mParameters["slave_sub_model_part_name"].GetString());
      const int bucket_size           = mParameters["bucket_size"].GetInt();
      NodesArrayType &r_nodes_master  = master_model_part.Nodes();
      NodesArrayType &r_nodes_slave   = slave_model_part.Nodes();

      NodeVector master_node_list(r_nodes_master.size());
      this->CreateListOfNodesOfMasterSubModelPart(master_node_list);

      KDTree::Pointer search_tree =
        Kratos::shared_ptr<KDTree>(new KDTree(master_node_list.begin(), master_node_list.end(), bucket_size));

      const int max_number_of_neighbors = 1;

      for(NodeType& node_i : r_nodes_slave)
      {
        double neighbor_search_radius = mParameters["neighbor_search_radius"].GetDouble();
        SizeType number_of_neighbors = 0;
        NodeVector neighbor_nodes( max_number_of_neighbors );
        DoubleVector resulting_squared_distances( max_number_of_neighbors );

        neighbor_nodes.clear();
        resulting_squared_distances.clear();
        //find nodal neighbors
        number_of_neighbors = search_tree->SearchInRadius( node_i,
                                                           neighbor_search_radius,
                                                           neighbor_nodes.begin(),
                                                           resulting_squared_distances.begin(),
                                                           max_number_of_neighbors );

        if (mParameters["must_find_neighbor"].GetBool())
        {
            if (number_of_neighbors<1)
                {
                    KRATOS_ERROR << "found no neighbor for slave node " << node_i.Id() << " " << node_i.Coordinates() << std::endl;
                }
        }

        if (number_of_neighbors>0)
        {
            if(mParameters["debug_info"].GetBool()) std::cout << "nr.ne.: " << number_of_neighbors << std::endl;
            //set all weights to 1 since we want to preserve the distance between the nodes
            DoubleVector list_of_weights( number_of_neighbors, 1.0 );

            //this->CalculateNodalWeights(resulting_squared_distances,list_of_weights,number_of_neighbors);
            this->CoupleSlaveToNeighborMasterNodes(node_i,neighbor_nodes,list_of_weights,number_of_neighbors);
        }

      }
      KRATOS_CATCH("");
    }

    void CoupleSlaveToNeighborMasterNodes(const NodeType& rCurrentSlaveNode,
      const NodeVector& rNeighborNodes, const DoubleVector& rNodalNeighborWeights,
      const SizeType& rNumberOfNeighbors)
    {
        ModelPart &master_model_part    = mrModelPart.GetSubModelPart(mParameters["master_sub_model_part_name"].GetString());
        ModelPart &slave_model_part     = mrModelPart.GetSubModelPart(mParameters["slave_sub_model_part_name"].GetString());
        NodesArrayType &r_nodes_master  = master_model_part.Nodes();
        NodesArrayType &r_nodes_slave   = slave_model_part.Nodes();

        for(SizeType dof_iterator=0;dof_iterator<mParameters["variable_names"].size();++dof_iterator)
        {
            VariableComponentType current_dof =
            KratosComponents<VariableComponentType>::Get(mParameters["variable_names"][dof_iterator].GetString());

            for(SizeType master_iterator =0;master_iterator<rNumberOfNeighbors;++master_iterator)
            {

                ModelPart::IndexType current_id = mrMainModelPart.NumberOfMasterSlaveConstraints()+1;

                const ConstraintPointer p_current_constraint = mrMainModelPart.CreateNewMasterSlaveConstraint(mParameters["constraint_set_name"].GetString(),current_id,
                r_nodes_master[rNeighborNodes[master_iterator]->Id()],current_dof,r_nodes_slave[rCurrentSlaveNode.Id()],
                current_dof,rNodalNeighborWeights[master_iterator],0.0);

                //p_current_constraint->Set(TO_ERASE);

                if(mParameters["debug_info"].GetBool()){
                    std::cout << rNeighborNodes[master_iterator]->Id()
                    << "-----" << rCurrentSlaveNode.Id() << "-----" << rNodalNeighborWeights[master_iterator]
                    << " ::: " << mParameters["variable_names"][dof_iterator].GetString() << std::endl;
                }
            } // each master node
        }  // each dof
    }
  */
  protected:

  private:

    ModelPart& mrModelPart;
    ModelPart& mrMainModelPart;
    Parameters mParameters;


}; //class

}; // namespace
#endif
