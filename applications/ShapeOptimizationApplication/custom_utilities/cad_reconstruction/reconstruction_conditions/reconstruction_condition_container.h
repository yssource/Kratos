// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
//
// ==============================================================================

#ifndef RECONSTRUCTION_CONDITION_CONTAINER_H
#define RECONSTRUCTION_CONDITION_CONTAINER_H

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <iostream>
#include <string>

// ------------------------------------------------------------------------------
// External includes
// ------------------------------------------------------------------------------

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "includes/define.h"
#include "processes/process.h"
#include "includes/node.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "includes/kratos_flags.h"
#include "../basic_nurbs_brep_handling/brep_gauss_point.h"
#include "../data_management/cad_projection_utility.h"
#include "../reconstruction_conditions/reconstruction_condition_displacement_mapping.h"
#include "../reconstruction_conditions/reconstruction_condition_distance_minimization.h"
#include "../reconstruction_conditions/reconstruction_constraint_displacement_coupling.h"
#include "../reconstruction_conditions/reconstruction_constraint_rotation_coupling.h"
#include "../reconstruction_conditions/reconstruction_constraint_tangent_enforcement.h"
#include "../reconstruction_conditions/regularization_condition_min_control_point_displacement.h"
#include "../reconstruction_conditions/regularization_condition_min_control_point_distance_to_surface.h"

// ==============================================================================

namespace Kratos
{

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

/// Short class definition.
/** Detail class definition.

*/

class ReconstructionConditionContainer
{
public:
    ///@name Type Definitions
    ///@{

    typedef Node<3> NodeType;
    typedef std::vector<NodeType::Pointer> NodeVector;   
    typedef std::vector<Patch> PatchVector;
    typedef std::vector<BREPElement> BREPElementVector;
    typedef std::vector<BREPGaussPoint> BREPGaussPointVector;    
    typedef Element::GeometryType::IntegrationMethod IntegrationMethodType;    
        
    /// Pointer definition of ReconstructionConditionContainer
    KRATOS_CLASS_POINTER_DEFINITION(ReconstructionConditionContainer);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ReconstructionConditionContainer( ReconstructionDataBase& reconstruction_data_base, Parameters& reconstruction_parameters )
    : mrReconstructionDataBase( reconstruction_data_base ),
      mrReconstructionParameters( reconstruction_parameters )    
    {
    }

    /// Destructor.
    virtual ~ReconstructionConditionContainer()
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    // ==============================================================================
    void CreateDisplacementMappingConditions()
    {
        std::cout << "\n> Starting to create displacement mapping conditions...";

        PatchVector& patch_vector = mrReconstructionDataBase.GetPatchVector();
        ModelPart& fe_model_part = mrReconstructionDataBase.GetFEModelPart();
        
        IntegrationMethodType fem_integration_method = GeometryData::GI_GAUSS_5;
        int integration_degree = mrReconstructionParameters["solution_parameters"]["strategy_specifc_parameters"]["fem_gauss_integration_degree"].GetInt();
        switch(integration_degree)
        {
            case 1 : fem_integration_method = GeometryData::GI_GAUSS_1; break;
            case 2 : fem_integration_method = GeometryData::GI_GAUSS_2; break;
            case 3 : fem_integration_method = GeometryData::GI_GAUSS_3; break;
            case 4 : fem_integration_method = GeometryData::GI_GAUSS_4; break;
            case 5 : fem_integration_method = GeometryData::GI_GAUSS_5; break;
        }
        
        CADProjectionUtility FE2CADProjector( patch_vector, mrReconstructionParameters["solution_parameters"]["projection_parameters"] );
        FE2CADProjector.Initialize();
        
        std::cout << "> Starting to create a condition for every integration point..." << std::endl;
        boost::timer timer;   
        
        ReconstructionCondition::Pointer NewCondition;
        
        for (auto & elem_i : fe_model_part.Elements())
        {
            Element::GeometryType& geom_i = elem_i.GetGeometry();
            const Element::GeometryType::IntegrationPointsArrayType& integration_points = geom_i.IntegrationPoints(fem_integration_method);

            for (auto & integration_point_i : integration_points)
            {
                int integration_point_number = &integration_point_i - &integration_points[0];
                NodeType::CoordinatesArrayType ip_coordinates = geom_i.GlobalCoordinates(ip_coordinates, integration_point_i.Coordinates());
                NodeType node_of_interest(1, ip_coordinates);

                array_1d<double,2> parameter_values_of_nearest_point;
                int patch_index_of_nearest_point = -1;

                FE2CADProjector.DetermineNearestCADPoint( node_of_interest, parameter_values_of_nearest_point, patch_index_of_nearest_point );

                NewCondition = ReconstructionCondition::Pointer( new DisplacementMappingCondition( geom_i,
                                                                                                   fem_integration_method,
                                                                                                   integration_point_number,
                                                                                                   patch_vector[patch_index_of_nearest_point],
                                                                                                   parameter_values_of_nearest_point ));
                mListOfReconstructionConditions.push_back( NewCondition );
            }
        }
        std::cout << "> Time needed for looping over integration points: " << timer.elapsed() << " s" << std::endl;            
        std::cout << "> Finished creating displacement mapping conditions: " << std::endl;      
    }

    // --------------------------------------------------------------------------
    void CreateDistanceMinimizationConditions()
    {
        std::cout << "\n> Starting to create distance minimization conditions...";

        PatchVector& patch_vector = mrReconstructionDataBase.GetPatchVector();
        ModelPart& fe_model_part = mrReconstructionDataBase.GetFEModelPart();
                
        CADProjectionUtility FE2CADProjector( patch_vector, mrReconstructionParameters["solution_parameters"]["projection_parameters"] );
        FE2CADProjector.Initialize();
        
        std::cout << "> Starting to create a condition for every fem point..." << std::endl;
        boost::timer timer;   
        
        ReconstructionCondition::Pointer NewCondition;
        
        for (auto & node_i : fe_model_part.Nodes())
        {
            array_1d<double,2> parameter_values_of_nearest_point;
            int patch_index_of_nearest_point = -1;

            FE2CADProjector.DetermineNearestCADPoint( node_i, parameter_values_of_nearest_point, patch_index_of_nearest_point );

            NewCondition = ReconstructionCondition::Pointer( new DistanceMinimizationCondition( node_i,
                                                                                                patch_vector[patch_index_of_nearest_point],
                                                                                                parameter_values_of_nearest_point ));
            mListOfReconstructionConditions.push_back( NewCondition );
        }
        std::cout << "> Time needed for looping over fem points: " << timer.elapsed() << " s" << std::endl;            
        std::cout << "> Finished creating distance minimization mapping conditions: " << std::endl;      
    }    

    // --------------------------------------------------------------------------
    void CreateDisplacementCouplingConstraintsOnAllCouplingPoints()
    {
        std::cout << "\n> Starting to create displacement coupling constraints..." << std::endl;
        boost::timer timer;   
        
        ReconstructionConstraint::Pointer NewConstraint;
        double penalty_factor = mrReconstructionParameters["solution_parameters"]["constraints"]["penalty_factor_for_displacement_coupling"].GetDouble();        

        BREPElementVector& brep_elements_vector = mrReconstructionDataBase.GetBREPElements();
        for(auto & brep_element_i : brep_elements_vector)
        {
            if(brep_element_i.HasCouplingCondition())
            {
                BREPGaussPointVector& coupling_gauss_points = brep_element_i.GetGaussPoints();
                for(auto & gauss_point_i : coupling_gauss_points)
                {
                        unsigned int master_patch_id = gauss_point_i.GetMasterPatchId();
                        Patch& master_patch = mrReconstructionDataBase.GetPatchFromPatchId( master_patch_id );
            
                        unsigned int slave_patch_id = gauss_point_i.GetSlavePatchId();
                        Patch& slave_patch = mrReconstructionDataBase.GetPatchFromPatchId( slave_patch_id );                
                        
                        NewConstraint = ReconstructionConstraint::Pointer( new DisplacementCouplingConstraint( gauss_point_i,
                                                                                                               master_patch,
                                                                                                               slave_patch,
                                                                                                               penalty_factor ));
                        mListOfReconstructionConstraints.push_back( NewConstraint );
                }
            }
        }
        std::cout << "> Time needed for creating displacement constraints: " << timer.elapsed() << " s" << std::endl;              
    }

    // --------------------------------------------------------------------------
    void CreateRotationCouplingConstraintsOnAllCouplingPoints( Parameters& rListOfEdgeIdsWithEnforcedTangent )
    {
        std::cout << "\n> Starting to create rotation coupling constraints..." << std::endl;
        boost::timer timer;

        std::vector<int> list_of_ids_of_edges_with_enforced_tangent;
        for(int index=0; index<rListOfEdgeIdsWithEnforcedTangent.size(); index++)
            list_of_ids_of_edges_with_enforced_tangent.push_back(rListOfEdgeIdsWithEnforcedTangent.GetArrayItem(index).GetInt());  
                   
        double penalty_factor = mrReconstructionParameters["solution_parameters"]["constraints"]["penalty_factor_for_rotation_coupling"].GetDouble();
        ReconstructionConstraint::Pointer NewConstraint;

        BREPElementVector& brep_elements_vector = mrReconstructionDataBase.GetBREPElements();
        for(auto & brep_element_i : brep_elements_vector)
        {
            if(brep_element_i.HasCouplingCondition())
            {
                unsigned int edge_id_corresponding_to_brep_element_i = brep_element_i.GetEdgeId();
                bool has_edge_enforced_tangent_constraint = (std::find(list_of_ids_of_edges_with_enforced_tangent.begin(), list_of_ids_of_edges_with_enforced_tangent.end(), edge_id_corresponding_to_brep_element_i) != list_of_ids_of_edges_with_enforced_tangent.end());

                // Skip current rotation coupling if for this edge tangent continuity shall be enforced
                if(has_edge_enforced_tangent_constraint)
                {
                    std::cout << "> Skipping brep_element on the following edge in the rotation coupling since tangent continuity shall be enforced here: " << edge_id_corresponding_to_brep_element_i << std::endl;
                    continue;
                }

                BREPGaussPointVector& coupling_gauss_points = brep_element_i.GetGaussPoints();
                for(auto & gauss_point_i : coupling_gauss_points)
                {
                        unsigned int master_patch_id = gauss_point_i.GetMasterPatchId();
                        Patch& master_patch = mrReconstructionDataBase.GetPatchFromPatchId( master_patch_id );
            
                        unsigned int slave_patch_id = gauss_point_i.GetSlavePatchId();
                        Patch& slave_patch = mrReconstructionDataBase.GetPatchFromPatchId( slave_patch_id );                
                        
                        NewConstraint = ReconstructionConstraint::Pointer( new RotationCouplingConstraint( gauss_point_i,
                                                                                                           master_patch,
                                                                                                           slave_patch,
                                                                                                           penalty_factor ));
                        mListOfReconstructionConstraints.push_back( NewConstraint );
                }
            }
        }
        std::cout << "> Time needed for creating rotation constraints: " << timer.elapsed() << " s" << std::endl;              
    }    

    // --------------------------------------------------------------------------
    void CreateDirichletConditions( Parameters& rListOfEdgeIds )
    {
    }

    // --------------------------------------------------------------------------
    void CreateTangentContinuityConditions( Parameters& rListOfEdgeIds )
    {
        std::cout << "\n> Starting to create tangenent enforcement constraints..." << std::endl;
        boost::timer timer;   
        
        double penalty_factor = mrReconstructionParameters["solution_parameters"]["constraints"]["penalty_factor_for_tangent_continuity_constraints"].GetDouble();
        ReconstructionConstraint::Pointer NewConstraint;

        std::vector<int> list_of_edge_ids;
        for(int index=0; index<rListOfEdgeIds.size(); index++)
            list_of_edge_ids.push_back(rListOfEdgeIds.GetArrayItem(index).GetInt());

        BREPElementVector& brep_elements_vector = mrReconstructionDataBase.GetBREPElements();
        for(auto & brep_element_i : brep_elements_vector)
        {
            unsigned int edge_id_corresponding_to_brep_element_i = brep_element_i.GetEdgeId();

            if( std::find(list_of_edge_ids.begin(), list_of_edge_ids.end(), edge_id_corresponding_to_brep_element_i) != list_of_edge_ids.end())
            {
                BREPGaussPointVector& coupling_gauss_points = brep_element_i.GetGaussPoints();
                for(auto & gauss_point_i : coupling_gauss_points)
                {
                        unsigned int master_patch_id = gauss_point_i.GetMasterPatchId();
                        Patch& master_patch = mrReconstructionDataBase.GetPatchFromPatchId( master_patch_id );
            
                        unsigned int slave_patch_id = gauss_point_i.GetSlavePatchId();
                        Patch& slave_patch = mrReconstructionDataBase.GetPatchFromPatchId( slave_patch_id );                
                        
                        NewConstraint = ReconstructionConstraint::Pointer( new TangentEnforcementConstraint( gauss_point_i,
                                                                                                             master_patch,
                                                                                                             slave_patch,
                                                                                                             penalty_factor ));
                        mListOfReconstructionConstraints.push_back( NewConstraint );
                }
            }
        }
        std::cout << "> Time needed for creating tangenent enforcement constraints: " << timer.elapsed() << " s" << std::endl;             
    }    
    
    // --------------------------------------------------------------------------
    void CreateMinimalControlPointDisplacementCondition()
    {
        std::cout << "\n> Starting to create beta-regularization condition (minimal control point displacement)..." << std::endl;
        
        RegularizationCondition::Pointer NewCondition = RegularizationCondition::Pointer( new MinimalControlPointDisplacementCondition( mrReconstructionDataBase, mrReconstructionParameters ) );
        mListOfRegularizationConditions.push_back( NewCondition );
        std::cout << "> Finished creating beta-regularization condition (minimal control point displacement)." << std::endl;                      
    }     

    // --------------------------------------------------------------------------
    void CreateMinimalControlPointDistanceToSurfaceCondition()
    {
        std::cout << "\n> Starting to create alpha-regularization condition (minimal control point distance to surface)..." << std::endl;
        
        RegularizationCondition::Pointer NewCondition = RegularizationCondition::Pointer( new MinimalControlPointDistanceToSurfaceCondition( mrReconstructionDataBase, mrReconstructionParameters ) );
        mListOfRegularizationConditions.push_back( NewCondition );
        std::cout << "> Finished creating alpha-regularization condition (minimal control point distance to surface)." << std::endl;                  
    }     


    // --------------------------------------------------------------------------
    std::vector<ReconstructionCondition::Pointer>& GetReconstructionConditions()
    {
        return mListOfReconstructionConditions;
    }

    // --------------------------------------------------------------------------
    std::vector<ReconstructionConstraint::Pointer>& GetReconstructionConstraints()
    {
        return mListOfReconstructionConstraints;
    }    

    // --------------------------------------------------------------------------
    std::vector<RegularizationCondition::Pointer>& GetRegularizationConditions()
    {
        return mListOfRegularizationConditions;
    }      

    // ==============================================================================

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
    virtual std::string Info() const
    {
        return "ReconstructionConditionContainer";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream &rOStream) const
    {
        rOStream << "ReconstructionConditionContainer";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream &rOStream) const
    {
    }

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

    // ==============================================================================
    // Initialized by class constructor
    // ==============================================================================
    ReconstructionDataBase& mrReconstructionDataBase;
    Parameters& mrReconstructionParameters;
    
    // ==============================================================================
    // Additional variables
    // ==============================================================================
    std::vector<ReconstructionCondition::Pointer> mListOfReconstructionConditions;
    std::vector<ReconstructionConstraint::Pointer> mListOfReconstructionConstraints;
    std::vector<RegularizationCondition::Pointer> mListOfRegularizationConditions;    
    

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

    /// Assignment operator.
    //      ReconstructionConditionContainer& operator=(ReconstructionConditionContainer const& rOther);

    /// Copy constructor.
    //      ReconstructionConditionContainer(ReconstructionConditionContainer const& rOther);

    ///@}

}; // Class ReconstructionConditionContainer

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.

#endif // RECONSTRUCTION_CONDITION_CONTAINER_H