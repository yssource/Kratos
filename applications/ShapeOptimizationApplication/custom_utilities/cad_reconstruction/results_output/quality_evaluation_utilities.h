// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
//
// ==============================================================================

#ifndef QUALITY_EVALUATION_UTILITY_H
#define QUALITY_EVALUATION_UTILITY_H

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <string>
#include <cmath>
#include <iostream>
#include <numeric>

// ------------------------------------------------------------------------------
// External includes
// ------------------------------------------------------------------------------

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "includes/define.h"
#include "includes/node.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "../basic_nurbs_brep_handling/brep_gauss_point.h"
#include "../data_management/cad_projection_utility.h"
#include "../data_management/reconstruction_data_base.h"
#include "quality_evaluation_vtk_output_utilities.h"

// ==============================================================================

namespace Kratos
{
class QualityEvaluationUtility
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
    
    /// Pointer definition of QualityEvaluationUtility
    KRATOS_CLASS_POINTER_DEFINITION(QualityEvaluationUtility);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    QualityEvaluationUtility( ReconstructionDataBase& reconstruction_data_base, ReconstructionConditionContainer& condition_container, Parameters& reconstruction_parameters)
    :mrReconstructionDataBase( reconstruction_data_base ),
     mrReconstructionConditions( condition_container.GetReconstructionConditions() ),
     mrReconstructionParameters( reconstruction_parameters ) 
    {      
    }

    /// Destructor.
    virtual ~QualityEvaluationUtility()
    {
    }

    // --------------------------------------------------------------------------
    void EvaluateSurfaceReconstruction()
    {
        std::cout << "\n> Starting to evaluate surface reconstruction..." << std::endl;
        boost::timer timer;  

        // Variable to define deformation
        std::string shape_change_variable_name = mrReconstructionParameters["inpute_parameters"]["shape_change_variable_name"].GetString();
        Variable<array_1d<double,3>> shape_change_variable = KratosComponents< Variable<array_1d<double,3>> >::Get(shape_change_variable_name);      

        // Quality measures
        std::vector<double> distances_by_projection;
        std::vector<double> distances_to_projected_cad_point_after_reconstruction;
        std::vector<double> distances_to_nearest_cad_point_after_reconstruction;
        std::vector<unsigned int> is_projected_cad_point_inside;

        // Quality evaluation
        bool is_vtk_output_required = mrReconstructionParameters["output_parameters"]["quality_evaluation_parameters"]["compile_results_in_vtk_file"].GetBool();
        if(is_vtk_output_required)
        {
            ModelPart::Pointer mdpa_to_evaluate_surface_reconstruction = ModelPart::Pointer( new ModelPart("EvaluationPart",1) );
            EvaluateQualityWithVTKOutput( shape_change_variable,
                                          distances_by_projection, 
                                          distances_to_projected_cad_point_after_reconstruction,
                                          distances_to_nearest_cad_point_after_reconstruction, 
                                          is_projected_cad_point_inside,
                                          mdpa_to_evaluate_surface_reconstruction );
            WriteVTKFileWithQualityResults( mdpa_to_evaluate_surface_reconstruction ); 
        }                                         
        else
        {
            EvaluateQualityWithoutVTKOutput( shape_change_variable,
                                             distances_by_projection, 
                                             distances_to_projected_cad_point_after_reconstruction,
                                             distances_to_nearest_cad_point_after_reconstruction, 
                                             is_projected_cad_point_inside );
        }

        WriteStatisticsToConsole( distances_by_projection,
                                  distances_to_projected_cad_point_after_reconstruction,
                                  distances_to_nearest_cad_point_after_reconstruction, 
                                  is_projected_cad_point_inside );

        std::cout << "> Finished evaluating surface reconstruction in " << timer.elapsed() << " s." << std::endl;                        
    }

    // --------------------------------------------------------------------------
    void EvaluateDisplacementCoupling()
    {
        std::cout << "\n> Starting to evaluate displacement coupling..." << std::endl;
        boost::timer timer;  

        double integral_d = 0;
        double length = 0;
        DoubleVector distances;

        BREPElementVector& brep_elements_vector = mrReconstructionDataBase.GetBREPElements();
        for(auto & brep_element_i : brep_elements_vector)
        {
            if(brep_element_i.HasCouplingCondition())
            {
                BREPGaussPointVector& coupling_gauss_points = brep_element_i.GetGaussPoints();
                for(auto & gauss_point_i : coupling_gauss_points)
                {
                    // Read information from Gauss point
                    unsigned int master_patch_id = gauss_point_i.GetMasterPatchId();
                    Patch& master_patch = mrReconstructionDataBase.GetPatchFromPatchId( master_patch_id );
                    unsigned int slave_patch_id = gauss_point_i.GetSlavePatchId();
                    Patch& slave_patch = mrReconstructionDataBase.GetPatchFromPatchId( slave_patch_id );                
                    double gp_i_weight = gauss_point_i.GetWeight();
                    array_1d<double,2> location_on_master_patch = gauss_point_i.GetLocationOnMasterInParameterSpace();
                    array_1d<double,2> location_on_slave_patch = gauss_point_i.GetLocationOnSlaveInParameterSpace();
                    array_1d<double,2> tangent_on_master_patch = gauss_point_i.GetTangentOnMasterInParameterSpace();
                    array_1d<double,2> tangent_on_slave_patch = gauss_point_i.GetTangentOnSlaveInParameterSpace();
                    array_1d<int, 2> knotspans_on_master_patch = master_patch.ComputeSurfaceKnotSpans(location_on_master_patch);
                    array_1d<int, 2> knotspans_on_slave_patch = slave_patch.ComputeSurfaceKnotSpans(location_on_slave_patch);

                    // Compute Jacobian J1
                    matrix<double> g_master = master_patch.ComputeBaseVectors(knotspans_on_master_patch, location_on_master_patch);
                    Vector g1 = ZeroVector(3);
                    g1(0) = g_master(0,0);
                    g1(1) = g_master(1,0);
                    g1(2) = g_master(2,0);
                    Vector g2 = ZeroVector(3);
                    g2(0) = g_master(0,1);
                    g2(1) = g_master(1,1);
                    g2(2) = g_master(2,1);
                    double J1 = norm_2( g1* tangent_on_master_patch(0) + g2* tangent_on_master_patch(1) );

                    // Compute distance between master and slave Gauss point
                    Point<3> slave_point;
                    Point<3> master_point;
                    slave_patch.EvaluateSurfacePoint(location_on_slave_patch, slave_point);
                    master_patch.EvaluateSurfacePoint(location_on_master_patch, master_point);
                    double distance;
                    EvaluateDistanceBetweenPoints(master_point, slave_point, distance);
                    distances.push_back(distance);

                    integral_d += distance * J1 * gp_i_weight;
                    length += J1 * gp_i_weight;
                }
            }
        }
      
        // analyze the global quality of reconstruction
        double average_dist = std::accumulate(distances.begin(), distances.end(), 0.0) / distances.size();
        double max = *std::max_element(distances.begin(), distances.end());
        std::cout << "\n> Quality of reconstruction: G0 (displacement coupling)" << std::endl;
        std::cout << "\t average distance = " << std::fixed << std::scientific << std::setprecision(6) << average_dist << std::endl;
        std::cout << "\t max distance = " << std::fixed << std::setprecision(6) << max << std::endl;
        std::cout << "\t integral mean = " << std::fixed << std::setprecision(6) << integral_d/length << std::endl;
      
        std::cout << "\n> Finished evaluating displacement coupling in " << timer.elapsed() << " s." << std::endl;  
    }

    // --------------------------------------------------------------------------
    void EvaluateRotationCoupling()
    {
        std::cout << "\n> Starting to evaluate rotation coupling..." << std::endl;
        boost::timer timer;  

        double integral_angle = 0;
        double length = 0;
        DoubleVector angles;

        BREPElementVector& brep_elements_vector = mrReconstructionDataBase.GetBREPElements();
        for(auto & brep_element_i : brep_elements_vector)
        {
            if(brep_element_i.HasCouplingCondition())
            {
                BREPGaussPointVector& coupling_gauss_points = brep_element_i.GetGaussPoints();
                for(auto & gauss_point_i : coupling_gauss_points)
                {
                    // Read information from Gauss point
                    unsigned int master_patch_id = gauss_point_i.GetMasterPatchId();
                    Patch& master_patch = mrReconstructionDataBase.GetPatchFromPatchId( master_patch_id );
                    unsigned int slave_patch_id = gauss_point_i.GetSlavePatchId();
                    Patch& slave_patch = mrReconstructionDataBase.GetPatchFromPatchId( slave_patch_id );                
                    double gp_i_weight = gauss_point_i.GetWeight();
                    array_1d<double,2> location_on_master_patch = gauss_point_i.GetLocationOnMasterInParameterSpace();
                    array_1d<double,2> location_on_slave_patch = gauss_point_i.GetLocationOnSlaveInParameterSpace();
                    array_1d<double,2> tangent_on_master_patch = gauss_point_i.GetTangentOnMasterInParameterSpace();
                    array_1d<double,2> tangent_on_slave_patch = gauss_point_i.GetTangentOnSlaveInParameterSpace();
                    array_1d<int, 2> knotspans_on_master_patch = master_patch.ComputeSurfaceKnotSpans(location_on_master_patch);
                    array_1d<int, 2> knotspans_on_slave_patch = slave_patch.ComputeSurfaceKnotSpans(location_on_slave_patch);

                    // Compute Jacobian J1
                    matrix<double> g_master = master_patch.ComputeBaseVectors(knotspans_on_master_patch, location_on_master_patch);
                    Vector g1 = ZeroVector(3);
                    g1(0) = g_master(0,0);
                    g1(1) = g_master(1,0);
                    g1(2) = g_master(2,0);
                    Vector g2 = ZeroVector(3);
                    g2(0) = g_master(0,1);
                    g2(1) = g_master(1,1);
                    g2(2) = g_master(2,1);
                    double J1 = norm_2( g1* tangent_on_master_patch(0) + g2* tangent_on_master_patch(1) );

                    // Compute normals on master and slave
                    Vector g3_m = ZeroVector(3);
                    g3_m(0) = g_master(0,2);
                    g3_m(1) = g_master(1,2);
                    g3_m(2) = g_master(2,2);
                    matrix<double> g_slave = slave_patch.ComputeBaseVectors(knotspans_on_slave_patch, location_on_slave_patch);
                    Vector g3_s = ZeroVector(3);
                    g3_s(0) = g_slave(0,2);
                    g3_s(1) = g_slave(1,2);
                    g3_s(2) = g_slave(2,2);

                    // Compute angle between normals on master and slave Gauss point
                    double cosine = inner_prod(g3_m, g3_s);
                    double angle_rad = acos(std::abs(cosine)); // abs() => ( 0 < angle_rad < PI/2 )
                    double angle_deg = angle_rad * 180 / 3.1415926535;
                    angles.push_back(angle_deg);

                    // Sum contribution
                    integral_angle += angle_deg * J1 * gp_i_weight;
                    length += J1 * gp_i_weight;
                }
            }
        }
      
        // analyze the global quality of reconstruction
        double average_angle = std::accumulate(angles.begin(), angles.end(), 0.0) / angles.size();
        double max = *std::max_element(angles.begin(), angles.end());
        std::cout << "\n> Quality of reconstruction: G1 (rotation coupling)" << std::endl;
        std::cout << "\t average angle = " << std::fixed << std::scientific << std::setprecision(6) << average_angle << " deg" << std::endl;
        std::cout << "\t max angle = " << std::fixed << std::setprecision(6) << max << " deg" << std::endl;
        std::cout << "\t integral mean = " << std::fixed << std::setprecision(6) << integral_angle/length << " deg" << std::endl;

        std::cout << "\n> Finished evaluating rotation coupling in " << timer.elapsed() << " s." << std::endl;        
    }

    /// Turn back information as a string.
    virtual std::string Info() const
    {
		return "QualityEvaluationUtility";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream &rOStream) const
    {
		rOStream << "QualityEvaluationUtility";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream &rOStream) const
    {
    }

private:

    ///@}
    ///@name Private member variables
    ///@{

    // Initialized by constructor 
    ReconstructionDataBase& mrReconstructionDataBase;
    std::vector<ReconstructionCondition::Pointer>& mrReconstructionConditions;
    Parameters& mrReconstructionParameters;
    
    // Additional variables    
    int points_counter=0;
    int outside_points_counter=0;  

    ///@}
    ///@name Private member functions
    ///@{    
    void EvaluateQualityWithVTKOutput( Variable<array_1d<double,3>> shape_change_variable,
                                       std::vector<double>& distances_by_projection, 
                                       std::vector<double>& distances_to_projected_cad_point_after_reconstruction, 
                                       std::vector<double>& distances_to_nearest_cad_point_after_reconstruction, 
                                       std::vector<unsigned int>& is_projected_cad_point_inside,
                                       ModelPart::Pointer mdpa_to_evaluate_surface_reconstruction )
    {    
        mdpa_to_evaluate_surface_reconstruction->AddNodalSolutionStepVariable(DISTANCE_BY_PROJECTION);
        mdpa_to_evaluate_surface_reconstruction->AddNodalSolutionStepVariable(PATCH_ID);
        mdpa_to_evaluate_surface_reconstruction->AddNodalSolutionStepVariable(IS_PROJECTED_CAD_POINT_INSIDE_VISIBLE_PATCH_REGION);
        mdpa_to_evaluate_surface_reconstruction->AddNodalSolutionStepVariable(SHAPE_CHANGE_ABSOLUTE);
        mdpa_to_evaluate_surface_reconstruction->AddNodalSolutionStepVariable(DISTANCE_TO_PROJECTED_CAD_POINT_AFTER_RECONSTRUCTION);
        mdpa_to_evaluate_surface_reconstruction->AddNodalSolutionStepVariable(DISTANCE_TO_NEAREST_CAD_POINT_AFTER_RECONSTRUCTION);   

        CADProjectionUtility FE2CADProjector( mrReconstructionDataBase.GetPatchVector(), mrReconstructionParameters["output_parameters"]["quality_evaluation_parameters"]["projection_parameters"] );
        FE2CADProjector.Initialize();

        for(auto & condition_i : mrReconstructionConditions)
        {
            int condition_index = &condition_i-&mrReconstructionConditions[0];

            // Get condition data
            array_1d<double,3> fe_coordinates_undeformed;
            array_1d<double,3> fe_coordinates_deformed;
            array_1d<double,3> cad_coordinates_undeformed;
            array_1d<double,3> cad_coordinates_deformed;
            array_1d<double,3> nearest_cad_coordinates_deformed; 
            condition_i->DetermineFECoordinatesInUndeformedConfiguration(fe_coordinates_undeformed);
            condition_i->DetermineFECoordinatesInDeformedConfiguration(shape_change_variable, fe_coordinates_deformed);
            condition_i->DetermineCADCoordinatesInUndeformedConfiguration(cad_coordinates_undeformed);
            condition_i->DetermineCADCoordinatesInDeformedConfiguration(cad_coordinates_deformed);
            
            Node<3> fe_point(1, fe_coordinates_deformed);
            FE2CADProjector.DetermineNearestCADPointInGeometrySpace( fe_point, nearest_cad_coordinates_deformed );
            
            bool projected_cad_point_is_inside = condition_i->IsProjectedCADPointInsideVisiblePatchRegion();

            // compute distances
            double distance_by_projection =  EvaluateDistanceBetweenCoordinates(cad_coordinates_undeformed, fe_coordinates_undeformed);
            double distance_to_projected_cad_point_after_reconstruction =  EvaluateDistanceBetweenCoordinates(cad_coordinates_deformed, fe_coordinates_deformed);
            double distance_to_nearest_cad_point_after_reconstruction =  EvaluateDistanceBetweenCoordinates(nearest_cad_coordinates_deformed, fe_coordinates_deformed);

            // push_back results
            if(projected_cad_point_is_inside)
            {
                distances_by_projection.push_back(distance_by_projection);
                distances_to_projected_cad_point_after_reconstruction.push_back(distance_to_projected_cad_point_after_reconstruction);
                is_projected_cad_point_inside.push_back(0);
            }
            else
            {
                is_projected_cad_point_inside.push_back(1);
            }
            distances_to_nearest_cad_point_after_reconstruction.push_back(distance_to_nearest_cad_point_after_reconstruction);

            // create node at FE coordinates in undeformed configuration
            NodeType::Pointer new_node = Node <3>::Pointer(new Node<3>(condition_index, fe_coordinates_undeformed));
            new_node->SetSolutionStepVariablesList(&(mdpa_to_evaluate_surface_reconstruction->GetNodalSolutionStepVariablesList()));
            mdpa_to_evaluate_surface_reconstruction->AddNode(new_node);

            // Add nodal results
            new_node->FastGetSolutionStepValue(DISTANCE_BY_PROJECTION) = cad_coordinates_undeformed - fe_coordinates_undeformed;
            new_node->FastGetSolutionStepValue(PATCH_ID) = condition_i->GetAffectedPatch().GetId();
            new_node->FastGetSolutionStepValue(IS_PROJECTED_CAD_POINT_INSIDE_VISIBLE_PATCH_REGION) = projected_cad_point_is_inside;
            new_node->FastGetSolutionStepValue(SHAPE_CHANGE_ABSOLUTE) = fe_coordinates_deformed - fe_coordinates_undeformed;
            new_node->FastGetSolutionStepValue(DISTANCE_TO_PROJECTED_CAD_POINT_AFTER_RECONSTRUCTION) = cad_coordinates_deformed - fe_coordinates_deformed;
            new_node->FastGetSolutionStepValue(DISTANCE_TO_NEAREST_CAD_POINT_AFTER_RECONSTRUCTION) = nearest_cad_coordinates_deformed - fe_coordinates_deformed;                 
        }        
    }

    // --------------------------------------------------------------------------
    void EvaluateQualityWithoutVTKOutput( Variable<array_1d<double,3>>& shape_change_variable,
                                          std::vector<double>& distances_by_projection, 
                                          std::vector<double>& distances_to_projected_cad_point_after_reconstruction, 
                                          std::vector<double>& distances_to_nearest_cad_point_after_reconstruction, 
                                          std::vector<unsigned int>& is_projected_cad_point_inside )
    {         
        const Parameters& projection_parameters =  mrReconstructionParameters["output_parameters"]["quality_evaluation_parameters"]["projection_parameters"];
        CADProjectionUtility FE2CADProjector( mrReconstructionDataBase.GetPatchVector(), projection_parameters );
        FE2CADProjector.Initialize();

        for(auto & condition_i : mrReconstructionConditions)
        {
            int condition_index = &condition_i-&mrReconstructionConditions[0];

            // Get condition data
            array_1d<double,3> fe_coordinates_undeformed;
            array_1d<double,3> fe_coordinates_deformed;
            array_1d<double,3> cad_coordinates_undeformed;
            array_1d<double,3> cad_coordinates_deformed;
            array_1d<double,3> nearest_cad_coordinates_deformed; 
            condition_i->DetermineFECoordinatesInUndeformedConfiguration(fe_coordinates_undeformed);
            condition_i->DetermineFECoordinatesInDeformedConfiguration(shape_change_variable, fe_coordinates_deformed);
            condition_i->DetermineCADCoordinatesInUndeformedConfiguration(cad_coordinates_undeformed);
            condition_i->DetermineCADCoordinatesInDeformedConfiguration(cad_coordinates_deformed);

            Node<3> fe_point(1, fe_coordinates_deformed);
            FE2CADProjector.DetermineNearestCADPointInGeometrySpace( fe_point, nearest_cad_coordinates_deformed );
            
            bool projected_cad_point_is_inside = condition_i->IsProjectedCADPointInsideVisiblePatchRegion();

            // compute distances
            double distance_by_projection =  EvaluateDistanceBetweenCoordinates(cad_coordinates_undeformed, fe_coordinates_undeformed);
            double distance_to_projected_cad_point_after_reconstruction =  EvaluateDistanceBetweenCoordinates(cad_coordinates_deformed, fe_coordinates_deformed);
            double distance_to_nearest_cad_point_after_reconstruction =  EvaluateDistanceBetweenCoordinates(nearest_cad_coordinates_deformed, fe_coordinates_deformed);

            // push_back results
            if(projected_cad_point_is_inside)
            {
                distances_by_projection.push_back(distance_by_projection);
                distances_to_projected_cad_point_after_reconstruction.push_back(distance_to_projected_cad_point_after_reconstruction);
                is_projected_cad_point_inside.push_back(0);
            }
            else
            {
                is_projected_cad_point_inside.push_back(1);
            }
            distances_to_nearest_cad_point_after_reconstruction.push_back(distance_to_nearest_cad_point_after_reconstruction);              
        }        
    }

    // --------------------------------------------------------------------------
    void WriteVTKFileWithQualityResults( ModelPart::Pointer mdpa_to_evaluate_surface_reconstruction )
    {
        Parameters vtk_output_parameters(R"({
            "nodal_results": [ "DISTANCE_BY_PROJECTION",
                               "PATCH_ID", 
                               "IS_PROJECTED_CAD_POINT_INSIDE_VISIBLE_PATCH_REGION",
                               "SHAPE_CHANGE_ABSOLUTE",
                               "DISTANCE_TO_PROJECTED_CAD_POINT_AFTER_RECONSTRUCTION", 
                               "DISTANCE_TO_NEAREST_CAD_POINT_AFTER_RECONSTRUCTION" ],
            "vtk_type"     : "Ascii"
        })");
        vtk_output_parameters.AddValue("output_folder", mrReconstructionParameters["output_parameters"]["output_folder"]);
        vtk_output_parameters.AddValue("output_filename", mrReconstructionParameters["output_parameters"]["quality_evaluation_parameters"]["vtk_filename"]);
        
        QualityEvaluationVTKOutputUtilities vtk_io( *mdpa_to_evaluate_surface_reconstruction, vtk_output_parameters );
        vtk_io.PrintOutput();        
    }

    // --------------------------------------------------------------------------
    void WriteStatisticsToConsole( std::vector<double>& distances_by_projection, 
                                   std::vector<double>& distances_to_projected_cad_point_after_reconstruction,
                                   std::vector<double>& distances_to_nearest_cad_point_after_reconstruction, 
                                   std::vector<unsigned int>& is_projected_cad_point_inside )
    {
        double average_dist = std::accumulate( distances_by_projection.begin(), distances_by_projection.end(), 0.0 ) / distances_by_projection.size();
        double max = *std::max_element( distances_by_projection.begin(), distances_by_projection.end());
        int number_of_qk_outside = std::accumulate( is_projected_cad_point_inside.begin(), is_projected_cad_point_inside.end(), 0 );
        double percentage_of_qk_outside = (100.0 * number_of_qk_outside) / is_projected_cad_point_inside.size();
        std::cout << "\n> Quality of projection:" << std::endl;
        std::cout << "\t average distance = " << std::fixed << std::scientific << std::setprecision(6) << average_dist << std::endl;
        std::cout << "\t max distance = " << std::fixed << std::setprecision(6) << max << std::endl;
        std::cout << "\t " << number_of_qk_outside << " out of " << is_projected_cad_point_inside.size() << " points are outside of trimming boundaries (" << std::fixed << std::setprecision(2) << percentage_of_qk_outside << "%)" << std::endl;
        
        average_dist = std::accumulate( distances_to_projected_cad_point_after_reconstruction.begin(), distances_to_projected_cad_point_after_reconstruction.end(), 0.0 ) / distances_to_projected_cad_point_after_reconstruction.size();
        max = *std::max_element( distances_to_projected_cad_point_after_reconstruction.begin(), distances_to_projected_cad_point_after_reconstruction.end());
        std::cout << "\n> Quality of reconstruction: | P - Q_k |" << std::endl;
        std::cout << "\t average distance = " << std::fixed << std::scientific << std::setprecision(6) << average_dist << std::endl;
        std::cout << "\t max distance = " << std::fixed << std::setprecision(6) << max << std::endl;
        
        average_dist = std::accumulate( distances_to_nearest_cad_point_after_reconstruction.begin(), distances_to_nearest_cad_point_after_reconstruction.end(), 0.0 ) / distances_to_nearest_cad_point_after_reconstruction.size();
        max = *std::max_element( distances_to_nearest_cad_point_after_reconstruction.begin(), distances_to_nearest_cad_point_after_reconstruction.end());
        std::cout << "\n> Quality of reconstruction: | P - Q |" << std::endl;
        std::cout << "\t average distance = " << std::fixed << std::scientific << std::setprecision(6) << average_dist << std::endl;
        std::cout << "\t max distance = " << std::fixed << std::setprecision(6) << max << std::endl; 
    }    

    // --------------------------------------------------------------------------
    void EvaluateDistanceBetweenPoints(Point<3>& first, Point<3>& second, double& rdistance)
    {
      Vector distance_vector = ZeroVector(3);
      distance_vector(0) = first.X() - second.X();
      distance_vector(1) = first.Y() - second.Y();
      distance_vector(2) = first.Z() - second.Z();
      rdistance = norm_2(distance_vector);
    }

    // --------------------------------------------------------------------------
    double EvaluateDistanceBetweenCoordinates(array_1d<double,3>& first, array_1d<double,3>& second)
    {
      array_1d<double,3> d  = first - second;
      return std::sqrt( d[0]*d[0] + d[1]*d[1] + d[2]*d[2] );
    }

    /// Assignment operator.
    //      QualityEvaluationUtility& operator=(QualityEvaluationUtility const& rOther);

    /// Copy constructor.
    //      QualityEvaluationUtility(QualityEvaluationUtility const& rOther);

}; // Class QualityEvaluationUtility
} // namespace Kratos.

#endif // QUALITY_EVALUATION_UTILITY_H