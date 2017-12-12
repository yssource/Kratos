// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
//
// ==============================================================================

#ifndef RECONSTRUCTION_OUTPUT_UTILITIES_H
#define RECONSTRUCTION_OUTPUT_UTILITIES_H

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
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
#include "../data_management/reconstruction_data_base.h"

// ==============================================================================

namespace Kratos
{
class ReconstructionOutputUtilities
{
public:
    ///@name Type Definitions
    ///@{

    typedef Node<3> NodeType;
    typedef Element::GeometryType::IntegrationMethod IntegrationMethodType;
    
    /// Pointer definition of ReconstructionOutputUtilities
    KRATOS_CLASS_POINTER_DEFINITION(ReconstructionOutputUtilities);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ReconstructionOutputUtilities( ReconstructionDataBase& reconstruction_data_base, Parameters& parameters )
    : mrReconstructionDataBase( reconstruction_data_base ),
      mOutputFolder( parameters["output_parameters"]["output_folder"].GetString() ),
      mOriginalGeorhinoFile( parameters["output_parameters"]["original_georhino_filename"].GetString() ),
      mRhinoResultsFile( parameters["output_parameters"]["rhino_results_filename"].GetString() )
    {      
    }

    /// Destructor.
    virtual ~ReconstructionOutputUtilities()
    {
    }

    // --------------------------------------------------------------------------
    void OutputCADSurfacePoints( std::string output_filename, Parameters parameters )
    {
	    std::cout << "\n> Start writing surface points of given CAD geometry to file..." << std::endl;
        std::ofstream  output_file( mOutputFolder + "/" + output_filename );

        int u_resolution = parameters["output_parameters"]["parameter_resolution_for_output_of_surface_points"][0].GetInt();
        int v_resolution = parameters["output_parameters"]["parameter_resolution_for_output_of_surface_points"][1].GetInt();      
	
        // Set a max value for the coordinte output to avoid clutter through points that moving uncontrolled
        double max_coordinate = 10000;
        std::cout << "> Max value for control point coordinates set to " << max_coordinate << "." << std::endl;

        for(auto & patch_i : mrReconstructionDataBase.GetPatchVector()) 
        {
            std::vector<double>& knot_vec_u_i = patch_i.GetSurfaceKnotVectorU();
            std::vector<double>& knot_vec_v_i = patch_i.GetSurfaceKnotVectorV();

            double u_min = knot_vec_u_i[0];
            double u_max = knot_vec_u_i[knot_vec_u_i.size()-1];
            double v_min = knot_vec_v_i[0];
            double v_max = knot_vec_v_i[knot_vec_v_i.size()-1];
            double delta_u = (u_max-u_min) / u_resolution;
            double delta_v = (v_max-v_min) / v_resolution;

            // Loop over all u & v according to specified resolution
            array_1d<double,2> point_in_parameter_space;
            for(
                int i=0; i<=u_resolution; i++)
            {
                point_in_parameter_space[0] = u_min + i*delta_u;

                for(
                    int j=0; j<=v_resolution; j++)
                {
                    point_in_parameter_space[1] = v_min + j*delta_v;

                    bool point_is_inside = patch_i.IsPointInside(point_in_parameter_space);
                    if(point_is_inside)
                    {
                        Point<3> cad_point;
                        patch_i.EvaluateSurfacePoint( point_in_parameter_space, cad_point );

                        if(std::abs(cad_point.X())>max_coordinate)
                            cad_point.X() = MathUtils<int>::Sign(cad_point.X()) * max_coordinate;
                        if(std::abs(cad_point.Y())>max_coordinate)
                            cad_point.Y() = MathUtils<int>::Sign(cad_point.Y()) * max_coordinate;
                        if(std::abs(cad_point.Z())>max_coordinate)
                            cad_point.Z() = MathUtils<int>::Sign(cad_point.Z()) * max_coordinate;														

                        output_file << cad_point.X() << " " << cad_point.Y() << " " << cad_point.Z() << std::endl;
                    }
                }
            }
        }
        output_file.close();
        std::cout << "> Finished writing surface points of given CAD geometry to file." << std::endl;
    }

    // --------------------------------------------------------------------------
    void OutputGaussPointsOfFEMesh( std::string output_filename, int integration_degree )
    {
	    std::cout << "\n> Start writing Gauss points of given FE-mesh..." << std::endl;
        std::ofstream output_file( mOutputFolder + "/" + output_filename );
  
        ModelPart& r_fe_model_part = mrReconstructionDataBase.GetFEModelPart();

        IntegrationMethodType fem_integration_method = GeometryData::GI_GAUSS_5;
        switch(integration_degree)
        {
            case 1 : fem_integration_method = GeometryData::GI_GAUSS_1; break;
            case 2 : fem_integration_method = GeometryData::GI_GAUSS_2; break;
            case 3 : fem_integration_method = GeometryData::GI_GAUSS_3; break;
            case 4 : fem_integration_method = GeometryData::GI_GAUSS_4; break;
            case 5 : fem_integration_method = GeometryData::GI_GAUSS_5; break;
        }   
      
        for (auto & elem_i : r_fe_model_part.Elements())
        {
            Element::GeometryType& geom_i = elem_i.GetGeometry();
            const Element::GeometryType::IntegrationPointsArrayType& integration_points = geom_i.IntegrationPoints( fem_integration_method );

            for (auto & integration_point_i : integration_points)
            {
                NodeType::CoordinatesArrayType ip_coordinates = geom_i.GlobalCoordinates(ip_coordinates, integration_point_i.Coordinates());
                output_file << ip_coordinates[0] << " " << ip_coordinates[1] << " " << ip_coordinates[2] << std::endl;
            }
        }
        output_file.close();
        std::cout << "> Finished writing Gauss points of given FE-mesh." << std::endl;
    }

    // --------------------------------------------------------------------------
    void OutputResultsInRhinoFormat( int iteration )
    {
        if(iteration==1)
        {
            std::cout << "\n> Starting to write modified georhino file to visualize displacements of all control points relevant for reconstruction..." << std::endl;
            OutputModifiedGeoRhinoFileToIncludeAllRelevantControlPoints();
            std::cout << "> Finished writing modified georhino file." << std::endl;
        }

        std::cout << "\n> Adding displacements of control points into Rhino results file..." << std::endl;
        OutputResultingControlPointDisplacementsInRhinoFormat( iteration );
        std::cout << "> Finished adding displacements of control points into Rhino results file." << std::endl;
    }

    // --------------------------------------------------------------------------
    void OutputModifiedGeoRhinoFileToIncludeAllRelevantControlPoints()
    {
        std::ifstream input_file(mOriginalGeorhinoFile);
        std::ofstream georhino_output_file( mOutputFolder + "/" + mOriginalGeorhinoFile );
        
        std::string line; // last line to be read
        
        // copy all lines from the old file before the node section
        while(std::getline(input_file, line))
        {
            // if the line contains the word "NODE"
            if(find_substring(line, "NODE"))
            break;
            // else copy line
            georhino_output_file << line << std::endl;
        }
        
        // read and ignore node section
        while(std::getline(input_file, line))
        {
            // stop as soon as the line doesn't contain the word "NODE" anymore
            if(!find_substring(line, "NODE"))
                break;
        }
        
        // write node section
        
        int control_point_iterator = 0;
        for(auto & patch_i : mrReconstructionDataBase.GetPatchVector()) 
        {
            for(auto & control_point_i : patch_i.GetSurfaceControlPoints())
            {
                control_point_iterator++;
                if(control_point_i.IsRelevantForReconstruction())
                    georhino_output_file << "  NODE  " << control_point_iterator << "  X " << control_point_i.GetX0() << "  Y "<< control_point_i.GetY0() << "  Z " << control_point_i.GetZ0() << std::endl;
            }
        }
        georhino_output_file << line << std::endl;
        
        // substitute all the lines with "NODE_ID", ignore those with "GP_POINT_GEO", copy the others
        std::vector<ControlPoint*> control_points_vector;
        for(auto & patch_i : mrReconstructionDataBase.GetPatchVector()) 
        {
            for(auto & control_point_i : patch_i.GetSurfaceControlPoints())
                control_points_vector.push_back(&control_point_i);
        }
        
        control_point_iterator = 0;
        while(std::getline(input_file, line))
        {
            if(find_substring(line, "NODE_ID"))
            {
                control_point_iterator++;
                ControlPoint* cp = control_points_vector[control_point_iterator-1];
                
                if(cp->IsRelevantForReconstruction())
                    georhino_output_file << "  NODE_ID  " << control_point_iterator << "  W  " << cp->GetWeight() << std::endl;
                else
                    georhino_output_file << "  NODE_ID  0  X  " << cp->GetX0() << "  Y  " << cp->GetY0() << "  Z  " << cp->GetZ0() << "  W  " << cp->GetWeight() << std::endl;
            }
            else if(find_substring(line, "GP_POINT_GEO"))
                ;// ignore line
            else
                georhino_output_file << line << std::endl; // copy line
        }
        georhino_output_file.close();
    }

    // --------------------------------------------------------------------------
    void OutputResultingControlPointDisplacementsInRhinoFormat( int iteration )
    {
        std::ofstream restult_output_file;

        if(iteration==1)
        {
            restult_output_file.open( mOutputFolder + "/" + mRhinoResultsFile, std::ofstream::out | std::ofstream::trunc );        
            restult_output_file << "Rhino Post Results File 1.0" << std::endl;
        }
        else
            restult_output_file.open( mOutputFolder + "/" + mRhinoResultsFile, std::ofstream::out | std::ofstream::app );        
        
        restult_output_file << "Result \"Displacement\" \"Load Case\" " << std::to_string(iteration) << " Vector OnNodes" << std::endl;
        restult_output_file << "Values" << std::endl;

        
        int control_point_iterator = 0;
        for(auto & patch_i : mrReconstructionDataBase.GetPatchVector()) 
        {
            for(auto & control_point_i : patch_i.GetSurfaceControlPoints())
            {
                control_point_iterator++;
                if(control_point_i.IsRelevantForReconstruction())                
                    restult_output_file << control_point_iterator << " " << control_point_i.GetdX() << " " << control_point_i.GetdY() << " " << control_point_i.GetdZ() << std::endl;
            }
        }
        restult_output_file << "End Values" << std::endl;
        restult_output_file.close();
    }      

    // --------------------------------------------------------------------------
    bool find_substring(std::string str, std::string substring)
    {
        std::size_t found = str.find(substring);
        if(found != std::string::npos) 
            return true;
        return false;
    }

    // ==============================================================================

    /// Turn back information as a string.
    virtual std::string Info() const
    {
		return "ReconstructionOutputUtilities";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream &rOStream) const
    {
		rOStream << "ReconstructionOutputUtilities";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream &rOStream) const
    {
    }

private:

    ReconstructionDataBase& mrReconstructionDataBase;
    std::string mOutputFolder;
    std::string mOriginalGeorhinoFile;
    std::string mRhinoResultsFile;

    /// Assignment operator.
    //      ReconstructionOutputUtilities& operator=(ReconstructionOutputUtilities const& rOther);

    /// Copy constructor.
    //      ReconstructionOutputUtilities(ReconstructionOutputUtilities const& rOther);

}; // Class ReconstructionOutputUtilities
} // namespace Kratos.

#endif // RECONSTRUCTION_OUTPUT_UTILITIES_H
