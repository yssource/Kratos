// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
//
// ==============================================================================

#ifndef CAD_RECONSTRUCTOR_H
#define CAD_RECONSTRUCTOR_H

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------

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
#include "reconstruction_data_base.h"
#include "cad_projection_utility.h"
#include "reconstruction_conditions/surface_displacement_mapping_condition.h"

// ==============================================================================

namespace Kratos
{
class CADReconstructor
{
  public:
    ///@name Type Definitions
    ///@{

    typedef boost::python::extract<int> ExtractInt;
    typedef Node<3> NodeType;
    typedef std::vector<NodeType::Pointer> NodeVector;   
    typedef std::vector<Patch> PatchVector; 
    typedef Element::GeometryType::IntegrationMethod IntegrationMethodType;

    /// Pointer definition of CADReconstructor
    KRATOS_CLASS_POINTER_DEFINITION(CADReconstructor);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    CADReconstructor( ReconstructionDataBase& reconstruction_data_base )
    : mrReconstructionDataBase( reconstruction_data_base )
    {      
    }

    /// Destructor.
    virtual ~CADReconstructor()
    {
    }

    // --------------------------------------------------------------------------
    void CreateSurfaceDisplacementMappingConditions( boost::python::list rParameterResolutionForCoarseNeighborSearch, 
                                                     int integration_degree, 
                                                     int max_iterations, 
                                                     double projection_tolerance )
    {
      std::cout << "\n> Starting to CreateSurfaceDisplacementMappingConditions...";
      
      PatchVector& patch_vector = mrReconstructionDataBase.GetPatchVector();
      ModelPart& fe_model_part = mrReconstructionDataBase.GetFEModelPart();
      IntegrationMethodType fem_integration_method = DefineIntegrationMethod( integration_degree );
      
      CADProjectionUtility FE2CADProjector( patch_vector, max_iterations, projection_tolerance );
      FE2CADProjector.Initialize( rParameterResolutionForCoarseNeighborSearch );
       
      std::cout << "> Starting to loop over integration points..." << std::endl;
      boost::timer timer;   
         
      for (ModelPart::ElementsContainerType::iterator elem_i = fe_model_part.ElementsBegin(); elem_i != fe_model_part.ElementsEnd(); ++elem_i)
      {
        Element::GeometryType& geom_i = elem_i->GetGeometry();
        const Element::GeometryType::IntegrationPointsArrayType& integration_pionts = geom_i.IntegrationPoints(fem_integration_method);

        for (auto & integration_point_i : integration_pionts)
        {
          int integration_point_number = &integration_point_i - &integration_pionts[0];
          NodeType::CoordinatesArrayType ip_coordinates = geom_i.GlobalCoordinates(ip_coordinates, integration_point_i.Coordinates());
          NodeType::Pointer node_of_interest = Node <3>::Pointer(new Node<3>(1, ip_coordinates));

          array_1d<double,2> parameter_values_of_nearest_point;
          array_1d<int,2> parameter_spans_of_nearest_point;
          int patch_index_of_nearest_point = -1;

          FE2CADProjector.DetermineNearestCADPoint( node_of_interest, 
                                                    parameter_values_of_nearest_point, 
                                                    parameter_spans_of_nearest_point, 
                                                    patch_index_of_nearest_point );
                                                                                                        
          mListOfConditions.push_back( SurfaceDisplacementMappingCondition::Pointer( new SurfaceDisplacementMappingCondition( geom_i,
                                                                                                                              fem_integration_method,
                                                                                                                              integration_point_number,
                                                                                                                              patch_vector[patch_index_of_nearest_point],
                                                                                                                              parameter_values_of_nearest_point,
                                                                                                                              parameter_spans_of_nearest_point )));
          }
      }
      std::cout << "> Time needed for looping over integration points: " << timer.elapsed() << " s" << std::endl;            
		  std::cout << "> Finished creating surface mapping conditions: " << std::endl;      
    }

    // --------------------------------------------------------------------------
    IntegrationMethodType DefineIntegrationMethod(int integration_degree )
    {
      IntegrationMethodType fem_integration_method;
      switch(integration_degree)
      {
        case 1 : fem_integration_method = GeometryData::GI_GAUSS_1; break;
        case 2 : fem_integration_method = GeometryData::GI_GAUSS_2; break;
        case 3 : fem_integration_method = GeometryData::GI_GAUSS_3; break;
        case 4 : fem_integration_method = GeometryData::GI_GAUSS_4; break;
        case 5 : fem_integration_method = GeometryData::GI_GAUSS_5; break;
      }
      return fem_integration_method;
    }

    // --------------------------------------------------------------------------
    void CreateCouplingConditionsOnAllCouplingPoints()
    {
    }
    
    // --------------------------------------------------------------------------
    void CreateDirichletConditions()
    {

    }   
    
    // --------------------------------------------------------------------------
    void IdentifyControlPointsRelevantForReconstruction()
    {

    }      

    // ==============================================================================

    /// Turn back information as a string.
    virtual std::string Info() const
    {
		return "CADReconstructor";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream &rOStream) const
    {
		rOStream << "CADReconstructor";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream &rOStream) const
    {
    }

  private:

    // ==============================================================================
    // Initialized by class constructor
    // ==============================================================================
    ReconstructionDataBase& mrReconstructionDataBase;

    // ==============================================================================
    // Additional variables
    // ==============================================================================
    std::vector<ReconstructionCondition::Pointer> mListOfConditions;

    /// Assignment operator.
    //      CADReconstructor& operator=(CADReconstructor const& rOther);

    /// Copy constructor.
    //      CADReconstructor(CADReconstructor const& rOther);

}; // Class CADReconstructor
} // namespace Kratos.

#endif // CAD_RECONSTRUCTOR_H
