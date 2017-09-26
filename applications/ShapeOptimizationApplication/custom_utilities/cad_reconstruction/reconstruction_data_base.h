// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
//
// ==============================================================================

#ifndef RECONSTRUCTION_DATA_BASE_H
#define RECONSTRUCTION_DATA_BASE_H

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "basic_nurbs_brep_handling\patch.h"
#include "basic_nurbs_brep_handling\brep_element.h"
#include "cad_model_reader.h"

// ==============================================================================

namespace Kratos
{
class ReconstructionDataBase
{
public:

    /// Pointer definition of ReconstructionDataBase
    KRATOS_CLASS_POINTER_DEFINITION(ReconstructionDataBase);

    /// Default constructor.
    ReconstructionDataBase(ModelPart& fe_model_part, boost::python::dict cad_geometry, boost::python::dict cad_integration_data)
	: mpFEModelPart(fe_model_part),
      mrCADGeometry(cad_geometry),
      mrCADIntegrationData(cad_integration_data)
    {     
    }

    /// Destructor.
    virtual ~ReconstructionDataBase()
    {
    }

    // ------------------------------------------------------------------------------
    void Create()
    {
        CADModelReader cad_reader(mrCADGeometry,mrCADIntegrationData);
        
        cad_reader.ReadGeometry(mPatches);
            if(len(mrCADIntegrationData.keys())>0)
        cad_reader.ReadIntegrationData(mBrepElements);      
    }

    // // --------------------------------------------------------------------------
    // void UpdateControlPoints(PatchVector& patches, Vector& ds)
    // {
    //     /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //     // 1. Step: Update C++ data base
    //     /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //     // Map to identify control point for given global id (needed for python update later)
    //     std::map<unsigned int, ControlPoint*> control_point_corresponding_to_global_id;

    //     for (PatchVector::iterator patch_i = patches.begin(); patch_i != patches.end(); ++patch_i)
    //     {
    //         for (ControlPointVector::iterator cp_i = patch_i->GetSurface().GetControlPoints().begin(); cp_i != patch_i->GetSurface().GetControlPoints().end(); ++cp_i)
    //         {
    //             if(cp_i->IsRelevantForMapping())
    //             {
    //                 // Updating c++ data base
    //                 unsigned int cp_mapping_matrix_id = cp_i->GetMappingMatrixId();
    //                 cp_i->SetdX( ds[3*cp_mapping_matrix_id+0] );
    //                 cp_i->SetdY( ds[3*cp_mapping_matrix_id+1] );
    //                 cp_i->SetdZ( ds[3*cp_mapping_matrix_id+2] );

    //             }
    //             // Filling map to be used later
    //             unsigned int cp_global_id = cp_i->GetGlobalId();
    //             control_point_corresponding_to_global_id[cp_global_id] = &(*cp_i);
    //         }
    //     }

    //     /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //     // 2. Step: Update pyhon data base
    //     /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        
    //     // loop over patches / faces in cad geometry file
    //     for (int i = 0; i < boost::python::len(mr_cad_geometry_in_json["faces"]); i++)
    //     {
    //         for (int cp_idx = 0; cp_idx < boost::python::len(mr_cad_geometry_in_json["faces"][i]["surface"][0]["control_points"]); cp_idx++)
    //         {
    //             unsigned int global_id = extractInt(mr_cad_geometry_in_json["faces"][i]["surface"][0]["control_points"][cp_idx][0]);

    //             ControlPoint* cp_j = control_point_corresponding_to_global_id[global_id];

    //             double sx = cp_j->GetX();
    //             double sy = cp_j->GetY();
    //             double sz = cp_j->GetZ();

    //             // Update python data base which is store as a reference (so it is also updated on python level)
    //             mr_cad_geometry_in_json["faces"][i]["surface"][0]["control_points"][cp_idx][1][0] = sx;
    //             mr_cad_geometry_in_json["faces"][i]["surface"][0]["control_points"][cp_idx][1][1] = sy;
    //             mr_cad_geometry_in_json["faces"][i]["surface"][0]["control_points"][cp_idx][1][2] = sz;
    //         }
    //     }
    // }

    // ------------------------------------------------------------------------------
    ModelPart& GetFEModelPart()
    {
        return mpFEModelPart;   
    }    

    // ------------------------------------------------------------------------------
    std::vector<Patch>& GetPatchVector()
    {
        return mPatches;   
    }    

    // ------------------------------------------------------------------------------
    std::vector<BREPElement>& GetBREPElementsPatches()
    {
        return mBrepElements;   
    } 

    // ==============================================================================


    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "ReconstructionDataBase";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "ReconstructionDataBase";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
    }


private:
    
    // ==============================================================================
    // Initialized by class constructor
    // ==============================================================================
    ModelPart& mpFEModelPart;
    boost::python::dict mrCADGeometry;
	boost::python::dict mrCADIntegrationData;
    std::vector<Patch> mPatches;
	std::vector<BREPElement> mBrepElements;

}; // Class ReconstructionDataBase

}  // namespace Kratos.

#endif // RECONSTRUCTION_DATA_BASE_H
