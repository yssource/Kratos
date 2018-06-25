//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Peter Wilson
//                   Philipp Bucher
//


// System includes


// External includes


// Project includes
#include "assign_material_orientation_utility.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{

void AssignMaterialOrientationUtility::Execute(Parameters MethodParameters)
{
    // Parsing the parameters and selecting the method


    // Computing the orientation and assigning it to the element
    for (auto& r_elem : mrModelPart.Elements()) // TODO omp?
    {
        // TODO compute theta
        double theta = 0.0;





        // set required rotation in element
        r_elem.SetValue(MATERIAL_ORIENTATION_ANGLE, theta);

        KRATOS_INFO_IF("AssignMaterialOrientationUtility", mEchoLevel > 1)
            << "Element " << element.GetId() << "; orientation = " << theta << std::endl;
    }
}

void AssignMaterialOrientationUtility::WriteFiberAngles(const std::string& rFileName)
{
    std::ofstream element_data_file(rFileName);

    if (element_data_file.is_open())
    {
        element_data_file << "Begin ElementalData MATERIAL_ORIENTATION_ANGLE" << std::endl;
        for (auto& element : mrModelPart.Elements())
        {
            element_data_file << "\t" << element.Id()
                              << " " << element.GetValue(MATERIAL_ORIENTATION_ANGLE)
                              << std::endl;
        }
        element_data_file << "End ElementalData" << std::endl;

        element_data_file.close();
    }
    else
        KRATOS_WARNING("AssignMaterialOrientationUtility")
            << "Unable to open file " << rFileName << std::endl;

    std::cout << "Use \"cat " << rFileName << " >> your_mdpa_file.mdpa\" to append "
                << "the MATERIAL_ORIENTATION_ANGLE values to the mdpa-file" << std::endl;
}

}  // namespace Kratos.


