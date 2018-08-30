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
#include <math.h>
namespace Kratos
{

void AssignMaterialOrientationUtility::Execute(Parameters MethodParameters)
{
	std::cout << std::endl
			  << "Assigning material orientation angle to elements using method 11" << std::endl;
	auto specific_parameters = MethodParameters["method_specific_settings"];
	Vector3 global_vector;

	//global_vector is the global fiber direction given by user
	//read this global direction
	CheckAndReadVectors(specific_parameters, "global_fiber_direction", global_vector);

	//Normalize global_vector
	global_vector /= std::sqrt(inner_prod(global_vector, global_vector));

	auto current_process_info = mrModelPart.GetProcessInfo();

	// Declare working variables
	Matrix LCSOrientation;
	Vector3 localAxis1 = ZeroVector(3);
	Vector3 localAxis2 = ZeroVector(3);
	Vector3 localAxis3 = ZeroVector(3);

	// declaration of copy of global_vector
	Vector3 global_vector_copy;

	// Loop over all elements in part
	for (auto &element : mrModelPart.Elements())
	{

		// make a copy of global_vector
		global_vector_copy = global_vector;

		// get local axis in cartesian coordinates
		element.Calculate(LOCAL_ELEMENT_ORIENTATION, LCSOrientation, current_process_info);

		// get element local axis vectors (global cartesian)
		for (size_t i = 0; i < 3; i++)
		{
			localAxis1[i] = LCSOrientation(0, i);
			localAxis2[i] = LCSOrientation(1, i);
			localAxis3[i] = LCSOrientation(2, i);
		}

		// normalise local axis vectors (global cartesian)
		localAxis1 /= std::sqrt(inner_prod(localAxis1, localAxis1));
		localAxis2 /= std::sqrt(inner_prod(localAxis2, localAxis2));
		localAxis3 /= std::sqrt(inner_prod(localAxis3, localAxis3));

		// (Abaqus default projection)
		// http://130.149.89.49:2080/v6.8/books/gsa/default.htm?startat=ch05s03.html
		// Shell local axis 1 is the projection of Global X vector onto the shell surface.
		// If the Global X vector is normal to the shell surface,
		// the shell local 1-direction is the projection of the
		// Global Z vector onto the shell surface

		// First, check if specified global_vector is normal to the shell surface
		if (std::abs(inner_prod(global_vector_copy, localAxis1)) < 1E-6 && std::abs(inner_prod(global_vector_copy, localAxis2)) < 1E-6)
		{
			std::cout << std::endl
					  << "Global direction is perpendicular to element " << element.GetId()
					  << " so skipped assigning direction to this element " << std::endl;
			element.SetValue(MATERIAL_ORIENTATION_ANGLE,0);
		}
		else
		{

			// Second, project the global vector onto the shell surface
			// http://www.euclideanspace.com/maths/geometry/elements/plane/lineOnPlane/index.htm
			// vector to be projected = A
			// Surface normal = B
			Vector3 A = global_vector_copy;
			Vector3 B = localAxis3;

			Vector3 ACrossB;
			MathUtils<double>::CrossProduct(ACrossB, A, B);
			Vector3 projectedResult;
			MathUtils<double>::CrossProduct(projectedResult, B, ACrossB);
			//noramlize projected result
			projectedResult /= std::sqrt(inner_prod(projectedResult, projectedResult));

			// Third, find the angle between our projected direction and the
			// current shell localAxis1
			double cosTheta = inner_prod(localAxis1, projectedResult);
			double theta = std::acos(cosTheta);
			// make sure the angle is positively defined according to right
			// hand rule
			double dotCheck = inner_prod(localAxis2, projectedResult);
			if (dotCheck < 0.0)
			{
				// theta is currently negative, flip to positive definition
				theta *= -1.0;
			}

			// set required rotation in element
			element.SetValue(MATERIAL_ORIENTATION_ANGLE, theta);
		}
	}
	std::cout << std::endl
			  << ".......done assigning direction for all elements......." << std::endl;
}

void AssignMaterialOrientationUtility::CheckAndReadVectors(Parameters ThisParameters, const std::string KeyName, Vector3 &rVector)
{
	

	if (ThisParameters[KeyName].size() != 3)
	{
		KRATOS_ERROR << "\" " << KeyName << "\" is not of size 3" << std::endl;
	}

	
	rVector[0] = ThisParameters[KeyName][0].GetDouble();
	rVector[1] = ThisParameters[KeyName][1].GetDouble();
	rVector[2] = ThisParameters[KeyName][2].GetDouble();

	if (inner_prod(rVector, rVector) < 1E-6)
	{
		KRATOS_ERROR << "Vector \" " << KeyName << "\" has zero length" << std::endl;
	}
}
void AssignMaterialOrientationUtility::WriteFiberAngles(const std::string &rFileName)
{
	std::ofstream element_data_file(rFileName);

	if (element_data_file.is_open())
	{
		element_data_file << "Begin ElementalData MATERIAL_ORIENTATION_ANGLE" << std::endl;
		for (auto &r_elem : mrModelPart.Elements())
		{
			element_data_file << "\t" << r_elem.Id()
							  << " " << r_elem.GetValue(MATERIAL_ORIENTATION_ANGLE)
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

} // namespace Kratos.
