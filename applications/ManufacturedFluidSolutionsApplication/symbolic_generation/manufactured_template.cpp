//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    AuthorName
//


// System includes


// External includes


// Project includes
#include "manufactured_template.h"
#include "includes/variables.h"
#include "includes/cfd_variables.h"


namespace Kratos
{

ManufacturedTemplate::ManufacturedTemplate(
    Properties::Pointer pProperties,
    Parameters::Pointer pParameters)
     : ManufacturedSolution(pProperties, pParameters)
{
    // Getting the fluid properties
    // Custom definition of the dynamics from the properties

    // Getting the manufactured settings
    Parameters default_parameters( R"(
    {
    })");

    mpParameters->ValidateAndAssignDefaults(default_parameters);

    // custom definition of the kinematics from the parameters
}


bool ManufacturedTemplate::IsInsideDomain(array_1d<double, 3>& rCoords)
{
    bool is_inside = false;
    // perform the check
    return is_inside;
}


// Velocity components definition
double ManufacturedTemplate::U1(array_1d<double, 3>& rCoords, double& rTime)
{
    double x = rCoords[0];
    double y = rCoords[1];
    return // substitute U1;
}

double ManufacturedTemplate::U2(array_1d<double, 3>& rCoords, double& rTime)
{
    double x = rCoords[0];
    double y = rCoords[1];
    return // substitute U2;
}


// Velocity time derivative definition
double ManufacturedTemplate::DU1_DT(array_1d<double, 3>& rCoords, double& rTime)
{
    double x = rCoords[0];
    double y = rCoords[1];
    return // substitute DU1_DT;
}

double ManufacturedTemplate::DU2_DT(array_1d<double, 3>& rCoords, double& rTime)
{
    double x = rCoords[0];
    double y = rCoords[1];
    return // substitute DU2_DT;
}


// Velocity first derivatives definition
double ManufacturedTemplate::DU1_DX1(array_1d<double, 3>& rCoords, double& rTime)
{
    double x = rCoords[0];
    double y = rCoords[1];
    return // substitute DU1_DX1;
}

double ManufacturedTemplate::DU1_DX2(array_1d<double, 3>& rCoords, double& rTime)
{
    double x = rCoords[0];
    double y = rCoords[1];
    return // substitute DU1_DX2;
}

double ManufacturedTemplate::DU2_DX1(array_1d<double, 3>& rCoords, double& rTime)
{
    double x = rCoords[0];
    double y = rCoords[1];
    return // substitute DU2_DX1;
}

double ManufacturedTemplate::DU2_DX2(array_1d<double, 3>& rCoords, double& rTime)
{
    double x = rCoords[0];
    double y = rCoords[1];
    return // substitute DU2_DX2;
}


// Velocity second derivatives definition
double ManufacturedTemplate::DDU1_DX11(array_1d<double, 3>& rCoords, double& rTime)
{
    double x = rCoords[0];
    double y = rCoords[1];
    return // substitute DDU1_DX11;
}

double ManufacturedTemplate::DDU1_DX22(array_1d<double, 3>& rCoords, double& rTime)
{
    double x = rCoords[0];
    double y = rCoords[1];
    return // substitute DDU1_DX22;
}

double ManufacturedTemplate::DDU2_DX11(array_1d<double, 3>& rCoords, double& rTime)
{
    double x = rCoords[0];
    double y = rCoords[1];
    return // substitute DDU2_DX11;
}

double ManufacturedTemplate::DDU2_DX22(array_1d<double, 3>& rCoords, double& rTime)
{
    double x = rCoords[0];
    double y = rCoords[1];
    return // substitute DDU2_DX22;
}


// Pressure and derivatives definition
double ManufacturedTemplate::P(array_1d<double, 3>& rCoords, double& rTime)
{
    double x = rCoords[0];
    double y = rCoords[1];
    return // substitute P;
}

double ManufacturedTemplate::DP_DX1(array_1d<double, 3>& rCoords, double& rTime)
{
    double x = rCoords[0];
    double y = rCoords[1];
    return // substitute DP_DX1;
}

double ManufacturedTemplate::DP_DX2(array_1d<double, 3>& rCoords, double& rTime)
{
    double x = rCoords[0];
    double y = rCoords[1];
    return // substitute DP_DX2;
}


}  // namespace Kratos.


