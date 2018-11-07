/*
//  KRATOS  _____________
//         /  _/ ____/   |
//         / // / __/ /| |
//       _/ // /_/ / ___ |
//      /___/\____/_/  |_| Application
//
//  Main authors:   Lukas Rauch
*/

// System includes
#include "includes/define.h"
#include "includes/variables.h"
#include "utilities/math_utils.h"
#include <typeinfo>

// #include <math.h>
// #include "Tools.h"
#include <fstream>

// External includes

// Project includes
#include "iga_base_element.h"
#include "iga_beam_element.h"
#include "iga_application_variables.h"
#include "custom_elements/base_discrete_element.h"
#include "custom_utilities/iga_debug.h"

namespace Kratos {

Element::Pointer IgaBeamElement::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    auto geometry = GetGeometry().Create(ThisNodes);

    return Kratos::make_shared<IgaBeamElement>(NewId, geometry,
        pProperties);
}

void IgaBeamElement::GetDofList(
    DofsVectorType& rElementalDofList,
    ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    rElementalDofList.resize(NumberOfDofs());

    for (std::size_t i = 0; i < NumberOfNodes(); i++) {
        SetElementDof(rElementalDofList, i, 0, DISPLACEMENT_X);
        SetElementDof(rElementalDofList, i, 1, DISPLACEMENT_Y);
        SetElementDof(rElementalDofList, i, 2, DISPLACEMENT_Z);
        SetElementDof(rElementalDofList, i, 3, DISPLACEMENT_ROTATION);
    }

    KRATOS_CATCH("")
}

void IgaBeamElement::EquationIdVector(
    EquationIdVectorType& rResult,
    ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    rResult.resize(NumberOfDofs());

    for (std::size_t i = 0; i < NumberOfNodes(); i++) {
        SetElementEquationId(rResult, i, 0, DISPLACEMENT_X);
        SetElementEquationId(rResult, i, 1, DISPLACEMENT_Y);
        SetElementEquationId(rResult, i, 2, DISPLACEMENT_Z);
        SetElementEquationId(rResult, i, 3, DISPLACEMENT_ROTATION);
    }

    KRATOS_CATCH("")
}

void IgaBeamElement::Initialize()
{
}

void IgaBeamElement::CalculateAll(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo,
    const bool ComputeLeftHandSide,
    const bool ComputeRightHandSide)
{
    KRATOS_TRY;

    // tmporary debug data
    auto expected_data = Parameters(GetValue(DEBUG_EXPECTED_DATA));
   

    // get integration data
    
    const double integration_weight = GetValue(INTEGRATION_WEIGHT);
    Vector& shape_function_values = GetValue(SHAPE_FUNCTION_VALUES);
    Matrix& shape_derivatives = GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);
    /// Testprit
    IgaDebug::CheckVector(expected_data,"shape_0", shape_function_values);

    // get properties
    const auto& properties = GetProperties();

    const double _emod = properties[YOUNG_MODULUS];
    const double _gmod = properties[SHEAR_MODULUS];
    const double _area = properties[CROSS_AREA];
    const double prestress = properties[PRESTRESS_CAUCHY];
    const double Phi = properties[PHI]; 
    const double Phi_der_0 = properties[PHI_0_DER];
    const double _m_inert_y = properties[MOMENT_OF_INERTIA_Y];
    const double _m_inert_z = properties[MOMENT_OF_INERTIA_Z];
    const double _mt_inert = properties[MOMENT_OF_INERTIA_T];

 
    // TODO: Add stiffness andd force

    std::cout << expected_data["t"] << std::endl;

    // Reference and Actual Configuration
    Vector3 R1;
    Vector3 R2;
    double A;
    double B;
    Vector3 r1;
    Vector3 r2;
    double a;
    double b;

    ComputeGeometryReference(R1, R2, A, B);
    ComputeGeometryActual(r1, r2, a, b);

    /// Debug Check
    IgaDebug::CheckVector(expected_data, "R_1", R1);
    IgaDebug::CheckVector(expected_data, "R_2", R2); 
    IgaDebug::CheckDouble(expected_data, "A", A);
    IgaDebug::CheckDouble(expected_data, "B", B);    
    IgaDebug::CheckVector(expected_data, "r_1", r1);
    IgaDebug::CheckVector(expected_data, "r_2", r2);
    IgaDebug::CheckDouble(expected_data, "a", a);
    IgaDebug::CheckDouble(expected_data, "b", b);

    // Create empty Stiffness Matrix
    MatrixType _gke; 
    VectorType _gfie; 

    ElementStiffnessMatrixNonlinear(_emod, _gmod, _area, _m_inert_y, _m_inert_z, _mt_inert, integration_weight, _gke, _gfie);
    // IgaDebug::CheckVector(expected_data,"stiffness", _gke);
    // IgaDebug::CheckVector(expected_data,"external_force", _gkfi);







//# Testarea ###
    // Vector3 coords;     // Coordinates of the Nodes
    // Vector coords2;
    // for (size_t i = 0; i < NumberOfNodes(); i++){
    //         coords = GetGeometry()[i].GetInitialPosition();
    //         coords2 = GetGeometry()[i].GetInitialPosition();
    //         std::cout << "GetGeometry: " << GetGeometry()[i] << std::endl;
    //         std::cout << "coords: " << coords << std::endl;
    //         std::cout << "coords2: " << coords2 << std::endl;
    // }

    // std::cout << "GetGeometry().size() = " << GetGeometry().size() << std::endl;


    KRATOS_CATCH("");
}


// //#################################################################################
// //#################################################################################
// //#
// //#                         +++ Iga_Beam_Element +++
// //#
// //#################################################################################
// //#################################################################################
// //#
// /** 
//      * 
//      * @return
//      * 
//      * @param[out]   
//      * 
//      * @author L.Rauch (10/2018)
//      * 
//      * @note by A.Bauer (10/2014)
//      */
// //#--------------------------------------------------------------------------------
// void IgaBeamElement::IgaBeamElement(
//         int _ID, 
//         Part* _Part_Prt,
//         std::vector<Node*> _Node_Vec,
//         MaterialBasis* _Mat_Prt,
//         PropertyIgaBeamElement* _Prop_Prt,
//         Nurbs1D* _Nurbs_Curve,
//         int _knotspan_index, 
//         std::vector<std::vector<double,2>> _phi_axis_n,
//         std::vector<double,3> _t0,
//         std::vector<EvalPtBasis*> _eval_points, 
//         KnotSpan_Belonging _elementbelongig, 
//         PropertyNURBS_Curve_Stab* _prop_stab )
//     :ElementCurve_Nurbs(
//         _ID, _Part_Prt, _Node_Vec, _Mat_Prt, 
//         _Nurbs_Curve, _knotspan_index, _Prop_Prt,
//         _eval_points, _elementbelongig, _prop_stab),
//         Prop_Prt(_Prop_Prt),
//         phi_axis_n(_phi_axis_n)

// {
// // Set disretisation type
// discretisation_type = nurbs1d; 

// //# TODO ### function [ GetPDeg() ] implementieren
// p_deg = _Nurbs_Curve -> GetPDeg();

// dof_node = 4

// N_dof = dof_node * (p_deg+1); 

// // Activate the DOSs at each Node that needed for this Element Type (ElementMembrane1)
// std::vector<int> act_dofs; 

// // get Actual DOFs
// GetDofTypesPerNode(act_dofs);

// //# TODO ### funktion [Set_Active_DOFs() ] implementieren
// for(size_t i=0; i < Node_Vec.size(); i++ ){
//     Node_Vec[i] -> SetActiveDOFs(act_dofs);
// }

// // Set Tollereance Level

// int tol = 1.0e-8;

// }  




// #################################################################################
// #################################################################################
// #
// #                       +++ Get_Dof_Types_PerNode +++
// #
// #################################################################################
// #################################################################################
// #
/** gets the degree of freedom types of this element at each node
     * 
     * @return
     * 
     * @param[out]   _act_dofs       Acting Degrees of Freedom per Node
     * 
     * @author L.Rauch (10/2018)
     * 
     * @note by A.Bauer (10/2014)
     */
// #--------------------------------------------------------------------------------

void IgaBeamElement::GetDofTypesPerNode(
    std::vector<int>& _act_dofs){
    
//# TODO ### die Anzahl der Freiheitsgrade pro Node festlegen. Irgendwo ganz am Anfang ( in caratt++ in erster Funtion)
    int dof_node = 4;

    _act_dofs.resize(dof_node);         // _act_dof Vektor auf richtige Länge einstellen
    _act_dofs[0] = DISPLACEMENT_X;
    _act_dofs[1] = DISPLACEMENT_Y;
    _act_dofs[2] = DISPLACEMENT_Z; 
    if(dof_node == 4 ) _act_dofs[3] = DISPLACEMENT_ROTATION;
}

// #################################################################################
// #################################################################################
// #
// #                          +++ Iga_Beam_Element +++
// #
// #################################################################################
// #################################################################################
// #
/** 
     * @author L.Rauch (10/2018)
     * 
     * @note by A.Bauer (10/2014)
     */
// #--------------------------------------------------------------------------------

IgaBeamElement::~IgaBeamElement(void)
{ 
}



//#################################################################################
//#################################################################################
//#
//#                 +++ Stiffnes Matrix Element Linear +++
//#
//#################################################################################
//#################################################################################
//#
/** Evaluates the Stiffnes Matrix and Residual Vector at the Integration Point
     * 
     * @param[in]   _u_act        NURBS coordinates of integration pont u-direction
     * @param[in]   _emod         Youns Modulus
     * @param[in]   _area         Area of Cross Section
     * @param[in]   _m_inert      Moment of Intertia
     * @param[out]  _gke          Elemental Stiffness Matrix per Gauss Point
     * 
     * @author  L.Rauch (10/2018)
     * 
     * @note    By M.Breitenberger(08/2012)
     */
//#--------------------------------------------------------------------------------
void IgaBeamElement::StffnessMatrixElementLinear(
    double _emod,
    double _gmod,
    double _area,
    double _m_inert_n,
    double _m_inert_v,
    double _mt_inert,
    double _dl,
        MatrixType& gke)
{
    gke.resize(NumberOfDofs(), NumberOfDofs());

    Vector& shape_function = GetValue(SHAPE_FUNCTION_VALUES);
    Matrix& shape_derivatives = GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);

    // Declaration 
    Vector3 R_1;   // 1st Derivative of the Curve Undeformed Config
    Vector3 R_2;   // 2nd Derivative of the Curve Undeformed Config
    double A;               // Length of the Base Vector
    double B;               // Curvature Including the Metric
    double B_n;
    double B_v;
    double C_12;
    double C_13;
    double Phi;
    double Phi_der;

    Vector r_1;   // 1st Derivative of the Curve Deformed Config. 
    Vector r_2;   // 2nd DErivative of the Curve Deformed Config.
    double a;               // Length of the Base Vector in Deformed Config.
    double b;               // Curvatur Including the Metric in Defomred Config. 
    double b_n;
    double b_v;
    double c_12;
    double c_13;
    double phi;
    double phi_der;
    
    // Declaration of Pricipal Axis
    Vector3 N;     // Principal Axis 1 of Cross Section in UNdeformed Config.
    Vector3 V;   // Prinzipal Axis 2 of Cross Section in UNdeformed Config. 
    Vector3 n;   // Principal Axis 1 Of Cross Section in Defomred Config.
    Vector3 v;  // Principal Axis 2 of Cross Section in Deformed Config. 
    Vector3 N0;    // Principal Axis 1 of Cross Section in Reference Config. 
    Vector3 V0;    // Principal Axis 2 of Cross Section in Reference Config. 


    // Material and Cross Section
    double emod_A = _emod * _area;
    double emod_I_n = _emod * _m_inert_n;
    double emod_I_v = _emod * _m_inert_v;
    double emod_It = _gmod * _mt_inert; 

    // Membrane Stiffness
    Matrix kem(NumberOfDofs(), NumberOfDofs());

    // Bending Stiffness (Principal Axis 1)
    Matrix keb_n(NumberOfDofs(), NumberOfDofs());

    // Bending Stiffness (Principal Axis 2) 
    Matrix keb_v(NumberOfDofs(), NumberOfDofs());

    // Torsional Stiffness
    Matrix ket_n(NumberOfDofs(), NumberOfDofs());
    Matrix ket_v(NumberOfDofs(), NumberOfDofs());


    // Get Previous Results
    Vector3 tmp_ini_disp;
    tmp_ini_disp.resize(4);
    tmp_ini_disp.clear();

    Phi = 0;
    Phi_der = 0;
    phi = 0;
    phi_der = 0;

    std::vector<int> act_dofs;
    GetDofTypesPerNode(act_dofs); 

    for (size_t i = 0; i < NumberOfNodes(); i++)
    {
        tmp_ini_disp = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT); 
        if (size(act_dofs) > 3)
        {
            phi += shape_function(i) * tmp_ini_disp[3];
            phi_der += shape_derivatives(0,i) * tmp_ini_disp[3];
        }
    }

    // Compute Geometry Reference
    ComputeGeometryReference(R_1, R_2, A, B);
    ComputeGeometryInitial(r_1, r_2, a, b);

    // ComputeCrossSectionGeometryReference(R_1, R_2, t0_0, N, V, N0, V0, B_n, B_v, C_12, C_13, Phi, Phi_der); 



}

//#################################################################################
//#################################################################################
//#
//#                 +++ Stiffness Matrix Element Nonlinear +++
//#
//#################################################################################
//#################################################################################
//#
/** Calculates the non-linear element stiffness matrix
     * 
     * @param[in]   _u_act        NURBS coordinates of integration pont u-direction
     * @param[in]   _emod         Youns Modulus
     * @param[in]   _area         Area of Cross Section
     * @param[in]   _m_inert      Moment of Intertia
     * @param[out]  _gke          Elemental Stiffness Matrix per Gauss Point
     * @param[out]  _gfie         Internal Force Vector per Gauss Point 
     * 
     * @author L.Rauch (10/2018)
     * 
     * @note
     */
//#--------------------------------------------------------------------------------
void IgaBeamElement::ElementStiffnessMatrixNonlinear(
    double _emod,
    double _gmod,
    double _area,
    double _m_inert_y,
    double _m_inert_z,
    double _mt_inert,
    double _dl,
        MatrixType& _gke,
        VectorType& _gfie){
        

    // SetUp empty stiffness Matrix
    // Matrix _gke;
    _gke.resize(NumberOfDofs(), NumberOfDofs());
    _gke.clear(); 

    // Vector _gfie;
    _gfie.resize(NumberOfDofs());
    _gfie.clear(); 
        
    Vector func;    // Vector for Basis function values
    Vector deriv1;  // Vector for Basis function Values
    Vector deriv2;  // Vector for Basis function values 

// //# TODO ############
// // compute_shape_Function ( output: func )
// // Compute_Shape_Function_Deriv1 ( output: deriv1 )
// // Compute_Shape_Function_Deriv2 ( output: deriv2 )
// //###################

    // declaration
    Vector R_1;   // 1st Derivative of the Curve Undeformed Config
    Vector R_2;   // 2nd Derivative of the Curve Undeformed Config
    double A;               // Length of the Base Vector
    double B;               // Curvature Including the Metric
    Vector r_1;   // 1st Derivative of the Curve Deformed Config. 
    Vector r_2;   // 2nd DErivative of the Curve Deformed Config.
    double a;               // Length of the Base Vector in Deformed Config.
    double b;               // Curvatur Including the Metric in Defomred Config. 

    // Further Declarations
    double B_n;
    double B_v;
    double C_12;
    double C_13;
    double Phi;
    double Phi_der;
    double b_n;
    double b_v;
    double c_12;
    double c_13;
    double phi;
    double phi_der;

// Declaration of Pricipal Axis
//# TODO ### How to use Vector3. Do i need to #include vector.h ?

    Vector3 N;     // Principal Axis 1 of Cross Section in UNdeformed Config.
    Vector3 V;   // Prinzipal Axis 2 of Cross Section in UNdeformed Config. 
    Vector3 n;   // Principal Axis 1 Of Cross Section in Defomred Config.
    Vector3 v;  // Principal Axis 2 of Cross Section in Deformed Config. 
    Vector3 N0;    // Principal Axis 1 of Cross Section in Reference Config. 
    Vector3 V0;    // Principal Axis 2 of Cross Section in Reference Config. 


    // Material And Cross Section
    double emod_A = _emod * _area; 
    double emod_I_n = _emod * _m_inert_z;
    double emod_I_v = _emod * _m_inert_y; 
    double gmod_It = _gmod * _mt_inert; 


    // Prestress
    //# TODO ### Introduce Presstress 
    // double prestress = [var properties] -> [ get_act_Presstress() ];    // Prestress in Normal Direction
    // double prestress_bend_n = [var propperties] -> [ get_act_Presstress_bend_n() bzw. get_act_Presstress_bend1() ];      // Prestress arround Z Axis
    // double prestress_bend_v = [ var propperties ] -> [ get_act_Presstress_bend_v()];     // Prestress arrund the Y Axis
    // double prestress_tor = [ var propperties ] -> [ get_act_Presstress_Tor() ];      // Torsional Presstress
    // bool prestress_bend1_auto = prop_prt -> GetPrestressBend1Auto(); 
    // bool prestress_bend2_auto = prop_prt -> GetPrestressBend2Auto();
    // bool prestress_tor_auto = prop_prt -> GetPrestressTorAuto(); 


    // Declarate Empty Stiffness Matixes

    // Membrane Stiffness Matrix
    Matrix kem(NumberOfDofs(), NumberOfDofs());

    // Bendign Stiffness Matrix (Principal Axis N)
    Matrix keb_n(NumberOfDofs(), NumberOfDofs());

    // Bendign Stiffness Matrix (Principal Axis V)
    Matrix keb_v(NumberOfDofs(), NumberOfDofs());

    // Bending Membrane Interaction Stiffness Matrix (Pricipal Axis N) 
    Matrix keb_ekn(NumberOfDofs(), NumberOfDofs());

    // Bending Membrane Interaction Stiffness Matrix (Principal Axis V)
    Matrix keb_ekv(NumberOfDofs(), NumberOfDofs());

    // Ineraction Bending Stiffness
    Matrix keb_knkv(NumberOfDofs(), NumberOfDofs());

    // Torsional Stiffness
    Matrix ket_n(NumberOfDofs(), NumberOfDofs()); 

    // Torsional Stiffness 
    Matrix ket_v(NumberOfDofs(), NumberOfDofs());

    // Torsional Stiffness
    Matrix ket(NumberOfDofs(), NumberOfDofs()); 


//     // Compute the Vectors R1 R2 and the Lentghts A and B in undeformed and deformed state
//     ComputeGeometryReference(_u_act, deriv1, deriv2, R_1, R_2, A, B);
//     ComputeGeometryActual(_u_act, deriv1, deriv2, r_1, r_2, a, b); 


    // get Previous Results
    Vector tmp_dof; 
    tmp_dof.resize(4,false);
    tmp_dof.clear();
    Vector tmp_ini_dof;
    tmp_ini_dof.resize(4,false);
    tmp_ini_dof.clear();

    std::vector<int> act_dofs;
    GetDofTypesPerNode(act_dofs); 
    
    phi = 0;
    phi_der = 0;
    Phi = 0;
    Phi_der = 0;

    









//     // std::cout << expected_data["t"].GetDouble() << std::endl;
//     // std::cout << expected_data["external_forces"].GetVector() << std::endl;

    // KRATOS_CATCH("")
}

//#################################################################################
//#################################################################################
//#
//#                      +++ Compute Geometry Initial +++
//#
//#################################################################################
//#################################################################################
//#
/** Computes the base Vector
     * 
     * @return
     * 
     * @param[in]   _deriv1          1st Derivative of the Curve in Referens Configuration
     * @param[in]   _deriv2          1st Derivative of the Curve in Referens Configuration 
     * @param[out]  r1              
     * @param[out]  r2              
     * @param[out]  a_ini           Length of the initial tangent Vector
     * @param[out]  b_ini           
     * 
     * @author L.Rauch (10/2018)
     * 
     * @note by A.Bauer (10/2017)
     */
//#--------------------------------------------------------------------------------
void IgaBeamElement::ComputeGeometryInitial(
        Vector& r1,
        Vector& r2,
        double& a_ini,
        double& b_ini)
{
    r1.resize(3);
    r1.clear();
    r2.resize(3);
    r2.clear();

    Matrix& shape_derivatives = GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);

    // Get previous Results
    Vector3 coords;         // Create the Coordinats Vector of the Nodes
    Vector3 tmp_disp_ini;       // Create the initial Displacement Vector of the Nodes
    

    for (size_t i = 0; i < NumberOfNodes(); i++){
        coords = GetGeometry()[i].GetInitialPosition();
        tmp_disp_ini = 0; //# TODO ### Get.Geometry INITIAL DISPLACEMENTS + GetGeometry()[i].FastGetSolutionStepValue(PREDISPLACEMENT);
    
        coords[0] += tmp_disp_ini[0];
        coords[1] += tmp_disp_ini[1];
        coords[2] += tmp_disp_ini[2];      
        
        r1[0] += shape_derivatives(0,i) * coords[0];
        r1[1] += shape_derivatives(0,i) * coords[1];
        r1[2] += shape_derivatives(0,i) * coords[2];

        r2[0] += shape_derivatives(1,i) * coords[0];
        r2[1] += shape_derivatives(1,i) * coords[1];
        r2[2] += shape_derivatives(1,i) * coords[2]; 
    }
    // Set Tollerance if not done bevore
    float tol = 1.0e-8;

    a_ini = norm_2(r1);     // Length of the base Vector

    double tmp = inner_prod(r2, r2) - std::pow( inner_prod(r1, r2), 2) /std::pow(a_ini, 2); 

    // Bending
    if (fabs(tmp) > tol){
         b_ini = std::sqrt(tmp);
    }else{
        b_ini = 0;
    }
}

//#################################################################################
//#################################################################################
//#
//#                      +++ Compute Geometry Initial +++
//#
//#################################################################################
//#################################################################################
//#
/** Computes the base Vector
     * 
     * @return
     * 
     * @param[in]   _deriv1         1st Derivative of the Curve in Referens Configuration
     * @param[in]   _deriv2         2nd Derivative of the Curve in Referens Configuration 
     * @param[in]   _deriv3         3rd Derivative of the Curve in Referens Configuration 
     * @param[out]  r1              
     * @param[out]  r2              
     * @param[out]  r3              
     * @param[out]  a_ini           Length of the initial tangent Vector
     * @param[out]  b_ini           
     * 
     * @author L.Rauch (10/2018)
     * 
     * @note by A.Bauer (02/2017)
     */
//#--------------------------------------------------------------------------------
void IgaBeamElement::ComputeGeometryInitial(
        Vector& r1,
        Vector& r2,
        Vector& r3,
        double& a_ini,
        double& b_ini)
{
    r1.resize(3);
    r1.clear();
    r2.resize(3);
    r2.clear();
    r3.resize(3);
    r3.clear();

    Matrix& shape_derivatives = GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);    
    
    // Get previous Results
    Vector3 coords;         // Create the Coordinats Vector of the Nodes
    Vector3 tmp_disp_ini;       // Create the initial Displacement Vector of the Nodes
    

    for (size_t i = 0; i < NumberOfNodes(); i++){
        coords = GetGeometry()[i].GetInitialPosition();
        tmp_disp_ini = 0; //# TODO ### Get.Geometry INITIAL DISPLACEMENTS + GetGeometry()[i].FastGetSolutionStepValue(PREDISPLACEMENT);
    
        coords[0] += tmp_disp_ini[0];
        coords[1] += tmp_disp_ini[1];
        coords[2] += tmp_disp_ini[2];      
        
        r1[0] += shape_derivatives(0,i) * coords[0];
        r1[1] += shape_derivatives(0,i) * coords[1];
        r1[2] += shape_derivatives(0,i) * coords[2];

        r2[0] += shape_derivatives(1,i) * coords[0];
        r2[1] += shape_derivatives(1,i) * coords[1];
        r2[2] += shape_derivatives(1,i) * coords[2]; 
        
        r3[0] += shape_derivatives(2,i) * coords[0];
        r3[1] += shape_derivatives(2,i) * coords[1];
        r3[2] += shape_derivatives(2,i) * coords[2]; 
    }
    // Set Tollerance if not done bevore
    float tol = 1.0e-8;

    a_ini = norm_2(r1);     // Length of the base Vector

    double tmp = inner_prod(r2, r2) - std::pow( inner_prod(r1, r2), 2) /std::pow(a_ini, 2); 

    // Bending
    if (fabs(tmp) > tol){
         b_ini = std::sqrt(tmp);
    }else{
        b_ini = 0;
    }
} 

//#################################################################################
//#################################################################################
//#
//#                      +++ Compute_Geometry_Reference +++
//#
//#################################################################################
//#################################################################################
//#
/** Computes the base Vector
     * 
     * @return
     * 
     * @param[in]   _u_act          NURBS coordinates of integration pont u-direction
     * @param[in]   _deriv1         1st derivatives of the Basis functions at u
     * @param[in]   _deriv2         2nd derivatives of the Basis functions at u
     * @param[out]  R1              1st Derivative of the Curve in Referens Configuration
     * @param[out]  R2              2nd Derivative ot the Curve in Referens Configuration
     * @param[out]  A           Length of the tangent Vector
     * @param[out]  B           
     * 
     * @author L.Rauch (10/2018)
     * 
     * @note by A.Bauer (10/2014)
     */
//#--------------------------------------------------------------------------------
void IgaBeamElement::ComputeGeometryReference(
        Vector3& R1,
        Vector3& R2, 
        double& A, 
        double& B)
{
    R1.resize(3);  // Resize 1st Derivative of the Curve
    R2.resize(3);
    R1.clear();     // Clear 1st Derivative of the Curve
    R2.clear();     // Clear 2nd Derivative of the Curve

    // Get the 1st and 2nd Shape Function Deriatides from the ModelPart
    Matrix& shape_derivatives = GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);      

    // Computation of the Basis Functions
    for (size_t i = 0; i < NumberOfNodes(); i++){
       
        const Node<3>& node = GetGeometry()[i];
        
        R1 += shape_derivatives(0,i) * node.GetInitialPosition();
        R2 += shape_derivatives(1,i) * node.GetInitialPosition();
    }

    // Set Tollerance if not done bevore
    float tol = 1.0e-8;

    A = norm_2(R1);     // Length of the base Vector

    double tmp = inner_prod(R2, R2) - std::pow( inner_prod(R1, R2), 2) /std::pow(A, 2); 

    if(fabs(tmp) > tol)
        B = std::sqrt(tmp);
    else 
        B = 0;  
}

//#################################################################################
//#################################################################################
//#
//#                      +++ Compute_Geometry_Reference +++
//#
//#################################################################################
//#################################################################################
//#
/** Computes the base Vector
     * 
     * @return
     * 
     * @param[in]   _u_act          NURBS coordinates of integration pont u-direction
     * @param[in]   _deriv1         1st derivatives of the Basis functions at u
     * @param[in]   _deriv2         2nd derivatives of the Basis functions at u
     * @param[in]   _deriv3         3rd derivatives of the Basis functions at u
     * @param[out]  R1              1st Derivative of the Curve in Referens Configuration
     * @param[out]  R2              2nd Derivative ot the Curve in Referens Configuration
     * @param[out]  R3              3rd Derivative ot the Curve in Referens Configuration
     * @param[out]  A_ref          Length of the tangent Vector
     * @param[out]  B          
     * 
     * @author L.Rauch (10/2018)
     * 
     * @note by A.Bauer (02/2018)
     */
//#--------------------------------------------------------------------------------

void IgaBeamElement::ComputeGeometryReference(
        Vector& R1,
        Vector& R2,
        Vector& R3,
        double& A,
        double& B){

    R1.resize(3);
    R1.clear();        // Clears 1st Derivative of the Curve
    R2.resize(3);
    R2.clear();        // Clears 2nd Derivative of the Curve
    R3.resize(3);
    R3.clear();        // Clears 3rd Derivative of the Curve

    Matrix& shape_derivatives = GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);      

    Vector3 coords;     // Coordinates of the Node 

    for (size_t i = 0; i < NumberOfNodes(); i++){
        const Node<3>& node = GetGeometry()[i];
        
        R1 += shape_derivatives(0,i) * node.GetInitialPosition(); 
        R2 += shape_derivatives(1,i) * node.GetInitialPosition(); 
        R3 += shape_derivatives(2,i) * node.GetInitialPosition(); 
    }

    A = norm_2(R1);       // Length of base Vector A

    double tmp = inner_prod(R2, R2) - std::pow(inner_prod(R1, R2), 2) / std::pow(A, 2);

    // Set Tollerance if not done bevore
    float tol = 1.0e-8;
    
    if(fabs(tmp) > tol)
        B = std::sqrt(tmp);
    else 
        B = 0;      
}

// //#################################################################################
// //#################################################################################
// //#
// //#                      +++ Compute_Geometry_Actual +++
// //#
// //#################################################################################
// //#################################################################################
// //#
// /** Computes the base Vector
//      * 
//      * @return
//      * 
//      * @param[in]   _u_act          NURBS coordinates of integration pont u-direction
//      * @param[in]   _deriv1         1st derivatives of the Basis functions at u
//      * @param[in]   _deriv2         2nd derivatives of the Basis functions at u
//      * @param[in]   _deriv3         3rd derivatives of the Basis functions at u
//      * @param[out]  r1             1st Derivative of the Curve in Referens Configuration
//      * @param[out]  r2             2nd Derivative ot the Curve in Referens Configuration
//      * @param[out]  r3             3rd Derivative ot the Curve in Referens Configuration
//      * @param[out]  a_ref          Length of the tangent Vector
//      * @param[out]  b_ref          
//      * 
//      * @author L.Rauch (10/2018)
//      * 
//      * @note by A.Bauer (10/2014)
//      */
// //#--------------------------------------------------------------------------------

void IgaBeamElement::ComputeGeometryActual(
        Vector3& r1,
        Vector3& r2,
        double& a,
        double& b)
{
    r1.resize(3);
    r1.clear();         // Clear the 1st Derivative of the Curve
    r2.resize(3);
    r2.clear();         // Clear the 2nd Derivative of the Curve

    // Compute the Basis Functions
    Matrix& shape_derivatives = GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);      
  
    // Compute Reference Geometry at each Node
    for (size_t i = 0; i < NumberOfNodes(); i++){
        const Node<3>& node = GetGeometry()[i];
        
        r1 += shape_derivatives(0,i) * node.Coordinates();
        r2 += shape_derivatives(1,i) * node.Coordinates();
    }

    a = norm_2(r1);     // Legth of the Base Vector

    double tmp = inner_prod(r2, r2) - std::pow(inner_prod(r1, r2), 2) / std::pow(a, 2);

    // Set Tollerance if not done bevore
    float tol = 1.0e-8;

    // Bending
    if (fabs(tmp) > tol)
    {
        b = std::sqrt(tmp);
    }else {
        b = 0;
    }
}

//#################################################################################
//#################################################################################
//#
//#                      +++ Compute_Geometry_Actual +++
//#
//#################################################################################
//#################################################################################
//#
/** Computes the base Vector
     * 
     * @return
     * 
     * @param[in]   _u_act          NURBS coordinates of integration pont u-direction
     * @param[in]   _deriv1         1st derivatives of the Basis functions at u
     * @param[in]   _deriv2         2nd derivatives of the Basis functions at u
     * @param[in]   _deriv3         3rd derivatives of the Basis functions at u
     * @param[out]  r1              1st Derivative of the Curve in Referens Configuration
     * @param[out]  r2              2nd Derivative ot the Curve in Referens Configuration
     * @param[out]  r3              3rd Derivative ot the Curve in Referens Configuration
     * @param[out]  a               Length of the tangent Vector
     * @param[out]  b          
     * 
     * @author L.Rauch (10/2018)
     * 
     * @note by A.Bauer (02/2018)
     */
//#--------------------------------------------------------------------------------

void IgaBeamElement::ComputeGeometryActual(
        Vector& r1,
        Vector& r2,
        Vector& r3,
        double& a,
        double& b)
{    
    r1.resize(3);
    r1.clear();        // Clears 1st Derivative of the Curve
    r2.resize(3);
    r2.clear();        // Clears 2nd Derivative of the Curve
    r3.resize(3);
    r3.clear();        // Clears 3rd Derivative of the Curve

    // Compute the Basis Functions
    Matrix& shape_derivatives = GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);      

    for (size_t i = 0; i < NumberOfNodes(); i++){
        const Node<3>& node = GetGeometry()[i];

        r1 += shape_derivatives(0,i) * node.Coordinates(); 
        r2 += shape_derivatives(1,i) * node.Coordinates();
        r3 += shape_derivatives(2,i) * node.Coordinates(); 
    }

    a = norm_2(r1);       // Length of the Base Vector in Deformed State

    double tmp = inner_prod(r2, r2) - pow(inner_prod(r1, r2), 2) / pow(a, 2);

    // Set Tollerance if not done bevore
    float tol = 1.0e-8;

    // Bendnding
    if (fabs(tmp) > tol)
        b = std::sqrt(tmp);
    else
        b = 0;
}

//#################################################################################
//#################################################################################
//#
//#              +++ Compute Cross Section Geometry Reference +++
//#
//#################################################################################
//#################################################################################
//#
/** Computes the base vector defined by the property. Is added to input!
     * 
     * @param[in]   _u_act        NURBS coordinates of integration pont u-direction
     * @param[out]  R1            1st Derivative of the Curve
     * @param[out]  R2            2nd Derivative of the Curve
     * @param[out]  R3            3rd Derivative of the Curve
     * @param[out]  A             Length of the Tangent Vector
     * 
     * @author L.Rauch (10/2018)
     * 
     * @note   A.Bauer (10/2014)
     */
//#--------------------------------------------------------------------------------
void IgaBeamElement::ComputeCrossSectionGeometryReference(
        BoundedVector<double,3>& R1,
        BoundedVector<double,3>& R2,
        BoundedVector<double,3>& T0_vec,
        Vector3 n_act,
        Vector3 v_act,
        Vector3 n0,
        Vector3 v0,
        Vector3 B_n,
        double C_12,
        double C_13,
        double Phi,
        double Phi_0_der)
{
    //# TODO ### ComputePhiReferencePRoperty implementirern
    // ComputePhiReferenceProperty(_Phi, _Phi_0_der);

    BoundedMatrix<double,3,3> mat_lamb;
    BoundedMatrix<double,3,3> mat_lamb_deriv;
    BoundedMatrix<double,3,3> mat_rod;
    BoundedMatrix<double,3,3> mat_rod_deriv;
    BoundedMatrix<double,3,3> mat_Ax1;
    mat_Ax1.clear();

    double R1_dl = norm_2(R1);

    BoundedVector<double,3> T_deriv;
    BoundedVector<double,3> T0_deriv;
    T0_deriv.clear(); 

    T_deriv = R2/R1_dl - inner_prod(R1,R2) / std::pow(R1_dl,3) * R1; 

    BoundedVector<double,3> _T_vec = R1/R1_dl; 

    //# TODO ### Continue with function [ mal_lambda]



}

// //#################################################################################
// //#################################################################################
// //#
// //#                          +++ Get Element Results +++
// //#
// //#################################################################################
// //#################################################################################
// //#
// /** Computes the Element results Stored at the Nodes
//      * 
//      * @param[in]   _typ          Nodal Result Tupe
//      * @param[return]   results   returns the Element results
//      * 
//      * @author L.Rauch (10/2018)
//      * 
//      * @note   A.Bauer (10/2014)
//      * @note   Nicht verwendete Funtion? 
//      */
// //#--------------------------------------------------------------------------------
// Vector IgaBeamElement::GetElementResult(
//     node_result_type _type)
// {
//     vector<double> results;
//     return results;
// }


//#################################################################################
//#################################################################################
//#
//#                          +++ Compute Epsilon DOF +++
//#
//#################################################################################
//#################################################################################
//#
/** Computes the Variation of the Axial Strain
     * 
     * @param[in]       _r1                         1st Derivative of the Curvature
     * @param[in]       _shape_function_deriv       1st Derivative of the Basis Functions
     * 
     * 
     * @param[return]   Vector with the Variation of the Normal Strain
     * 
     * @author L.Rauch (10/2018)
     * 
     * @note   A.Bauer (10/2014)
     */
//#--------------------------------------------------------------------------------
Vector IgaBeamElement::CoumputeEplsilonDof(Vector& _r1)
{
    int n_dof = NumberOfDofs() ; //# TODO ### n_dof = NumberOfDofs * P_Deg ### P_Deg = Nurbs_Curve->get_P_Deg();
    
    Vector epsilon_var;
    epsilon_var.resize(n_dof);
    double r1_L2 = norm_2(_r1);

    Matrix& shape_derivatives = GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);

    for (int i = 0; i < n_dof; i++)
    {
        int xyz_r = i%NumberOfDofs();
        int r = i/NumberOfDofs();

        if  (xyz_r > 2)
            epsilon_var[i] = 0;
        else
            epsilon_var[i] = _r1[xyz_r] * shape_derivatives(0,r);
    }
    return epsilon_var;
}

//#################################################################################
//#################################################################################
//#
//#                  +++ Compute Phi Reference Property +++
//#
//#################################################################################
//#################################################################################
//#
/** 
     * 
     * @param[in]   
     * @param[out]  Phi          // Phi Reference
     * @param[out]  Phi_0_der    // Phi Initial
     * 
     * @author L.Rauch (10/2018)
     * 
     * @note   A.Bauer (10/2014)
     */
//#--------------------------------------------------------------------------------
void IgaBeamElement::ComputePhiReferenceProperty(
        double Phi,
        double Phi_0_der)
{
    // Define Variables
    double phi_0;
    double phi_1;
    double diff_phi;
    double u_0;
    double u_1;
    int n_size;
    
    // Construct Vector of Temporarry dislplacements
    array_1d<double,4> tmp_ini_disp; 
    // Get Shape Function Derivatieves from the Element
    Vector& shape_function = GetValue(SHAPE_FUNCTION_VALUES);
    Matrix& shape_derivatives = GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);

    for (size_t i = 0; i < NumberOfNodes(); i++)
    {
        tmp_ini_disp = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT); 
      
        Phi += shape_function(i) * tmp_ini_disp[3];
        Phi_0_der += shape_derivatives(0,i) * tmp_ini_disp[3];
    }

}

//#################################################################################
//#################################################################################
//#
//#                             +++ Compute Phi DOF +++
//#
//#################################################################################
//#################################################################################
//#
/** Computes the Variation of the Rotation of the Curvature
     * 
     * @param[in]       _func         
     * @param[out]
     * 
     * @author L.Rauch (10/2018)
     * 
     * @note   A.Bauer (06/2014)
     */
//#--------------------------------------------------------------------------------
Vector IgaBeamElement::ComputePhieDof(Vector& _func)
{
    int n_dof = NumberOfDofs();  // * P_Deg; //# TODO ### P_Deg in Kratos
    Vector phi_var;     // Variation of the Axial Strain
    phi_var.resize(n_dof);  //# TODO ### Finktion [ N_Dof in Kratos ]
    phi_var.clear();

    for (int i = 0; i < n_dof; i++)
    {
        int xyz = i%NumberOfDofs(); 
        int r = i/NumberOfDofs(); 

        if (xyz < 2) 
            phi_var[i] = 0;
        else
            phi_var[i] = _func[r];
    }
    return phi_var;
}

//#################################################################################
//#################################################################################
//#
//#                        +++ Compute Matrix Lambda +++
//#
//#################################################################################
//#################################################################################
//#
/** Computes the Rotation Matrix Lambda (vec1 --> Vec2)
     * 
     * @param[in]   _matrix_lambda      (empty) Rotaition Matrix 
     * @param[in]   _vec1               initial Vector
     * @param[in]   _vec2               actual Vector
     * 
     * @param[return]   Matrix to Transform from vec1 to vec2
     * 
     * @author L.Rauch (10/2018)
     * 
     * @note   A.Bauer (11/2014)
     */
//#--------------------------------------------------------------------------------

void IgaBeamElement::ComputeMatrixLambda(
    BoundedMatrix<double,3,3>& _matrix_lambda,
    BoundedVector<double,3> _vec1,
    BoundedVector<double,3> _vec2)
{
    _matrix_lambda.clear();      // initialization by 0

    // Create Temporary Varables
    BoundedMatrix<double,3,3> matrix_lambda_tmp;
    matrix_lambda_tmp.clear();
    double tmp;
    BoundedMatrix<double,3,3> matrix_identity = IdentityMatrix(3);    // ZeroMatrix(3,3);

    // Compute CrossProduct between Vector 1 and Vector 2
    BoundedVector<double,3> vec1_x_vec2 = Cross(_vec1, _vec2);
    double l_vec1_x_vec2 = norm_2(vec1_x_vec2);
    BoundedVector<double,3> e_hat = vec1_x_vec2;

    double l_tol_crossprodunct = 1.0e-12;  // Tollerance for Crossproduct beeing numericaly Zero

    if (l_vec1_x_vec2 > l_tol_crossprodunct)
        e_hat = e_hat/l_vec1_x_vec2;
    
    double tol = 1.0e-8;        // Tollerance for Inner Product if not defined bevore
    if (inner_prod(_vec1, _vec2) > tol)
    {
        for (int i = 0; i < 3; i++) {_matrix_lambda(i,i) = inner_prod(_vec1, _vec2) ;}
        // Cross Product Vector - Matrix [hardcoded Implementation] /# TODO ### Mit Karat Function ersetzen
        // _matrix_lambda += CrossProductVectorMatrix(Cross(_vec1, _vec2), matrix_identity);
        for (int i = 0; i < 3; i++){
            for (int j = 0; j < 3; j++){
                matrix_lambda_tmp(i,j) = vec1_x_vec2(i) * vec1_x_vec2(j);
            }   // matrix_lambda_tmp = EinsteinSummation(vec1_x_vec2, vec1_vec2);
        }       

        tmp = 1.0 / (1.0 + inner_prod(_vec1, _vec2)); 
        matrix_lambda_tmp = matrix_lambda_tmp * tmp; 
        _matrix_lambda += matrix_lambda_tmp;
    }
    else
    {
        for (int i = 0; i < 3; i++){
            for (int j = 0; j < 3; j++){
                _matrix_lambda(i,j) = (e_hat(i) * e_hat(j)) * (1 - inner_prod(_vec1, _vec2)) + inner_prod(_vec1, _vec2) * matrix_identity(i,j); 
            }
        }
    }
    //# TODO function [ CroddProductVectorMatrix] implementrieren
    _matrix_lambda += CrossProductVectorMatrix(vec1_x_vec2, matrix_identity);
}

/**
 * Calculates the Cross Product between a Vector and a Matrix
 * @param[in]   vec     Input Vector [3,1]
 * @param[in]   mat     Input Matrix [3,3]
 * @return      vec_mat Output Matrix [3,3]
 */
BoundedMatrix<double,3,3> IgaBeamElement::CrossProductVectorMatrix(
    BoundedVector<double,3> vec,    
    BoundedMatrix<double,3,3> mat)
{
    BoundedMatrix<double,3,3> vec_mat;
    vec_mat.clear();

    int permutation[3][3][3];
    for (int i = 0; i < 3; i++){
        for (int j = 0; j < 3; j++){
            for (int k = 0; k < 3; k++){
                permutation[i][j][k] = 0; 
            }
        }
    }
    
    permutation[0][1][2] = 1;
    permutation[2][0][1] = 1;
    permutation[1][2][0] = 1;

    permutation[0][2][1] = -1;
    permutation[1][0][2] = -1;
    permutation[2][1][0] = -1; 

    for (int i = 0; i < 3; i++){
        for (int j = 0; j < 3; j++){
            for (int k = 0; k < 3; k++){
                for (int l = 0; l < 3; l++){
                    vec_mat(i,j) += permutation[i][k][l] * vec(k) * mat(l,j); 
                }
            }
        }
    }
    return vec_mat;
}




void IgaBeamElement::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "\"IgaBeamElement\" #" << Id();
}

//#################################################################################
//#################################################################################
//#
//#                          +++ Claculate Dead Load +++
//#
//#################################################################################
//#################################################################################
//#
/** Computes the dead Load Vecotor aof an Element
     * 
     * @param[in]   _dead_load      // Dead Load Object
     * @param[out]  lokal_load_vec  // Element Load Vector
     * 
     * @author L.Rauch (10/2018)
     * 
     * @note   A.Bauer (10/2014)
     */
//#--------------------------------------------------------------------------------
// void IgaBeamElement::CalculateDeadLoad(

// )


// }

} // namespace Kratos