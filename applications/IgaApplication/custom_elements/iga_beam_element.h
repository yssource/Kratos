/*
//  KRATOS  _____________
//         /  _/ ____/   |
//         / // / __/ /| |
//       _/ // /_/ / ___ |
//      /___/\____/_/  |_| Application
//
//  Main authors:   Lukas Rauch
*/

#if !defined(KRATOS_IGA_BEAM_ELEMENT_H_INCLUDED)
#define KRATOS_IGA_BEAM_ELEMENT_H_INCLUDED

// System includes
#include "includes/define.h"
#include "includes/element.h"

// External includes

// Project includes
#include "iga_base_element.h"


namespace Kratos
{

class IgaBeamElement
    : public IgaBaseElement<4>
{
public:
    KRATOS_CLASS_POINTER_DEFINITION( IgaBeamElement );

    using IgaBaseElementType::IgaBaseElementType;

    // IgaBaseElement();
    // IgaBeamElement(
    //     int _ID, 
    //     Part* _Part_Prt,
    //     std::vector<Node*> _Node_Vec,
    //     MaterialBasis* _Mat_Prt,
    //     PropertyIgaBeamElement* _Prop_Prt,
    //     Nurbs1D* _Nurbs_Curve,
    //     int _knotspan_index, 
    //     std::vector<std::vector<double,2>> _phi_axis_n,
    //     std::vector<double,3> _t0,
    //     std::vector<EvalPtBasis*> _eval_points, 
    //     KnotSpan_Belonging _elementbelongig=INS, 
    //     PropertyNURBS_Curve_Stab* _prop_stab=0 )

    ~IgaBeamElement() override;

    Element::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties) const override;

    void GetDofList(
        DofsVectorType& rElementalDofList,
        ProcessInfo& rCurrentProcessInfo) override;

    void EquationIdVector(
        EquationIdVectorType& rResult,
        ProcessInfo& rCurrentProcessInfo) override;

    void Initialize() override;

    void CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo,
        const bool ComputeLeftHandSide,
        const bool ComputeRightHandSide) override;

    void PrintInfo(std::ostream& rOStream) const override;

    void IgaBeamElement::GetDofTypesPerNode(
        std::vector<int>& _act_dofs); 

    void IgaBeamElement::ComputeGeometryInitial(
            Vector3& r1,
            Vector3& r2,
            double& a_ini,
            double& b_ini) ;

    void IgaBeamElement::ComputeGeometryInitial(
            Vector3& r1,
            Vector3& r2,
            Vector3& r3,
            double& a_ini,
            double& b_ini) ;

    void IgaBeamElement::ComputeGeometryReference(
            Vector3& R1,
            Vector3& R2, 
            double& A, 
            double& B);

    void IgaBeamElement::ComputeGeometryReference(
            Vector3& R1,
            Vector3& R2,
            Vector3& R3,
            double& A,
            double& B) ;
                
    void IgaBeamElement::ComputeGeometryActual(
            Vector3& r1,
            Vector3& r2,
            double& a,
            double& b) ; 

    void IgaBeamElement::ComputeGeometryActual(
            Vector3& r1,
            Vector3& r2,
            Vector3& r3,
            double& a,
            double& b) ;

    void IgaBeamElement::ComputeCrossSectionGeometryReference(
        Vector3 _R1,
        Vector3 _R2,
        Vector3& A2,
        Vector3& A3,
        Vector3 A01,
        Vector3& A02,
        Vector3& A03,
        double& B_2,
        double& B_3,
        double& C_12,
        double& C_13,
        double Phi,
        double Phi_der ) ;

    void IgaBeamElement::ComputeCrossSectionGeometryActual(
        Vector3 A1,
        Vector3 A1_1,
        Vector3 a1,
        Vector3 a1_1,
        Vector3 A02,
        Vector3 A03,
        Vector3& a2,
        Vector3& a3,
        double& b_2,
        double& b_3,
        double& c_12,
        double& c_13,
        double Phi,
        double Phi_der,
        double phi,
        double phi_der) ;

    void IgaBeamElement::ComputePhiReferenceProperty(
        double& phi,
        double& phi_0_der) ;

    void IgaBeamElement::GetDofsPerNode(
        std::vector<int>& _act_dofs) ;

    // Vector IgaBeamElement::GetElementResult(
    //     node_result_type _type); 

    Vector IgaBeamElement::ComputeEpsilonFirstDerivative(Vector3 _r1) ;

    Matrix IgaBeamElement::ComputeEpsilonSecondDerivative();        //change 17.11   (Vector3 _r1) ;

    void IgaBeamElement::ComputeRotationalDof(
    VectorType& _curv_var_t,
    VectorType& _curv_var_n,
    VectorType& _curv_var_v,
    Vector3 _R1,
    Vector3 _r1,
    Vector3 _N0,
    Vector3 _V0,
    Vector3 _N_ref,
    Vector3 _V_ref,
    Vector3 _n_act,
    Vector3 _v_act,
    double _Phi,
    double _phi ) ;

    // Vector IgaBeamElement::ComputePhiDof(Vector& _func) ;        //changes 17.11     funktion nicht benötigt

    void IgaBeamElement::ComputeMatrixLambda(
        BoundedMatrix<double,3,3>& _matrix_lambda,
        BoundedVector<double,3> _vec1,
        BoundedVector<double,3> _vec2) ;

    void IgaBeamElement::ComputeMatrixLambdaFirstDerivative(
        BoundedMatrix<double,3,3>& _matrix_lambda_der,
        BoundedVector<double,3> _vec1,
        BoundedVector<double,3> _vec2,
        BoundedVector<double,3> _vec1_der,
        BoundedVector<double,3> _vec2_der) ;

    void IgaBeamElement::ComputeMatrixLambdaSecondDerivative(
        BoundedMatrix<double,3,3>& _matrix_lambda_derder,
        BoundedVector<double,3> _vec1,
        BoundedVector<double,3> _vec2,
        BoundedVector<double,3> _vec1_der1,
        BoundedVector<double,3> _vec2_der1,
        BoundedVector<double,3> _vec1_der2,
        BoundedVector<double,3> _vec2_der2) ;

    void IgaBeamElement::ComputeMatrixLambdaVariation(
        Matrix& _matrix_lambda_var,
        BoundedVector<double,3> _vec1,
        BoundedVector<double,3> _vec2,
        BoundedVector<double,3> _vec2_var) ;

    void IgaBeamElement::ComputeMatrixLambdaFirstDerivativeVariation(
        Matrix _matrix_lambda_der1var,
        BoundedVector<double,3> _vec1,
        BoundedVector<double,3> _vec2,
        BoundedVector<double,3> _vec1der1,
        BoundedVector<double,3> _vec2der1,
        Vector _vec2var,
        Vector _vec2der1var) ;

    void IgaBeamElement::ComputeMatrixLambdaSecondDerivativeVariation(
        Matrix& _matrix_lambda_der2var,
        BoundedVector<double,3> _vec1,
        BoundedVector<double,3> _vec2,
        BoundedVector<double,3> _vec1der1,
        BoundedVector<double,3> _vec2der1,
        BoundedVector<double,3> _vec1der2,
        BoundedVector<double,3> _vec2der2,
        Vector _vec2var,
        Vector _vec2der1var,
        Vector _vec2der2var) ; 

    void IgaBeamElement::ComputeMatrixLambdaSecondVariation(
        Matrix& _matrix_lambda_second_variation,
        BoundedVector<double,3> _vec1,
        BoundedVector<double,3> _vec2,
        Vector _vec2var1,
        Matrix _vec2var2) ;

    void IgaBeamElement::ComputeMatrixLambdaFirstDerivativeSecondVariation(
        Matrix _matrix_lambda_der1var2,
        BoundedVector<double,3> _vec1,
        BoundedVector<double,3> _vec2,
        BoundedVector<double,3> _vec1der1,
        BoundedVector<double,3> _vec2der1,
        Vector _vec2var1,
        Matrix _vec2var2,
        Vector _vec2der1var1,
        Matrix _vec2der1var2) ;

    void IgaBeamElement::ComputeMatrixLambdaAll(
        Matrix& _matrix_lambda_var1,
        Matrix& _matrix_lambda_var2,
        Matrix& _matrix_lambda_der1_var1,
        Matrix& _matrix_lambda_der1_var2,
        Vector3 _vec1,
        Vector3 _vec2,
        Vector3 _vec1_der1,
        Vector3 _vec2_der1,
        Vector  _vec2_var1,
        Vector  _vec2_der1_var1,
        Matrix  _vec2_var2,
        Matrix  _vec2_der1_var2) ;

    void IgaBeamElement::StffnessMatrixElementLinear(
        double _emod,
        double _gmod,
        double _area,
        double _m_inert_n,
        double _m_inert_v,
        double _mt_inert,
        double _dl,
            MatrixType& gke) ; 

    void IgaBeamElement::ElementStiffnessMatrixNonlinear(
    double _emod,
    double _gmod,
    double _area,
    double _m_inert_y,
    double _m_inert_z,
    double _mt_inert,
        MatrixType& _gke,
        VectorType& _gfie,
        double& _dL) ; 

    void IgaBeamElement::ComputeMatrixRodriuez(
        BoundedMatrix<double,3,3>& _matrix_rodrigues,
        BoundedVector<double,3> _vec,
        double phi) ;

    void IgaBeamElement::ComputeMatrixRodriuezFirstDerivative(
        BoundedMatrix<double,3,3>& _matrix_rodrigues_der,
        BoundedVector<double,3> _vec,
        BoundedVector<double,3> _vec_der,
        double phi,
        double phi_der) ;

    void IgaBeamElement::ComputeTVariation(
        Vector& _t_var,
        Vector  _r1,
        Matrix  _deriv);

    void IgaBeamElement::ComputeDofNonlinear(
    Vector& _curve_var_n,
    Vector& _curve_var_v,
    Vector& _tor_var_n,
    Vector& _tor_var_v,
    Matrix& _curve_var_var_n,
    Matrix& _curve_var_var_v,
    Matrix& _tor_var_var_n,
    Matrix& _tor_var_var_v,
    Vector3 _R1,
    Vector3 _R2,
    Vector3 _r1,
    Vector3 _r2,
    Vector3 _N0,
    Vector3 _V0,
    double  _Phi,
    double  _Phi_der,
    double  _phi,
    double  _phi_der) ;


    void ComputeTVar(
        const Vector3& r1,
        std::vector<Vector3>& t_var);

    void ComputeTVarVar(
        const Vector3& r1,
        std::vector<Vector3>& t_var_var);

    void ComputeTDerVar(
        const Vector3& r1,
        const Vector3& r2,
        std::vector<Vector3>& t_der_var);

    void ComputeTDerVarVar(
        const Vector3& r1,
        const Vector3& r2,
        std::vector<Vector3>& t_der_var_var);

    void ComputeRodrigues(
        const Vector3& v,
        const double& phi,
        BoundedMatrix<double, 3, 3>& rodrigues);

    void ComputeRodriguesDer(
        const Vector3& v,
        const Vector3& v_der,
        const double& phi,
        const double& phi_der,
        BoundedMatrix<double, 3, 3>& rodrigues_der);

    void ComputeRodriguesDerVar(
        const Vector3& v,
        const std::vector<Vector3>& v_var,
        const Vector3& v_der,
        const std::vector<Vector3>& v_der_var,
        const double& phi,
        const double& phi_der,
        std::vector<BoundedMatrix<double, 3, 3>>& rodrigues_der_var);

    void ComputeRodriguesDerVarVar(
        const Vector3& v, 
        const std::vector<Vector3>& v_var, 
        const std::vector<Vector3>& v_var_var, 
        const Vector3& v_der, 
        const std::vector<Vector3>& v_der_var, 
        const std::vector<Vector3>& v_der_var_var, 
        const double& phi, 
        const double& phi_der, 
        std::vector<BoundedMatrix<double, 3, 3>>& rodrigues_der_var_var);

    void ComputeRodriguesVar(
        const Vector3& v,
        const std::vector<Vector3>& v_var, 
        const double& phi, 
        std::vector<BoundedMatrix<double, 3, 3>>& rodrigues_var);

    void ComputeRodriguesVarVar(
        const Vector3& v, 
        const std::vector<Vector3>& v_var, 
        const std::vector<Vector3>& v_var_var, 
        const double& phi, 
        std::vector<BoundedMatrix<double, 3, 3>>& rodrigues_var_var);

    void ComputeLambda(
        const Vector3& v1,
        const Vector3& v2,
        BoundedMatrix<double, 3, 3>& lambda);

    void ComputeLambdaVar(
        const Vector3& v1,
        const Vector3& v2,
        const std::vector<Vector3>& v2_var,
        std::vector<BoundedMatrix<double, 3, 3>>& lambda_var);

    void ComputeLambdaVarVar(
        const Vector3& v1,
        const Vector3& v2,
        const std::vector<Vector3>& v2_var,
        const std::vector<Vector3>& v2_var_var,
        std::vector<BoundedMatrix<double, 3, 3>>& lambda_var_var);

    void ComputeLambdaDer(
        const Vector3& v1,
        const Vector3& v1_der,
        const Vector3& v2,
        const Vector3& v2_der,
        BoundedMatrix<double, 3, 3>& lambda_der);

    void ComputeLambdaDerVar(
        const Vector3& v1,
        const Vector3& v1_der,
        const Vector3& v2,
        const std::vector<Vector3>& v2_var,
        const Vector3& v2_der,
        const std::vector<Vector3>& v2_der_var,
        std::vector<BoundedMatrix<double, 3, 3>>& lambda_der_var);

    void ComputeLambdaDerVarVar(
        const Vector3& v1,
        const Vector3& v1_der,
        const Vector3& v2,
        const std::vector<Vector3>& v2_var,
        const std::vector<Vector3>& v2_var_var,
        const Vector3& v2_der,
        const std::vector<Vector3>& v2_der_var,
        const std::vector<Vector3>& v2_der_var_var,
        std::vector<BoundedMatrix<double, 3, 3>>& lambda_der_var_var);


// Adition Math_utilities
BoundedMatrix<double,3,3> CrossProductVectorMatrix(
    BoundedVector<double,3> vec,
    BoundedMatrix<double,3,3> mat) ;

BoundedMatrix<double, 3, 3> CrossVectorIdentity(
    const Vector3& v);

Matrix IgaBeamElement::CrossProductVectorMatrix(
    Vector vec,
    Matrix mat) ;

BoundedMatrix<double,3,3> IgaBeamElement::CrossProductMatrixVector(
    BoundedMatrix<double,3,3> mat,
    BoundedVector<double,3> vec) ;



void IgaBeamElement::CalculateLoadMoment( ) ;

void IgaBeamElement::ComputeStressNonlinear( ) ;

void IgaBeamElement::ComputeShearForcesNonlinear( 
    Vector shear_forces_n, 
    Vector shear_forces_v, 
        Vector3 _N0,
        Vector3 _V0,
        Vector3 _R1,
        Vector3 _R2,
        Vector3 _R3,
        Vector3 _r1,
        Vector3 _r2,
        Vector3 _r3, 
        double  _Phi,
        double  _Phi_der1,
        double  _Phi_der2,
        double  _phi_der,
        double  _phi_der1,
        double  _phi_der2,
        bool _prestress_bend1_auto,
        bool _prestress_bend2_auto ) ;

void IgaBeamElement::ComputeInitialCrossRotation(Vector3& A1, Vector3& A2, Vector3& A3, double Theta, double dL) ;
    
}; // class IgaBeamElement

} // namespace Kratos

#endif // !defined(KRATOS_IGA_BEAM_ELEMENT_H_INCLUDED)
