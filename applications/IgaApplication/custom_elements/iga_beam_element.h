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
        Vector3 R1,
        Vector3 R2,
        Vector3& n_act,
        Vector3& v_act,
        Vector3 T0,
        Vector3& N0,
        Vector3& V0,
        double& B_n,
        double& B_v,
        double& C_12,
        double& C_13,
        double Phi,
        double Phi_der) ;

void IgaBeamElement::ComputeCrossSectionGeometryActual(
        Vector3 _R1,
        Vector3 _R2,
        Vector3 _r1,
        Vector3 _r2,
        Vector3 _N0,
        Vector3 _V0,
        Vector3& _n_act,
        Vector3& _v_act,
        double& _b_n,
        double& _b_v,
        double& _c_12,
        double& _c_13,
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

    Vector IgaBeamElement::CoumputeEplsilonDof(Vector& _r1) ;

    Vector IgaBeamElement::ComputePhieDof(Vector& _func) ;

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
    double _dl,
        MatrixType& _gke,
        VectorType& _gfie) ; 

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





// Adition Math_utilities
BoundedMatrix<double,3,3> CrossProductVectorMatrix(
    BoundedVector<double,3> vec,
    BoundedMatrix<double,3,3> mat) ;

Matrix IgaBeamElement::CrossProductVectorMatrix(
    Vector vec,
    Matrix mat) ;

BoundedMatrix<double,3,3> IgaBeamElement::CrossProductMatrixVector(
    BoundedMatrix<double,3,3> mat,
    BoundedVector<double,3> vec) ;
    
}; // class IgaBeamElement

} // namespace Kratos

#endif // !defined(KRATOS_IGA_BEAM_ELEMENT_H_INCLUDED)
