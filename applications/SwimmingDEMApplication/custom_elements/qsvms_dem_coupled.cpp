//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//

// Project includes
#include "swimming_DEM_application.h"
#include "qsvms_dem_coupled.h"
#include "includes/cfd_variables.h"
#include "includes/checks.h"

#include "custom_utilities/qsvms_data.h"
//#include "../FluidDynamicsApplication/custom_utilities/fluid_element_data.h"
//#include "../FluidDynamicsApplication/custom_utilities/qsvms_data.h"
//#include "../FluidDynamicsApplication/custom_elements/qs_vms.h"
//#include "../FluidDynamicsApplication/custom_elements/fluid_element.h"
#include "../FluidDynamicsApplication/custom_utilities/fluid_element_utilities.h"

namespace Kratos
{

//////////////////////////Life cycle

template< class TElementData >
QSVMSDEMCoupled<TElementData>::QSVMSDEMCoupled(IndexType NewId):
    QSVMS<TElementData>(NewId)
{}

template< class TElementData >
QSVMSDEMCoupled<TElementData>::QSVMSDEMCoupled(IndexType NewId, const NodesArrayType& ThisNodes):
    QSVMS<TElementData>(NewId,ThisNodes)
{}


template< class TElementData >
QSVMSDEMCoupled<TElementData>::QSVMSDEMCoupled(IndexType NewId, GeometryType::Pointer pGeometry):
    QSVMS<TElementData>(NewId,pGeometry)
{}


template< class TElementData >
QSVMSDEMCoupled<TElementData>::QSVMSDEMCoupled(IndexType NewId, GeometryType::Pointer pGeometry, Properties::Pointer pProperties):
     QSVMS<TElementData>(NewId,pGeometry,pProperties)
{}

///////////Destructor

template< class TElementData >
QSVMSDEMCoupled<TElementData>::~QSVMSDEMCoupled()
{}

template< class TElementData >
Element::Pointer QSVMSDEMCoupled<TElementData>::Create(IndexType NewId,NodesArrayType const& ThisNodes,Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<QSVMSDEMCoupled>(NewId, this->GetGeometry().Create(ThisNodes), pProperties);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @see QSVMSDEMCoupled::EquationIdVector
 **/
template < class TElementData >
void QSVMSDEMCoupled<TElementData>::EquationIdVector(EquationIdVectorType& rResult,
                                         ProcessInfo& rCurrentProcessInfo)
{
    if (rCurrentProcessInfo[FRACTIONAL_STEP] == 1) {

        unsigned int LocalIndex = 0;

        if (rResult.size() != LocalSize)
            rResult.resize(LocalSize, false);

        unsigned int vpos = this->GetGeometry()[0].GetDofPosition(VELOCITY_X);
        unsigned int ppos = this->GetGeometry()[0].GetDofPosition(PRESSURE);

        for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
        {
            rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_X,vpos).EquationId();
            rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_Y,vpos+1).EquationId();
            rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(PRESSURE,ppos).EquationId();
        }
    }

    else {

        unsigned int LocalIndex = 0;

        unsigned int lappos = this->GetGeometry()[0].GetDofPosition(VELOCITY_LAPLACIAN_X);

        if (rResult.size() != LocalSize)
            rResult.resize(LocalSize, false);

        for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
        {
            rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_LAPLACIAN_X,lappos).EquationId();
            rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_LAPLACIAN_Y,lappos+1).EquationId();
        }
    }

}

/**
 * @see QSVMSDEMCoupled::GetDofList
 */
template < class TElementData >
void QSVMSDEMCoupled<TElementData>::GetDofList(DofsVectorType& rElementalDofList,
                        ProcessInfo& rCurrentProcessInfo)
{
    if (rCurrentProcessInfo[FRACTIONAL_STEP] == 1) {

        if (rElementalDofList.size() != LocalSize)
            rElementalDofList.resize(LocalSize);

        unsigned int LocalIndex = 0;

        for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
        {
            rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_X);
            rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_Y);
            rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(PRESSURE);
        }
    }

    else {

        if (rElementalDofList.size() != LocalSize)
            rElementalDofList.resize(LocalSize);

        unsigned int LocalIndex = 0;

        for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
        {
            rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_LAPLACIAN_X);
            rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_LAPLACIAN_Y);
        }
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Inquiry

template< class TElementData >
int QSVMSDEMCoupled<TElementData>::Check(const ProcessInfo &rCurrentProcessInfo)
{
    int out = QSVMS<TElementData>::Check(rCurrentProcessInfo);
    KRATOS_ERROR_IF_NOT(out == 0)
        << "Error in base class Check for Element " << this->Info() << std::endl
        << "Error code is " << out << std::endl;

    // Extra variables
    KRATOS_CHECK_VARIABLE_KEY(ACCELERATION);
    KRATOS_CHECK_VARIABLE_KEY(NODAL_AREA);

    // Output variables (for Calculate() functions)
    KRATOS_CHECK_VARIABLE_KEY(SUBSCALE_VELOCITY);
    KRATOS_CHECK_VARIABLE_KEY(SUBSCALE_PRESSURE);

    for(unsigned int i=0; i<NumNodes; ++i)
    {
        Node<3>& rNode = this->GetGeometry()[i];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ACCELERATION,rNode);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(NODAL_AREA,rNode);
    }

    return out;
}
///////////////////////////////////////////////////////////////////////////////////////////////////

template< class TElementData >
void QSVMSDEMCoupled<TElementData>::GetValueOnIntegrationPoints(Variable<array_1d<double, 3 > > const& rVariable,
                                                                std::vector<array_1d<double, 3 > >& rOutput,
                                                                ProcessInfo const& rCurrentProcessInfo)
{

    if (rVariable == VORTICITY)
    {
        // Set output vector (for a single integration point)
        rOutput.resize(1);
        array_1d<double, 3 > & rVorticity = rOutput[0];
        rVorticity[0] = 0.0;
        rVorticity[1] = 0.0;
        rVorticity[2] = 0.0;

        double Area;
        array_1d<double, NumNodes> N;
        BoundedMatrix<double, NumNodes, Dim> DN_DX;
        GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DX, N, Area);

        for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
        {
            const array_1d<double, 3 > & rVelocity = this->GetGeometry()[iNode].FastGetSolutionStepValue(VELOCITY);
            rVorticity[2] += DN_DX(iNode,0)*rVelocity[1] - DN_DX(iNode,1)*rVelocity[0];
        }
    }
    else if (rVariable == SUBSCALE_VELOCITY)
    {
        if(0) //this->GetValue(TRACK_SUBSCALES) == 1 )
        {
            rOutput.resize(1);
            const QSVMSDEMCoupled<TElementData>* const_this = static_cast< QSVMSDEMCoupled<TElementData>* >(this);
            rOutput[0] = const_this->GetValue(rVariable);
        }
        else
        {
            double Area;
            array_1d<double, NumNodes> N;
            BoundedMatrix<double, NumNodes, Dim> DN_DX;
            GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DX, N, Area);

            array_1d<double,3> AdvVel;
            this->GetAdvectiveVel(AdvVel,N);

            double Density,KinViscosity,Viscosity;
            this->EvaluateInPoint(Density,DENSITY,N);
            this->EvaluateInPoint(KinViscosity,VISCOSITY,N);
            this->GetEffectiveViscosity(Density,KinViscosity,N,DN_DX,Viscosity,rCurrentProcessInfo);

            double TauOne,TauTwo;
            this->CalculateTau(TauOne,TauTwo,AdvVel,Area,Density,Viscosity,rCurrentProcessInfo);

            // Set output vector (for a single integration point)
            rOutput.resize(1);
            array_1d<double,3> MomError(3,0.0);
            if (rCurrentProcessInfo[OSS_SWITCH]==1)
            {
                this->OSSMomResidual(AdvVel,Density,MomError,N,DN_DX,1.0);
            }
            else
            {
                this->ASGSMomResidual(AdvVel,Density,MomError,N,DN_DX,1.0);
            }
            MomError *= TauOne;
            array_1d<double,3>& rSubscale = rOutput[0];
            rSubscale[0] = MomError[0];
            rSubscale[1] = MomError[1];
            rSubscale[2] = 0.0;
        }
    }
    else // Default behaviour (returns elemental data)
    {
        rOutput.resize(1);
        /*
         The cast is done to avoid modification of the element's data. Data modification
         would happen if rVariable is not stored now (would initialize a pointer to &rVariable
         with associated value of 0.0). This is catastrophic if the variable referenced
         goes out of scope.
         */
        const QSVMSDEMCoupled<TElementData>* const_this = static_cast< QSVMSDEMCoupled<TElementData>* >(this);
        rOutput[0] = const_this->GetValue(rVariable);
    }
}

template < class TElementData >
void QSVMSDEMCoupled< TElementData >::CalculateWeights(ShapeFunctionDerivativesArrayType& rDN_DX,
                                                       Matrix& rNContainer,
                                                       Vector& rGaussWeights)
{

  const GeometryType& rGeom = this->GetGeometry();
  Vector DetJ;
  rGeom.ShapeFunctionsIntegrationPointsGradients(rDN_DX, DetJ, GeometryData::GI_GAUSS_2);
  rNContainer = rGeom.ShapeFunctionsValues(GeometryData::GI_GAUSS_2);
  const GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints(GeometryData::GI_GAUSS_2);

  rGaussWeights.resize(rGeom.IntegrationPointsNumber(GeometryData::GI_GAUSS_2), false);

  for (unsigned int g = 0; g < rGeom.IntegrationPointsNumber(GeometryData::GI_GAUSS_2); g++)
      rGaussWeights[g] = DetJ[g] * IntegrationPoints[g].Weight();
}

/**
 * @see QSVMSDEMCoupled::ElementSize
 */
template < class TElementData >
double QSVMSDEMCoupled< TElementData >::ElementSize(const double Variable)
{
    if (Dim == 2){
        return 1.128379167 * sqrt(Variable); //Diameter of circumference of given Area
    } else{
        return 0.60046878 * pow(Variable, 0.333333333333333333333);
    }
}

/**
 * Returns the squared element size, estimated as h^2 = 2*Area
 * Returns the squared element size, estimated from the assumption V = (1/6) * h^3
 * @see QSVMSDEMCoupled::FilterWidth
 */
template <class TElementData>
double QSVMSDEMCoupled<TElementData>::FilterWidth()
{
    if (Dim == 2){
        double FilterWidth = GeometryUtils::CalculateVolume2D(this->GetGeometry());
        return 2.0 * FilterWidth;
    } else {
        const double TwoThirds = 2.0 / 3.0;
        double FilterWidth = GeometryUtils::CalculateVolume3D(this->GetGeometry());
        FilterWidth *= 6.0;
        return pow(FilterWidth, TwoThirds);
    }
}
/**
 * Returns the square of the minimum element height, to be used as filter width in the Smagorinsky model
 * @see QSVMSDEMCoupled::FilterWidth
 */
template <class TElementData>
double QSVMSDEMCoupled<TElementData>::FilterWidth(const BoundedMatrix<double, NumNodes, Dim>& DN_DX)
{
    double inv_h_max = 0.0;
    for(unsigned int i=0; i<NumNodes; i++)
    {
        double inv_h = 0.0;
        for(unsigned int d=0; d<Dim; d++)
            inv_h += DN_DX(i,d)*DN_DX(i,d);

        if(inv_h > inv_h_max) inv_h_max = inv_h;
    }

    double DeltaSquared = 1.0/inv_h_max;

    return DeltaSquared ;
}
/**
 * See QSVMSDEMCoupled::CalculateB
 */
template < class TElementData >
void QSVMSDEMCoupled<TElementData>::CalculateB(BoundedMatrix<double, (Dim * NumNodes) / 2, Dim * NumNodes >& rB,
                                               const BoundedMatrix<double, NumNodes, Dim >& rShapeDeriv)
{
    KRATOS_TRY

    if (Dim == 2){
        for (unsigned int i = 0; i < NumNodes; i++)
        {
            unsigned int index = 2 * i;

            rB(0, index) = rShapeDeriv(i, 0);
            rB(0, index + 1) = 0.0;
            rB(1, index) = 0.0;
            rB(1, index + 1) = rShapeDeriv(i, 1);
            rB(2, index) = rShapeDeriv(i, 1);
            rB(2, index + 1) = rShapeDeriv(i, 0);
        }
    } else if (Dim == 3){

        for (unsigned int i = 0; i < NumNodes; ++i)
        {
            unsigned int index = Dim*i;

            rB(0, index) = rShapeDeriv(i, 0);
            rB(0, index + 1) = 0.0;
            rB(0, index + 2) = 0.0;
            rB(1, index) = 0.0;
            rB(1, index + 1) = rShapeDeriv(i, 1);
            rB(1, index + 2) =0.0;
            rB(2, index) = 0.0;
            rB(2, index + 1) = 0.0;
            rB(2, index + 2) = rShapeDeriv(i, 2);
            rB(3, index) = rShapeDeriv(i, 1);
            rB(3, index + 1) = rShapeDeriv(i, 0);
            rB(3, index + 2) = 0.0;
            rB(4, index) = 0.0;
            rB(4, index + 1) = rShapeDeriv(i, 2);
            rB(4, index + 2) = rShapeDeriv(i, 1);
            rB(5, index) = rShapeDeriv(i, 2);
            rB(5, index + 1) = 0.0;
            rB(5, index + 2) = rShapeDeriv(i, 0);
        }
    }
    KRATOS_CATCH("")
}

/**
 * See QSVMSDEMCoupled::CalculateC
 */
template < class TElementData >
void QSVMSDEMCoupled<TElementData>::CalculateC(BoundedMatrix<double, (Dim * NumNodes)/2, (Dim*NumNodes)/2> & rC,
                                               const double Viscosity)
{
    if (Dim == 2){
        noalias(rC) = ZeroMatrix(3,3);
    }else{
        noalias(rC) = ZeroMatrix(6,6);
    }
    rC(0, 0) =  Viscosity*(1.3333333333333333333333333333333);
    rC(0, 1) = -Viscosity*(0.666666666666666666666666666667);
    rC(1, 0) = -Viscosity*(0.666666666666666666666666666667);
    rC(1, 1) =  Viscosity*(1.3333333333333333333333333333);
    if (Dim == 2){
        rC(0, 2) = 0.0;
        rC(1, 2) = 0.0;
        rC(2, 0) = 0.0;
        rC(2, 1) = 0.0;
        rC(2, 2) = Viscosity;
    }else{
        rC(0, 2) = -Viscosity*(0.666666666666666666666666666667);
        rC(1, 2) = -Viscosity*(0.666666666666666666666666666667);
        rC(2, 0) = -Viscosity*(0.666666666666666666666666666667);
        rC(2, 1) = -Viscosity*(0.666666666666666666666666666667);
        rC(2, 2) =  Viscosity*(1.3333333333333333333333333333);
        rC(3, 3) = Viscosity;
        rC(4, 4) = Viscosity;
        rC(5, 5) = Viscosity;
    }
}
/**
 * @see QSVMSDEMCoupled::AddViscousTerm
 */
template < class TElementData >
void QSVMSDEMCoupled<TElementData>::AddViscousTerm(MatrixType& rDampMatrix,
                                                   const BoundedMatrix<double,NumNodes,Dim>& rShapeDeriv,
                                                   const double Weight)
{
    const double OneThird = 1.0 / 3.0;
    const double FourThirds = 4.0 / 3.0;
    const double nTwoThirds = -2.0 / 3.0;

    unsigned int FirstRow(0),FirstCol(0);
    if (Dim == 2){
        for (unsigned int j = 0; j < NumNodes; ++j)
        {
            for (unsigned int i = 0; i < NumNodes; ++i)
            {
                // First Row
                rDampMatrix(FirstRow,FirstCol) += Weight * ( FourThirds * rShapeDeriv(i,0) * rShapeDeriv(j,0) + rShapeDeriv(i,1) * rShapeDeriv(j,1) );
                rDampMatrix(FirstRow,FirstCol+1) += Weight * ( nTwoThirds * rShapeDeriv(i,0) * rShapeDeriv(j,1) + rShapeDeriv(i,1) * rShapeDeriv(j,0) );

                // Second Row
                rDampMatrix(FirstRow+1,FirstCol) += Weight * ( nTwoThirds * rShapeDeriv(i,1) * rShapeDeriv(j,0) + rShapeDeriv(i,0) * rShapeDeriv(j,1) );
                rDampMatrix(FirstRow+1,FirstCol+1) += Weight * ( FourThirds * rShapeDeriv(i,1) * rShapeDeriv(j,1) + rShapeDeriv(i,0) * rShapeDeriv(j,0) );

                // Update Counter
                FirstRow += NumNodes;
        }
            FirstRow = 0;
            FirstCol += NumNodes;
    }
    } else{
        for (unsigned int j = 0; j < NumNodes; ++j)
        {
            for (unsigned int i = 0; i < NumNodes; ++i)
            {
                // (dN_i/dx_k dN_j/dx_k)
                const double Diag =  rShapeDeriv(i,0) * rShapeDeriv(j,0) + rShapeDeriv(i,1) * rShapeDeriv(j,1) + rShapeDeriv(i,2) * rShapeDeriv(j,2);

                // First Row
                rDampMatrix(FirstRow,FirstCol)   += Weight * ( OneThird * rShapeDeriv(i,0) * rShapeDeriv(j,0) + Diag );
                rDampMatrix(FirstRow,FirstCol+1) += Weight * ( nTwoThirds * rShapeDeriv(i,0) * rShapeDeriv(j,1) + rShapeDeriv(i,1) * rShapeDeriv(j,0) );
                rDampMatrix(FirstRow,FirstCol+2) += Weight * ( nTwoThirds * rShapeDeriv(i,0) * rShapeDeriv(j,2) + rShapeDeriv(i,2) * rShapeDeriv(j,0) );

                // Second Row
                rDampMatrix(FirstRow+1,FirstCol)   += Weight * ( nTwoThirds * rShapeDeriv(i,1) * rShapeDeriv(j,0) + rShapeDeriv(i,0) * rShapeDeriv(j,1) );
                rDampMatrix(FirstRow+1,FirstCol+1) += Weight * ( OneThird * rShapeDeriv(i,1) * rShapeDeriv(j,1) + Diag );
                rDampMatrix(FirstRow+1,FirstCol+2) += Weight * ( nTwoThirds * rShapeDeriv(i,1) * rShapeDeriv(j,2) + rShapeDeriv(i,2) * rShapeDeriv(j,1) );

                // Third Row
                rDampMatrix(FirstRow+2,FirstCol)   += Weight * ( nTwoThirds * rShapeDeriv(i,2) * rShapeDeriv(j,0) + rShapeDeriv(i,0) * rShapeDeriv(j,2) );
                rDampMatrix(FirstRow+2,FirstCol+1) += Weight * ( nTwoThirds * rShapeDeriv(i,2) * rShapeDeriv(j,1) + rShapeDeriv(i,1) * rShapeDeriv(j,2) );
                rDampMatrix(FirstRow+2,FirstCol+2) += Weight * ( OneThird * rShapeDeriv(i,2) * rShapeDeriv(j,2) + Diag );

                // Update Counter
                FirstRow += NumNodes;
            }
            FirstRow = 0;
            FirstCol += NumNodes;
        }
    }
}

template<class TElementData>
double QSVMSDEMCoupled<TElementData>::ConsistentMassCoef(const double Variable)
{
    if(Dim == 2){
        const double Coef = 1.0/12.0;
        return Variable * Coef;
    }else{
        const double Coef = 1.0/20.0;
        return Variable * Coef;
    }

}

template<class TElementData>
void QSVMSDEMCoupled<TElementData>::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                          VectorType& rRightHandSideVector,
                          ProcessInfo& rCurrentProcessInfo)
    {
        if (rCurrentProcessInfo[FRACTIONAL_STEP] == 1) {
            const unsigned int LocalSize = (Dim + 1) * NumNodes;

            // Check sizes and initialize matrix
            if (rLeftHandSideMatrix.size1() != LocalSize)
                rLeftHandSideMatrix.resize(LocalSize, LocalSize, false);

            noalias(rLeftHandSideMatrix) = ZeroMatrix(LocalSize, LocalSize);

            // Calculate RHS
            this->CalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo);
        }

        else {
            const unsigned int LocalSize = Dim * NumNodes;

            if (rLeftHandSideMatrix.size1() != LocalSize)
                rLeftHandSideMatrix.resize(LocalSize, LocalSize, false);

            noalias(rLeftHandSideMatrix) = ZeroMatrix(LocalSize, LocalSize);
            CalculateLaplacianMassMatrix(rLeftHandSideMatrix, rCurrentProcessInfo);
            noalias(rRightHandSideVector) = ZeroVector(LocalSize);
            this->CalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo);
        }
    }

template<class TElementData>
void QSVMSDEMCoupled<TElementData>::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                                                          ProcessInfo& rCurrentProcessInfo)
    {
        if (rCurrentProcessInfo[FRACTIONAL_STEP] == 1) {

            const unsigned int LocalSize = (Dim + 1) * NumNodes;

            if (rLeftHandSideMatrix.size1() != LocalSize) rLeftHandSideMatrix.resize(LocalSize, LocalSize, false);

            noalias(rLeftHandSideMatrix) = ZeroMatrix(LocalSize, LocalSize);
        }

        else {
            const unsigned int LocalSize = Dim * NumNodes;

            if (rLeftHandSideMatrix.size1() != LocalSize) rLeftHandSideMatrix.resize(LocalSize, LocalSize, false);

            noalias(rLeftHandSideMatrix) = ZeroMatrix(LocalSize, LocalSize);
            CalculateLaplacianMassMatrix(rLeftHandSideMatrix, rCurrentProcessInfo);
        }
    }

template<class TElementData>
void QSVMSDEMCoupled<TElementData>::CalculateRightHandSide(VectorType& rRightHandSideVector,
                                                           ProcessInfo& rCurrentProcessInfo)
    {
        // Calculate this element's geometric parameters
        double Area;
        array_1d<double, NumNodes> N;
        BoundedMatrix<double, NumNodes, Dim> DN_DX;
        GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DX, N, Area);

        // Calculate this element's fluid properties
        double Density;
        this->EvaluateInPoint(Density, DENSITY, N);

        if (rCurrentProcessInfo[FRACTIONAL_STEP] == 1) {

            const unsigned int LocalSize = (Dim + 1) * NumNodes;

            // Check sizes and initialize
            if (rRightHandSideVector.size() != LocalSize)
                rRightHandSideVector.resize(LocalSize, false);

            noalias(rRightHandSideVector) = ZeroVector(LocalSize);

            // Calculate Momentum RHS contribution
            this->AddMomentumRHS(rRightHandSideVector, Density, N, Area);
    //G
            const double& DeltaTime = rCurrentProcessInfo[DELTA_TIME];
            static const double arr[] = {1.0,-1.0};
            std::vector<double> SchemeWeights (arr, arr + sizeof(arr) / sizeof(arr[0]));
            this->AddMassRHS(rRightHandSideVector, Density, N, Area, SchemeWeights, DeltaTime);
        }

        else {
            const unsigned int LocalSize = Dim * NumNodes;

            // Check sizes and initialize
            if (rRightHandSideVector.size() != LocalSize)
                rRightHandSideVector.resize(LocalSize, false);

            noalias(rRightHandSideVector) = ZeroVector(LocalSize);

            this->AddRHSLaplacian(rRightHandSideVector, DN_DX, Area);
        }
        if (rCurrentProcessInfo[OSS_SWITCH] == 1)
        {
            array_1d<double, 3 > AdvVel;
            this->GetAdvectiveVel(AdvVel, N);

            double KinViscosity;
            this->EvaluateInPoint(KinViscosity, VISCOSITY, N);

            double Viscosity;
            this->GetEffectiveViscosity(Density,KinViscosity, N, DN_DX, Viscosity, rCurrentProcessInfo);

            // Calculate stabilization parameters
            double TauOne, TauTwo;
            this->CalculateTau(TauOne, TauTwo, AdvVel, Area,Density, Viscosity, rCurrentProcessInfo);

            this->AddProjectionToRHS(rRightHandSideVector, AdvVel, Density, TauOne, TauTwo, N, DN_DX, Area, rCurrentProcessInfo[DELTA_TIME]);
        }
    }

template<class TElementData>
void QSVMSDEMCoupled<TElementData>::CalculateTau(double& TauOne, double& TauTwo, const array_1d< double, 3 > & rAdvVel,
                                                 const double Area,
                                                 const double Density,
                                                 const double KinViscosity,
                                                 const ProcessInfo& rCurrentProcessInfo)
    {
        // Compute mean advective velocity norm
        double AdvVelNorm = 0.0;
        for (unsigned int d = 0; d <  Dim; ++d)
            AdvVelNorm += rAdvVel[d] * rAdvVel[d];

        AdvVelNorm = sqrt(AdvVelNorm);

        const double Element_Size = this->ElementSize(Area);
//G
      // TauOne = 1.0 / (Density * ( rCurrentProcessInfo[DYNAMIC_TAU] / rCurrentProcessInfo[DELTA_TIME] + 4 * KinViscosity / (Element_Size * Element_Size) + 2.0 * AdvVelNorm / Element_Size) );
         TauOne = 1.0 / (Density * ( rCurrentProcessInfo[DYNAMIC_TAU] / rCurrentProcessInfo[DELTA_TIME] + 5.6666666666 * KinViscosity / (Element_Size * Element_Size) + 2.0 * AdvVelNorm / Element_Size) );
//Z
        TauTwo = Density * (KinViscosity + 0.5 * Element_Size * AdvVelNorm);
        //TauOne = 1.0 / (Density * ( rCurrentProcessInfo[DYNAMIC_TAU] / rCurrentProcessInfo[DELTA_TIME] + 12.0 * KinViscosity / (Element_Size * Element_Size) + 2.0 * AdvVelNorm / Element_Size) );
        //TauTwo = Density * (KinViscosity + Element_Size * AdvVelNorm / 6.0);


    }

template<class TElementData>
void QSVMSDEMCoupled<TElementData>::CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo)
    {
        const unsigned int LocalSize = (Dim + 1) * NumNodes;

        // Resize and set to zero
        if (rMassMatrix.size1() != LocalSize)
            rMassMatrix.resize(LocalSize, LocalSize, false);

        rMassMatrix = ZeroMatrix(LocalSize, LocalSize);

        // Get the element's geometric parameters
        double Area;
        array_1d<double, NumNodes> N;
        BoundedMatrix<double, NumNodes, Dim> DN_DX;
        GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DX, N, Area);

        // Calculate this element's fluid properties
        double Density;
        this->EvaluateInPoint(Density, DENSITY, N);

        // Add 'classical' mass matrix (lumped)
        double Coeff = Density * Area / NumNodes; //Optimize!
//G
        this->CalculateLumpedMassMatrix(rMassMatrix, Coeff);

        /* For ASGS: add dynamic stabilization terms.
         These terms are not used in OSS, as they belong to the finite element
         space and cancel out with their projections.
         */
        if (rCurrentProcessInfo[OSS_SWITCH] != 1)
        {
            double KinViscosity;
            this->EvaluateInPoint(KinViscosity, VISCOSITY, N);

            double Viscosity;
            this->GetEffectiveViscosity(Density,KinViscosity, N, DN_DX, Viscosity, rCurrentProcessInfo);

            // Get Advective velocity
            array_1d<double, 3 > AdvVel;
            this->GetAdvectiveVel(AdvVel, N);

            // Calculate stabilization parameters
            double TauOne, TauTwo;
            this->CalculateTau(TauOne, TauTwo, AdvVel, Area, Density, Viscosity, rCurrentProcessInfo);

            this->AddMassStabTerms(rMassMatrix, Density, AdvVel, TauOne, N, DN_DX, Area);

//Z
        }
    }

template<class TElementData>
void QSVMSDEMCoupled<TElementData>::AddMassStabTerms(MatrixType& rLHSMatrix,
                                                     const double Density,
                                                     const array_1d<double, 3 > & rAdvVel,
                                                     const double TauOne,
                                                     const array_1d<double, NumNodes>& rShapeFunc,
                                                     const BoundedMatrix<double, NumNodes, Dim>& rShapeDeriv,
                                                     const double Weight)
    {
        const unsigned int BlockSize = Dim + 1;

        double Coef = Weight * TauOne;
        unsigned int FirstRow(0), FirstCol(0);
        double K; // Temporary results

        // If we want to use more than one Gauss point to integrate the convective term, this has to be evaluated once per integration point
        array_1d<double, NumNodes> AGradN, AGradNMod;
        this->GetConvectionOperator(AGradN, rAdvVel, rShapeDeriv); // Get a * grad(Ni)
//G
        double AdvVelDiv = 0.0;
        this->GetAdvectiveVelDivergence(AdvVelDiv, rShapeDeriv);
        //this->GetModifiedConvectionOperator(AGradNMod, rAdvVel, AdvVelDiv, rShapeFunc, rShapeDeriv); // Get a * grad(Ni) + div(a) * Ni
        double FluidFraction;
        this->EvaluateInPoint(FluidFraction, FLUID_FRACTION, rShapeFunc);
//Z

        // Note: Dof order is (vx,vy,[vz,]p) for each node
        for (unsigned int i = 0; i < NumNodes; ++i)
        {
            // Loop over columns
            for (unsigned int j = 0; j < NumNodes; ++j)
            {
                // Delta(u) * TauOne * [ AdvVel * Grad(v) ] in velocity block
//G
                K = Coef * Density * AGradN[i] * Density * rShapeFunc[j];
                //K = Coef * Density * AGradNMod[i] * Density * rShapeFunc[j];
//Z

                for (unsigned int d = 0; d < Dim; ++d) // iterate over dimensions for velocity Dofs in this node combination
                {
                    rLHSMatrix(FirstRow + d, FirstCol + d) += K;
                    // Delta(u) * TauOne * Grad(q) in q * Div(u) block
//G
    //                rLHSMatrix(FirstRow + Dim, FirstCol + d) += Coef * Density * rShapeDeriv(i, d) * rShapeFunc[j];
                      rLHSMatrix(FirstRow + Dim, FirstCol + d) += Coef * Density * FluidFraction * rShapeDeriv(i, d) * rShapeFunc[j]; // Delta(u) * TauOne * alpha * Grad(q)
//Z
                }
                // Update column index
                FirstCol += BlockSize;
            }
            // Update matrix indices
            FirstRow += BlockSize;
            FirstCol = 0;
        }
    }

template<class TElementData>
void QSVMSDEMCoupled<TElementData>::CalculateLaplacianMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo)
    {
        const unsigned int LocalSize = Dim * NumNodes;

        // Resize and set to zero
        if (rMassMatrix.size1() != LocalSize)
            rMassMatrix.resize(LocalSize, LocalSize, false);

        rMassMatrix = ZeroMatrix(LocalSize, LocalSize);

        // Get the element's geometric parameters
        double Area;
        array_1d<double, NumNodes> N;
        BoundedMatrix<double, NumNodes, Dim> DN_DX;
        GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DX, N, Area);

        // Calculate this element's fluid properties

        // Add 'classical' mass matrix (lumped)
        double Coeff = Area / NumNodes; //Optimize!

        this->CalculateLaplacianLumpedMassMatrix(rMassMatrix, Coeff);
    }

template<class TElementData>
void QSVMSDEMCoupled<TElementData>::CalculateLocalVelocityContribution(MatrixType& rDampingMatrix,
                                                            VectorType& rRightHandSideVector,
                                                            ProcessInfo& rCurrentProcessInfo)
    {
        const unsigned int LocalSize = (Dim + 1) * NumNodes;

        // Resize and set to zero the matrix
        // Note that we don't clean the RHS because it will already contain body force (and stabilization) contributions
        if (rDampingMatrix.size1() != LocalSize)
            rDampingMatrix.resize(LocalSize, LocalSize, false);

        noalias(rDampingMatrix) = ZeroMatrix(LocalSize, LocalSize);

        // Get this element's geometric properties
        double Area;
        array_1d<double, NumNodes> N;
        BoundedMatrix<double, NumNodes, Dim> DN_DX;
        GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DX, N, Area);

        // Calculate this element's fluid properties
        double Density, KinViscosity;
        this->EvaluateInPoint(Density, DENSITY, N);
        this->EvaluateInPoint(KinViscosity, VISCOSITY, N);

        double Viscosity;
        this->GetEffectiveViscosity(Density,KinViscosity, N, DN_DX, Viscosity, rCurrentProcessInfo);

        // Get Advective velocity
        array_1d<double, 3 > AdvVel;
        this->GetAdvectiveVel(AdvVel, N);

        // Calculate stabilization parameters
        double TauOne, TauTwo;
        this->CalculateTau(TauOne, TauTwo, AdvVel, Area, Density, Viscosity, rCurrentProcessInfo);
//G
        this->AddIntegrationPointVelocityContribution(rDampingMatrix, rRightHandSideVector, Density, Viscosity, AdvVel, TauOne, TauTwo, N, DN_DX, Area);

        VectorType U = ZeroVector(LocalSize);
        int LocalIndex = 0;

        for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
        {
            array_1d< double, 3 > & rVel = this->GetGeometry()[iNode].FastGetSolutionStepValue(VELOCITY);
            for (unsigned int d = 0; d < Dim; ++d) // Velocity Dofs
            {
                U[LocalIndex] = rVel[d];
                ++LocalIndex;
            }
            U[LocalIndex] = this->GetGeometry()[iNode].FastGetSolutionStepValue(PRESSURE); // Pressure Dof
            ++LocalIndex;
        }

        noalias(rRightHandSideVector) -= prod(rDampingMatrix, U);
    }

template<class TElementData>
void QSVMSDEMCoupled<TElementData>::FinalizeSolutionStep(ProcessInfo &rCurrentProcessInfo)
    {
        for(unsigned int iNode = 0; iNode < NumNodes; ++iNode)
        {
            this->GetGeometry()[iNode].SetLock();
            this->GetGeometry()[iNode].FastGetSolutionStepValue(FLUID_FRACTION_OLD) =  this->GetGeometry()[iNode].FastGetSolutionStepValue(FLUID_FRACTION);
            this->GetGeometry()[iNode].UnSetLock();
        }
    }

    /// Write the divergence of the advective velocity evaluated at this point to an array
    /**
     * Writes the value of the divergence of the advective velocity evaluated at a point inside
     * the element to a double
     * @param rAdvVelDiv: Output array
     * @param rShapeDeriv: Derivatives of shape functions evaluated at the integration point
     * @param Step: The time Step
     */
template<class TElementData>
void QSVMSDEMCoupled<TElementData>::GetAdvectiveVelDivergence(double & rAdvVelDiv,
                                                              const BoundedMatrix<double, NumNodes, Dim >& rShapeDeriv)
    {
        rAdvVelDiv = 0.0;

        for (unsigned int iNode = 1; iNode < NumNodes; ++iNode){ // loop over nodes
            const array_1d< double, 3 > vel_at_nodes = this->GetGeometry()[iNode].FastGetSolutionStepValue(VELOCITY) - this->GetGeometry()[iNode].FastGetSolutionStepValue(MESH_VELOCITY);

            for (unsigned int d = 1; d < Dim; ++d){
                // loop over components
                rAdvVelDiv += vel_at_nodes[d] * rShapeDeriv(iNode, d);
              }

          }

    }

template<class TElementData>
void QSVMSDEMCoupled<TElementData>::GetAdvectiveVelDivergence(double & rAdvVelDiv,
                                                              const BoundedMatrix<double, NumNodes, Dim >& rShapeDeriv,
                                                              const std::size_t Step)
    {
        rAdvVelDiv = 0.0;

        for (unsigned int iNode = 1; iNode < NumNodes; ++iNode){ // loop over nodes
            array_1d< double, 3 > vel_at_nodes = this->GetGeometry()[iNode].FastGetSolutionStepValue(VELOCITY, Step) - this->GetGeometry()[iNode].FastGetSolutionStepValue(MESH_VELOCITY, Step);

            for (unsigned int d = 1; d < Dim; ++d){
                // loop over components
                rAdvVelDiv += vel_at_nodes[d] * rShapeDeriv(iNode, d);
              }

          }

    }

template<class TElementData>
void QSVMSDEMCoupled<TElementData>::Calculate(const Variable<array_1d<double, 3 > >& rVariable,
                           array_1d<double, 3 > & rOutput,
                           const ProcessInfo& rCurrentProcessInfo)
    {
        if (rVariable == ADVPROJ) // Compute residual projections for OSS
        {
            // Get the element's geometric parameters
            double Area;
            array_1d<double, NumNodes> N;
            BoundedMatrix<double, NumNodes, Dim> DN_DX;
            GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DX, N, Area);

            // Calculate this element's fluid properties
            double Density;
            this->EvaluateInPoint(Density, DENSITY, N);

            // Get Advective velocity
            array_1d<double, 3 > AdvVel;
            this->GetAdvectiveVel(AdvVel, N);

            // Output containers
            array_1d< double, 3 > ElementalMomRes(3, 0.0);
            double ElementalMassRes(0);

            this->AddProjectionResidualContribution(AdvVel, Density, ElementalMomRes, ElementalMassRes, N, DN_DX, Area);

            if (rCurrentProcessInfo[OSS_SWITCH] == 1)
            {
                // Carefully write results to nodal variables, to avoid parallelism problems
                for (unsigned int i = 0; i < NumNodes; ++i)
                {
                    this->GetGeometry()[i].SetLock(); // So it is safe to write in the node in OpenMP
                    array_1d< double, 3 > & rAdvProj = this->GetGeometry()[i].FastGetSolutionStepValue(ADVPROJ);
                    for (unsigned int d = 0; d < Dim; ++d)
                        rAdvProj[d] += N[i] * ElementalMomRes[d];

                    this->GetGeometry()[i].FastGetSolutionStepValue(DIVPROJ) += N[i] * ElementalMassRes;
                    this->GetGeometry()[i].FastGetSolutionStepValue(NODAL_AREA) += Area * N[i];
                    this->GetGeometry()[i].UnSetLock(); // Free the node for other threads
                }
            }

            /// Return output
            rOutput = ElementalMomRes;
        }
        else if (rVariable == SUBSCALE_VELOCITY)
        {
            // Get the element's geometric parameters
            double Area;
            array_1d<double, NumNodes> N;
            BoundedMatrix<double, NumNodes, Dim> DN_DX;
            GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DX, N, Area);

            // Calculate this element's fluid properties
            double Density;
            this->EvaluateInPoint(Density, DENSITY, N);

            // Get Advective velocity
            array_1d<double, 3 > AdvVel;
            this->GetAdvectiveVel(AdvVel, N);

            // Output containers
            array_1d< double, 3 > ElementalMomRes(3,0.0);
            double ElementalMassRes(0.0);

            this->AddProjectionResidualContribution(AdvVel, Density, ElementalMomRes, ElementalMassRes, N, DN_DX, Area);

            if (rCurrentProcessInfo[OSS_SWITCH] == 1)
            {
                /* Projections of the elemental residual are computed with
                 * Newton-Raphson iterations of type M(lumped) dx = ElemRes - M(consistent) * x
                 */
                const double Weight = ConsistentMassCoef(Area); // Consistent mass matrix is Weight * ( Ones(NumNodes,NumNodes) + Identity(NumNodes,NumNodes) )
                // Carefully write results to nodal variables, to avoid parallelism problems
                for (unsigned int i = 0; i < NumNodes; ++i)
                {
                    this->GetGeometry()[i].SetLock(); // So it is safe to write in the node in OpenMP

                    // Add elemental residual to RHS
                    array_1d< double, 3 > & rMomRHS = this->GetGeometry()[i].GetValue(ADVPROJ);
                    double& rMassRHS = this->GetGeometry()[i].GetValue(DIVPROJ);
                    for (unsigned int d = 0; d < Dim; ++d)
                        rMomRHS[d] += N[i] * ElementalMomRes[d];

                    rMassRHS += N[i] * ElementalMassRes;

                    // Write nodal area
                    this->GetGeometry()[i].FastGetSolutionStepValue(NODAL_AREA) += Area * N[i];

                    // Substract M(consistent)*x(i-1) from RHS
                    for(unsigned int j = 0; j < NumNodes; ++j) // RHS -= Weigth * Ones(NumNodes,NumNodes) * x(i-1)
                    {
                        for(unsigned int d = 0; d < Dim; ++d)
                            rMomRHS[d] -= Weight * this->GetGeometry()[j].FastGetSolutionStepValue(ADVPROJ)[d];
                        rMassRHS -= Weight * this->GetGeometry()[j].FastGetSolutionStepValue(DIVPROJ);
                    }
                    for(unsigned int d = 0; d < Dim; ++d) // RHS -= Weigth * Identity(NumNodes,NumNodes) * x(i-1)
                        rMomRHS[d] -= Weight * this->GetGeometry()[i].FastGetSolutionStepValue(ADVPROJ)[d];
                    rMassRHS -= Weight * this->GetGeometry()[i].FastGetSolutionStepValue(DIVPROJ);

                    this->GetGeometry()[i].UnSetLock(); // Free the node for other threads
                }
            }

            /// Return output
            rOutput = ElementalMomRes;
        }
    }

template<class TElementData>
void QSVMSDEMCoupled<TElementData>::CalculateLaplacianLumpedMassMatrix(MatrixType& rLHSMatrix,
                                                                       const double Mass)
    {
        unsigned int DofIndex = 0;
        for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
        {
            for (unsigned int d = 0; d < Dim; ++d)
            {
                rLHSMatrix(DofIndex, DofIndex) += Mass;
                ++DofIndex;
            }
        }
    }

template<class TElementData>
void QSVMSDEMCoupled<TElementData>::CalculateStaticTau(double& TauOne,
                                                       const array_1d< double, 3 > & rAdvVel,
                                                       const double Area,
                                                       const double Density,
                                                       const double KinViscosity)
    {
        // Compute mean advective velocity norm
        double AdvVelNorm = 0.0;
        for (unsigned int d = 0; d < Dim; ++d)
            AdvVelNorm += rAdvVel[d] * rAdvVel[d];

        AdvVelNorm = sqrt(AdvVelNorm);

        const double Element_Size = this->ElementSize(Area);

        TauOne = 1.0 / (Density*(4.0 * KinViscosity / (Element_Size * Element_Size) + 2.0 * AdvVelNorm / Element_Size));
    }

    /// Add the momentum equation contribution to the RHS (body forces)
template<class TElementData>
void QSVMSDEMCoupled<TElementData>::AddMomentumRHS(VectorType& F,
                                                   const double Density,
                                                   const array_1d<double, NumNodes>& rShapeFunc,
                                                   const double Weight)
    {
        double Coef = Density * Weight;

        array_1d<double, 3 > BodyForce(3, 0.0);
        this->EvaluateInPoint(BodyForce, BODY_FORCE, rShapeFunc);

        // Add the results to the velocity components (Local Dofs are vx, vy, [vz,] p for each node)
        int LocalIndex = 0;

        for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
        {
            for (unsigned int d = 0; d < Dim; ++d)
            {
                F[LocalIndex++] += Coef * rShapeFunc[iNode] * BodyForce[d];
            }
            ++LocalIndex; // Skip pressure Dof
        }
    }

template<class TElementData>
void QSVMSDEMCoupled<TElementData>::AddRHSLaplacian(VectorType& F,
                                        const BoundedMatrix<double, NumNodes, Dim>& rShapeDeriv,
                                        const double Weight)
    {
        double Coef = Weight;
        array_1d<double, 3 > Velocity(3, 0.0);

        int LocalIndex = 0;
        int LocalNodalIndex = 0;

        for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
        {
            Velocity = this->GetGeometry()[iNode].FastGetSolutionStepValue(VELOCITY);
            for (unsigned int d = 0; d < Dim; ++d)
            {
                F[LocalIndex++] -= Coef * rShapeDeriv(LocalNodalIndex, d) * Velocity[d] * rShapeDeriv(iNode, d);
            }
            LocalNodalIndex++;
        }
    }

template<class TElementData>
void QSVMSDEMCoupled<TElementData>::AddMassRHS(VectorType& F,
                                               const double Density,
                                               const array_1d<double, NumNodes>& rShapeFunc,
                                               const double Weight,
                                               const std::vector<double>& TimeSchemeWeights,
                                               const double& DeltaTime)
    {
      double FluidFractionRate = 0.0;
      this->EvaluateTimeDerivativeInPoint(FluidFractionRate, FLUID_FRACTION_RATE, rShapeFunc, DeltaTime, TimeSchemeWeights);
      //this->EvaluateInPoint(FluidFractionRate, FLUID_FRACTION_RATE, rShapeFunc);
      // Add the results to the pressure components (Local Dofs are vx, vy, [vz,] p for each node)
      int LocalIndex = Dim;

      for (unsigned int iNode = 0; iNode < NumNodes; ++iNode){
          F[LocalIndex] -= Weight * rShapeFunc[iNode] * FluidFractionRate;
          LocalIndex += Dim + 1;
      }

    }

template<class TElementData>
void QSVMSDEMCoupled<TElementData>::AddProjectionToRHS(VectorType& RHS,
                                                       const array_1d<double, 3 > & rAdvVel,
                                                       const double Density,
                                                       const double TauOne,
                                                       const double TauTwo,
                                                       const array_1d<double, NumNodes>& rShapeFunc,
                                                       const BoundedMatrix<double, NumNodes, Dim>& rShapeDeriv,
                                                       const double Weight,
                                                       const double DeltaTime)
    {
        const unsigned int BlockSize = Dim + 1;

        array_1d<double, NumNodes> AGradN;
        this->GetConvectionOperator(AGradN, rAdvVel, rShapeDeriv); // Get a * grad(Ni)
        array_1d<double,3> MomProj(3,0.0);
        double DivProj = 0.0;
        this->EvaluateInPoint(MomProj,ADVPROJ,rShapeFunc);
        this->EvaluateInPoint(DivProj,DIVPROJ,rShapeFunc);

        MomProj *= TauOne;
        DivProj *= TauTwo;

        unsigned int FirstRow = 0;

        for (unsigned int i = 0; i < NumNodes; i++)
        {
//G
            double FluidFraction = this->GetGeometry()[i].FastGetSolutionStepValue(FLUID_FRACTION);
            array_1d<double,3> FluidFractionGradient(3,0.0);
            FluidFractionGradient[0] += rShapeDeriv(i, 0) * FluidFraction;
            FluidFractionGradient[1] += rShapeDeriv(i, 1) * FluidFraction;
            FluidFractionGradient[2] += rShapeDeriv(i, 2) * FluidFraction;
//Z

            for (unsigned int d = 0; d < Dim; d++)
            {
//G
//              RHS[FirstRow+d] -= Weight * (Density * AGradN[i] * MomProj[d] + rShapeDeriv(i,d) * DivProj); // TauOne * ( a * Grad(v) ) * MomProjection + TauTwo * Div(v) * MassProjection
                RHS[FirstRow+d] -= Weight * (Density * AGradN[i] * MomProj[d] + (FluidFraction * rShapeDeriv(i,d) + FluidFractionGradient[d] * rShapeFunc[i]) * DivProj); // TauOne * (a * Grad(v)) * MomProjection + TauTwo * Div(v) * MassProjection
//Z
                RHS[FirstRow+Dim] -= Weight * rShapeDeriv(i,d) * MomProj[d]; // TauOne * Grad(q) * MomProjection
            }
            FirstRow += BlockSize;
        }
    }

template<class TElementData>
void QSVMSDEMCoupled<TElementData>::EvaluateTimeDerivativeInPoint(double& rResult,
                                                                  const Variable< double >& rVariable,
                                                                  const array_1d< double,  NumNodes >& rShapeFunc,
                                                                  const double& DeltaTime,
                                                                  const std::vector<double>& rSchemeWeigths)
    {
        // Compute the time derivative of a nodal variable as a liner contribution of weighted value of the nodal variable in the (Gauss) Point

      if (rVariable == FLUID_FRACTION_RATE){
//          int index = 0;
          double delta_time_inv = 1.0 / DeltaTime;

//          for (unsigned int iWeight = 0; iWeight < rSchemeWeigths.size(); ++iWeight){

//              for (unsigned int iNode = 0; iNode <  NumNodes; ++iNode){

//                  rResult += rSchemeWeigths[iWeight] * rShapeFunc[iNode] * this->GetGeometry()[iNode].FastGetSolutionStepValue(FLUID_FRACTION, index);
//              }

//              ++index;

//          }


           for (unsigned int iNode = 0; iNode <  NumNodes; ++iNode){
              double rate = delta_time_inv * (this->GetGeometry()[iNode].FastGetSolutionStepValue(FLUID_FRACTION) - this->GetGeometry()[iNode].FastGetSolutionStepValue(FLUID_FRACTION_OLD));
                this->GetGeometry()[iNode].SetLock();
                this->GetGeometry()[iNode].FastGetSolutionStepValue(FLUID_FRACTION_RATE) = rate;
                this->GetGeometry()[iNode].UnSetLock();
              rResult += rShapeFunc[iNode] * rate;
             }

          }

    }

template<class TElementData>
void QSVMSDEMCoupled<TElementData>::CalculateLumpedMassMatrix(MatrixType& rLHSMatrix,
                                                              const double Mass)
    {
        unsigned int DofIndex = 0;
        for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
        {
            for (unsigned int d = 0; d < Dim; ++d)
            {
                rLHSMatrix(DofIndex, DofIndex) += Mass;
                ++DofIndex;
            }
            ++DofIndex; // Skip pressure Dof
        }
    }

    /// Add a the contribution from a single integration point to the velocity contribution
template<class TElementData>
void QSVMSDEMCoupled<TElementData>::AddIntegrationPointVelocityContribution(MatrixType& rDampingMatrix,
                                                                            VectorType& rDampRHS,
                                                                            const double Density,
                                                                            const double Viscosity,
                                                                            const array_1d< double, 3 > & rAdvVel,
                                                                            const double TauOne,
                                                                            const double TauTwo,
                                                                            const array_1d< double, NumNodes >& rShapeFunc,
                                                                            const BoundedMatrix<double, NumNodes, Dim >& rShapeDeriv,
                                                                            const double Weight)
    {
        const unsigned int BlockSize = Dim + 1;

        // If we want to use more than one Gauss point to integrate the convective term, this has to be evaluated once per integration point
        array_1d<double, NumNodes> AGradN, AGradNMod;
        this->GetConvectionOperator(AGradN, rAdvVel, rShapeDeriv); // Get a * grad(Ni)
//G
        double AdvVelDiv = 0.0;
        this->GetAdvectiveVelDivergence(AdvVelDiv, rShapeDeriv);
        this->GetModifiedConvectionOperator(AGradNMod, rAdvVel, AdvVelDiv, rShapeFunc, rShapeDeriv); // Get a * grad(Ni) + div(a) * Ni
//Z
        // Build the local matrix and RHS
        unsigned int FirstRow(0), FirstCol(0); // position of the first term of the local matrix that corresponds to each node combination
        double K, G, PDivV, L, qF; // Temporary results

        array_1d<double,3> BodyForce(3,0.0);
        this->EvaluateInPoint(BodyForce,BODY_FORCE,rShapeFunc);
        BodyForce *= Density;
//G
        double PAlphaDivV, GAlpha, FluidFraction, FluidFractionRate;
        array_1d<double,3> FluidFractionGradient(3,0.0);
        this->EvaluateInPoint(FluidFraction, FLUID_FRACTION, rShapeFunc);
        this->EvaluateGradientOfScalarInPoint(FluidFractionGradient, FLUID_FRACTION, rShapeDeriv);

        for (unsigned int i = 0; i < NumNodes; ++i) {
            this->GetGeometry()[i].FastGetSolutionStepValue(FLUID_FRACTION_GRADIENT) = FluidFractionGradient;
        }

        this->EvaluateInPoint(FluidFractionRate,FLUID_FRACTION_RATE,rShapeFunc);

        const double EpsilonInside = false;

        if (EpsilonInside){
            array_1d<double,3> rGradEpsOverEps;
            rGradEpsOverEps = 1.0 / FluidFraction * FluidFractionGradient ;
        }

//Z
        for (unsigned int i = 0; i < NumNodes; ++i) // iterate over rows
        {
            for (unsigned int j = 0; j < NumNodes; ++j) // iterate over columns
            {
                // Calculate the part of the contributions that is constant for each node combination

                // Velocity block
                K = Density * rShapeFunc[i] * AGradN[j]; // Convective term: v * ( a * Grad(u) )
                //K = 0.5 * Density * (rShapeFunc[i] * AGradN[j] - AGradN[i] * rShapeFunc[j]); // Skew-symmetric convective term 1/2( v*grad(u)*u - grad(v) uu )
//G
                K += TauOne * Density * AGradN[i] * Density * AGradN[j]; // Stabilization: (a * Grad(v)) * TauOne * (a * Grad(u))
                //K += TauOne * Density * AGradNMod[i] * Density * AGradN[j]; // Stabilization: (a * Grad(v) + Div(a) * v) * TauOne * (a * Grad(u))
//Z
                K *= Weight;

                // q-p stabilization block (reset result)
                L = 0;

                for (unsigned int m = 0; m < Dim; ++m) // iterate over v components (vx,vy[,vz])
                {
                    // Velocity block
                    //K += Weight * Density * Viscosity * rShapeDeriv(i, m) * rShapeDeriv(j, m); // Diffusive term: Viscosity * Grad(v) * Grad(u)

                    // v * Grad(p) block
//G
                    G = TauOne * Density * AGradN[i] * rShapeDeriv(j, m); // Stabilization: (a * Grad(v)) * TauOne * Grad(p)
                    //G = TauOne * Density * AGradNMod[i] * rShapeDeriv(j, m); // Stabilization: (a * Grad(v) + Div(a) * v) * TauOne * Grad(p)

                    GAlpha = TauOne * Density * AGradN[i] * (FluidFraction * rShapeDeriv(j, m)); // Stabilization: (a * Grad(u)) * TauOne * (alpha * Grad(q))

                    PAlphaDivV = (FluidFraction * rShapeDeriv(i, m) + FluidFractionGradient[m] * rShapeFunc[i]) * rShapeFunc[j]; // alpha * q * Div(u) + q * Grad(alpha) * u
//Z
                    PDivV = rShapeDeriv(i, m) * rShapeFunc[j]; // Div(v) * p
                    // Write v * Grad(p) component
                    rDampingMatrix(FirstRow + m, FirstCol + Dim) += Weight * (G - PDivV);
                    // Use symmetry to write the q * Div(u) component
//G
    //              rDampingMatrix(FirstCol + Dim, FirstRow + m) += Weight * (G + PDivV);
                    rDampingMatrix(FirstCol + Dim, FirstRow + m) += Weight * (GAlpha + PAlphaDivV); // note that here PAlphaDivV includes G; do not look for it in GAlpha!

    //              q-p stabilization block
    //              L += rShapeDeriv(i, m) * rShapeDeriv(j, m); // Stabilization: Grad(q) * TauOne * Grad(p)
                    L += FluidFraction * rShapeDeriv(i, m) * rShapeDeriv(j, m); // Stabilization: alpha * Grad(q) * TauOne * Grad(p)
//Z
                    for (unsigned int n = 0; n < Dim; ++n) // iterate over u components (ux,uy[,uz])
                    {
                        // Velocity block
//G
    //                  rDampingMatrix(FirstRow + m, FirstCol + n) += Weight * TauTwo * rShapeDeriv(i, m) * rShapeDeriv(j, n); // Stabilization: Div(v) * TauTwo * Div(u)
                        rDampingMatrix(FirstRow + m, FirstCol + n) += Weight * TauTwo * rShapeDeriv(i, m) * (FluidFraction * rShapeDeriv(j, n) + FluidFractionGradient[n] * rShapeFunc[j]); // Stabilization: Div(v) * TauTwo * (alpha * Div(u) + Grad(alpha) * u)
//Z
                    }

                }

                // Write remaining terms to velocity block
                for (unsigned int d = 0; d < Dim; ++d)
                    rDampingMatrix(FirstRow + d, FirstCol + d) += K;

                // Write q-p stabilization block
                rDampingMatrix(FirstRow + Dim, FirstCol + Dim) += Weight * TauOne * L;


                // Update reference column index for next iteration
                FirstCol += BlockSize;
            }

            // Operate on RHS
            qF = 0.0;

            for (unsigned int d = 0; d < Dim; ++d)
            {
                //rDampRHS[FirstRow + d] += Weight * TauOne * Density * AGradN[i] * BodyForce[d]; // ( a * Grad(v) ) * TauOne * (Density * BodyForce)
//G
//GG                rDampRHS[FirstRow + d] += Weight * (TauOne * Density * AGradNMod[i] * BodyForce[d] - TauTwo * rShapeDeriv(i, d) * FluidFractionRate); // ( a * Grad(v) ) * TauOne * (Density * BodyForce) + Div(v) * TauTwo * (- DAlphaDt)
                rDampRHS[FirstRow + d] += Weight * (TauOne * Density * AGradN[i] * BodyForce[d] - TauTwo * rShapeDeriv(i, d) * FluidFractionRate); // ( a * Grad(v) ) * TauOne * (Density * BodyForce) + Div(v) * TauTwo * (- DAlphaDt)

    //          qF += rShapeDeriv(i, d) * BodyForce[d];
                qF += (FluidFraction * rShapeDeriv(i, d)) * BodyForce[d];
//Z
            }
            rDampRHS[FirstRow + Dim] += Weight * TauOne * qF; // Grad(q) * TauOne * (Density * BodyForce)

            // Update reference indices
            FirstRow += BlockSize;
            FirstCol = 0;
        }

//            this->AddBTransCB(rDampingMatrix,rShapeDeriv,Viscosity*Coef);
        this->AddViscousTerm(rDampingMatrix,rShapeDeriv,Viscosity*Density*Weight);
    }

template<class TElementData>
void QSVMSDEMCoupled<TElementData>::Calculate(const Variable<double>& rVariable,
                                              double& rOutput,
                                              const ProcessInfo& rCurrentProcessInfo)
    {
        if (rVariable == ERROR_RATIO)
        {
            // Get the element's geometric parameters
            double Area;
            array_1d<double, NumNodes> N;
            BoundedMatrix<double, NumNodes, Dim> DN_DX;
            GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DX, N, Area);

            // Calculate this element's fluid properties
            double Density, KinViscosity;
            this->EvaluateInPoint(Density, DENSITY, N);
            this->EvaluateInPoint(KinViscosity, VISCOSITY, N);

            double Viscosity;
            this->GetEffectiveViscosity(Density,KinViscosity, N, DN_DX, Viscosity, rCurrentProcessInfo);

            // Get Advective velocity
            array_1d<double, 3 > AdvVel;
            this->GetAdvectiveVel(AdvVel, N);

            // Output container
            array_1d< double, 3 > ElementalMomRes(3, 0.0);

            // Calculate stabilization parameter. Note that to estimate the subscale velocity, the dynamic coefficient in TauOne is assumed zero.
            double TauOne;
            this->CalculateStaticTau(TauOne, AdvVel, Area,Density, Viscosity);

            if ( rCurrentProcessInfo[OSS_SWITCH] != 1 ) // ASGS
            {
//G
            this->ASGSMomResidual(AdvVel,Density,ElementalMomRes,N,DN_DX,1.0);


/*
                MatrixType NContainer;
                ShapeFunctionDerivativesArrayType DN_DXContainer;
                VectorType GaussWeights;
                this->CalculateWeights(DN_DXContainer, NContainer, GaussWeights);
                const SizeType NumGauss = NContainer.size1();

                for (SizeType g = 0; g < NumGauss; g++){
                    const double GaussWeight = GaussWeights[g];
                    const ShapeFunctionsType& Ng = row(NContainer, g);
                    this->GetAdvectiveVel(AdvVel, Ng);
                    this->ASGSMomResidual(AdvVel, Density, ElementalMomRes, Ng, DN_DX, GaussWeight);
                  }

*/
//Z
                ElementalMomRes *= TauOne;
            }
            else // OSS
            {
                this->OSSMomResidual(AdvVel,Density,ElementalMomRes,N,DN_DX,1.0);;
                ElementalMomRes *= TauOne;
            }

            // Error estimation ( ||U'|| / ||Uh_gauss|| ), taking ||U'|| = TauOne ||MomRes||
            double ErrorRatio(0.0);//, UNorm(0.0);
//                array_1d< double, 3 > UGauss(3, 0.0);
//                this->AddPointContribution(UGauss, VELOCITY, N);

            for (unsigned int i = 0; i < Dim; ++i)
            {
                ErrorRatio += ElementalMomRes[i] * ElementalMomRes[i];
//                    UNorm += UGauss[i] * UGauss[i];
            }
            ErrorRatio = sqrt(ErrorRatio); // / UNorm);
            ErrorRatio /= Density;
            this->SetValue(ERROR_RATIO, ErrorRatio);
            rOutput = ErrorRatio;
        }
        else if (rVariable == NODAL_AREA)
        {
            // Get the element's geometric parameters
            double Area;
            array_1d<double, NumNodes> N;
            BoundedMatrix<double, NumNodes, Dim> DN_DX;
            GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DX, N, Area);

            // Carefully write results to nodal variables, to avoid parallelism problems
            for (unsigned int i = 0; i < NumNodes; ++i)
            {
                this->GetGeometry()[i].SetLock(); // So it is safe to write in the node in OpenMP
                this->GetGeometry()[i].FastGetSolutionStepValue(NODAL_AREA) += Area * N[i];
                this->GetGeometry()[i].UnSetLock(); // Free the node for other threads
            }
        }
    }

template<class TElementData>
void QSVMSDEMCoupled<TElementData>::GetModifiedConvectionOperator(array_1d< double,  NumNodes >& rResult,
                                                                  const array_1d< double, 3 > & rVelocity,
                                                                  const double & rVelocityDiv,
                                                                  const array_1d< double,  NumNodes >& rShapeFunc,
                                                                  const BoundedMatrix<double,  NumNodes,  Dim >& rShapeDeriv)
    {
        // Evaluate (and weight) the a * Grad(Ni) + div(a) * Ni operator in the integration point, for each node i
        for (unsigned int iNode = 0; iNode <  NumNodes; ++iNode){ // Loop over nodes{
            // Initialize result
            rResult[iNode] = rVelocityDiv * rShapeFunc[iNode];

            for (unsigned int d = 0; d <  Dim; ++d){ // loop over components
                rResult[iNode] += rVelocity[d] * rShapeDeriv(iNode, d);
              }

        }
    }

template<class TElementData>
void QSVMSDEMCoupled<TElementData>::GetValueOnIntegrationPoints(const Variable<double>& rVariable,
                                                                std::vector<double>& rValues,
                                                                const ProcessInfo& rCurrentProcessInfo)
    {
        if (rVariable == TAUONE || rVariable == TAUTWO || rVariable == MU)
        {
            double TauOne, TauTwo;
            double Area;
            array_1d<double, NumNodes> N;
            BoundedMatrix<double, NumNodes, Dim> DN_DX;
            GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DX, N, Area);

            array_1d<double, 3 > AdvVel;
            this->GetAdvectiveVel(AdvVel, N);

            double Density,KinViscosity;
            this->EvaluateInPoint(Density, DENSITY, N);
            this->EvaluateInPoint(KinViscosity, VISCOSITY, N);

            double Viscosity;
            this->GetEffectiveViscosity(Density,KinViscosity, N, DN_DX, Viscosity, rCurrentProcessInfo);

            this->CalculateTau(TauOne, TauTwo, AdvVel, Area, Density, Viscosity, rCurrentProcessInfo);

            rValues.resize(1, false);
            if (rVariable == TAUONE)
            {
                rValues[0] = TauOne;
            }
            else if (rVariable == TAUTWO)
            {
                rValues[0] = TauTwo;
            }
            else if (rVariable == MU)
            {
                rValues[0] = Density * Viscosity;
            }
        }
        else if(rVariable == SUBSCALE_PRESSURE)
        {
            double TauOne, TauTwo;
            double Area;
            array_1d<double, NumNodes> N;
            BoundedMatrix<double, NumNodes, Dim> DN_DX;
            GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DX, N, Area);

            array_1d<double, 3 > AdvVel;
            this->GetAdvectiveVel(AdvVel, N);

            double Density,KinViscosity;
            this->EvaluateInPoint(Density, DENSITY, N);
            this->EvaluateInPoint(KinViscosity, VISCOSITY, N);

            double Viscosity;
            this->GetEffectiveViscosity(Density,KinViscosity, N, DN_DX, Viscosity, rCurrentProcessInfo);

            this->CalculateTau(TauOne, TauTwo, AdvVel, Area, Density, Viscosity, rCurrentProcessInfo);

            double DivU = 0.0;
            for(unsigned int i=0; i < NumNodes; i++)
            {
                for(unsigned int d = 0; d < Dim; d++)
                    DivU -= DN_DX(i,d) * this->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY)[d];
            }

            rValues.resize(1, false);

            rValues[0] = TauTwo * DivU;// *Density?? decide on criteria and use the same for SUBSCALE_VELOCITY

            if(rCurrentProcessInfo[OSS_SWITCH]==1)
            {
                double Proj = 0.0;
                for(unsigned int i=0; i < NumNodes; i++)
                {
                    Proj += N[i]*this->GetGeometry()[i].FastGetSolutionStepValue(DIVPROJ);
                }
                rValues[0] -= TauTwo*Proj;
            }
        }
        else if (rVariable == NODAL_AREA && Dim == 3)
        {
            MatrixType J = ZeroMatrix(3,3);
            const array_1d<double,3>& X0 = this->GetGeometry()[0].Coordinates();
            const array_1d<double,3>& X1 = this->GetGeometry()[1].Coordinates();
            const array_1d<double,3>& X2 = this->GetGeometry()[2].Coordinates();
            const array_1d<double,3>& X3 = this->GetGeometry()[3].Coordinates();

            J(0,0) = X1[0]-X0[0];
            J(0,1) = X2[0]-X0[0];
            J(0,2) = X3[0]-X0[0];
            J(1,0) = X1[1]-X0[1];
            J(1,1) = X2[1]-X0[1];
            J(1,2) = X3[1]-X0[1];
            J(2,0) = X1[2]-X0[2];
            J(2,1) = X2[2]-X0[2];
            J(2,2) = X3[2]-X0[2];

            double DetJ = J(0,0)*( J(1,1)*J(2,2) - J(1,2)*J(2,1) ) + J(0,1)*( J(1,2)*J(2,0) - J(1,0)*J(2,2) ) + J(0,2)*( J(1,0)*J(2,1) - J(1,1)*J(2,0) );
            rValues.resize(1, false);
            rValues[0] = DetJ;
        }
        else // Default behaviour (returns elemental data)
        {
            rValues.resize(1, false);
            /*
             The cast is done to avoid modification of the element's data. Data modification
             would happen if rVariable is not stored now (would initialize a pointer to &rVariable
             with associated value of 0.0). This is catastrophic if the variable referenced
             goes out of scope.
             */
            const QSVMSDEMCoupled<TElementData>* const_this = static_cast<QSVMSDEMCoupled<TElementData>* > (this);
            rValues[0] = const_this->GetValue(rVariable);
        }
    }

template<class TElementData>
void QSVMSDEMCoupled<TElementData>::ASGSMomResidual(const array_1d< double, 3 > & rAdvVel,
                         const double Density,
                         array_1d< double, 3 > & rElementalMomRes,
                         const array_1d< double, NumNodes >& rShapeFunc,
                         const BoundedMatrix<double, NumNodes, Dim >& rShapeDeriv,
                         const double Weight)
    {
        // If we want to use more than one Gauss point to integrate the convective term, this has to be evaluated once per integration point
        array_1d<double, NumNodes> AGradN;
        this->GetConvectionOperator(AGradN, rAdvVel, rShapeDeriv); // Get a * grad(Ni)
        // Compute contribution to Kij * Uj, with Kij = Ni * Residual(Nj); Uj = (v,p)Node_j (column vector)
        for (unsigned int i = 0; i < NumNodes; ++i) // Iterate over element nodes
        {

            // Variable references
            const array_1d< double, 3 > & rVelocity = this->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
            const array_1d< double, 3 > & rAcceleration = this->GetGeometry()[i].FastGetSolutionStepValue(ACCELERATION);
            const array_1d< double, 3 > & rBodyForce = this->GetGeometry()[i].FastGetSolutionStepValue(BODY_FORCE);
            const double& rPressure = this->GetGeometry()[i].FastGetSolutionStepValue(PRESSURE);

            // Compute this node's contribution to the residual (evaluated at inegration point)
            for (unsigned int d = 0; d < Dim; ++d)
            {
                rElementalMomRes[d] += Weight * (Density * (rShapeFunc[i] * (rBodyForce[d] - rAcceleration[d]) - AGradN[i] * rVelocity[d]) - rShapeDeriv(i, d) * rPressure);
            }
        }
    }

    /// Assemble the contribution from an integration point to the element's residual.
    /**
     * OSS version. Note that rElementalMomRes should be initialized before calling this.
     * @param rAdvVel Convection velocity (not including subscale)
     * @param Density Fluid density evaluated at integration point
     * @param rElementalMomRes Result
     * @param rShapeFunc Shape functions evaluated at integration point
     * @param rShapeDeriv Shape function derivatives evaluated at integration point
     * @param Weight Integration point weight (as a fraction of area or volume)
     */

template<class TElementData>
void QSVMSDEMCoupled<TElementData>::OSSMomResidual(const array_1d< double, 3 > & rAdvVel,
                        const double Density,
                        array_1d< double, 3 > & rElementalMomRes,
                        const array_1d< double, NumNodes >& rShapeFunc,
                        const BoundedMatrix<double, NumNodes, Dim >& rShapeDeriv,
                        const double Weight)
    {
        // If we want to use more than one Gauss point to integrate the convective term, this has to be evaluated once per integration point
        array_1d<double, NumNodes> AGradN;
        this->GetConvectionOperator(AGradN, rAdvVel, rShapeDeriv); // Get a * grad(Ni)
        // Compute contribution to Kij * Uj, with Kij = Ni * Residual(Nj); Uj = (v,p)Node_j (column vector)
        for (unsigned int i = 0; i < NumNodes; ++i) // Iterate over element nodes
        {

            // Variable references
            const array_1d< double, 3 > & rVelocity = this->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
            const array_1d< double, 3 > & rBodyForce = this->GetGeometry()[i].FastGetSolutionStepValue(BODY_FORCE);
            const array_1d< double, 3 > & rProjection = this->GetGeometry()[i].FastGetSolutionStepValue(ADVPROJ);
            const double& rPressure = this->GetGeometry()[i].FastGetSolutionStepValue(PRESSURE);

            // Compute this node's contribution to the residual (evaluated at inegration point)
            for (unsigned int d = 0; d < Dim; ++d)
            {
                rElementalMomRes[d] += Weight * (Density * (rShapeFunc[i] * rBodyForce[d] - AGradN[i] * rVelocity[d]) - rShapeDeriv(i, d) * rPressure);
                rElementalMomRes[d] -= Weight * rShapeFunc[i] * rProjection[d];

            }
        }
    }

template<class TElementData>
void QSVMSDEMCoupled<TElementData>::AddProjectionResidualContribution(const array_1d< double, 3 > & rAdvVel,
                                                                      const double Density,
                                                                      array_1d< double, 3 > & rElementalMomRes,
                                                                      double& rElementalMassRes,
                                                                      const array_1d< double, NumNodes >& rShapeFunc,
                                                                      const BoundedMatrix<double, NumNodes, Dim >& rShapeDeriv,
                                                                      const double Weight)
    {
        // If we want to use more than one Gauss point to integrate the convective term, this has to be evaluated once per integration point
        array_1d<double, NumNodes> AGradN;
        this->GetConvectionOperator(AGradN, rAdvVel, rShapeDeriv); // Get a * grad(Ni)
//G
        double FluidFraction, FluidFractionRate;
        array_1d<double,3> FluidFractionGradient(3,0.0);
        this->EvaluateInPoint(FluidFraction, FLUID_FRACTION, rShapeFunc);
        this->EvaluateGradientOfScalarInPoint(FluidFractionGradient, FLUID_FRACTION, rShapeDeriv);
        this->EvaluateInPoint(FluidFractionRate,FLUID_FRACTION_RATE,rShapeFunc);
//Z
        // Compute contribution to Kij * Uj, with Kij = Ni * Residual(Nj); Uj = (v,p)Node_j (column vector)
        for (unsigned int i = 0; i < NumNodes; ++i) // Iterate over element nodes
        {

            // Variable references
            const array_1d< double, 3 > & rVelocity = this->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
            const array_1d< double, 3 > & rBodyForce = this->GetGeometry()[i].FastGetSolutionStepValue(BODY_FORCE);
            const double& rPressure = this->GetGeometry()[i].FastGetSolutionStepValue(PRESSURE);

            // Compute this node's contribution to the residual (evaluated at inegration point)
            for (unsigned int d = 0; d < Dim; ++d)
            {

                rElementalMomRes[d] += Weight * (Density * (rShapeFunc[i] * rBodyForce[d] - AGradN[i] * rVelocity[d]) - rShapeDeriv(i, d) * rPressure);
//G
    //          rElementalMassRes -= Weight * rShapeDeriv(i, d) * rVelocity[d];
                rElementalMassRes -= Weight * (FluidFraction * rShapeDeriv(i, d) * rVelocity[d] + FluidFractionGradient[d] * rShapeFunc[i] * rVelocity[d]);
            }

        }

        rElementalMassRes -= Weight * FluidFractionRate;
//Z
    }

template<class TElementData>
void QSVMSDEMCoupled<TElementData>::EvaluateGradientOfScalarInPoint(array_1d< double, 3 >& rResult,
                                                                    const Variable< double >& rVariable,
                                                                    const BoundedMatrix<double, NumNodes, Dim >& rShapeDeriv)
    {

        for (unsigned int i = 0; i < NumNodes; ++i) {
            double& scalar = this->GetGeometry()[i].FastGetSolutionStepValue(rVariable);

            for (unsigned int d = 0; d < Dim; ++d){
                rResult[d] += rShapeDeriv(i, d) * scalar;
            }
        }
    }

template<class TElementData>
void QSVMSDEMCoupled<TElementData>::AddBTransCB(MatrixType& rDampingMatrix,
                                                const BoundedMatrix<double, NumNodes, Dim >& rShapeDeriv,
                                                const double Weight)
    {
        BoundedMatrix<double, (Dim * NumNodes)/2, Dim*NumNodes > B;
        BoundedMatrix<double, (Dim * NumNodes)/2, (Dim*NumNodes)/2 > C;
        this->CalculateB(B, rShapeDeriv);
        this->CalculateC(C, Weight);

        const unsigned int BlockSize = Dim + 1;
        const unsigned int StrainSize = (Dim*NumNodes)/2;

        DenseVector<unsigned int> aux(Dim*NumNodes);
        for(unsigned int i=0; i<NumNodes; i++)
        {
            int base_index = Dim*i;
            int aux_index = BlockSize*i;
            for(unsigned int j=0; j<Dim; j++)
            {
                aux[base_index+j] = aux_index+j;
            }
        }

        for(unsigned int k=0; k< StrainSize; k++)
        {
            for(unsigned int l=0; l< StrainSize; l++)
            {
                const double Ckl = C(k,l);
                for (unsigned int i = 0; i < Dim*NumNodes; ++i) // iterate over v components (vx,vy[,vz])
                {
                    const double Bki=B(k,i);
                    for (unsigned int j = 0; j < Dim*NumNodes; ++j) // iterate over u components (ux,uy[,uz])
                    {
                        rDampingMatrix(aux[i],aux[j]) += Bki*Ckl*B(l,j);
                    }

                }

            }
        }
    }

template<class TElementData>
void QSVMSDEMCoupled<TElementData>::EvaluateInPoint(double& rResult,
                                                    const Variable< double >& rVariable,
                                                    const array_1d< double, NumNodes >& rShapeFunc)
    {
        // Compute the weighted value of the nodal variable in the (Gauss) Point
        rResult = rShapeFunc[0] * this->GetGeometry()[0].FastGetSolutionStepValue(rVariable);

        for (unsigned int iNode = 1; iNode < NumNodes; ++iNode)
            rResult += rShapeFunc[iNode] * this->GetGeometry()[iNode].FastGetSolutionStepValue(rVariable);
    }

template<class TElementData>
void QSVMSDEMCoupled<TElementData>::EvaluateInPoint(array_1d< double, 3 > & rResult,
                                                    const Variable< array_1d< double, 3 > >& rVariable,
                                                    const array_1d< double, NumNodes >& rShapeFunc)
    {
        // Compute the weighted value of the nodal variable in the (Gauss) Point
        rResult = rShapeFunc[0] * this->GetGeometry()[0].FastGetSolutionStepValue(rVariable);
        for (unsigned int iNode = 1; iNode < NumNodes; ++iNode)
            rResult += rShapeFunc[iNode] * this->GetGeometry()[iNode].FastGetSolutionStepValue(rVariable);
    }

template<class TElementData>
void QSVMSDEMCoupled<TElementData>::ModulatedGradientDiffusion(MatrixType& rDampingMatrix,
            const BoundedMatrix<double, NumNodes, Dim >& rDN_DX,
            const double Weight)
    {
        const GeometryType& rGeom = this->GetGeometry();

        // Velocity gradient
        MatrixType GradU = ZeroMatrix(Dim,Dim);
        for (unsigned int n = 0; n < NumNodes; n++)
        {
            const array_1d<double,3>& rVel = this->GetGeometry()[n].FastGetSolutionStepValue(VELOCITY);
            for (unsigned int i = 0; i < Dim; i++)
                for (unsigned int j = 0; j < Dim; j++)
                    GradU(i,j) += rDN_DX(n,j)*rVel[i];
        }

        // Element lengths
        array_1d<double,3> Delta(3,0.0);
        Delta[0] = fabs(rGeom[NumNodes-1].X()-rGeom[0].X());
        Delta[1] = fabs(rGeom[NumNodes-1].Y()-rGeom[0].Y());
        Delta[2] = fabs(rGeom[NumNodes-1].Z()-rGeom[0].Z());

        for (unsigned int n = 1; n < NumNodes; n++)
        {
            double hx = fabs(rGeom[n].X()-rGeom[n-1].X());
            if (hx > Delta[0]) Delta[0] = hx;
            double hy = fabs(rGeom[n].Y()-rGeom[n-1].Y());
            if (hy > Delta[1]) Delta[1] = hy;
            double hz = fabs(rGeom[n].Z()-rGeom[n-1].Z());
            if (hz > Delta[2]) Delta[2] = hz;
        }

        double AvgDeltaSq = Delta[0];
        for (unsigned int d = 1; d < Dim; d++)
            AvgDeltaSq *= Delta[d];
        AvgDeltaSq = std::pow(AvgDeltaSq,2./Dim);

        Delta[0] = Delta[0]*Delta[0]/12.0;
        Delta[1] = Delta[1]*Delta[1]/12.0;
        Delta[2] = Delta[2]*Delta[2]/12.0;

        // Gij
        MatrixType G = ZeroMatrix(Dim,Dim);
        for (unsigned int i = 0; i < Dim; i++)
            for (unsigned int j = 0; j < Dim; j++)
                for (unsigned int d = 0; d < Dim; d++)
                    G(i,j) += Delta[d]*GradU(i,d)*GradU(j,d);

        // Gij:Sij
        double GijSij = 0.0;
        for (unsigned int i = 0; i < Dim; i++)
            for (unsigned int j = 0; j < Dim; j++)
                GijSij += 0.5*G(i,j)*( GradU(i,j) + GradU(j,i) );

        if (GijSij < 0.0) // Otherwise model term is clipped
        {
            // Gkk
            double Gkk = G(0,0);
            for (unsigned int d = 1; d < Dim; d++)
                Gkk += G(d,d);

            // C_epsilon
            const double Ce = 1.0;

            // ksgs
            double ksgs = -4*AvgDeltaSq*GijSij/(Ce*Ce*Gkk);

            // Assembly of model term
            unsigned int RowIndex = 0;
            unsigned int ColIndex = 0;

            for (unsigned int i = 0; i < NumNodes; i++)
            {
                for (unsigned int j = 0; j < NumNodes; j++)
                {
                    for (unsigned int d = 0; d < Dim; d++)
                    {
                        double Aux = rDN_DX(i,d) * Delta[0] * G(d,0)*rDN_DX(j,0);
                        for (unsigned int k = 1; k < Dim; k++)
                            Aux += rDN_DX(i,d) *Delta[k] * G(d,k)*rDN_DX(j,k);
                        rDampingMatrix(RowIndex+d,ColIndex+d) += Weight * 2.0*ksgs *  Aux;
                    }

                    ColIndex += Dim;
                }
                RowIndex += Dim;
                ColIndex = 0;
            }
        }

    }

template<class TElementData>
void QSVMSDEMCoupled<TElementData>::GetEffectiveViscosity(const double Density, const double MolecularViscosity,
                                                          const array_1d<double, NumNodes>& rShapeFunc,
                                                          const BoundedMatrix<double, NumNodes, Dim >& rShapeDeriv,
                                                          double& TotalViscosity, const ProcessInfo& rCurrentProcessInfo)
    {
        const double C = this->GetValue(C_SMAGORINSKY);

        TotalViscosity = MolecularViscosity;
        if (C != 0.0 )
        {
            // The filter width in Smagorinsky is typically the element size h. We will store the square of h, as the final formula involves the squared filter width
            const double FilterWidth = this->FilterWidth(rShapeDeriv);

            const double NormS = this->SymmetricGradientNorm(rShapeDeriv);

            // Total Viscosity
            TotalViscosity += 2.0 * C * C * FilterWidth * NormS;
        }
    }

template<class TElementData>
double QSVMSDEMCoupled<TElementData>::SymmetricGradientNorm(const BoundedMatrix<double, NumNodes, Dim >& rShapeDeriv)
    {
        const unsigned int GradientSize = (Dim*(Dim+1))/2; // Number of different terms in the symmetric gradient matrix
        array_1d<double,GradientSize> GradientVector( GradientSize, 0.0 );
        unsigned int Index;

        // Compute Symmetric Grad(u). Note that only the lower half of the matrix is calculated
        for (unsigned int k = 0; k < NumNodes; ++k)
        {
            const array_1d< double, 3 > & rNodeVel = this->GetGeometry()[k].FastGetSolutionStepValue(VELOCITY);
            Index = 0;
            for (unsigned int i = 0; i < Dim; ++i)
            {
                for (unsigned int j = 0; j < i; ++j) // Off-diagonal
                    GradientVector[Index++] += 0.5 * (rShapeDeriv(k, j) * rNodeVel[i] + rShapeDeriv(k, i) * rNodeVel[j]);
                GradientVector[Index++] += rShapeDeriv(k, i) * rNodeVel[i]; // Diagonal
            }
        }

        // Norm[ Symmetric Grad(u) ] = ( 2 * Sij * Sij )^(1/2)
        Index = 0;
        double NormS(0.0);
        for (unsigned int i = 0; i < Dim; ++i)
        {
            for (unsigned int j = 0; j < i; ++j)
            {
                NormS += 2.0 * GradientVector[Index] * GradientVector[Index]; // Using symmetry, lower half terms of the matrix are added twice
                ++Index;
            }
            NormS += GradientVector[Index] * GradientVector[Index]; // Diagonal terms
            ++Index; // Diagonal terms
        }

        NormS = sqrt( 2.0 * NormS );
        return NormS;
    }

    /// Write the advective velocity evaluated at this point to an array
    /**
     * Writes the value of the advective velocity evaluated at a point inside
     * the element to an array_1d
     * @param rAdvVel: Output array
     * @param rShapeFunc: Shape functions evaluated at the point of interest
     */
template<class TElementData>
void QSVMSDEMCoupled<TElementData>::GetAdvectiveVel(array_1d< double, 3 > & rAdvVel,
                                 const array_1d< double, NumNodes >& rShapeFunc)
    {
        // Compute the weighted value of the advective velocity in the (Gauss) Point
        rAdvVel = rShapeFunc[0] * (this->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY) - this->GetGeometry()[0].FastGetSolutionStepValue(MESH_VELOCITY));
        for (unsigned int iNode = 1; iNode < NumNodes; ++iNode)
            rAdvVel += rShapeFunc[iNode] * (this->GetGeometry()[iNode].FastGetSolutionStepValue(VELOCITY) - this->GetGeometry()[iNode].FastGetSolutionStepValue(MESH_VELOCITY));
    }

    /// Write the advective velocity evaluated at this point to an array
    /**
     * Writes the value of the advective velocity evaluated at a point inside
     * the element to an array_1d
     * @param rAdvVel: Output array
     * @param rShapeFunc: Shape functions evaluated at the point of interest
     * @param Step: The time Step
     */

template<class TElementData>
void QSVMSDEMCoupled<TElementData>::GetAdvectiveVel(array_1d< double, 3 > & rAdvVel,
                                 const array_1d< double, NumNodes >& rShapeFunc,
                                 const std::size_t Step)
    {
        // Compute the weighted value of the advective velocity in the (Gauss) Point
        rAdvVel = rShapeFunc[0] * (this->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY, Step) - this->GetGeometry()[0].FastGetSolutionStepValue(MESH_VELOCITY, Step));
        for (unsigned int iNode = 1; iNode < NumNodes; ++iNode)
            rAdvVel += rShapeFunc[iNode] * (this->GetGeometry()[iNode].FastGetSolutionStepValue(VELOCITY, Step) - this->GetGeometry()[iNode].FastGetSolutionStepValue(MESH_VELOCITY, Step));
    }

template<class TElementData>
void QSVMSDEMCoupled<TElementData>::GetConvectionOperator(array_1d< double, NumNodes >& rResult,
                               const array_1d< double, 3 > & rVelocity,
                               const BoundedMatrix<double, NumNodes, Dim >& rShapeDeriv)
    {
        // Evaluate (and weight) the a * Grad(Ni) operator in the integration point, for each node i
        for (unsigned int iNode = 0; iNode < NumNodes; ++iNode) // Loop over nodes
        {
            // Initialize result
            rResult[iNode] = rVelocity[0] * rShapeDeriv(iNode, 0);
            for (unsigned int d = 1; d < Dim; ++d) // loop over components
                rResult[iNode] += rVelocity[d] * rShapeDeriv(iNode, d);
        }
    }
///////////////////////////////////////////////////////////////////////////////////////////////////
// Private functions
///////////////////////////////////////////////////////////////////////////////////////////////////

// serializer

template< class TElementData >
void QSVMSDEMCoupled<TElementData>::save(Serializer& rSerializer) const
{
    typedef QSVMS<TElementData> BaseElement;
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseElement );
}


template< class TElementData >
void QSVMSDEMCoupled<TElementData>::load(Serializer& rSerializer)
{
    typedef QSVMS<TElementData> BaseElement;
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseElement);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation

template class QSVMSDEMCoupled<QSVMSData< 2, 3 >>;
template class QSVMSDEMCoupled<QSVMSData< 3, 4 >>;

}

