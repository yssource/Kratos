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

#include "../FluidDynamicsApplication/fluid_dynamics_application_variables.h"
#include "../FluidDynamicsApplication/fluid_dynamics_application.h"
#include "../FluidDynamicsApplication/custom_utilities/qsvms_data.h"
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

template< class TElementData >
Element::Pointer QSVMSDEMCoupled<TElementData>::Create(IndexType NewId,GeometryType::Pointer pGeom,Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<QSVMSDEMCoupled>(NewId, pGeom, pProperties);
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

template <class TElementData>
void QSVMSDEMCoupled<TElementData>::Calculate(const Variable<double>& rVariable,
    double& rOutput, const ProcessInfo& rCurrentProcessInfo) {}

template <class TElementData>
void QSVMSDEMCoupled<TElementData>::Calculate(const Variable<array_1d<double, 3>>& rVariable,
                                              array_1d<double, 3>& rOutput,
                                              const ProcessInfo& rCurrentProcessInfo) {

    QSVMS<TElementData>::Calculate(rVariable, rOutput, rCurrentProcessInfo);

}

template <class TElementData>
void QSVMSDEMCoupled<TElementData>::Calculate(const Variable<Vector>& rVariable,
    Vector& rOutput, const ProcessInfo& rCurrentProcessInfo) {}

template <class TElementData>
void QSVMSDEMCoupled<TElementData>::Calculate(const Variable<Matrix>& rVariable,
    Matrix& rOutput, const ProcessInfo& rCurrentProcessInfo) {}

///////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @see QSVMSDEMCoupled::EquationIdVector
 **/
template < class TElementData >
void QSVMSDEMCoupled<TElementData>::EquationIdVector(EquationIdVectorType& rResult,
                                                     ProcessInfo& rCurrentProcessInfo)
{
    if (rCurrentProcessInfo[FRACTIONAL_STEP] == 1) {
        QSVMS<TElementData>::EquationIdVector(rResult, rCurrentProcessInfo);
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
                if(Dim == 3) rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_LAPLACIAN_Z,lappos+2).EquationId();
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
            QSVMS<TElementData>::GetDofList(rElementalDofList, rCurrentProcessInfo);

        }
        else {

            if (rElementalDofList.size() != LocalSize)
                rElementalDofList.resize(LocalSize);

            unsigned int LocalIndex = 0;

            for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
            {
                rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_LAPLACIAN_X);
                rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_LAPLACIAN_Y);
                if (Dim == 3) rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_LAPLACIAN_Z);
            }
        }
}

template< class TElementData >
std::string QSVMSDEMCoupled<TElementData>::Info() const
{
    std::stringstream buffer;
    buffer << "QSVMSDEMCoupled #" << this->Id();
    return buffer.str();
}


template< class TElementData >
void QSVMSDEMCoupled<TElementData>::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "QSVMSDEMCoupled" << Dim << "D";
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

    return DeltaSquared;
}

template <class TElementData>
void QSVMSDEMCoupled<TElementData>::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                                      VectorType& rRightHandSideVector,
                                                      ProcessInfo& rCurrentProcessInfo)
    {
        QSVMS<TElementData>::CalculateLocalSystem(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo);
        if (rCurrentProcessInfo[FRACTIONAL_STEP] == 1) {
            this->CalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo);
        }else{
            noalias(rLeftHandSideMatrix) = ZeroMatrix(LocalSize, LocalSize);
            CalculateLaplacianMassMatrix(rLeftHandSideMatrix, rCurrentProcessInfo);
            this->CalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo);
        }
    }

template <class TElementData>
void QSVMSDEMCoupled<TElementData>::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                                                          ProcessInfo& rCurrentProcessInfo)
{
    // Resize and intialize output
    QSVMS<TElementData>::CalculateLeftHandSide(rLeftHandSideMatrix, rCurrentProcessInfo);
    if (rCurrentProcessInfo[FRACTIONAL_STEP] != 1) {
        noalias(rLeftHandSideMatrix) = ZeroMatrix(LocalSize, LocalSize);
        CalculateLaplacianMassMatrix(rLeftHandSideMatrix, rCurrentProcessInfo);
    }
}

template<class TElementData>
void QSVMSDEMCoupled<TElementData>::CalculateRightHandSide(VectorType& rRightHandSideVector,
                                                           ProcessInfo& rCurrentProcessInfo)
    {
        TElementData data;
        data.Initialize(*this, rCurrentProcessInfo);
        // Calculate this element's geometric parameters
        double Area = data.Weight;
        array_1d<double, 3> convective_velocity;

        // Calculate this element's fluid properties
        QSVMS<TElementData>::CalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo);

        if (rCurrentProcessInfo[FRACTIONAL_STEP] == 1) {

            this->AddMomentumRHS(rRightHandSideVector, Area, data);
            const double& DeltaTime = rCurrentProcessInfo[DELTA_TIME];
            static const double arr[] = {1.0,-1.0};
            std::vector<double> SchemeWeights (arr, arr + sizeof(arr) / sizeof(arr[0]));

            this->AddMassRHS(rRightHandSideVector, data.Density, data.N, Area, SchemeWeights, DeltaTime, data);
        }
        else {
            const unsigned int LocalSize = Dim * NumNodes;

            // Check sizes and initialize
            if (rRightHandSideVector.size() != LocalSize)
                rRightHandSideVector.resize(LocalSize, false);

            noalias(rRightHandSideVector) = ZeroVector(LocalSize);
            this->AddRHSLaplacian(rRightHandSideVector, data.DN_DX, Area);
        }
        if (data.UseOSS == 1.0)

            this->GetAdvectiveVel(convective_velocity, data.N);

            double tau_one, tau_two;
            this->CalculateTau(data, convective_velocity, tau_one, tau_two);

            this->AddProjectionToRHS(rRightHandSideVector, convective_velocity, data, tau_one, tau_two, Area, rCurrentProcessInfo[DELTA_TIME]);
    }

template<class TElementData>
void QSVMSDEMCoupled<TElementData>::AddMassStabilization(TElementData& rData,
                                                         MatrixType &rMassMatrix)
    {

        const double density = this->GetAtCoordinate(rData.Density, rData.N);

        double tau_one;
        double tau_two;
        array_1d<double, 3> convective_velocity=
        this->GetAtCoordinate(rData.Velocity, rData.N) -
        this->GetAtCoordinate(rData.MeshVelocity, rData.N);

        this->CalculateTau(rData, convective_velocity, tau_one, tau_two);

        double K; // Temporary results
        double weight = rData.Weight * tau_one * density; // This density is for the dynamic term in the residual (rho*Du/Dt)
        // If we want to use more than one Gauss point to integrate the convective term, this has to be evaluated once per integration point

        Vector AGradN, AGradNMod;
        this->ConvectionOperator(AGradN, convective_velocity, rData.DN_DX); // Get a * grad(Ni)

        AGradN *= density;

        const double fluid_fraction = this->GetAtCoordinate(rData.FluidFraction, rData.N);

        // Note: Dof order is (vx,vy,[vz,]p) for each node
        for (unsigned int i = 0; i < NumNodes; ++i)
        {
            unsigned int row = i*BlockSize;
            // Loop over columns
            for (unsigned int j = 0; j < NumNodes; ++j)
            {
                unsigned int col = j*BlockSize;
                K = weight * AGradN[i] * rData.N[j];

                for (unsigned int d = 0; d < Dim; ++d) // iterate over dimensions for velocity Dofs in this node combination
                {
                    rMassMatrix(row+d, col+d) += K;
                    rMassMatrix(row+Dim,col+d) += weight * fluid_fraction * rData.DN_DX(i, d) * rData.N[j]; // Delta(u) * TauOne * alpha * Grad(q)

                }
            }
        }
    }

template <class TElementData>
void QSVMSDEMCoupled<TElementData>::GetAdvectiveVel(array_1d< double, 3 > & rAdvVel,
                                                    const typename TElementData::ShapeFunctionsType& rShapeFunc)
    {
        // Compute the weighted value of the advective velocity in the (Gauss) Point
        rAdvVel = rShapeFunc[0] * (this->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY) - this->GetGeometry()[0].FastGetSolutionStepValue(MESH_VELOCITY));

        for (unsigned int iNode = 1; iNode < NumNodes; ++iNode){
            rAdvVel += rShapeFunc[iNode] * (this->GetGeometry()[iNode].FastGetSolutionStepValue(VELOCITY) - this->GetGeometry()[iNode].FastGetSolutionStepValue(MESH_VELOCITY));
        }
    }

template <class TElementData>
void QSVMSDEMCoupled<TElementData>::GetEffectiveViscosity(
    TElementData& rData, double& rViscosity) {

    double c_s = rData.CSmagorinsky;
    //this->GetAtCoordinate(rData.DynamicViscosity, rData.N);
    const double ElementSize = rData.ElementSize;
    rViscosity = rData.DynamicViscosity;

    if (c_s != 0.0) {
        const double density = this->GetAtCoordinate(rData.Density,rData.N);
        const auto& r_velocities = rData.Velocity;
        const auto& r_dndx = rData.DN_DX;
        const double FilterWidth = this->FilterWidth(r_dndx);

        // Calculate Symetric gradient
        MatrixType strain_rate = ZeroMatrix(Dim, Dim);
        for (unsigned int n = 0; n < NumNodes; ++n) {
            for (unsigned int i = 0; i < Dim; ++i)
                for (unsigned int j = 0; j < Dim; ++j)
                    strain_rate(i, j) +=
                        0.5 * (r_dndx(n, j) * r_velocities(n, i) +
                                  r_dndx(n, i) * r_velocities(n, j));
        }

        // Norm of symetric gradient
        double strain_rate_norm = 0.0;
        for (unsigned int i = 0; i < Dim; ++i)
            for (unsigned int j = 0; j < Dim; ++j)
                strain_rate_norm += strain_rate(i, j) * strain_rate(i, j);
        strain_rate_norm = sqrt(2.0 * strain_rate_norm);

        // Nu_sgs = (c_s * Delta)^2 * (2*Sij*Sij)^(1/2)
        rViscosity +=
            2.0 * c_s * c_s * ElementSize * ElementSize * FilterWidth * strain_rate_norm;
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
        TElementData data;
        data.Initialize(*this, rCurrentProcessInfo);
        // Get the element's geometric parameters
        double Area = data.Weight;

        // Calculate this element's fluid properties

        // Add 'classical' mass matrix (lumped)
        double Coeff = Area / NumNodes; //Optimize!

        this->CalculateLaplacianLumpedMassMatrix(rMassMatrix, Coeff);
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
void QSVMSDEMCoupled<TElementData>::AddMomentumRHS(VectorType& F,
                                                   const double Weight,
                                                   TElementData& rData)
    {
        double Coef = rData.Density * Weight;

        array_1d<double,3> body_force = this->GetAtCoordinate(rData.BodyForce, rData.N);

        // Add the results to the velocity components (Local Dofs are vx, vy, [vz,] p for each node)
        int LocalIndex = 0;

        for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
        {
            for (unsigned int d = 0; d < Dim; ++d)
            {
                F[LocalIndex++] += Coef * rData.N[iNode] * body_force[d];
            }
            ++LocalIndex; // Skip pressure Dof
        }
    }

template<class TElementData>
void QSVMSDEMCoupled<TElementData>::AddMassRHS(VectorType& F,
                                               const double Density,
                                               const array_1d<double, NumNodes>& rShapeFunc,
                                               const double Weight,
                                               const std::vector<double>& TimeSchemeWeights,
                                               const double& DeltaTime,
                                               TElementData& rData)
    {
      double FluidFractionRate = this->GetAtCoordinate(rData.FluidFractionRate, rData.N);
      //this->EvaluateTimeDerivativeInPoint(FluidFractionRate, FLUID_FRACTION_RATE, rShapeFunc, DeltaTime, TimeSchemeWeights);
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
                                                       TElementData& rData,
                                                       const double TauOne,
                                                       const double TauTwo,
                                                       const double Weight,
                                                       const double DeltaTime)
    {
        const unsigned int BlockSize = Dim + 1;
        const double density = this->GetAtCoordinate(rData.Density, rData.N);
        Vector AGradN;
        this->ConvectionOperator(AGradN, rAdvVel, rData.DN_DX); // Get a * grad(Ni)
        array_1d<double,Dim> MomProj(Dim,0.0);
        MomProj = this->GetAtCoordinate(rData.MomentumProjection, rData.N);
        double DivProj = 0.0;
        DivProj = this->GetAtCoordinate(rData.MassProjection, rData.N);
        const double fluid_fraction = this->GetAtCoordinate(rData.FluidFraction, rData.N);
        unsigned int FirstRow = 0;

        for (unsigned int i = 0; i < NumNodes; i++)
        {
            array_1d<double,Dim> FluidFractionGradient(Dim,0.0);
            for (unsigned int d = 0; d < Dim; d++)
            {

                FluidFractionGradient[d] += rData.DN_DX(i, d) * rData.FluidFraction[d];
                RHS[FirstRow+d] -= Weight * (density * AGradN[i] * TauOne * MomProj[d] + (fluid_fraction * rData.DN_DX(i,d) + FluidFractionGradient[d] * rData.N[i]) * TauTwo * DivProj); // TauOne * (a * Grad(v)) * MomProjection + TauTwo * Div(v) * MassProjection
                RHS[FirstRow+Dim] -= Weight * rData.DN_DX(i,d) * TauOne * MomProj[d]; // TauOne * Grad(q) * MomProjection
            }
            FirstRow += BlockSize;
        }
    }

template< class TElementData >
void QSVMSDEMCoupled<TElementData>::CalculateTau(
    TElementData& rData,
    const array_1d<double,3> &Velocity,
    double &TauOne,
    double &TauTwo)
{
    constexpr double c1 = 8.0;
    constexpr double c2 = 2.0;

    const double h = rData.ElementSize;
    const double density = this->GetAtCoordinate(rData.Density,rData.N);
    this->GetEffectiveViscosity(rData, rData.EffectiveViscosity);
    const double viscosity = this->GetAtCoordinate(rData.EffectiveViscosity,rData.N);
    double velocity_norm = Velocity[0]*Velocity[0];
    for (unsigned int d = 1; d < Dim; d++)
        velocity_norm += Velocity[d]*Velocity[d];
    velocity_norm = std::sqrt(velocity_norm);

    double inv_tau = c1 * viscosity / (h*h) + density * (rData.DynamicTau/rData.DeltaTime + c2 * velocity_norm / h);
    TauOne = 1.0/inv_tau;

    TauTwo = viscosity + c2 * density * velocity_norm * h / c1;

}

// Add a the contribution from a single integration point to the velocity contribution
template< class TElementData >
void QSVMSDEMCoupled<TElementData>::AddVelocitySystem(TElementData& rData,
                                                      MatrixType &rLocalLHS,
                                                      VectorType &rLocalRHS)
{
    auto& LHS = rData.LHS;
    LHS.clear();

    // Interpolate nodal data on the integration point
    const double density = this->GetAtCoordinate(rData.Density, rData.N);
    array_1d<double,3> body_force = this->GetAtCoordinate(rData.BodyForce,rData.N);
    array_1d<double,3> momentum_projection = this->GetAtCoordinate(rData.MomentumProjection, rData.N);
    double mass_projection = this->GetAtCoordinate(rData.MassProjection, rData.N);

    double tau_one;
    double tau_two;
    array_1d<double, 3> convective_velocity =
        this->GetAtCoordinate(rData.Velocity, rData.N) -
        this->GetAtCoordinate(rData.MeshVelocity, rData.N);

    this->CalculateTau(rData, convective_velocity, tau_one, tau_two);

    Vector AGradN;
    this->ConvectionOperator(AGradN, convective_velocity, rData.DN_DX);

    // Multiplying some quantities by density to have correct units
    body_force *= density; // Force per unit of volume
    AGradN *= density; // Convective term is always multiplied by density

    double PAlphaDivV, GAlpha;
    array_1d<double,3> FluidFractionGradient(3,0.0);
    const double fluid_fraction = this->GetAtCoordinate(rData.FluidFraction, rData.N);

    this->EvaluateGradientOfScalarInPoint(FluidFractionGradient, rData.FluidFraction, rData.DN_DX);

    for (unsigned int i = 0; i < NumNodes; ++i) {
        this->GetGeometry()[i].FastGetSolutionStepValue(FLUID_FRACTION_GRADIENT) = FluidFractionGradient;
    }

    const double fluid_fraction_rate = this->GetAtCoordinate(rData.FluidFractionRate, rData.N);

    // Temporary containers
    double K, G, PDivV;

    // Note: Dof order is (u,v,[w,]p) for each node
    for (unsigned int i = 0; i < NumNodes; i++)
    {
        unsigned int row = i*BlockSize;

        // LHS terms
        for (unsigned int j = 0; j < NumNodes; j++)
        {
            unsigned int col = j*BlockSize;

            // Some terms are the same for all velocity components, calculate them once for each i,j
            K = rData.N[i] * AGradN[j];
            K += AGradN[i] * tau_one * (AGradN[j]); // Stabilization: u*grad(v) * tau_one * u*grad(u)
            K *= rData.Weight;

            // q-p stabilization block (initialize result)
            double laplacian = 0;

            for (unsigned int d = 0; d < Dim; d++)
            {
                LHS(row+d,col+d) += K;

                // v * Grad(p) block
                G = tau_one * AGradN[i] * rData.DN_DX(j,d); // Stabilization: (a * Grad(v)) * tau_one * Grad(p)
                GAlpha = tau_one * AGradN[i] * (fluid_fraction * rData.DN_DX(j, d));
                PAlphaDivV = (fluid_fraction * rData.DN_DX(i, d) + FluidFractionGradient[d] * rData.N[i]) *  rData.N[j]; // alpha * q * Div(u) + q * Grad(alpha) * u

                PDivV = rData.DN_DX(i,d) * rData.N[j]; // Div(v) * p

                // Write v * Grad(p) component
                LHS(row+d,col+Dim) += rData.Weight * (G - PDivV);
                // Use symmetry to write the q * Div(u) component

                LHS(col+Dim,row+d) += rData.Weight * (GAlpha + PAlphaDivV); // note that here PAlphaDivV includes G; do not look for it in GAlpha!

                // q-p stabilization block
                laplacian += fluid_fraction * rData.DN_DX(i,d) * rData.DN_DX(j,d); // Stabilization: Grad(q) * tau_one * Grad(p)

                for (unsigned int e = 0; e < Dim; e++) // Stabilization: Div(v) * tau_two * Div(u)
                    LHS(row+d,col+e) += rData.Weight * tau_two * rData.DN_DX(i,d) * (fluid_fraction * rData.DN_DX(j,e) + FluidFractionGradient[e] * rData.N[j]);
            }

            // Write q-p term
            LHS(row+Dim,col+Dim) += rData.Weight * tau_one * laplacian;
        }

        // RHS terms
        double forcing = 0.0;
        for (unsigned int d = 0; d < Dim; ++d)
        {
            rLocalRHS[row+d] += rData.Weight * rData.N[i] * body_force[d]; // v*BodyForce
            rLocalRHS[row+d] += rData.Weight * tau_one * AGradN[i] * ( body_force[d] - momentum_projection[d]); // ( a * Grad(v) ) * tau_one * (Density * BodyForce)
            rLocalRHS[row+d] -= rData.Weight * tau_two * rData.DN_DX(i,d) * mass_projection * fluid_fraction_rate;
            forcing += fluid_fraction * rData.DN_DX(i,d) * (body_force[d] - momentum_projection[d]);
        }
        rLocalRHS[row + Dim] += rData.Weight * tau_one * forcing; // Grad(q) * tau_one * (Density * BodyForce)
    }

    // Write (the linearized part of the) local contribution into residual form (A*dx = b - A*x)
    array_1d<double,LocalSize> values;
    this->GetCurrentValuesVector(rData,values);
    noalias(rLocalRHS) -= prod(LHS, values);
    /* Viscous contribution (with symmetric gradient 2*nu*{E(u) - 1/3 Tr(E)} )
     * For a generic (potentially non-linear) constitutive law, one cannot assume that RHS = F - LHS*current_values.
     * Because of this, the AddViscousTerm function manages both the LHS and the RHS.
     */
    this->AddViscousTerm(rData,LHS,rLocalRHS);

    noalias(rLocalLHS) += LHS;

}

template<class TElementData>
void QSVMSDEMCoupled<TElementData>::GetModifiedConvectionOperator(array_1d< double,  NumNodes >& rResult,
                                                                  array_1d< double, 3 > & rVelocity,
                                                                  const double & rVelocityDiv,
                                                                  const typename TElementData::ShapeFunctionsType& rShapeFunc,
                                                                  const typename TElementData::ShapeDerivativesType& rShapeDeriv)
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
void QSVMSDEMCoupled<TElementData>::MassProjTerm(const TElementData& rData,
                                                 double &rMassRHS) const
    {
        const auto& r_velocities = rData.Velocity;
        array_1d<double,3> FluidFractionGradient(3,0.0);

        const double fluid_fraction = this->GetAtCoordinate(rData.FluidFraction, rData.N);

        this->EvaluateGradientOfScalarInPoint(FluidFractionGradient, rData.FluidFraction, rData.DN_DX);

        const double fluid_fraction_rate = this->GetAtCoordinate(rData.FluidFractionRate, rData.N);
        // Compute this node's contribution to the residual (evaluated at inegration point)
        for (unsigned int i = 0; i < NumNodes; i++) {
            for (unsigned int d = 0; d < Dim; ++d)
            {
                rMassRHS -= (fluid_fraction * rData.DN_DX(i, d) * r_velocities(i, d)) + FluidFractionGradient[d] * rData.N[i] * r_velocities(i, d);
            }
        }
        rMassRHS -= fluid_fraction_rate;
    }

template<class TElementData>
void QSVMSDEMCoupled<TElementData>::EvaluateGradientOfScalarInPoint(array_1d< double, 3 >& rResult,
                                                                    const typename TElementData::NodalScalarData& variable,
                                                                    const typename TElementData::ShapeDerivativesType& rShapeDeriv) const
    {

        for (unsigned int i = 0; i < NumNodes; ++i) {
            const double& scalar = variable[i];

            for (unsigned int d = 0; d < Dim; ++d){
                rResult[d] += rShapeDeriv(i, d) * scalar;
            }
        }
    }

    /// Write the advective velocity evaluated at this point to an array
    /**
     * Writes the value of the advective velocity evaluated at a point inside
     * the element to an array_1d
     * @param rAdvVel: Output array
     * @param rShapeFunc: Shape functions evaluated at the point of interest
     */
    /// Write the advective velocity evaluated at this point to an array
    /**
     * Writes the value of the advective velocity evaluated at a point inside
     * the element to an array_1d
     * @param rAdvVel: Output array
     * @param rShapeFunc: Shape functions evaluated at the point of interest
     * @param Step: The time Step
     */

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

