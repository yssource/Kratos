void ComputeVelocityGradientTensor(const ProcessInfo& rCurrentProcessInfo)
{
    ShapeFunctionDerivativesType DN_DX;
    array_1d< double, TNumNodes > N;
    double Volume;

    GeometryUtils::CalculateGeometryData(this->GetGeometry(),DN_DX,N,Volume);

    BoundedMatrix< double, TDim, TDim > DensityGradVel;
    this->CalculateVelocityGradient(DensityGradVel,DN_DX);

    this->SetValue(MATRIX_VELOCITY_GRADIENT_TENSOR, DensityGradVel);
}

void ComputeVelocityDidacticProduct(const ProcessInfo& rCurrentProcessInfo)
{
    ShapeFunctionDerivativesType DN_DX;
    array_1d< double, TNumNodes > N;
    double Volume;

    GeometryUtils::CalculateGeometryData(this->GetGeometry(),DN_DX,N,Volume);

    array_1d< double, TDim > Velocity;
    this->EvaluateInPoint(Velocity,VELOCITY,N);

    BoundedMatrix< double, TDim, TDim > VelocityDidacticProduct;
    for (unsigned int i=0; i < TDim; i++)
        for (unsigned int j=0; j < TDim; j++)
            VelocityDidacticProduct(i,j) = Velocity[i]*Velocity[j];

    this->SetValue(MATRIX_VELOCITY_DIDACTIC_PRODUCT, VelocityDidacticProduct);
}

void ComputePressureGradientTensor(const ProcessInfo& rCurrentProcessInfo)
{
    ShapeFunctionDerivativesType DN_DX;
    array_1d< double, TNumNodes > N;
    double Volume;

    GeometryUtils::CalculateGeometryData(this->GetGeometry(),DN_DX,N,Volume);

    // Grad(p)
    array_1d< double, TDim > GradP;
    this->CalculatePressureGradient(GradP,DN_DX);

    this->SetValue(VECTOR_PRESSURE_GRADIENT, GradP);
}

void ProcessMatrices(const ProcessInfo& rCurrentProcessInfo)
{
    ComputeVelocityGradientTensor(rCurrentProcessInfo);
    ComputeVelocityDidacticProduct(rCurrentProcessInfo);
    ComputePressureGradientTensor(rCurrentProcessInfo);
}

// void CalculateReynoldsStressTensor(MatrixType& rReynoldsStressTensor,
//                                    const ProcessInfo& rCurrentProcessInfo)
// {
//     CalculateReynoldsStressTensorLinearEddyViscosity(rReynoldsStressTensor, rCurrentProcessInfo);
// }

// void CalculateReynoldsStressTensorFirstDerivativeLHS(MatrixType& rAdjointMatrix,
//                                                      const ProcessInfo& rCurrentProcessInfo)
// {
//     CalculateReynoldsStressTensorFirstDerivativeLHSLinearEddyViscosity(
//         rAdjointMatrix, rCurrentProcessInfo);
// }

// void CalculateReynoldsStressTensorShapeGradient(MatrixType& rAdjointMatrix,
//                                                 const ProcessInfo& rCurrentProcessInfo)
// {
//     CalculateReynoldsStressTensorShapeGradientLinearEddyViscosity(
//         rAdjointMatrix, rCurrentProcessInfo);
// }

// void CalculateReynoldsStressTensorLinearEddyViscosity(MatrixType& rReynoldsStressTensor,
//                                                       const ProcessInfo& rCurrentProcessInfo)
// {
//     KRATOS_TRY

//     if (rReynoldsStressTensor.size1() != TDim || rReynoldsStressTensor.size2() != TDim)
//         rReynoldsStressTensor.resize(TDim,TDim, false);


//     const double nu_t = this->GetValue(TURBULENT_KINEMATIC_VISCOSITY);
//     const double turbulent_kinetic_energy = this->GetValue(TURBULENT_KINETIC_ENERGY);

//     ShapeFunctionDerivativesType DN_DX;
//     array_1d< double, TNumNodes > N;
//     double Volume;

//     GeometryUtils::CalculateGeometryData(this->GetGeometry(),DN_DX,N,Volume);

//     BoundedMatrix< double, TDim, TDim > GradVel;
//     this->CalculateVelocityGradient(GradVel,DN_DX);

//     const BoundedMatrix< double, TDim, TDim > mSymGradVel = 0.5*(GradVel) + 0.5*trans(GradVel);
//     const IdentityMatrix<TDim> mI(TDim);
//     double SymGradVelTrace = 0.0;

//     for (unsigned int i=0; i < TDim; i++)
//         SymGradVelTrace += GradVel(i,i);

//     noalias(rReynoldsStressTensor) =
//         -2 * nu_t * (mSymGradVel - (1 / TDim) * SymGradVelTrace * mI) +
//         (2 / TDim) * turbulent_kinetic_energy * mI;

//     KRATOS_CATCH("")
// }

// void CalculateReynoldsStressTensorFirstDerivativeLHSLinearEddyViscosity(
//     MatrixType& rAdjointMatrix, const ProcessInfo& rCurrentProcessInfo)
// {
//     KRATOS_TRY

//     if (rAdjointMatrix.size1() != TFluidLocalSize || rAdjointMatrix.size2() != TFluidLocalSize)
//         rAdjointMatrix.resize(TFluidLocalSize, TFluidLocalSize, false);

//     rAdjointMatrix.clear();

//     const double nu_t = this->GetValue(TURBULENT_KINEMATIC_VISCOSITY);
//     const double tke = this->GetValue(TURBULENT_KINETIC_ENERGY);

//     ShapeFunctionDerivativesType DN_DX;
//     array_1d<double, TNumNodes> N;
//     double Volume;

//     GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DX, N, Volume);

//     // Density
//     double Density;
//     this->EvaluateInPoint(Density, DENSITY, N);

//     // u
//     array_1d< double, TDim > Velocity;
//     this->EvaluateInPoint(Velocity,VELOCITY,N);
//     const double turbulent_intensity = tke/(3.0*norm_2(Velocity)^2);

//     BoundedMatrix<double, TNumNodes, TNumNodes> NodalShapeGradients;
//     NodalShapeGradients.clear();

//     for (IndexType a; a < TNumNodes; ++a)
//         for (IndexType c; c < TNumNodes; ++c)
//             for (IndexType l; l < TDim; ++l)
//                 NodalShapeGradients(a,c) += DN_DX(a,l)*DN_DX(c,l);


//     IndexType FirstRow(0), FirstCol(0);
//     // Loop over nodes
//     for (IndexType a = 0; a < TNumNodes; ++a)
//     {
//         for (IndexType c = 0; c < TNumNodes; ++c)
//         {
//             for (IndexType i = 0; i < TDim; ++i)
//             {
//                 double diag = -(nu_t * NodalShapeGradients(a,c));

//                 for (IndexType k = 0; k < TDim; ++k)
//                 {
//                     double valik = 0.0;

//                     valik -= nu_t*DN_DX(a,k)*DN_DX(c,i);
//                     valik += (2/TDim)*nu_t*DN_DX(c,k)*DN_DX(a,i);
//                     valik += (6/TDim)*(turbulent_intensity^2)*N[c]*Velocity[k]*DN_DX(a,i);
//                     rAdjointMatrix(FirstRow + i, FirstCol + k) += Volume * Density * valik;
//                 }
//                 rAdjointMatrix(FirstRow + i, FirstCol + i) += Volume * Density * diag;
//             }
//             FirstCol += TBlockSize;
//         } // Node block columns
//         FirstRow += TBlockSize;
//         FirstCol = 0;
//     } // Node block rows

//     KRATOS_CATCH("")
// }

// void CalculateReynoldsStressTensorShapeGradientLinearEddyViscosity(MatrixType& rShapeDerivativesMatrix,
//                                                                    const ProcessInfo& rCurrentProcessInfo)
// {
//     KRATOS_TRY

//         if (rShapeDerivativesMatrix.size1() != TCoordLocalSize
//                 || rShapeDerivativesMatrix.size2() != TFluidLocalSize)
//         {
//             rShapeDerivativesMatrix.resize(TCoordLocalSize,TFluidLocalSize,false);
//         }

//         // Get shape functions, shape function gradients and element volume (area in
//         // 2D). Only one integration point is used so the volume is its weight.
//         ShapeFunctionDerivativesType DN_DX;
//         array_1d< double, TNumNodes > N;
//         double Volume;

//         GeometryUtils::CalculateGeometryData(this->GetGeometry(),DN_DX,N,Volume);

//         // Density
//         double Density;
//         this->EvaluateInPoint(Density,DENSITY,N);

//         // Dynamic viscosity
//         double Viscosity;
//         this->EvaluateInPoint(Viscosity,VISCOSITY,N);
//         Viscosity *= Density;

//         // u
//         array_1d< double, TDim > Velocity;
//         this->EvaluateInPoint(Velocity,VELOCITY,N);

//         // u * Grad(N)
//         array_1d< double, TNumNodes > DensityVelGradN;
//         noalias(DensityVelGradN) = Density * prod(DN_DX,Velocity);

//         // Det(J)
//         const double InvDetJ = 1.0 / this->GetGeometry().DeterminantOfJacobian(0);
//         array_1d< double, TCoordLocalSize > DetJDerivatives;
//         this->CalculateDeterminantOfJacobianDerivatives(DetJDerivatives);

//         // Stabilization parameters TauOne, TauTwo
//         double VelNorm = norm_2(Velocity);
//         double ElemSize = this->CalculateElementSize(Volume);
//         double TauOne, TauTwo;
//         this->CalculateStabilizationParameters(TauOne,TauTwo,VelNorm,ElemSize,
//                 Density,Viscosity,rCurrentProcessInfo);

//         // External body force
//         array_1d< double, TDim > BodyForce;
//         this->EvaluateInPoint(BodyForce,BODY_FORCE,N);
//         BodyForce *= Density;

//         array_1d< double, TFluidLocalSize > FluidValues;

//         IndexType DofIndex = 0;
//         for (IndexType iNode = 0; iNode < TNumNodes; iNode++)
//         {
//             array_1d< double, 3 >& rVelocity = this->GetGeometry()[iNode].FastGetSolutionStepValue(VELOCITY);
//             for (IndexType d = 0; d < TDim; d++)
//                 FluidValues[DofIndex++] = rVelocity[d];
//             FluidValues[DofIndex++] = this->GetGeometry()[iNode].FastGetSolutionStepValue(PRESSURE);
//         }

//         // We compute the derivative of the residual w.r.t each coordinate of each
//         // node and assign it to the corresponding row of the shape derivatives
//         // matrix.
//         for (IndexType iCoord = 0; iCoord < TCoordLocalSize; ++iCoord)
//         {
//             // Det(J)'
//             double DetJDeriv = DetJDerivatives[iCoord];

//             // DN_DX'
//             BoundedMatrix< double, TNumNodes, TDim > DN_DX_Deriv;
//             for (IndexType i = 0; i < TNumNodes; ++i)
//                 for (IndexType d = 0; d < TDim; ++d)
//                     DN_DX_Deriv(i,d) = -DN_DX(iCoord / TDim,d) * DN_DX(i,iCoord % TDim);

//             // Volume'
//             double VolumeDeriv = Volume * InvDetJ * DetJDeriv;

//             // u * Grad(N)'
//             array_1d< double, TNumNodes > DensityVelGradNDeriv;
//             noalias(DensityVelGradNDeriv) = Density * prod(DN_DX_Deriv,Velocity);

//             // TauOne', TauTwo'
//             double TauOneDeriv, TauTwoDeriv;
//             this->CalculateStabilizationParametersDerivative(
//                     TauOneDeriv,TauTwoDeriv,TauOne,TauTwo,VelNorm,ElemSize,Density,
//                     Viscosity,DetJDeriv);

//             BoundedMatrix<double, TFluidLocalSize, TFluidLocalSize> LHS;
//             array_1d< double, TFluidLocalSize > RHS;
//             for (IndexType i = 0; i < TFluidLocalSize; i++)
//             {
//                 RHS[i] = 0.0;
//                 for (IndexType j = 0; j < TFluidLocalSize; j++)
//                     LHS(i,j) = 0.0;
//             }

//             for (IndexType i = 0; i < TNumNodes; ++i)
//             {
//                 for (IndexType j = 0; j < TNumNodes; ++j)
//                 {
//                     // Left-hand side matrix
//                     double diag = 0.0;
//                     double ddiag = 0.0;

//                     // Convective term, v * (u * Grad(u))
//                     diag += N[i] * DensityVelGradN[j];
//                     ddiag += N[i] * DensityVelGradNDeriv[j];

//                     // Stabilization, lsq convection
//                     // (u * Grad(v)) * TauOne * (u * Grad(u))
//                     diag += DensityVelGradN[i] * TauOne * DensityVelGradN[j];
//                     ddiag += DensityVelGradNDeriv[i] * TauOne * DensityVelGradN[j]
//                     + DensityVelGradN[i] * TauOneDeriv * DensityVelGradN[j]
//                     + DensityVelGradN[i] * TauOne * DensityVelGradNDeriv[j];

//                     for (IndexType m = 0; m < TDim; ++m)
//                     {
//                         for (IndexType n = 0; n < TDim; ++n)
//                         {
//                             // Stabilization, lsq divergence
//                             // Div(v) * TauTwo * Div(u)
//                             double valmn = DN_DX(i,m) * TauTwo * DN_DX(j,n);
//                             double dvalmn = DN_DX_Deriv(i,m) * TauTwo * DN_DX(j,n)
//                             + DN_DX(i,m) * TauTwoDeriv * DN_DX(j,n)
//                             + DN_DX(i,m) * TauTwo * DN_DX_Deriv(j,n);

//                             LHS(i*TBlockSize+m,j*TBlockSize+n) += VolumeDeriv*valmn
//                             + Volume*dvalmn;
//                         }
//                         LHS(i*TBlockSize+m,j*TBlockSize+m) += VolumeDeriv * diag
//                         + Volume * ddiag;

//                         double valmp = 0.0;
//                         double dvalmp = 0.0;
//                         // Pressure term
//                         // Div(v) * p
//                         valmp -= DN_DX(i,m) * N[j];
//                         dvalmp -= DN_DX_Deriv(i,m) * N[j];

//                         // Stabilization, convection-pressure
//                         // (u * Grad(v)) * TauOne * Grad(p)
//                         valmp += TauOne * DensityVelGradN[i] * DN_DX(j,m);
//                         dvalmp += TauOneDeriv * DensityVelGradN[i] * DN_DX(j,m)
//                         + TauOne * DensityVelGradNDeriv[i] * DN_DX(j,m)
//                         + TauOne * DensityVelGradN[i] * DN_DX_Deriv(j,m);

//                         double valpn = 0.0;
//                         double dvalpn = 0.0;
//                         // Divergence term
//                         // q * Div(u)
//                         valpn += N[i] * DN_DX(j,m);
//                         dvalpn += N[i] * DN_DX_Deriv(j,m);

//                         // Stabilization, pressure-convection
//                         // Grad(q) * TauOne * (u * Grad(u))
//                         valpn += TauOne * DensityVelGradN[j] * DN_DX(i,m);
//                         dvalpn += TauOneDeriv * DensityVelGradN[j] * DN_DX(i,m)
//                         + TauOne * DensityVelGradNDeriv[j] * DN_DX(i,m)
//                         + TauOne * DensityVelGradN[j] * DN_DX_Deriv(i,m);

//                         LHS(i*TBlockSize+m,j*TBlockSize+TDim) += VolumeDeriv * valmp
//                         + Volume * dvalmp;
//                         LHS(i*TBlockSize+TDim,j*TBlockSize+m) += VolumeDeriv * valpn
//                         + Volume * dvalpn;
//                     }

//                     double valpp = 0.0;
//                     double dvalpp = 0.0;
//                     // Stabilization, lsq pressure
//                     // TauOne * Grad(q) * Grad(p)
//                     for (IndexType d = 0; d < TDim; ++d)
//                     {
//                         valpp += DN_DX(i,d) * DN_DX(j,d) * TauOne;
//                         dvalpp += DN_DX_Deriv(i,d) * DN_DX(j,d) * TauOne
//                         + DN_DX(i,d) * DN_DX_Deriv(j,d) * TauOne
//                         + DN_DX(i,d) * DN_DX(j,d) * TauOneDeriv;
//                     }

//                     LHS(i*TBlockSize+TDim,j*TBlockSize+TDim) += VolumeDeriv * valpp
//                     + Volume * dvalpp;
//                 } // Node block columns

//                 // Right-hand side vector
//                 double DN_DX_BodyForce = 0.0;
//                 double DN_DX_BodyForceDeriv = 0.0;
//                 for (IndexType d = 0; d < TDim; ++d)
//                 {
//                     DN_DX_BodyForce += DN_DX(i,d) * BodyForce[d];
//                     DN_DX_BodyForceDeriv += DN_DX_Deriv(i,d) * BodyForce[d];
//                 }

//                 for (IndexType m = 0; m < TDim; ++m)
//                 {
//                     double valm = 0.0;
//                     double dvalm = 0.0;

//                     // External body force
//                     valm += N[i] * BodyForce[m];

//                     // Stabilization, convection-BodyForce
//                     // (u * Grad(v)) * TauOne * f
//                     valm += TauOne * DensityVelGradN[i] * BodyForce[m];
//                     dvalm += TauOneDeriv * DensityVelGradN[i] * BodyForce[m]
//                     + TauOne * DensityVelGradNDeriv[i] * BodyForce[m];

//                     RHS[i*TBlockSize+m] += VolumeDeriv * valm + Volume * dvalm;
//                 }

//                 double valp = TauOne * DN_DX_BodyForce;
//                 double dvalp = TauOneDeriv * DN_DX_BodyForce
//                 + TauOne * DN_DX_BodyForceDeriv;

//                 RHS[i*TBlockSize+TDim] += VolumeDeriv * valp + Volume * dvalp;
//             } // Node block rows

//             this->AddViscousTermDerivative(LHS,DN_DX,DN_DX_Deriv,Viscosity * Volume,
//                     Viscosity * VolumeDeriv);

//             // Assign the derivative of the residual w.r.t this coordinate to the
//             // shape derivatives matrix.
//             array_1d< double, TFluidLocalSize > ResidualDerivative;
//             noalias(ResidualDerivative) = RHS - prod(LHS,FluidValues);
//             for (IndexType k = 0; k < TFluidLocalSize; ++k)
//                 rShapeDerivativesMatrix(iCoord,k) = ResidualDerivative[k];
//         }

//         KRATOS_CATCH("")
// }
