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

void CalculateReynoldsStressTensor(MatrixType& rReynoldsStressTensor,
                                   const ProcessInfo& rCurrentProcessInfo)
{


}

void CalculateReynoldsStressTensorLinearEddyViscosity(MatrixType& rReynoldsStressTensor,
                                                      const ProcessInfo& rCurrentProcessInfo)
{
    if (rReynoldsStressTensor.size1() != TDim || rReynoldsStressTensor.size2() != TDim)
        rReynoldsStressTensor.resize(TDim,TDim, false);


    const double nu_t = this->GetValue(TURBULENT_KINEMATIC_VISCOSITY);
    const double turbulent_kinetic_energy = this->GetValue(TURBULENT_KINETIC_ENERGY);

    ShapeFunctionDerivativesType DN_DX;
    array_1d< double, TNumNodes > N;
    double Volume;

    GeometryUtils::CalculateGeometryData(this->GetGeometry(),DN_DX,N,Volume);

    BoundedMatrix< double, TDim, TDim > GradVel;
    this->CalculateVelocityGradient(GradVel,DN_DX);

    const BoundedMatrix< double, TDim, TDim > mSymGradVel = 0.5*(GradVel) + 0.5*trans(GradVel);
    const IdentityMatrix<TDim> mI(TDim);
    double SymGradVelTrace = 0.0;

    for (unsigned int i=0; i < TDim; i++)
        SymGradVelTrace += GradVel(i,i);

    noalias(rReynoldsStressTensor) = -2*nu_t*(mSymGradVel - (1/3)*mSymGradVel*mI) + (2/3)*turbulent_kinetic_energy*mI;
}

void CalculateReynoldsStressTensorLeftHandSideContributionLinearEddyViscosity(
    MatrixType& rAdjointMatrix, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if (rAdjointMatrix.size1() != TFluidLocalSize || rAdjointMatrix.size2() != TFluidLocalSize)
        rAdjointMatrix.resize(TFluidLocalSize,TFluidLocalSize,false);

    rAdjointMatrix.clear();


    const double nu_t = this->GetValue(TURBULENT_KINEMATIC_VISCOSITY);
    const double tke = this->GetValue(TURBULENT_KINETIC_ENERGY);

    ShapeFunctionDerivativesType DN_DX;
    array_1d< double, TNumNodes > N;
    double Volume;

    GeometryUtils::CalculateGeometryData(this->GetGeometry(),DN_DX,N,Volume);

    // Density
    double Density;
    this->EvaluateInPoint(Density,DENSITY,N);

    IndexType FirstRow(0), FirstCol(0);
    // Loop over nodes
    for (IndexType i = 0; i < TNumNodes; ++i)
    {
        for (IndexType j = 0; j < TNumNodes; ++j)
        {
            for (IndexType m = 0; m < TDim; ++m)
            {
                for (IndexType n = 0; n < TDim; ++n)
                {
                    rAdjointMatrix(FirstRow+m,FirstCol+n) += Volume * valmn;
                }

                rAdjointMatrix(FirstRow+m,FirstCol+m) += Volume * diag;
                rAdjointMatrix(FirstRow+m,FirstCol+TDim) += Volume * valmp;
                rAdjointMatrix(FirstRow+TDim,FirstCol+m) += Volume * valpn;
            }

            FirstCol += TBlockSize;
        }  // Node block columns

        FirstRow += TBlockSize;
        FirstCol = 0;
    }  // Node block rows

    KRATOS_CATCH("")
}

