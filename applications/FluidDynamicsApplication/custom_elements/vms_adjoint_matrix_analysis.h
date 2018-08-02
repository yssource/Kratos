// Warning: Alpha Bossak value is hardcoded in this file

double CalculateMatrixEnergy( const Vector& rVector1,
                              const Element::MatrixType& rMatrix,
                              const Vector& rVector2 )
{
    KRATOS_ERROR_IF(rVector1.size() != rMatrix.size1())<<"Error calculating matrix energy. Vector 1 size [" << rVector1.size() << "] and Matrix size [" << rMatrix.size1() << "] doesn't match."<<std::endl;
    KRATOS_ERROR_IF(rVector2.size() != rMatrix.size2())<<"Error calculating matrix energy. Vector 2 size [" << rVector2.size() << "] and Matrix size [" << rMatrix.size2() << "] doesn't match."<<std::endl;

    const Vector& temp = prod(trans(rMatrix), rVector1);
    double energy = inner_prod(temp, rVector2);

    return energy;
}

void CalculateEigenValues( const Element::MatrixType& rMatrix,
                           double& eigen_min,
                           double& eigen_max )
{
    if (rMatrix.size1() != rMatrix.size2())
        return;

    const unsigned int matrix_size = rMatrix.size1();

    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>  symmetric_matrix;
    symmetric_matrix.resize(matrix_size, matrix_size);

    for (unsigned int i = 0; i < matrix_size; i++)
        for (unsigned int j = i; j < matrix_size; j++)
        {
            double value = 0.5*rMatrix(i,j) + 0.5*rMatrix(j,i);
            symmetric_matrix(i,j) = value;
            symmetric_matrix(j,i) = value;
        }

    Eigen::EigenSolver<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> es(symmetric_matrix, false);
    const auto& values = es.eigenvalues().real();

    eigen_min = values[0];
    eigen_max = values[0];

    double temp;

    for (unsigned int i = 0; i < values.size(); i++)
    {
        temp = values[i];
        if (eigen_min > temp)
            eigen_min = temp;

        if (eigen_max < temp)
            eigen_max = temp;
    }
}

void CalculateTau( double& rTauOne,
                   double& rTauTwo,
                   const ProcessInfo& rCurrentProcessInfo )
{
    BoundedMatrix<double, TNumNodes, TDim> shape_function_derivatives;
    array_1d< double, TNumNodes > shape_function;
    double weight;

    GeometryUtils::CalculateGeometryData(this->GetGeometry(), shape_function_derivatives, shape_function, weight);

    // u
    array_1d< double, TDim > velocity;
    this->EvaluateInPoint(velocity, VELOCITY, shape_function);

    // Density
    double density;
    this->EvaluateInPoint(density, DENSITY, shape_function);

    // Dynamic viscosity
    double viscosity;
    this->EvaluateInPoint(viscosity, VISCOSITY, shape_function);
    viscosity *= density;

    double elem_size = this->CalculateElementSize(weight);

    // Compute mean advective velocity norm
    double adv_vel_norm = 0.0;
    for (unsigned int d = 0; d < TDim; ++d)
        adv_vel_norm += velocity[d] * velocity[d];

    adv_vel_norm = sqrt(adv_vel_norm);

    this->CalculateStabilizationParameters(rTauOne, rTauTwo, adv_vel_norm, elem_size,density, viscosity,rCurrentProcessInfo);
}


void CalculateVelocityNormDerivatives( BoundedMatrix<double, TNumNodes, TDim>& rVelocityNormDerivs )
{
    BoundedMatrix<double, TNumNodes, TDim> shape_function_derivatives;
    array_1d< double, TNumNodes > shape_function;
    double weight;

    GeometryUtils::CalculateGeometryData(this->GetGeometry(), shape_function_derivatives, shape_function, weight);

    // u
    array_1d< double, TDim > velocity;
    this->EvaluateInPoint(velocity,VELOCITY, shape_function);

    const double velocity_norm = norm_2(velocity);

    if (velocity_norm > 0.0)
        for (unsigned int i=0; i < TNumNodes; ++i)
            for (unsigned int j=0; j < TDim; ++j)
                rVelocityNormDerivs(i,j) = velocity[j]*shape_function[i]/velocity_norm;
}

void CalculateTau1Derivatives( BoundedMatrix<double, TNumNodes, TDim>& rTau1Derivs,
                               const ProcessInfo& rCurrentProcessInfo )
{
    BoundedMatrix<double, TNumNodes, TDim> shape_function_derivatives;
    array_1d< double, TNumNodes > shape_function;
    double weight;

    GeometryUtils::CalculateGeometryData(this->GetGeometry(), shape_function_derivatives, shape_function, weight);

    // Density
    double density;
    this->EvaluateInPoint(density, DENSITY, shape_function);

    double elem_size = this->CalculateElementSize(weight);

    double tau_one(0.0), tau_two(0.0);

    CalculateTau( tau_one, tau_two, rCurrentProcessInfo);

    BoundedMatrix<double, TNumNodes, TDim> velocity_norm_derivs;
    CalculateVelocityNormDerivatives( velocity_norm_derivs);
    double coeff = -2*density*tau_one*tau_one/elem_size;

    for (unsigned int i=0; i < TNumNodes; ++i)
    {
        for (unsigned int j=0; j<TDim; ++j)
        {
            rTau1Derivs(i,j) = coeff*velocity_norm_derivs(i,j);
        }
    }
}

void CalculateTau2Derivatives( BoundedMatrix<double, TNumNodes, TDim>& rTau2Derivs,
                               const ProcessInfo& rCurrentProcessInfo )
{
    BoundedMatrix<double, TNumNodes, TDim> shape_function_derivatives;
    array_1d< double, TNumNodes > shape_function;
    double weight;

    GeometryUtils::CalculateGeometryData(this->GetGeometry(), shape_function_derivatives, shape_function, weight);

    // Density
    double density;
    this->EvaluateInPoint(density, DENSITY, shape_function);

    double elem_size = this->CalculateElementSize(weight);

    BoundedMatrix<double, TNumNodes, TDim> velocity_norm_derivs;
    CalculateVelocityNormDerivatives( velocity_norm_derivs);

    double coeff = 0.5*density*elem_size;

    for (unsigned int i=0; i < TNumNodes; ++i)
    {
        for (unsigned int j=0; j<TDim; ++j)
        {
            rTau2Derivs(i,j) = coeff*velocity_norm_derivs(i,j);
        }
    }
}

void CalculateBossakSchemeConstants( double& rAlphaBossak,
                                     double& rBetaNewmark,
                                     double& rGammaNewmark,
                                     double& rDeltaTime,
                                     const ProcessInfo& rCurrentProcessInfo )
{
    rAlphaBossak = rCurrentProcessInfo[BOSSAK_ALPHA];
    KRATOS_ERROR_IF(rAlphaBossak == 0.0)<<" BOSSAK_ALPHA is zero.";
    rBetaNewmark = 0.25 * (1.0 - rAlphaBossak) * (1.0 - rAlphaBossak);
    rGammaNewmark = 0.5 - rAlphaBossak;
    rDeltaTime = -rCurrentProcessInfo[DELTA_TIME];
}


void GetConvectionOperator( array_1d< double, TNumNodes >& rResult,
                            const array_1d< double, TDim > & rVelocity,
                            const BoundedMatrix<double, TNumNodes, TDim >& rShapeDeriv)
{
    // Evaluate (and weight) the a * Grad(Ni) operator in the integration point, for each node i
    for (unsigned int iNode = 0; iNode < TNumNodes; ++iNode) // Loop over nodes
    {
        // Initialize result
        rResult[iNode] = rVelocity[0] * rShapeDeriv(iNode, 0);
        for (unsigned int d = 1; d < TDim; ++d) // loop over components
            rResult[iNode] += rVelocity[d] * rShapeDeriv(iNode, d);
    }
}

void GetPostVelocityGradient( BoundedMatrix<double, TDim, TDim >& rResult,
                              const BoundedMatrix<double, TNumNodes, TDim>& rShapeFunctionDerivs )
{
    rResult.clear();

    for (unsigned int i=0; i < TDim; ++i)
        for (unsigned int j=0; j<TDim; ++j)
            for (unsigned b=0; b<TNumNodes; ++b)
            {
                const array_1d<double, 3>& velocity = this->GetGeometry()[b].FastGetSolutionStepValue(VELOCITY,0);
                rResult(i,j) += rShapeFunctionDerivs(b,i)*velocity[j];
            }
}

void GetPostPressureGradient( BoundedVector<double, TDim >& rResult,
                              const BoundedMatrix<double, TNumNodes, TDim>& rShapeFunctionDerivs )
{
    rResult.clear();

    for (unsigned int i=0; i < TDim; ++i)
        for (unsigned b=0; b<TNumNodes; ++b)
        {
            const double pressure = this->GetGeometry()[b].FastGetSolutionStepValue(PRESSURE,0);
            rResult[i] += rShapeFunctionDerivs(b,i)*pressure;
        }
}

void Matrix_1( Element::MatrixType& rMatrix,
               const ProcessInfo& rCurrentProcessInfo )
{


    double alpha, beta, gamma, delta_time;
    CalculateBossakSchemeConstants(alpha, beta, gamma, delta_time, rCurrentProcessInfo);

    BoundedMatrix<double, TNumNodes, TDim> shape_function_derivatives;
    array_1d< double, TNumNodes > shape_function;
    double weight;

    const double coeff = gamma/(beta*delta_time);

    GeometryUtils::CalculateGeometryData(this->GetGeometry(), shape_function_derivatives, shape_function, weight);

    // u
    array_1d< double, TDim > velocity;
    this->EvaluateInPoint(velocity, VELOCITY, shape_function);

    double tau_one, tau_two;
    CalculateTau( tau_one, tau_two, rCurrentProcessInfo);

    BoundedMatrix<double, TNumNodes, TDim> tau_one_derivatives;
    CalculateTau1Derivatives( tau_one_derivatives, rCurrentProcessInfo);

    // Density
    double density;
    this->EvaluateInPoint(density, DENSITY, shape_function);

    // BodyForce
    array_1d< double, TDim > body_force;
    this->EvaluateInPoint(body_force, BODY_FORCE,shape_function);
    body_force *= density;

    array_1d< double, TNumNodes > convection_operator;
    GetConvectionOperator(convection_operator, velocity, shape_function_derivatives);

    rMatrix.resize(TCoordLocalSize, TCoordLocalSize, false);
    rMatrix.clear();

    for (unsigned int a=0; a < TNumNodes; ++a)
    {
        for (unsigned int i=0; i<TDim; ++i)
        {
            for (unsigned int c=0; c < TNumNodes; ++c)
            {
                for (unsigned int k=0; k < TDim; ++k)
                {
                    double val = density * convection_operator[a]*body_force[i]*tau_one_derivatives(c,k);
                    val += shape_function[c]*density*shape_function_derivatives(a,k)*tau_one*body_force[i];
                    rMatrix(a*TDim+i, c*TDim+k) = coeff*val*weight;
                }
            }
        }
    }

}

void Matrix_2( Element::MatrixType& rMatrix,
               const ProcessInfo& rCurrentProcessInfo )
{


    double alpha, beta, gamma, delta_time;
    CalculateBossakSchemeConstants(alpha, beta, gamma, delta_time, rCurrentProcessInfo);

    BoundedMatrix<double, TNumNodes, TDim> shape_function_derivatives;
    array_1d< double, TNumNodes > shape_function;
    double weight;

    const double coeff = -gamma/(beta*delta_time);

    GeometryUtils::CalculateGeometryData(this->GetGeometry(), shape_function_derivatives, shape_function, weight);

    // Density
    double density;
    this->EvaluateInPoint(density, DENSITY, shape_function);

    BoundedMatrix<double, TDim, TDim> velocity_gradient;
    GetPostVelocityGradient( velocity_gradient, shape_function_derivatives);

    rMatrix.resize(TCoordLocalSize, TCoordLocalSize, false);
    rMatrix.clear();

    for (unsigned int a=0; a < TNumNodes; ++a)
    {
        for (unsigned int i=0; i<TDim; ++i)
        {
            for (unsigned int c=0; c < TNumNodes; ++c)
            {
                for (unsigned int k=0; k < TDim; ++k)
                {
                    double val = density * shape_function[a] * shape_function[c] * velocity_gradient(k,i);
                    rMatrix(a*TDim+i, c*TDim+k) = coeff*val*weight;
                }
            }
        }
    }

}

void Matrix_3( Element::MatrixType& rMatrix,
               const ProcessInfo& rCurrentProcessInfo )
{


    double alpha, beta, gamma, delta_time;
    CalculateBossakSchemeConstants(alpha, beta, gamma, delta_time, rCurrentProcessInfo);

    BoundedMatrix<double, TNumNodes, TDim> shape_function_derivatives;
    array_1d< double, TNumNodes > shape_function;
    double weight;

    const double coeff = -gamma/(beta*delta_time);

    GeometryUtils::CalculateGeometryData(this->GetGeometry(), shape_function_derivatives, shape_function, weight);

    // u
    array_1d< double, TDim > velocity;
    this->EvaluateInPoint(velocity, VELOCITY, shape_function);

    double tau_one, tau_two;
    CalculateTau( tau_one, tau_two, rCurrentProcessInfo);

    BoundedMatrix<double, TNumNodes, TDim> tau_one_derivatives;
    CalculateTau1Derivatives( tau_one_derivatives, rCurrentProcessInfo);

    // Density
    double density;
    this->EvaluateInPoint(density, DENSITY, shape_function);

    array_1d< double, TNumNodes > convection_operator;
    GetConvectionOperator(convection_operator, velocity, shape_function_derivatives);

    BoundedMatrix<double, TDim, TDim> velocity_gradient;
    GetPostVelocityGradient( velocity_gradient, shape_function_derivatives);

    BoundedVector<double, TDim> convective_velocity;
    convective_velocity.clear();
    for (unsigned int i=0;i<TDim;++i)
        for (unsigned int j=0;j<TDim;j++)
            convective_velocity[i] += velocity[j]*velocity_gradient(j,i);

    rMatrix.resize(TCoordLocalSize, TCoordLocalSize, false);
    rMatrix.clear();

    for (unsigned int a=0; a < TNumNodes; ++a)
    {
        for (unsigned int i=0; i<TDim; ++i)
        {
            for (unsigned int c=0; c < TNumNodes; ++c)
            {
                for (unsigned int k=0; k < TDim; ++k)
                {
                    double val = density*density*tau_one_derivatives(c,k)*convection_operator[a]*convective_velocity[i];
                    val += density*density*tau_one*shape_function[c]*shape_function_derivatives(a,k)*convective_velocity[i];
                    val += density*density*tau_one*convection_operator[a]*shape_function[c]*velocity_gradient(k,i);
                    rMatrix(a*TDim+i, c*TDim+k) = coeff*val*weight;
                }
            }
        }
    }

}

void Matrix_4( Element::MatrixType& rMatrix,
               const ProcessInfo& rCurrentProcessInfo )
{


    double alpha, beta, gamma, delta_time;
    CalculateBossakSchemeConstants(alpha, beta, gamma, delta_time, rCurrentProcessInfo);

    BoundedMatrix<double, TNumNodes, TDim> shape_function_derivatives;
    array_1d< double, TNumNodes > shape_function;
    double weight;

    const double coeff = -gamma/(beta*delta_time);

    GeometryUtils::CalculateGeometryData(this->GetGeometry(), shape_function_derivatives, shape_function, weight);

    BoundedMatrix<double, TNumNodes, TDim> tau_two_derivatives;
    CalculateTau2Derivatives( tau_two_derivatives, rCurrentProcessInfo);

    BoundedMatrix<double, TDim, TDim> velocity_gradient;
    GetPostVelocityGradient( velocity_gradient, shape_function_derivatives);

    double divergence = 0.0;
    for (unsigned int i=0;i<TDim;++i)
        divergence += velocity_gradient(i,i);

    rMatrix.resize(TCoordLocalSize, TCoordLocalSize, false);
    rMatrix.clear();

    for (unsigned int a=0; a < TNumNodes; ++a)
    {
        for (unsigned int i=0; i<TDim; ++i)
        {
            for (unsigned int c=0; c < TNumNodes; ++c)
            {
                for (unsigned int k=0; k < TDim; ++k)
                {
                    double val = tau_two_derivatives(c,k)*shape_function_derivatives(a,i)*divergence;
                    rMatrix(a*TDim+i, c*TDim+k) = coeff*val*weight;
                }
            }
        }
    }

}

void Matrix_5( Element::MatrixType& rMatrix,
               const ProcessInfo& rCurrentProcessInfo )
{


    double alpha, beta, gamma, delta_time;
    CalculateBossakSchemeConstants(alpha, beta, gamma, delta_time, rCurrentProcessInfo);

    BoundedMatrix<double, TNumNodes, TDim> shape_function_derivatives;
    array_1d< double, TNumNodes > shape_function;
    double weight;

    const double coeff = -gamma/(beta*delta_time);

    GeometryUtils::CalculateGeometryData(this->GetGeometry(), shape_function_derivatives, shape_function, weight);

    // u
    array_1d< double, TDim > velocity;
    this->EvaluateInPoint(velocity, VELOCITY, shape_function);

    double tau_one, tau_two;
    CalculateTau( tau_one, tau_two, rCurrentProcessInfo);

    BoundedMatrix<double, TNumNodes, TDim> tau_one_derivatives;
    CalculateTau1Derivatives( tau_one_derivatives, rCurrentProcessInfo);

    // Density
    double density;
    this->EvaluateInPoint(density, DENSITY, shape_function);

    array_1d< double, TNumNodes > convection_operator;
    GetConvectionOperator(convection_operator, velocity, shape_function_derivatives);

    BoundedVector<double, TDim> pressure_gradient;
    GetPostPressureGradient( pressure_gradient, shape_function_derivatives);

    rMatrix.resize(TCoordLocalSize, TCoordLocalSize, false);
    rMatrix.clear();

    for (unsigned int a=0; a < TNumNodes; ++a)
    {
        for (unsigned int i=0; i<TDim; ++i)
        {
            for (unsigned int c=0; c < TNumNodes; ++c)
            {
                for (unsigned int k=0; k < TDim; ++k)
                {
                    double val = density*convection_operator[a]*tau_one_derivatives(c,k)*pressure_gradient[i];
                    val += density*shape_function[c]*shape_function_derivatives(a,k)*tau_one*pressure_gradient[i];
                    rMatrix(a*TDim+i, c*TDim+k) = coeff*val*weight;
                }
            }
        }
    }

}

void Matrix_6( Element::MatrixType& rMatrix,
               const ProcessInfo& rCurrentProcessInfo )
{


    double alpha, beta, gamma, delta_time;
    CalculateBossakSchemeConstants(alpha, beta, gamma, delta_time, rCurrentProcessInfo);

    BoundedMatrix<double, TNumNodes, TDim> shape_function_derivatives;
    array_1d< double, TNumNodes > shape_function;
    double weight;

    const double coeff = -gamma/(beta*delta_time);

    GeometryUtils::CalculateGeometryData(this->GetGeometry(), shape_function_derivatives, shape_function, weight);

    // u
    array_1d< double, TDim > velocity;
    this->EvaluateInPoint(velocity, VELOCITY, shape_function);

    // Density
    double density;
    this->EvaluateInPoint(density, DENSITY, shape_function);

    array_1d< double, TNumNodes > convection_operator;
    GetConvectionOperator(convection_operator, velocity, shape_function_derivatives);

    rMatrix.resize(TCoordLocalSize, TCoordLocalSize, false);
    rMatrix.clear();

    for (unsigned int a=0; a < TNumNodes; ++a)
    {
        for (unsigned int i=0; i<TDim; ++i)
        {
            for (unsigned int c=0; c < TNumNodes; ++c)
            {
                for (unsigned int k=0; k < TDim; ++k)
                {
                    double val = density*shape_function[a]*convection_operator[c];
                    if (i!=k)
                        val = 0.0;
                    rMatrix(a*TDim+i, c*TDim+k) = coeff*val*weight;
                }
            }
        }
    }
}

void Matrix_7( Element::MatrixType& rMatrix,
               const ProcessInfo& rCurrentProcessInfo )
{


    double alpha, beta, gamma, delta_time;
    CalculateBossakSchemeConstants(alpha, beta, gamma, delta_time, rCurrentProcessInfo);

    BoundedMatrix<double, TNumNodes, TDim> shape_function_derivatives;
    array_1d< double, TNumNodes > shape_function;
    double weight;

    const double coeff = -gamma/(beta*delta_time);

    GeometryUtils::CalculateGeometryData(this->GetGeometry(), shape_function_derivatives, shape_function, weight);

    // Density
    double density;
    this->EvaluateInPoint(density, DENSITY, shape_function);

    // Viscosity
    double viscosity;
    this->EvaluateInPoint(viscosity, VISCOSITY, shape_function);
    viscosity *= density;

    rMatrix.resize(TCoordLocalSize, TCoordLocalSize, false);
    rMatrix.clear();

    for (unsigned int a=0; a < TNumNodes; ++a)
    {
        for (unsigned int i=0; i<TDim; ++i)
        {
            for (unsigned int c=0; c < TNumNodes; ++c){
                double diag = 0.0;
                for (unsigned int l=0; l < TDim; ++l)
                {
                    diag += shape_function_derivatives(a,l)*shape_function_derivatives(c,l);
                }
                for (unsigned int k=0; k < TDim; ++k)
                {
                    double val = shape_function_derivatives(a,k)*shape_function_derivatives(c,i);
                    val -= (2.0/3.0)*shape_function_derivatives(a,i)*shape_function_derivatives(c,k);

                    if (i==k)
                        val += diag;
                    rMatrix(a*TDim+i, c*TDim+k) = coeff*viscosity*val*weight;
                }
            }
        }
    }
}

void Matrix_8( Element::MatrixType& rMatrix,
               const ProcessInfo& rCurrentProcessInfo )
{


    double alpha, beta, gamma, delta_time;
    CalculateBossakSchemeConstants(alpha, beta, gamma, delta_time, rCurrentProcessInfo);

    BoundedMatrix<double, TNumNodes, TDim> shape_function_derivatives;
    array_1d< double, TNumNodes > shape_function;
    double weight;

    const double coeff = -gamma/(beta*delta_time);

    GeometryUtils::CalculateGeometryData(this->GetGeometry(), shape_function_derivatives, shape_function, weight);

    // u
    array_1d< double, TDim > velocity;
    this->EvaluateInPoint(velocity, VELOCITY, shape_function);

    double tau_one, tau_two;
    CalculateTau( tau_one, tau_two, rCurrentProcessInfo);

    // Density
    double density;
    this->EvaluateInPoint(density, DENSITY, shape_function);

    array_1d< double, TNumNodes > convection_operator;
    GetConvectionOperator(convection_operator, velocity, shape_function_derivatives);

    rMatrix.resize(TCoordLocalSize, TCoordLocalSize, false);
    rMatrix.clear();

    for (unsigned int a=0; a < TNumNodes; ++a)
    {
        for (unsigned int i=0; i<TDim; ++i)
        {
            for (unsigned int c=0; c < TNumNodes; ++c)
            {
                for (unsigned int k=0; k < TDim; ++k)
                {
                    double val = 0.0;
                    if (i==k)
                        val = density*density*tau_one*convection_operator[a]*convection_operator[c];
                    rMatrix(a*TDim+i, c*TDim+k) = coeff*val*weight;
                }
            }
        }
    }
}

void Matrix_9( Element::MatrixType& rMatrix,
               const ProcessInfo& rCurrentProcessInfo )
{


    double alpha, beta, gamma, delta_time;
    CalculateBossakSchemeConstants(alpha, beta, gamma, delta_time, rCurrentProcessInfo);

    BoundedMatrix<double, TNumNodes, TDim> shape_function_derivatives;
    array_1d< double, TNumNodes > shape_function;
    double weight;

    const double coeff = -gamma/(beta*delta_time);

    GeometryUtils::CalculateGeometryData(this->GetGeometry(), shape_function_derivatives, shape_function, weight);

    double tau_one, tau_two;
    CalculateTau( tau_one, tau_two, rCurrentProcessInfo);

    rMatrix.resize(TCoordLocalSize, TCoordLocalSize, false);
    rMatrix.clear();

    for (unsigned int a=0; a < TNumNodes; ++a)
    {
        for (unsigned int i=0; i<TDim; ++i)
        {
            for (unsigned int c=0; c < TNumNodes; ++c)
            {
                for (unsigned int k=0; k < TDim; ++k)
                {
                    double val = tau_two*shape_function_derivatives(a,i)*shape_function_derivatives(c,k);
                    rMatrix(a*TDim+i, c*TDim+k) = coeff*val*weight;
                }
            }
        }
    }
}

void Matrix_10( Element::MatrixType& rMatrix,
                const ProcessInfo& rCurrentProcessInfo )
{


    double alpha, beta, gamma, delta_time;
    CalculateBossakSchemeConstants(alpha, beta, gamma, delta_time, rCurrentProcessInfo);

    BoundedMatrix<double, TNumNodes, TDim> shape_function_derivatives;
    array_1d< double, TNumNodes > shape_function;
    double weight;

    const double coeff = -gamma/(beta*delta_time);

    GeometryUtils::CalculateGeometryData(this->GetGeometry(), shape_function_derivatives, shape_function, weight);

    // u
    array_1d< double, TDim > velocity;
    this->EvaluateInPoint(velocity, VELOCITY, shape_function);

    // a
    array_1d< double, TDim > acceleration;
    this->EvaluateInPoint(acceleration, ACCELERATION, shape_function);

    double tau_one, tau_two;
    CalculateTau( tau_one, tau_two, rCurrentProcessInfo);

    BoundedMatrix<double, TNumNodes, TDim> tau_one_derivatives;
    CalculateTau1Derivatives( tau_one_derivatives, rCurrentProcessInfo);

    // Density
    double density;
    this->EvaluateInPoint(density, DENSITY, shape_function);

    array_1d< double, TNumNodes > convection_operator;
    GetConvectionOperator(convection_operator, velocity, shape_function_derivatives);

    rMatrix.resize(TCoordLocalSize, TCoordLocalSize, false);
    rMatrix.clear();

    for (unsigned int a=0; a < TNumNodes; ++a)
    {
        for (unsigned int i=0; i<TDim; ++i)
        {
            for (unsigned int c=0; c < TNumNodes; ++c)
            {
                for (unsigned int k=0; k < TDim; ++k)
                {
                    double val = density*density*convection_operator[a]*tau_one_derivatives(c,k)*acceleration[i];
                    val += density*density*shape_function[c]*shape_function_derivatives(a,k)*tau_one*acceleration[i];
                    rMatrix(a*TDim+i, c*TDim+k) = coeff*val*weight;
                }
            }
        }
    }
}

void Matrix_11( Element::MatrixType& rMatrix,
                const ProcessInfo& rCurrentProcessInfo )
{


    double alpha, beta, gamma, delta_time;
    CalculateBossakSchemeConstants(alpha, beta, gamma, delta_time, rCurrentProcessInfo);

    BoundedMatrix<double, TNumNodes, TDim> shape_function_derivatives;
    array_1d< double, TNumNodes > shape_function;
    double weight;

    const double coeff = -(1-alpha)/(beta*delta_time*delta_time);

    GeometryUtils::CalculateGeometryData(this->GetGeometry(), shape_function_derivatives, shape_function, weight);

    // Density
    double density;
    this->EvaluateInPoint(density, DENSITY, shape_function);

    rMatrix.resize(TCoordLocalSize, TCoordLocalSize, false);
    rMatrix.clear();

    double mass = density / static_cast<double>(TNumNodes);

    unsigned int DofIndex = 0;
    for (unsigned int a=0; a < TNumNodes; ++a)
    {
        for (unsigned int i=0; i<TDim; ++i)
        {
            rMatrix(DofIndex, DofIndex) = coeff*mass*weight;
            DofIndex++;
        }
    }
}

void Matrix_12( Element::MatrixType& rMatrix,
                const ProcessInfo& rCurrentProcessInfo )
{


    double alpha, beta, gamma, delta_time;
    CalculateBossakSchemeConstants(alpha, beta, gamma, delta_time, rCurrentProcessInfo);

    BoundedMatrix<double, TNumNodes, TDim> shape_function_derivatives;
    array_1d< double, TNumNodes > shape_function;
    double weight;

    const double coeff = -(1-alpha)/(beta*delta_time*delta_time);

    GeometryUtils::CalculateGeometryData(this->GetGeometry(), shape_function_derivatives, shape_function, weight);

    // u
    array_1d< double, TDim > velocity;
    this->EvaluateInPoint(velocity, VELOCITY, shape_function);

    double tau_one, tau_two;
    CalculateTau( tau_one, tau_two, rCurrentProcessInfo);

    // Density
    double density;
    this->EvaluateInPoint(density, DENSITY, shape_function);

    array_1d< double, TNumNodes > convection_operator;
    GetConvectionOperator(convection_operator, velocity, shape_function_derivatives);

    rMatrix.resize(TCoordLocalSize, TCoordLocalSize, false);
    rMatrix.clear();

    for (unsigned int a=0; a < TNumNodes; ++a)
    {
        for (unsigned int i=0; i<TDim; ++i)
        {
            for (unsigned int c=0; c < TNumNodes; ++c)
            {
                for (unsigned int k=0; k < TDim; ++k)
                {
                    double val = 0.0;
                    if (i==k)
                        val = density*density*convection_operator[a]*tau_one*shape_function[c];
                    rMatrix(a*TDim+i, c*TDim+k) = coeff*val*weight;
                }
            }
        }
    }
}

void Matrix_13( Element::MatrixType& rMatrix,
                const ProcessInfo& rCurrentProcessInfo )
{


    double alpha, beta, gamma, delta_time;
    CalculateBossakSchemeConstants(alpha, beta, gamma, delta_time, rCurrentProcessInfo);

    BoundedMatrix<double, TNumNodes, TDim> shape_function_derivatives;
    array_1d< double, TNumNodes > shape_function;
    double weight;

    const double coeff = gamma/(beta*delta_time);

    GeometryUtils::CalculateGeometryData(this->GetGeometry(), shape_function_derivatives, shape_function, weight);

    BoundedMatrix<double, TNumNodes, TDim> tau_one_derivatives;
    CalculateTau1Derivatives( tau_one_derivatives, rCurrentProcessInfo);

    // Density
    double density;
    this->EvaluateInPoint(density, DENSITY, shape_function);

    array_1d< double, TDim > body_force;
    this->EvaluateInPoint(body_force, BODY_FORCE, shape_function);
    body_force *= density;

    array_1d< double, TNumNodes > convection_operator;
    GetConvectionOperator(convection_operator, body_force, shape_function_derivatives);

    rMatrix.resize(TNumNodes, TCoordLocalSize, false);
    rMatrix.clear();

    for (unsigned int a=0; a < TNumNodes; ++a)
    {
        for (unsigned int c=0; c < TNumNodes; ++c)
        {
            for (unsigned int k=0; k < TDim; ++k)
            {
                double val = convection_operator[a]*tau_one_derivatives(c,k);
                rMatrix(a, c*TDim+k) = coeff*val*weight;
            }
        }
    }
}

void Matrix_14( Element::MatrixType& rMatrix,
                const ProcessInfo& rCurrentProcessInfo )
{


    double alpha, beta, gamma, delta_time;
    CalculateBossakSchemeConstants(alpha, beta, gamma, delta_time, rCurrentProcessInfo);

    BoundedMatrix<double, TNumNodes, TDim> shape_function_derivatives;
    array_1d< double, TNumNodes > shape_function;
    double weight;

    const double coeff = -gamma/(beta*delta_time);

    GeometryUtils::CalculateGeometryData(this->GetGeometry(), shape_function_derivatives, shape_function, weight);

    // u
    array_1d< double, TDim > velocity;
    this->EvaluateInPoint(velocity, VELOCITY, shape_function);

    double tau_one, tau_two;
    CalculateTau( tau_one, tau_two, rCurrentProcessInfo);

    BoundedMatrix<double, TNumNodes, TDim> tau_one_derivatives;
    CalculateTau1Derivatives( tau_one_derivatives, rCurrentProcessInfo);

    // Density
    double density;
    this->EvaluateInPoint(density, DENSITY, shape_function);

    BoundedMatrix<double, TDim, TDim> velocity_gradient;
    GetPostVelocityGradient( velocity_gradient, shape_function_derivatives);

    BoundedVector<double, TDim> convective_velocity;
    convective_velocity.clear();
    for (unsigned int i=0;i<TDim;++i)
        for (unsigned int j=0;j<TDim;j++)
            convective_velocity[i] += velocity[j]*velocity_gradient(j,i);

    rMatrix.resize(TNumNodes, TCoordLocalSize, false);
    rMatrix.clear();

    for (unsigned int a=0; a < TNumNodes; ++a)
    {
        for (unsigned int c=0; c < TNumNodes; ++c)
        {
            for (unsigned int k=0; k < TDim; ++k)
            {
                double val = 0.0;
                for (unsigned int i=0;i<TDim; ++i)
                {
                    val += shape_function_derivatives(a,i)*density*tau_one_derivatives(c,k)*convective_velocity[i];
                    val += shape_function_derivatives(a,i)*tau_one*density*shape_function[c]*velocity_gradient(k,i);
                }
                rMatrix(a, c*TDim+k) = coeff*val*weight;
            }
        }
    }
}

void Matrix_15( Element::MatrixType& rMatrix,
                const ProcessInfo& rCurrentProcessInfo )
{


    double alpha, beta, gamma, delta_time;
    CalculateBossakSchemeConstants(alpha, beta, gamma, delta_time, rCurrentProcessInfo);

    BoundedMatrix<double, TNumNodes, TDim> shape_function_derivatives;
    array_1d< double, TNumNodes > shape_function;
    double weight;

    const double coeff = -gamma/(beta*delta_time);

    GeometryUtils::CalculateGeometryData(this->GetGeometry(), shape_function_derivatives, shape_function, weight);

    BoundedMatrix<double, TNumNodes, TDim> tau_one_derivatives;
    CalculateTau1Derivatives( tau_one_derivatives, rCurrentProcessInfo);

    BoundedVector<double, TDim> pressure_gradient;
    GetPostPressureGradient( pressure_gradient, shape_function_derivatives);

    array_1d< double, TNumNodes > convection_operator;
    GetConvectionOperator(convection_operator, pressure_gradient, shape_function_derivatives);

    rMatrix.resize(TNumNodes, TCoordLocalSize, false);
    rMatrix.clear();

    for (unsigned int a=0; a < TNumNodes; ++a)
    {
        for (unsigned int c=0; c < TNumNodes; ++c)
        {
            for (unsigned int k=0; k < TDim; ++k)
            {
                double val = convection_operator[a]*tau_one_derivatives(c,k);
                rMatrix(a, c*TDim+k) = coeff*val*weight;
            }
        }
    }
}

void Matrix_16( Element::MatrixType& rMatrix,
                const ProcessInfo& rCurrentProcessInfo )
{


    double alpha, beta, gamma, delta_time;
    CalculateBossakSchemeConstants(alpha, beta, gamma, delta_time, rCurrentProcessInfo);

    BoundedMatrix<double, TNumNodes, TDim> shape_function_derivatives;
    array_1d< double, TNumNodes > shape_function;
    double weight;

    const double coeff = -gamma/(beta*delta_time);

    GeometryUtils::CalculateGeometryData(this->GetGeometry(), shape_function_derivatives, shape_function, weight);

    rMatrix.resize(TNumNodes, TCoordLocalSize, false);
    rMatrix.clear();

    for (unsigned int a=0; a < TNumNodes; ++a)
    {
        for (unsigned int c=0; c < TNumNodes; ++c)
        {
            for (unsigned int k=0; k < TDim; ++k)
            {
                double val = shape_function[a]*shape_function_derivatives(c,k);
                rMatrix(a, c*TDim+k) = coeff*val*weight;
            }
        }
    }
}

void Matrix_17( Element::MatrixType& rMatrix,
                const ProcessInfo& rCurrentProcessInfo )
{


    double alpha, beta, gamma, delta_time;
    CalculateBossakSchemeConstants(alpha, beta, gamma, delta_time, rCurrentProcessInfo);

    BoundedMatrix<double, TNumNodes, TDim> shape_function_derivatives;
    array_1d< double, TNumNodes > shape_function;
    double weight;

    const double coeff = -gamma/(beta*delta_time);

    GeometryUtils::CalculateGeometryData(this->GetGeometry(), shape_function_derivatives, shape_function, weight);

    // u
    array_1d< double, TDim > velocity;
    this->EvaluateInPoint(velocity, VELOCITY, shape_function);

    double tau_one, tau_two;
    CalculateTau( tau_one, tau_two, rCurrentProcessInfo);

    // Density
    double density;
    this->EvaluateInPoint(density, DENSITY, shape_function);

    array_1d< double, TNumNodes > convection_operator;
    GetConvectionOperator(convection_operator, velocity, shape_function_derivatives);

    rMatrix.resize(TNumNodes, TCoordLocalSize, false);
    rMatrix.clear();

    for (unsigned int a=0; a < TNumNodes; ++a)
    {
        for (unsigned int c=0; c < TNumNodes; ++c)
        {
            for (unsigned int k=0; k < TDim; ++k)
            {
                double val = shape_function_derivatives(a,k)*tau_one*density*convection_operator[c];
                rMatrix(a, c*TDim+k) = coeff*val*weight;
            }
        }
    }
}

void Matrix_18( Element::MatrixType& rMatrix,
                const ProcessInfo& rCurrentProcessInfo )
{


    double alpha, beta, gamma, delta_time;
    CalculateBossakSchemeConstants(alpha, beta, gamma, delta_time, rCurrentProcessInfo);

    BoundedMatrix<double, TNumNodes, TDim> shape_function_derivatives;
    array_1d< double, TNumNodes > shape_function;
    double weight;

    const double coeff = -gamma/(beta*delta_time);

    GeometryUtils::CalculateGeometryData(this->GetGeometry(), shape_function_derivatives, shape_function, weight);

    // a
    array_1d< double, TDim > acceleration;
    this->EvaluateInPoint(acceleration, ACCELERATION, shape_function);

    BoundedMatrix<double, TNumNodes, TDim> tau_one_derivatives;
    CalculateTau1Derivatives( tau_one_derivatives, rCurrentProcessInfo);

    // Density
    double density;
    this->EvaluateInPoint(density, DENSITY, shape_function);

    array_1d< double, TNumNodes > convection_operator;
    GetConvectionOperator(convection_operator, acceleration, shape_function_derivatives);

    rMatrix.resize(TNumNodes, TCoordLocalSize, false);
    rMatrix.clear();

    for (unsigned int a=0; a < TNumNodes; ++a)
    {
        for (unsigned int c=0; c < TNumNodes; ++c)
        {
            for (unsigned int k=0; k < TDim; ++k)
            {
                double val = density*convection_operator[a]*tau_one_derivatives(c,k);
                rMatrix(a, c*TDim+k) = coeff*val*weight;
            }
        }
    }
}

void Matrix_19( Element::MatrixType& rMatrix,
                const ProcessInfo& rCurrentProcessInfo )
{


    double alpha, beta, gamma, delta_time;
    CalculateBossakSchemeConstants(alpha, beta, gamma, delta_time, rCurrentProcessInfo);

    BoundedMatrix<double, TNumNodes, TDim> shape_function_derivatives;
    array_1d< double, TNumNodes > shape_function;
    double weight;

    const double coeff = -(1-alpha)/(beta*delta_time*delta_time);

    GeometryUtils::CalculateGeometryData(this->GetGeometry(), shape_function_derivatives, shape_function, weight);

    double tau_one, tau_two;
    CalculateTau( tau_one, tau_two, rCurrentProcessInfo);

    // Density
    double density;
    this->EvaluateInPoint(density, DENSITY, shape_function);

    rMatrix.resize(TNumNodes, TCoordLocalSize, false);
    rMatrix.clear();

    for (unsigned int a=0; a < TNumNodes; ++a)
    {
        for (unsigned int c=0; c < TNumNodes; ++c)
        {
            for (unsigned int k=0; k < TDim; ++k)
            {
                double val = density*shape_function_derivatives(a,k)*tau_one*shape_function[c];
                rMatrix(a, c*TDim+k) = coeff*val*weight;
            }
        }
    }
}

void Matrix_20( Element::MatrixType& rMatrix,
                const ProcessInfo& rCurrentProcessInfo )
{


    double alpha, beta, gamma, delta_time;
    CalculateBossakSchemeConstants(alpha, beta, gamma, delta_time, rCurrentProcessInfo);

    BoundedMatrix<double, TNumNodes, TDim> shape_function_derivatives;
    array_1d< double, TNumNodes > shape_function;
    double weight;

    const double coeff = gamma/(beta*delta_time);

    GeometryUtils::CalculateGeometryData(this->GetGeometry(), shape_function_derivatives, shape_function, weight);

    rMatrix.resize(TCoordLocalSize, TNumNodes, false);
    rMatrix.clear();

    for (unsigned int a=0; a < TNumNodes; ++a)
    {
        for (unsigned int i=0; i < TDim; ++i)
        {
            for (unsigned int c=0; c < TNumNodes; ++c)
            {
                double val = shape_function_derivatives(a,i)*shape_function[c];
                rMatrix(a*TDim+i, c) = coeff*val*weight;
            }
        }
    }
}

void Matrix_21( Element::MatrixType& rMatrix,
                const ProcessInfo& rCurrentProcessInfo )
{


    double alpha, beta, gamma, delta_time;
    CalculateBossakSchemeConstants(alpha, beta, gamma, delta_time, rCurrentProcessInfo);

    BoundedMatrix<double, TNumNodes, TDim> shape_function_derivatives;
    array_1d< double, TNumNodes > shape_function;
    double weight;

    const double coeff = -gamma/(beta*delta_time);

    GeometryUtils::CalculateGeometryData(this->GetGeometry(), shape_function_derivatives, shape_function, weight);

    // u
    array_1d< double, TDim > velocity;
    this->EvaluateInPoint(velocity, VELOCITY, shape_function);

    double tau_one, tau_two;
    CalculateTau( tau_one, tau_two, rCurrentProcessInfo);

    // Density
    double density;
    this->EvaluateInPoint(density, DENSITY, shape_function);

    array_1d< double, TNumNodes > convection_operator;
    GetConvectionOperator(convection_operator, velocity, shape_function_derivatives);

    rMatrix.resize(TCoordLocalSize, TNumNodes, false);
    rMatrix.clear();

    for (unsigned int a=0; a < TNumNodes; ++a)
    {
        for (unsigned int i=0; i < TDim; ++i)
        {
            for (unsigned int c=0; c < TNumNodes; ++c)
            {
                double val = density*convection_operator[a]*tau_one*shape_function_derivatives(c,i);
                rMatrix(a*TDim+i, c) = coeff*val*weight;
            }
        }
    }
}

void Matrix_22( Element::MatrixType& rMatrix,
                const ProcessInfo& rCurrentProcessInfo )
{


    double alpha, beta, gamma, delta_time;
    CalculateBossakSchemeConstants(alpha, beta, gamma, delta_time, rCurrentProcessInfo);

    BoundedMatrix<double, TNumNodes, TDim> shape_function_derivatives;
    array_1d< double, TNumNodes > shape_function;
    double weight;

    const double coeff = -gamma/(beta*delta_time);

    GeometryUtils::CalculateGeometryData(this->GetGeometry(), shape_function_derivatives, shape_function, weight);

    double tau_one, tau_two;
    CalculateTau( tau_one, tau_two, rCurrentProcessInfo);

    rMatrix.resize(TNumNodes, TNumNodes, false);
    rMatrix.clear();

    for (unsigned int a=0; a < TNumNodes; ++a)
    {
        for (unsigned int c=0; c < TNumNodes; ++c)
        {
            double val = 0.0;
            for (unsigned int i=0; i < TDim; ++i)
            {
                val += shape_function_derivatives(a,i)*shape_function_derivatives(c,i);
            }
            val *= tau_one;
            rMatrix(a, c) = coeff*val*weight;
        }
    }
}

void Matrix_25( Element::MatrixType& rMatrix,
                const ProcessInfo& rCurrentProcessInfo )
{
     KRATOS_TRY;

    if (rMatrix.size1() != TCoordLocalSize || rMatrix.size2() != TCoordLocalSize)
        rMatrix.resize(TCoordLocalSize,TCoordLocalSize,false);

    rMatrix.clear();

    // Get shape functions, shape function gradients and element volume (area in
    // 2D). Only one integration point is used so the volume is its weight.
    ShapeFunctionDerivativesType shape_function_derivatives;
    array_1d< double, TNumNodes > shape_functions;
    double volume;
    GeometryUtils::CalculateGeometryData(this->GetGeometry(), shape_function_derivatives, shape_functions, volume);

    // Density
    double density;
    this->EvaluateInPoint(density, DENSITY, shape_functions);

    // Dynamic viscosity
    double viscosity;
    this->EvaluateInPoint(viscosity, VISCOSITY, shape_functions);
    viscosity *= density;

    // u
    array_1d< double, TDim > velocity;
    this->EvaluateInPoint(velocity, VELOCITY, shape_functions);

    // a
    array_1d< double, TDim > acceleration;
    this->EvaluateInPoint(acceleration, AUX_ADJOINT_ACCELERATION, shape_functions);

    array_1d<double, TNumNodes> acceleration_grad_shape_function_derivatives;
    this->GetConvectionOperator(
        acceleration_grad_shape_function_derivatives,
        acceleration,
        shape_function_derivatives);

    array_1d<double, TNumNodes> velocity_grad_shape_function_derivatives;
    this->GetConvectionOperator(
        velocity_grad_shape_function_derivatives,
        velocity,
        shape_function_derivatives);

    const double velocity_norm = norm_2(velocity);
    const double elem_size = this->CalculateElementSize(volume);
    double tau_one, tau_two;

    this->CalculateStabilizationParameters(
        tau_one,
        tau_two,
        velocity_norm,
        elem_size,
        density,
        viscosity,
        rCurrentProcessInfo
    );

    double tau_one_time_gradient;
    this->CalculateTauOneTimeGradient(
        tau_one_time_gradient,
        tau_one,
        density,
        elem_size,
        velocity_norm,
        velocity,
        acceleration
    );

    const double coeff = density * density;

    for (unsigned int a = 0; a < TNumNodes; ++a)
    {
        for (unsigned int b = 0; b < TNumNodes; ++b)
        {
            for (unsigned i = 0; i < TDim; ++i)
            {
                double value = 0.0;
                value += acceleration_grad_shape_function_derivatives[a] * tau_one * shape_functions[b];
                value += velocity_grad_shape_function_derivatives[a] * tau_one_time_gradient * shape_functions[b];

                rMatrix(a*TDim +i, b*TDim +i) = coeff * value * volume;
            }
        }
    }

    KRATOS_CATCH("");
}

void Matrix_26( Element::MatrixType& rMatrix,
                const ProcessInfo& rCurrentProcessInfo )
{
     KRATOS_TRY;

    if (rMatrix.size1() != TNumNodes || rMatrix.size2() != TCoordLocalSize)
        rMatrix.resize(TNumNodes,TCoordLocalSize,false);

    rMatrix.clear();

    // Get shape functions, shape function gradients and element volume (area in
    // 2D). Only one integration point is used so the volume is its weight.
    ShapeFunctionDerivativesType shape_function_derivatives;
    array_1d< double, TNumNodes > shape_functions;
    double volume;
    GeometryUtils::CalculateGeometryData(this->GetGeometry(), shape_function_derivatives, shape_functions, volume);

    // Density
    double density;
    this->EvaluateInPoint(density, DENSITY, shape_functions);

    // Dynamic viscosity
    double viscosity;
    this->EvaluateInPoint(viscosity, VISCOSITY, shape_functions);
    viscosity *= density;

    // u
    array_1d< double, TDim > velocity;
    this->EvaluateInPoint(velocity, VELOCITY, shape_functions);

    // a
    array_1d< double, TDim > acceleration;
    this->EvaluateInPoint(acceleration, AUX_ADJOINT_ACCELERATION, shape_functions);

    const double velocity_norm = norm_2(velocity);
    const double elem_size = this->CalculateElementSize(volume);
    double tau_one, tau_two;

    this->CalculateStabilizationParameters(
        tau_one,
        tau_two,
        velocity_norm,
        elem_size,
        density,
        viscosity,
        rCurrentProcessInfo
    );

    double tau_one_time_gradient;
    this->CalculateTauOneTimeGradient(
        tau_one_time_gradient,
        tau_one,
        density,
        elem_size,
        velocity_norm,
        velocity,
        acceleration
    );

    const double coeff = density;

    for (unsigned int a = 0; a < TNumNodes; ++a)
    {
        for (unsigned int b = 0; b < TNumNodes; ++b)
        {
            for (unsigned j = 0; j < TDim; ++j)
            {
                double value = 0.0;

                value += shape_function_derivatives(a,j) * tau_one_time_gradient * shape_functions(b);

                rMatrix(a, b*TDim +j) = coeff * value * volume;
            }
        }
    }

    KRATOS_CATCH("");
}

void Matrix_27( Element::MatrixType& rMatrix,
                const ProcessInfo& rCurrentProcessInfo )
{
     KRATOS_TRY;

    if (rMatrix.size1() != TCoordLocalSize || rMatrix.size2() != TCoordLocalSize)
        rMatrix.resize(TCoordLocalSize,TCoordLocalSize,false);

    rMatrix.clear();

    // Get shape functions, shape function gradients and element volume (area in
    // 2D). Only one integration point is used so the volume is its weight.
    ShapeFunctionDerivativesType shape_function_derivatives;
    array_1d< double, TNumNodes > shape_functions;
    double volume;
    GeometryUtils::CalculateGeometryData(this->GetGeometry(), shape_function_derivatives, shape_functions, volume);

    // Density
    double density;
    this->EvaluateInPoint(density, DENSITY, shape_functions);

    // Dynamic viscosity
    double viscosity;
    this->EvaluateInPoint(viscosity, VISCOSITY, shape_functions);
    viscosity *= density;

    // u
    array_1d< double, TDim > velocity;
    this->EvaluateInPoint(velocity, VELOCITY, shape_functions);

    // a
    array_1d< double, TDim > acceleration;
    this->EvaluateInPoint(acceleration, AUX_ADJOINT_ACCELERATION, shape_functions);

    array_1d<double, TNumNodes> acceleration_grad_shape_function_derivatives;
    this->GetConvectionOperator(
        acceleration_grad_shape_function_derivatives,
        acceleration,
        shape_function_derivatives);

    array_1d<double, TNumNodes> velocity_grad_shape_function_derivatives;
    this->GetConvectionOperator(
        velocity_grad_shape_function_derivatives,
        velocity,
        shape_function_derivatives);

    const double velocity_norm = norm_2(velocity);
    const double elem_size = this->CalculateElementSize(volume);
    double tau_one, tau_two;

    this->CalculateStabilizationParameters(
        tau_one,
        tau_two,
        velocity_norm,
        elem_size,
        density,
        viscosity,
        rCurrentProcessInfo
    );

    double tau_one_time_gradient;
    this->CalculateTauOneTimeGradient(
        tau_one_time_gradient,
        tau_one,
        density,
        elem_size,
        velocity_norm,
        velocity,
        acceleration
    );

    BoundedMatrix<double, TNumNodes, TDim> tau_one_primal_derivatives;
    this->CalculateTauOnePrimalDerivatives(
        tau_one_primal_derivatives,
        tau_one,
        density,
        elem_size,
        velocity_norm,
        shape_functions,
        velocity
    );

    BoundedMatrix<double, TNumNodes, TDim> tau_one_time_gradient_primal_derivatives;
    this->CalculateTauOneTimeGradientPrimalDerivatives(
        tau_one_time_gradient_primal_derivatives,
        tau_one,
        density,
        elem_size,
        velocity_norm,
        shape_functions,
        velocity,
        acceleration
    );

    const double coeff = density * density;

    for (unsigned int a = 0; a < TNumNodes; ++a)
    {
        for (unsigned int c = 0; c < TNumNodes; ++c)
        {
            for (unsigned i = 0; i < TDim; ++i)
            {
                for (unsigned int k = 0; k < TDim; ++k)
                {
                    double value = 0.0;

                    value += acceleration_grad_shape_function_derivatives[a] * tau_one_primal_derivatives(c,k) * velocity[i];
                    value += shape_functions[c] * shape_function_derivatives(a,k) * tau_one_time_gradient * velocity[i];
                    value += velocity_grad_shape_function_derivatives[a] * tau_one_time_gradient_primal_derivatives(c,k) * velocity[i];

                    rMatrix(a*TDim+i, c*TDim + k) = coeff * value * volume;
                }
            }
        }
    }

    KRATOS_CATCH("");
}

void Matrix_28( Element::MatrixType& rMatrix,
                const ProcessInfo& rCurrentProcessInfo )
{
     KRATOS_TRY;

    if (rMatrix.size1() != TNumNodes || rMatrix.size2() != TCoordLocalSize)
        rMatrix.resize(TNumNodes,TCoordLocalSize,false);

    rMatrix.clear();

    // Get shape functions, shape function gradients and element volume (area in
    // 2D). Only one integration point is used so the volume is its weight.
    ShapeFunctionDerivativesType shape_function_derivatives;
    array_1d< double, TNumNodes > shape_functions;
    double volume;
    GeometryUtils::CalculateGeometryData(this->GetGeometry(), shape_function_derivatives, shape_functions, volume);

    // Density
    double density;
    this->EvaluateInPoint(density, DENSITY, shape_functions);

    // Dynamic viscosity
    double viscosity;
    this->EvaluateInPoint(viscosity, VISCOSITY, shape_functions);
    viscosity *= density;

    // u
    array_1d< double, TDim > velocity;
    this->EvaluateInPoint(velocity, VELOCITY, shape_functions);

    // a
    array_1d< double, TDim > acceleration;
    this->EvaluateInPoint(acceleration, AUX_ADJOINT_ACCELERATION, shape_functions);

    array_1d<double, TNumNodes> velocity_grad_shape_function_derivatives;
    this->GetConvectionOperator(
        velocity_grad_shape_function_derivatives,
        velocity,
        shape_function_derivatives);

    const double velocity_norm = norm_2(velocity);
    const double elem_size = this->CalculateElementSize(volume);
    double tau_one, tau_two;

    this->CalculateStabilizationParameters(
        tau_one,
        tau_two,
        velocity_norm,
        elem_size,
        density,
        viscosity,
        rCurrentProcessInfo
    );

    double tau_one_time_gradient;
    this->CalculateTauOneTimeGradient(
        tau_one_time_gradient,
        tau_one,
        density,
        elem_size,
        velocity_norm,
        velocity,
        acceleration
    );

    BoundedMatrix<double, TNumNodes, TDim> tau_one_primal_derivatives;
    this->CalculateTauOnePrimalDerivatives(
        tau_one_primal_derivatives,
        tau_one,
        density,
        elem_size,
        velocity_norm,
        shape_functions,
        velocity
    );

    BoundedMatrix<double, TNumNodes, TDim> tau_one_time_gradient_primal_derivatives;
    this->CalculateTauOneTimeGradientPrimalDerivatives(
        tau_one_time_gradient_primal_derivatives,
        tau_one,
        density,
        elem_size,
        velocity_norm,
        shape_functions,
        velocity,
        acceleration
    );

    const double coeff = density;

    for (unsigned int a = 0; a < TNumNodes; ++a)
    {
        for (unsigned int c = 0; c < TNumNodes; ++c)
        {
            for (unsigned int k = 0; k < TDim; ++k)
            {
                double value = 0.0;

                value += velocity_grad_shape_function_derivatives[a] * tau_one_time_gradient_primal_derivatives(c,k);

                rMatrix(a, c*TDim + k) = coeff * value * volume;
            }
        }
    }

    KRATOS_CATCH("");
}

void GetAdjointVelocity( BoundedVector<double, TCoordLocalSize>& rAdjointVelocity )
{
    rAdjointVelocity.clear();
    unsigned int local_index = 0;
    for (unsigned int iNode = 0; iNode < TNumNodes; ++iNode)
    {
        const array_1d<double, 3>& adjoint_vector = this->GetGeometry()[iNode].FastGetSolutionStepValue(ADJOINT_FLUID_VECTOR_1, 0);
        for (unsigned int i=0; i < TDim; ++i)
        {
            rAdjointVelocity[local_index] = adjoint_vector[i];
            local_index++;
        }
    }
}

void GetAdjointPressure( BoundedVector<double, TNumNodes>& rAdjointPressure )
{
    rAdjointPressure.clear();
    for (unsigned int iNode = 0; iNode < TNumNodes; ++iNode)
    {
        const double value = this->GetGeometry()[iNode].FastGetSolutionStepValue(ADJOINT_FLUID_SCALAR_1, 0);
        rAdjointPressure[iNode] = value;
    }
}

void ComputeMatrix( const unsigned int MatrixId,
                    const ProcessInfo& rCurrentProcessInfo )
{
    BoundedVector<double, TCoordLocalSize> adjoint_velocity;
    GetAdjointVelocity(adjoint_velocity);

    BoundedVector<double, TNumNodes> adjoint_pressure;
    GetAdjointPressure(adjoint_pressure);

    VectorType adjoint_vector;
    this->GetValuesVector(adjoint_vector);

    Element::MatrixType current_matrix;

    Variable<double> const* var_matrix_energy{};
    Variable<double> const* var_eigen_min{};
    Variable<double> const* var_eigen_max{};

    double matrix_energy(0.0), eigen_min(0.0), eigen_max(0.0);

    // Kuu matrices
    switch (MatrixId)
    {
        case 1:
        {
            Matrix_1(current_matrix, rCurrentProcessInfo);
            var_matrix_energy = &STABILIZATION_ANALYSIS_MATRIX_1_ENERGY;
            var_eigen_min = &STABILIZATION_ANALYSIS_MATRIX_1_EIGEN_MIN;
            var_eigen_max = &STABILIZATION_ANALYSIS_MATRIX_1_EIGEN_MAX;
            break;
        }
        case 2:
        {
            Matrix_2(current_matrix, rCurrentProcessInfo);
            var_matrix_energy = &STABILIZATION_ANALYSIS_MATRIX_2_ENERGY;
            var_eigen_min = &STABILIZATION_ANALYSIS_MATRIX_2_EIGEN_MIN;
            var_eigen_max = &STABILIZATION_ANALYSIS_MATRIX_2_EIGEN_MAX;
            break;
        }
        case 3:
        {
            Matrix_3(current_matrix, rCurrentProcessInfo);
            var_matrix_energy = &STABILIZATION_ANALYSIS_MATRIX_3_ENERGY;
            var_eigen_min = &STABILIZATION_ANALYSIS_MATRIX_3_EIGEN_MIN;
            var_eigen_max = &STABILIZATION_ANALYSIS_MATRIX_3_EIGEN_MAX;
            break;
        }
        case 4:
        {
            Matrix_4(current_matrix, rCurrentProcessInfo);
            var_matrix_energy = &STABILIZATION_ANALYSIS_MATRIX_4_ENERGY;
            var_eigen_min = &STABILIZATION_ANALYSIS_MATRIX_4_EIGEN_MIN;
            var_eigen_max = &STABILIZATION_ANALYSIS_MATRIX_4_EIGEN_MAX;

            break;
        }
        case 5:
        {
            Matrix_5(current_matrix, rCurrentProcessInfo);
            var_matrix_energy = &STABILIZATION_ANALYSIS_MATRIX_5_ENERGY;
            var_eigen_min = &STABILIZATION_ANALYSIS_MATRIX_5_EIGEN_MIN;
            var_eigen_max = &STABILIZATION_ANALYSIS_MATRIX_5_EIGEN_MAX;
            break;
        }
        case 6:
        {
            Matrix_6(current_matrix, rCurrentProcessInfo);
            var_matrix_energy = &STABILIZATION_ANALYSIS_MATRIX_6_ENERGY;
            var_eigen_min = &STABILIZATION_ANALYSIS_MATRIX_6_EIGEN_MIN;
            var_eigen_max = &STABILIZATION_ANALYSIS_MATRIX_6_EIGEN_MAX;
            break;
        }
        case 7:
        {
            Matrix_7(current_matrix, rCurrentProcessInfo);
            var_matrix_energy = &STABILIZATION_ANALYSIS_MATRIX_7_ENERGY;
            var_eigen_min = &STABILIZATION_ANALYSIS_MATRIX_7_EIGEN_MIN;
            var_eigen_max = &STABILIZATION_ANALYSIS_MATRIX_7_EIGEN_MAX;
            break;
        }
        case 8:
        {
            Matrix_8(current_matrix, rCurrentProcessInfo);
            var_matrix_energy = &STABILIZATION_ANALYSIS_MATRIX_8_ENERGY;
            var_eigen_min = &STABILIZATION_ANALYSIS_MATRIX_8_EIGEN_MIN;
            var_eigen_max = &STABILIZATION_ANALYSIS_MATRIX_8_EIGEN_MAX;
            break;
        }
        case 9:
        {
            Matrix_9(current_matrix, rCurrentProcessInfo);
            var_matrix_energy = &STABILIZATION_ANALYSIS_MATRIX_9_ENERGY;
            var_eigen_min = &STABILIZATION_ANALYSIS_MATRIX_9_EIGEN_MIN;
            var_eigen_max = &STABILIZATION_ANALYSIS_MATRIX_9_EIGEN_MAX;
            break;
        }
        case 10:
        {
            Matrix_10(current_matrix, rCurrentProcessInfo);
            var_matrix_energy = &STABILIZATION_ANALYSIS_MATRIX_10_ENERGY;
            var_eigen_min = &STABILIZATION_ANALYSIS_MATRIX_10_EIGEN_MIN;
            var_eigen_max = &STABILIZATION_ANALYSIS_MATRIX_10_EIGEN_MAX;
            break;
        }
        case 11:
        {
            Matrix_11(current_matrix, rCurrentProcessInfo);
            var_matrix_energy = &STABILIZATION_ANALYSIS_MATRIX_11_ENERGY;
            var_eigen_min = &STABILIZATION_ANALYSIS_MATRIX_11_EIGEN_MIN;
            var_eigen_max = &STABILIZATION_ANALYSIS_MATRIX_11_EIGEN_MAX;
            break;
        }
        case 12:
        {
            Matrix_12(current_matrix, rCurrentProcessInfo);
            var_matrix_energy = &STABILIZATION_ANALYSIS_MATRIX_12_ENERGY;
            var_eigen_min = &STABILIZATION_ANALYSIS_MATRIX_12_EIGEN_MIN;
            var_eigen_max = &STABILIZATION_ANALYSIS_MATRIX_12_EIGEN_MAX;
            break;
        }
    }

    if (MatrixId<13)
    {
        matrix_energy = CalculateMatrixEnergy(adjoint_velocity, current_matrix, adjoint_velocity);
    }

    // Kup matrices
    switch (MatrixId)
    {
        case 13:
        {
            Matrix_13(current_matrix, rCurrentProcessInfo);
            var_matrix_energy = &STABILIZATION_ANALYSIS_MATRIX_13_ENERGY;
            var_eigen_min = &STABILIZATION_ANALYSIS_MATRIX_13_EIGEN_MIN;
            var_eigen_max = &STABILIZATION_ANALYSIS_MATRIX_13_EIGEN_MAX;
            break;
        }
        case 14:
        {
            Matrix_14(current_matrix, rCurrentProcessInfo);
            var_matrix_energy = &STABILIZATION_ANALYSIS_MATRIX_14_ENERGY;
            var_eigen_min = &STABILIZATION_ANALYSIS_MATRIX_14_EIGEN_MIN;
            var_eigen_max = &STABILIZATION_ANALYSIS_MATRIX_14_EIGEN_MAX;
            break;
        }
        case 15:
        {
            Matrix_15(current_matrix, rCurrentProcessInfo);
            var_matrix_energy = &STABILIZATION_ANALYSIS_MATRIX_15_ENERGY;
            var_eigen_min = &STABILIZATION_ANALYSIS_MATRIX_15_EIGEN_MIN;
            var_eigen_max = &STABILIZATION_ANALYSIS_MATRIX_15_EIGEN_MAX;
            break;
        }
        case 16:
        {
            Matrix_16(current_matrix, rCurrentProcessInfo);
            var_matrix_energy = &STABILIZATION_ANALYSIS_MATRIX_16_ENERGY;
            var_eigen_min = &STABILIZATION_ANALYSIS_MATRIX_16_EIGEN_MIN;
            var_eigen_max = &STABILIZATION_ANALYSIS_MATRIX_16_EIGEN_MAX;
            break;
        }
        case 17:
        {
            Matrix_17(current_matrix, rCurrentProcessInfo);
            var_matrix_energy = &STABILIZATION_ANALYSIS_MATRIX_17_ENERGY;
            var_eigen_min = &STABILIZATION_ANALYSIS_MATRIX_17_EIGEN_MIN;
            var_eigen_max = &STABILIZATION_ANALYSIS_MATRIX_17_EIGEN_MAX;
            break;
        }
        case 18:
        {
            Matrix_18(current_matrix, rCurrentProcessInfo);
            var_matrix_energy = &STABILIZATION_ANALYSIS_MATRIX_18_ENERGY;
            var_eigen_min = &STABILIZATION_ANALYSIS_MATRIX_18_EIGEN_MIN;
            var_eigen_max = &STABILIZATION_ANALYSIS_MATRIX_18_EIGEN_MAX;
            break;
        }
        case 19:
        {
            Matrix_19(current_matrix, rCurrentProcessInfo);
            var_matrix_energy = &STABILIZATION_ANALYSIS_MATRIX_19_ENERGY;
            var_eigen_min = &STABILIZATION_ANALYSIS_MATRIX_19_EIGEN_MIN;
            var_eigen_max = &STABILIZATION_ANALYSIS_MATRIX_19_EIGEN_MAX;
            break;
        }
    }
    if (MatrixId>=13 and MatrixId < 20)
    {
        matrix_energy = CalculateMatrixEnergy(adjoint_pressure, current_matrix, adjoint_velocity);
    }

    // Kpu matrices
    switch (MatrixId)
    {
        case 20:
        {
            Matrix_20(current_matrix, rCurrentProcessInfo);
            var_matrix_energy = &STABILIZATION_ANALYSIS_MATRIX_20_ENERGY;
            var_eigen_min = &STABILIZATION_ANALYSIS_MATRIX_20_EIGEN_MIN;
            var_eigen_max = &STABILIZATION_ANALYSIS_MATRIX_20_EIGEN_MAX;
            break;
        }
        case 21:
        {
            Matrix_21(current_matrix, rCurrentProcessInfo);
            var_matrix_energy = &STABILIZATION_ANALYSIS_MATRIX_21_ENERGY;
            var_eigen_min = &STABILIZATION_ANALYSIS_MATRIX_21_EIGEN_MIN;
            var_eigen_max = &STABILIZATION_ANALYSIS_MATRIX_21_EIGEN_MAX;
            break;
        }
    }

    if (MatrixId>=20 and MatrixId < 22)
    {
        matrix_energy = CalculateMatrixEnergy(adjoint_velocity, current_matrix, adjoint_pressure);
    }

    if (MatrixId == 22)
    {
        Matrix_22(current_matrix, rCurrentProcessInfo);
        var_matrix_energy = &STABILIZATION_ANALYSIS_MATRIX_22_ENERGY;
        var_eigen_min = &STABILIZATION_ANALYSIS_MATRIX_22_EIGEN_MIN;
        var_eigen_max = &STABILIZATION_ANALYSIS_MATRIX_22_EIGEN_MAX;
        matrix_energy = CalculateMatrixEnergy(adjoint_pressure, current_matrix, adjoint_pressure);
    }

    double alpha, beta, gamma, delta_time;
    CalculateBossakSchemeConstants( alpha, beta, gamma, delta_time, rCurrentProcessInfo);

    // Validation matrices from the element code
    if (MatrixId == 23)
    {
        this->CalculatePrimalGradientOfVMSSteadyTerm(current_matrix, rCurrentProcessInfo);
        matrix_energy = (gamma/(beta*delta_time))*CalculateMatrixEnergy(adjoint_vector, current_matrix, adjoint_vector);
        var_matrix_energy = &STABILIZATION_ANALYSIS_MATRIX_23_ENERGY;
        var_eigen_min = &STABILIZATION_ANALYSIS_MATRIX_23_EIGEN_MIN;
        var_eigen_max = &STABILIZATION_ANALYSIS_MATRIX_23_EIGEN_MAX;
    }

    if (MatrixId == 24)
    {
        current_matrix = this->GetValue(STABILIZATION_ANALYSIS_MATRIX_24);
        matrix_energy = CalculateMatrixEnergy(adjoint_vector, current_matrix, adjoint_vector);
        var_matrix_energy = &STABILIZATION_ANALYSIS_MATRIX_24_ENERGY;
        var_eigen_min = &STABILIZATION_ANALYSIS_MATRIX_24_EIGEN_MIN;
        var_eigen_max = &STABILIZATION_ANALYSIS_MATRIX_24_EIGEN_MAX;
    }

    if (MatrixId == 25)
    {
        this->Matrix_25(current_matrix, rCurrentProcessInfo);
        matrix_energy = (gamma/(beta*delta_time))*CalculateMatrixEnergy(adjoint_velocity, current_matrix, adjoint_velocity);
        var_matrix_energy = &STABILIZATION_ANALYSIS_MATRIX_25_ENERGY;
        var_eigen_min = &STABILIZATION_ANALYSIS_MATRIX_25_EIGEN_MIN;
        var_eigen_max = &STABILIZATION_ANALYSIS_MATRIX_25_EIGEN_MAX;
    }

    if (MatrixId == 26)
    {
        this->Matrix_26(current_matrix, rCurrentProcessInfo);
        matrix_energy = (gamma/(beta*delta_time))*CalculateMatrixEnergy(adjoint_pressure, current_matrix, adjoint_velocity);
        var_matrix_energy = &STABILIZATION_ANALYSIS_MATRIX_26_ENERGY;
        var_eigen_min = &STABILIZATION_ANALYSIS_MATRIX_26_EIGEN_MIN;
        var_eigen_max = &STABILIZATION_ANALYSIS_MATRIX_26_EIGEN_MAX;
    }

    if (MatrixId == 27)
    {
        this->Matrix_27(current_matrix, rCurrentProcessInfo);
        matrix_energy = (gamma/(beta*delta_time))*CalculateMatrixEnergy(adjoint_velocity, current_matrix, adjoint_velocity);
        var_matrix_energy = &STABILIZATION_ANALYSIS_MATRIX_27_ENERGY;
        var_eigen_min = &STABILIZATION_ANALYSIS_MATRIX_27_EIGEN_MIN;
        var_eigen_max = &STABILIZATION_ANALYSIS_MATRIX_27_EIGEN_MAX;
    }

    if (MatrixId == 28)
    {
        this->Matrix_28(current_matrix, rCurrentProcessInfo);
        matrix_energy = (gamma/(beta*delta_time))*CalculateMatrixEnergy(adjoint_pressure, current_matrix, adjoint_velocity);
        var_matrix_energy = &STABILIZATION_ANALYSIS_MATRIX_28_ENERGY;
        var_eigen_min = &STABILIZATION_ANALYSIS_MATRIX_28_EIGEN_MIN;
        var_eigen_max = &STABILIZATION_ANALYSIS_MATRIX_28_EIGEN_MAX;
    }

    if (MatrixId == 29)
    {
        current_matrix = this->GetValue(STABILIZATION_ANALYSIS_MATRIX_29);
        matrix_energy = (gamma/(beta*delta_time))*CalculateMatrixEnergy(adjoint_vector, current_matrix, adjoint_vector);
        var_matrix_energy = &STABILIZATION_ANALYSIS_MATRIX_29_ENERGY;
        var_eigen_min = &STABILIZATION_ANALYSIS_MATRIX_29_EIGEN_MIN;
        var_eigen_max = &STABILIZATION_ANALYSIS_MATRIX_29_EIGEN_MAX;
    }

    CalculateEigenValues(current_matrix, eigen_min, eigen_max);

    this->SetValue(*var_eigen_min, eigen_min);
    this->SetValue(*var_eigen_max, eigen_max);
    this->SetValue(*var_matrix_energy, matrix_energy);
}

void ProcessSymmetricMatrices(const ProcessInfo& rCurrentProcessInfo)
{
    for (unsigned int iMatrix = 1; iMatrix < 30; ++iMatrix)
    {
        ComputeMatrix(iMatrix, rCurrentProcessInfo);
    }
}
