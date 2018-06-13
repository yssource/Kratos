#include <cmath>
#include <vector>

#include "includes/process_info.h"

template<unsigned int TDim>
double CalculateMatrixEnergy(const Element::MatrixType& rMatrix, const Vector& rVector )
{
    const Vector& temp = prod(rMatrix, rVector);
    double energy = inner_prod(temp, rVector);

    return energy;
}

template<unsigned int TDim>
void WriteMatrix(
    const unsigned int MatrixId,
    const Element::MatrixType& rMatrix,
    VMSAdjointElement<TDim>* pElement,
    const ProcessInfo& rCurrentProcessInfo
)
{
    const unsigned int matrix_size = rMatrix.size1();

    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>  characteristic_matrix;
    characteristic_matrix.resize(matrix_size, matrix_size);

    for (unsigned int i = 0; i < matrix_size; i++)
        for (unsigned int j = 0; j < matrix_size; j++)
            characteristic_matrix(i,j) = rMatrix(i,j);
   
    Eigen::EigenSolver<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> es(characteristic_matrix, false);
    const auto& values = es.eigenvalues().real();

    double eigen_min = values[0];
    double eigen_max = values[0];

    double temp;

    for (unsigned int i = 0; i < values.size(); i++)
    {
        temp = values[i];
        if (eigen_min > temp)
            eigen_min = temp;
        
        if (eigen_max < temp)
            eigen_max = temp;
    }

    double matrix_energy = 0.0;
    unsigned int time_step = 0;

    switch (MatrixId)
    {
        case 1:
        case 2:
        case 3:
        case 4:
        case 5:
        case 6:
        case 7:
        case 8:
        case 9:
        case 10:
        case 11:
        case 12:
        case 13:
        case 15:
        case 16:
        {
            BoundedVector<double, (TDim+1)*TDim> values;

            BoundedVector<array_1d<double, TDim>, TDim+1> nodal_velocity_vectors;
            for (unsigned int iNode = 0; iNode < TDim+1; iNode++)
            {
                nodal_velocity_vectors[iNode] = pElement->GetGeometry()[iNode].FastGetSolutionStepValue(ADJOINT_FLUID_VECTOR_1, time_step);
                for (unsigned int iDim = 0; iDim < TDim; iDim++)
                    values[iNode*TDim+iDim] = nodal_velocity_vectors[iNode][iDim];
            }
            matrix_energy = CalculateMatrixEnergy<TDim>(rMatrix, values);
            break;
        }
        case 14:
        {
            BoundedVector<double, TDim+1> values;

            for (unsigned int iNode = 0; iNode < TDim+1; iNode++)
            {
                values[iNode] = pElement->GetGeometry()[iNode].FastGetSolutionStepValue(ADJOINT_FLUID_SCALAR_1, time_step);
            }
            matrix_energy = CalculateMatrixEnergy<TDim>(rMatrix, values);
            break;
        }
        case 17:
        {
            Vector adjoint_values;
            pElement->GetValuesVector(adjoint_values, time_step);
            matrix_energy = CalculateMatrixEnergy<TDim>(rMatrix, adjoint_values);
            break;
        }
        case 18:
        {
            Vector adjoint_values;
            pElement->GetValuesVector(adjoint_values, time_step+3);
            matrix_energy = CalculateMatrixEnergy<TDim>(rMatrix, adjoint_values);
            break;
        }
    }

    switch (MatrixId)
    {
        case 1:
            pElement->SetValue(SYMMETRIC_MATRIX_1_EIGEN_MIN, eigen_min);
            pElement->SetValue(SYMMETRIC_MATRIX_1_EIGEN_MAX, eigen_max);
            pElement->SetValue(SYMMETRIC_MATRIX_1_ENERGY, matrix_energy);
            break;
        case 2:
            pElement->SetValue(SYMMETRIC_MATRIX_2_EIGEN_MIN, eigen_min);
            pElement->SetValue(SYMMETRIC_MATRIX_2_EIGEN_MAX, eigen_max);
            pElement->SetValue(SYMMETRIC_MATRIX_2_ENERGY, matrix_energy);
            break;
        case 3:
            pElement->SetValue(SYMMETRIC_MATRIX_3_EIGEN_MIN, eigen_min);
            pElement->SetValue(SYMMETRIC_MATRIX_3_EIGEN_MAX, eigen_max);
            pElement->SetValue(SYMMETRIC_MATRIX_3_ENERGY, matrix_energy);
            break;
        case 4:
            pElement->SetValue(SYMMETRIC_MATRIX_4_EIGEN_MIN, eigen_min);
            pElement->SetValue(SYMMETRIC_MATRIX_4_EIGEN_MAX, eigen_max);
            pElement->SetValue(SYMMETRIC_MATRIX_4_ENERGY, matrix_energy);
            break;
        case 5:
            pElement->SetValue(SYMMETRIC_MATRIX_5_EIGEN_MIN, eigen_min);
            pElement->SetValue(SYMMETRIC_MATRIX_5_EIGEN_MAX, eigen_max);
            pElement->SetValue(SYMMETRIC_MATRIX_5_ENERGY, matrix_energy);
            break;
        case 6:
            pElement->SetValue(SYMMETRIC_MATRIX_6_EIGEN_MIN, eigen_min);
            pElement->SetValue(SYMMETRIC_MATRIX_6_EIGEN_MAX, eigen_max);
            pElement->SetValue(SYMMETRIC_MATRIX_6_ENERGY, matrix_energy);
            break;
        case 7:
            pElement->SetValue(SYMMETRIC_MATRIX_7_EIGEN_MIN, eigen_min);
            pElement->SetValue(SYMMETRIC_MATRIX_7_EIGEN_MAX, eigen_max);
            pElement->SetValue(SYMMETRIC_MATRIX_7_ENERGY, matrix_energy);
            break;
        case 8:
            pElement->SetValue(SYMMETRIC_MATRIX_8_EIGEN_MIN, eigen_min);
            pElement->SetValue(SYMMETRIC_MATRIX_8_EIGEN_MAX, eigen_max);
            pElement->SetValue(SYMMETRIC_MATRIX_8_ENERGY, matrix_energy);
            break;
        case 9:
            pElement->SetValue(SYMMETRIC_MATRIX_9_EIGEN_MIN, eigen_min);
            pElement->SetValue(SYMMETRIC_MATRIX_9_EIGEN_MAX, eigen_max);
            pElement->SetValue(SYMMETRIC_MATRIX_9_ENERGY, matrix_energy);
            break;
        case 10:
            pElement->SetValue(SYMMETRIC_MATRIX_10_EIGEN_MIN, eigen_min);
            pElement->SetValue(SYMMETRIC_MATRIX_10_EIGEN_MAX, eigen_max);
            pElement->SetValue(SYMMETRIC_MATRIX_10_ENERGY, matrix_energy);
            break;
        case 11:
            pElement->SetValue(SYMMETRIC_MATRIX_11_EIGEN_MIN, eigen_min);
            pElement->SetValue(SYMMETRIC_MATRIX_11_EIGEN_MAX, eigen_max);
            pElement->SetValue(SYMMETRIC_MATRIX_11_ENERGY, matrix_energy);
            break;
        case 12:
            pElement->SetValue(SYMMETRIC_MATRIX_12_EIGEN_MIN, eigen_min);
            pElement->SetValue(SYMMETRIC_MATRIX_12_EIGEN_MAX, eigen_max);
            pElement->SetValue(SYMMETRIC_MATRIX_12_ENERGY, matrix_energy);
            break;
        case 13:
            pElement->SetValue(SYMMETRIC_MATRIX_13_EIGEN_MIN, eigen_min);
            pElement->SetValue(SYMMETRIC_MATRIX_13_EIGEN_MAX, eigen_max);
            pElement->SetValue(SYMMETRIC_MATRIX_13_ENERGY, matrix_energy);
            break;
        case 14:
            pElement->SetValue(SYMMETRIC_MATRIX_14_EIGEN_MIN, eigen_min);
            pElement->SetValue(SYMMETRIC_MATRIX_14_EIGEN_MAX, eigen_max);
            pElement->SetValue(SYMMETRIC_MATRIX_14_ENERGY, matrix_energy);
            break;
        case 15:
            pElement->SetValue(SYMMETRIC_MATRIX_15_EIGEN_MIN, eigen_min);
            pElement->SetValue(SYMMETRIC_MATRIX_15_EIGEN_MAX, eigen_max);
            pElement->SetValue(SYMMETRIC_MATRIX_15_ENERGY, matrix_energy);
            break;                                    
        case 16:
            pElement->SetValue(SYMMETRIC_MATRIX_16_EIGEN_MIN, eigen_min);
            pElement->SetValue(SYMMETRIC_MATRIX_16_EIGEN_MAX, eigen_max);
            pElement->SetValue(SYMMETRIC_MATRIX_16_ENERGY, matrix_energy);
            break;
        case 17:
            pElement->SetValue(SYMMETRIC_MATRIX_17_EIGEN_MIN, eigen_min);
            pElement->SetValue(SYMMETRIC_MATRIX_17_EIGEN_MAX, eigen_max);
            pElement->SetValue(SYMMETRIC_MATRIX_17_ENERGY, matrix_energy);
            break;
        case 18:
            pElement->SetValue(SYMMETRIC_MATRIX_18_EIGEN_MIN, eigen_min);
            pElement->SetValue(SYMMETRIC_MATRIX_18_EIGEN_MAX, eigen_max);
            pElement->SetValue(SYMMETRIC_MATRIX_18_ENERGY, matrix_energy);
            break;            
    }

}

template<unsigned int TDim>
double dot(
    const array_1d< double, TDim >& vec1,
    const array_1d< double, TDim >& vec2
) 
{
    double result = 0.0;

    for (unsigned int i=0; i < TDim; i++)
        result += vec1[i]*vec2[i];
    
    return result;
}

template<unsigned int TDim>
void Matrix_1(
            VMSAdjointElement<TDim>* pElement,
            array_1d< double, TDim >& Velocity,
            double Density,
            double DynamicViscosity, 
            ProcessInfo& rCurrentProcessInfo)
{
    constexpr unsigned int TNumNodes = TDim + 1;
    constexpr unsigned int TCoordLocalSize = TDim * TNumNodes;

    Element::MatrixType resultant_matrix;

    resultant_matrix.resize(TCoordLocalSize, TCoordLocalSize);
    resultant_matrix.clear();

    BoundedMatrix<double, TNumNodes, TDim> dn_dx;
    array_1d< double, TNumNodes > N;
    double volume;
    GeometryUtils::CalculateGeometryData(pElement->GetGeometry(),dn_dx,N,volume);

    unsigned int matrix_index_i, matrix_index_j;

    for (unsigned int a=0; a < TNumNodes; a++)
    {
        for (unsigned int r=0; r < TNumNodes; r++)
        {
            // delta_im, therefore always i = m
            for (unsigned int i=0; i < TDim; i++)
            {
                unsigned int m = i;

                matrix_index_i = a*TDim + i;
                matrix_index_j = r*TDim + m;

                double value = 0.0;

                array_1d< double, TNumNodes > DNr_DX;
                for (unsigned int k = 0; k < TDim; k++)
                    DNr_DX[k] = dn_dx(r,k);

                value += N[a]*dot<TDim>(Velocity, DNr_DX);

                array_1d< double, TNumNodes > DNa_DX;
                for (unsigned int k = 0; k < TDim; k++)
                    DNa_DX[k] = dn_dx(a,k);

                value += N[r]*dot<TDim>(Velocity, DNa_DX);

                value *= 0.5*Density;

                resultant_matrix(matrix_index_i, matrix_index_j) = value*volume;
            }
        }
    }

    WriteMatrix(1, resultant_matrix, pElement, rCurrentProcessInfo);
}

template<unsigned int TDim>
void Matrix_2(
    VMSAdjointElement<TDim>* pElement,
    array_1d< double, TDim >& Velocity,
    double Density,
    double DynamicViscosity, 
    ProcessInfo& rCurrentProcessInfo
)
{
    constexpr unsigned int TNumNodes = TDim + 1;
    constexpr unsigned int TCoordLocalSize = TDim * TNumNodes;

    Element::MatrixType resultant_matrix;

    resultant_matrix.resize(TCoordLocalSize, TCoordLocalSize);
    resultant_matrix.clear();

    BoundedMatrix<double, TNumNodes, TDim> dn_dx;
    array_1d< double, TNumNodes > N;
    double volume;
    GeometryUtils::CalculateGeometryData(pElement->GetGeometry(),dn_dx,N,volume);

    unsigned int matrix_index_i, matrix_index_j;

    BoundedVector<array_1d<double, TDim>, TNumNodes> nodal_velocity_vectors;
    for (unsigned int iNode = 0; iNode < TNumNodes; iNode++)
         nodal_velocity_vectors[iNode] = pElement->GetGeometry()[iNode].FastGetSolutionStepValue(VELOCITY);

    for (unsigned int a=0; a < TNumNodes; a++)
    {
        for (unsigned int r=0; r < TNumNodes; r++)
        {
            for (unsigned int i=0; i < TDim; i++)
            {
                for (unsigned int m=0; m < TDim; m++)
                {
                    matrix_index_i = a*TDim + i;
                    matrix_index_j = r*TDim + m;

                    double value = 0.0;

                    for (unsigned int k = 0; k < TNumNodes; k++)
                    {
                        value += N[a]*N[r]*dn_dx(k,m)*nodal_velocity_vectors[k][i];
                        value += N[r]*N[a]*dn_dx(k,i)*nodal_velocity_vectors[k][m];
                    }

                    value *= 0.5*Density;

                    resultant_matrix(matrix_index_i, matrix_index_j) = value*volume;
                }
            }
        }
    }

    WriteMatrix(2, resultant_matrix, pElement, rCurrentProcessInfo);
}

template<unsigned int TDim>
void Matrix_3(
    VMSAdjointElement<TDim>* pElement,
    array_1d< double, TDim >& Velocity,
    double Density,
    double DynamicViscosity, 
    ProcessInfo& rCurrentProcessInfo,
    double TauOne
)
{
    constexpr unsigned int TNumNodes = TDim + 1;
    constexpr unsigned int TCoordLocalSize = TDim * TNumNodes;

    Element::MatrixType resultant_matrix;

    resultant_matrix.resize(TCoordLocalSize, TCoordLocalSize);
    resultant_matrix.clear();

    BoundedMatrix<double, TNumNodes, TDim> dn_dx;
    array_1d< double, TNumNodes > N;
    double volume;
    GeometryUtils::CalculateGeometryData(pElement->GetGeometry(),dn_dx,N,volume);

    unsigned int matrix_index_i, matrix_index_j;

    BoundedVector<array_1d<double, TDim>, TNumNodes> nodal_velocity_vectors;
    for (unsigned int iNode = 0; iNode < TNumNodes; iNode++)
         nodal_velocity_vectors[iNode] = pElement->GetGeometry()[iNode].FastGetSolutionStepValue(VELOCITY);

         

    for (unsigned int a=0; a < TNumNodes; a++)
    {
        for (unsigned int r=0; r < TNumNodes; r++)
        {
            for (unsigned int i=0; i < TDim; i++)
            {
                for (unsigned int m=0; m < TDim; m++)
                {
                    matrix_index_i = a*TDim + i;
                    matrix_index_j = r*TDim + m;

                    double value = 0.0;

                    array_1d<double, TDim> del_u_i;
                    for (unsigned int k = 0; k < TDim; k++)
                    {
                        del_u_i[k] = 0.0;
                        for (unsigned int l = 0; l < TNumNodes; l++)
                            del_u_i[k] += dn_dx(l,k)*nodal_velocity_vectors[l][i];
                    }

                    array_1d<double, TDim> del_u_m;
                    for (unsigned int k = 0; k < TDim; k++)
                    {
                        del_u_m[k] = 0.0;
                        for (unsigned int l = 0; l < TNumNodes; l++)
                            del_u_m[k] += dn_dx(l,k)*nodal_velocity_vectors[l][m];
                    }

                    value += N[r]*dn_dx(a,m)*dot<TDim>(Velocity, del_u_i);
                    value += N[a]*dn_dx(r,i)*dot<TDim>(Velocity, del_u_m);

                    value *= Density*Density*TauOne;

                    resultant_matrix(matrix_index_i, matrix_index_j) = 0.5*value*volume;
                }
            }
        }
    }

    WriteMatrix(3, resultant_matrix, pElement, rCurrentProcessInfo);
}

template<unsigned int TDim>
void Matrix_4(
    VMSAdjointElement<TDim>* pElement,
    array_1d< double, TDim >& Velocity,
    double Density,
    double DynamicViscosity, 
    ProcessInfo& rCurrentProcessInfo,
    Element::MatrixType& rTauOneDerivatives
)
{
    constexpr unsigned int TNumNodes = TDim + 1;
    constexpr unsigned int TCoordLocalSize = TDim * TNumNodes;

    Element::MatrixType resultant_matrix;

    resultant_matrix.resize(TCoordLocalSize, TCoordLocalSize);
    resultant_matrix.clear();

    BoundedMatrix<double, TNumNodes, TDim> dn_dx;
    array_1d< double, TNumNodes > N;
    double volume;
    GeometryUtils::CalculateGeometryData(pElement->GetGeometry(),dn_dx,N,volume);

    unsigned int matrix_index_i, matrix_index_j;

    BoundedVector<array_1d<double, TDim>, TNumNodes> nodal_velocity_vectors;
    for (unsigned int iNode = 0; iNode < TNumNodes; iNode++)
         nodal_velocity_vectors[iNode] = pElement->GetGeometry()[iNode].FastGetSolutionStepValue(VELOCITY);


    for (unsigned int a=0; a < TNumNodes; a++)
    {
        for (unsigned int r=0; r < TNumNodes; r++)
        {
            for (unsigned int i=0; i < TDim; i++)
            {
                for (unsigned int m=0; m < TDim; m++)
                {
                    matrix_index_i = a*TDim + i;
                    matrix_index_j = r*TDim + m;

                    double value1 = 0.0, value2 = 0.0;

                    double temp1=0.0, temp2 = 0.0;
                    for (unsigned int k=0; k < TDim; k++)
                    {
                        temp1 += Velocity[k]*dn_dx(a,k);
                        temp2 += Velocity[k]*dn_dx(r,k);
                    }
                    
                    value1 += temp1;
                    value2 += temp2;

                    temp1 = 0.0;
                    temp2 = 0.0;
                    for (unsigned int k=0; k < TDim; k++)
                        for (unsigned int l=0; l < TNumNodes; l++)
                        {
                            temp1 += Velocity[k]*dn_dx(l,k)*nodal_velocity_vectors[l][i];
                            temp2 += Velocity[k]*dn_dx(l,k)*nodal_velocity_vectors[l][m];
                        }
                    
                    value1 *= temp1;
                    value2 *= temp2;

                    value1 *= rTauOneDerivatives(r,m);
                    value2 *= rTauOneDerivatives(a,i);

                    value1 *= Density*Density;
                    value2 *= Density*Density;

                    resultant_matrix(matrix_index_i, matrix_index_j) = 0.5*(value1 + value2)*volume;
                }
            }
        }
    }

    WriteMatrix(4, resultant_matrix, pElement, rCurrentProcessInfo);
}

template<unsigned int TDim>
void Matrix_5(
    VMSAdjointElement<TDim>* pElement,
    array_1d< double, TDim >& Velocity,
    double Density,
    double DynamicViscosity, 
    ProcessInfo& rCurrentProcessInfo,
    double TauOne
)
{
    constexpr unsigned int TNumNodes = TDim + 1;
    constexpr unsigned int TCoordLocalSize = TDim * TNumNodes;

    Element::MatrixType resultant_matrix;

    resultant_matrix.resize(TCoordLocalSize, TCoordLocalSize);
    resultant_matrix.clear();

    BoundedMatrix<double, TNumNodes, TDim> dn_dx;
    array_1d< double, TNumNodes > N;
    double volume;
    GeometryUtils::CalculateGeometryData(pElement->GetGeometry(),dn_dx,N,volume);

    unsigned int matrix_index_i, matrix_index_j;

    BoundedVector<array_1d<double, TDim>, TNumNodes> nodal_velocity_vectors;
    for (unsigned int iNode = 0; iNode < TNumNodes; iNode++)
         nodal_velocity_vectors[iNode] = pElement->GetGeometry()[iNode].FastGetSolutionStepValue(VELOCITY);

         

    for (unsigned int a=0; a < TNumNodes; a++)
    {
        for (unsigned int r=0; r < TNumNodes; r++)
        {
            for (unsigned int i=0; i < TDim; i++)
            {
                for (unsigned int m=0; m < TDim; m++)
                {
                    matrix_index_i = a*TDim + i;
                    matrix_index_j = r*TDim + m;

                    double value1 = 0.0, value2 = 0.0;

                    double temp1=0.0, temp2 = 0.0;
                    for (unsigned int k=0; k < TDim; k++)
                    {
                        temp1 += Velocity[k]*dn_dx(a,k);
                        temp2 += Velocity[k]*dn_dx(r,k);
                    }
                    
                    value1 += temp1;
                    value2 += temp2;

                    temp1 = 0.0;
                    temp2 = 0.0;
                    for (unsigned int k=0; k < TNumNodes; k++)
                    {
                        temp1 += dn_dx(k,m)*nodal_velocity_vectors[k][i];
                        temp2 += dn_dx(k,i)*nodal_velocity_vectors[k][m];
                    }                    
                    
                    value1 *= temp1;
                    value2 *= temp2;

                    value1 *= N[r] * Density * Density * TauOne;
                    value2 *= N[a] * Density * Density * TauOne;

                    resultant_matrix(matrix_index_i, matrix_index_j) = 0.5*(value1 + value2)*volume;
                }
            }
        }
    }

    WriteMatrix(5, resultant_matrix, pElement, rCurrentProcessInfo);
}

template<unsigned int TDim>
void Matrix_6(
    VMSAdjointElement<TDim>* pElement,
    array_1d< double, TDim >& Velocity,
    double Density,
    double DynamicViscosity, 
    ProcessInfo& rCurrentProcessInfo,
    double TauOne
)
{
    constexpr unsigned int TNumNodes = TDim + 1;
    constexpr unsigned int TCoordLocalSize = TDim * TNumNodes;

    Element::MatrixType resultant_matrix;

    resultant_matrix.resize(TCoordLocalSize, TCoordLocalSize);
    resultant_matrix.clear();

    BoundedMatrix<double, TNumNodes, TDim> dn_dx;
    array_1d< double, TNumNodes > N;
    double volume;
    GeometryUtils::CalculateGeometryData(pElement->GetGeometry(),dn_dx,N,volume);

    unsigned int matrix_index_i, matrix_index_j;

    for (unsigned int a=0; a < TNumNodes; a++)
    {
        for (unsigned int r=0; r < TNumNodes; r++)
        {
            for (unsigned int i=0; i < TDim; i++)
            {
                unsigned int m=i;

                matrix_index_i = a*TDim + i;
                matrix_index_j = r*TDim + m;

                double value = 0.0;
                
                for (unsigned int k=0; k < TDim; k++)
                    value += Velocity[k]*dn_dx(a,k);
                double temp = 0.0;
                for (unsigned int k=0; k < TDim; k++)
                    temp += Velocity[k]*dn_dx(r,k);
                
                value *= temp*TauOne*Density*Density;

                resultant_matrix(matrix_index_i, matrix_index_j) = value*volume;
            }
        }
    }
    
    WriteMatrix(6, resultant_matrix, pElement, rCurrentProcessInfo);
}

template<unsigned int TDim>
void Matrix_7(
    VMSAdjointElement<TDim>* pElement,
    array_1d< double, TDim >& Velocity,
    double Density,
    double DynamicViscosity, 
    ProcessInfo& rCurrentProcessInfo,
    double TauTwo
)
{
    constexpr unsigned int TNumNodes = TDim + 1;
    constexpr unsigned int TCoordLocalSize = TDim * TNumNodes;

    Element::MatrixType resultant_matrix;

    resultant_matrix.resize(TCoordLocalSize, TCoordLocalSize);
    resultant_matrix.clear();

    BoundedMatrix<double, TNumNodes, TDim> dn_dx;
    array_1d< double, TNumNodes > N;
    double volume;
    GeometryUtils::CalculateGeometryData(pElement->GetGeometry(),dn_dx,N,volume);

    unsigned int matrix_index_i, matrix_index_j;

    for (unsigned int a=0; a < TNumNodes; a++)
    {
        for (unsigned int r=0; r < TNumNodes; r++)
        {
            for (unsigned int i=0; i < TDim; i++)
            {
                for (unsigned int m=0; m < TDim; m++)
                {

                    matrix_index_i = a*TDim + i;
                    matrix_index_j = r*TDim + m;

                    double value = 0.0;
                    
                    value = dn_dx(a,i)*TauTwo*dn_dx(r,m);

                    resultant_matrix(matrix_index_i, matrix_index_j) = value*volume;
                }
            }
        }
    }
    
    WriteMatrix(7, resultant_matrix, pElement, rCurrentProcessInfo);
}

template<unsigned int TDim>
void Matrix_8(
    VMSAdjointElement<TDim>* pElement,
    array_1d< double, TDim >& Velocity,
    double Density,
    double DynamicViscosity, 
    ProcessInfo& rCurrentProcessInfo
)
{
    constexpr unsigned int TNumNodes = TDim + 1;
    constexpr unsigned int TCoordLocalSize = TDim * TNumNodes;

    Element::MatrixType resultant_matrix;

    resultant_matrix.resize(TCoordLocalSize, TCoordLocalSize);
    resultant_matrix.clear();

    BoundedMatrix<double, TNumNodes, TDim> dn_dx;
    array_1d< double, TNumNodes > N;
    double volume;
    GeometryUtils::CalculateGeometryData(pElement->GetGeometry(),dn_dx,N,volume);

    unsigned int matrix_index_i, matrix_index_j;

    for (unsigned int a=0; a < TNumNodes; a++)
    {
        for (unsigned int r=0; r < TNumNodes; r++)
        {
            for (unsigned int i=0; i < TDim; i++)
            {
                unsigned int m=i;

                matrix_index_i = a*TDim + i;
                matrix_index_j = r*TDim + m;

                double value = 0.0;

                for (unsigned int k=0; k < TDim; k++)
                    value += dn_dx(a,k)*dn_dx(r,k);
                
                resultant_matrix(matrix_index_i, matrix_index_j) = DynamicViscosity*value*volume;
            }
        }
    }
    
    WriteMatrix(8, resultant_matrix, pElement, rCurrentProcessInfo);
}

template<unsigned int TDim>
void Matrix_9(
    VMSAdjointElement<TDim>* pElement,
    array_1d< double, TDim >& Velocity,
    double Density,
    double DynamicViscosity, 
    ProcessInfo& rCurrentProcessInfo
)
{
    constexpr unsigned int TNumNodes = TDim + 1;
    constexpr unsigned int TCoordLocalSize = TDim * TNumNodes;

    Element::MatrixType resultant_matrix;

    resultant_matrix.resize(TCoordLocalSize, TCoordLocalSize);
    resultant_matrix.clear();

    BoundedMatrix<double, TNumNodes, TDim> dn_dx;
    array_1d< double, TNumNodes > N;
    double volume;
    GeometryUtils::CalculateGeometryData(pElement->GetGeometry(),dn_dx,N,volume);

    unsigned int matrix_index_i, matrix_index_j;

    for (unsigned int a=0; a < TNumNodes; a++)
    {
        for (unsigned int r=0; r < TNumNodes; r++)
        {
            for (unsigned int i=0; i < TDim; i++)
            {
                for (unsigned int m=0; m < TDim; m++)
                {
                    matrix_index_i = a*TDim + i;
                    matrix_index_j = r*TDim + m;

                    double value = dn_dx(a,m)*dn_dx(r,i);

                    resultant_matrix(matrix_index_i, matrix_index_j) = DynamicViscosity*value*volume;
                }
            }
        }
    }
    
    WriteMatrix(9, resultant_matrix, pElement, rCurrentProcessInfo);
}

template<unsigned int TDim>
void Matrix_10(
    VMSAdjointElement<TDim>* pElement,
    array_1d< double, TDim >& Velocity,
    double Density,
    double DynamicViscosity, 
    ProcessInfo& rCurrentProcessInfo
)
{
    constexpr unsigned int TNumNodes = TDim + 1;
    constexpr unsigned int TCoordLocalSize = TDim * TNumNodes;

    Element::MatrixType resultant_matrix;

    resultant_matrix.resize(TCoordLocalSize, TCoordLocalSize);
    resultant_matrix.clear();

    BoundedMatrix<double, TNumNodes, TDim> dn_dx;
    array_1d< double, TNumNodes > N;
    double volume;
    GeometryUtils::CalculateGeometryData(pElement->GetGeometry(),dn_dx,N,volume);

    unsigned int matrix_index_i, matrix_index_j;

    for (unsigned int a=0; a < TNumNodes; a++)
    {
        for (unsigned int r=0; r < TNumNodes; r++)
        {
            for (unsigned int i=0; i < TDim; i++)
            {
                for (unsigned int m=0; m < TDim; m++)
                {
                    matrix_index_i = a*TDim + i;
                    matrix_index_j = r*TDim + m;

                    double value = dn_dx(a,i)*dn_dx(r,m);

                    resultant_matrix(matrix_index_i, matrix_index_j) = -2.0*DynamicViscosity*value*volume/3.0;
                }
            }
        }
    }
    
    WriteMatrix(10, resultant_matrix, pElement, rCurrentProcessInfo);
}

template<unsigned int TDim>
void Matrix_11(
    VMSAdjointElement<TDim>* pElement,
    array_1d< double, TDim >& Velocity,
    double Density,
    double DynamicViscosity, 
    ProcessInfo& rCurrentProcessInfo,
    Element::MatrixType& TauTwoDerivatives
)
{
    constexpr unsigned int TNumNodes = TDim + 1;
    constexpr unsigned int TCoordLocalSize = TDim * TNumNodes;

    Element::MatrixType resultant_matrix;

    resultant_matrix.resize(TCoordLocalSize, TCoordLocalSize);
    resultant_matrix.clear();

    BoundedMatrix<double, TNumNodes, TDim> dn_dx;
    array_1d< double, TNumNodes > N;
    double volume;
    GeometryUtils::CalculateGeometryData(pElement->GetGeometry(),dn_dx,N,volume);

    unsigned int matrix_index_i, matrix_index_j;

    BoundedVector<array_1d<double, TDim>, TNumNodes> nodal_velocity_vectors;
    for (unsigned int iNode = 0; iNode < TNumNodes; iNode++)
         nodal_velocity_vectors[iNode] = pElement->GetGeometry()[iNode].FastGetSolutionStepValue(VELOCITY);

    double divU = 0.0;
    for (unsigned int k=0; k < TDim; k++)
        for (unsigned int l=0; l < TNumNodes; l++)
            divU += dn_dx(l,k)*nodal_velocity_vectors[l][k];

    for (unsigned int a=0; a < TNumNodes; a++)
    {
        for (unsigned int r=0; r < TNumNodes; r++)
        {
            for (unsigned int i=0; i < TDim; i++)
            {
                for (unsigned int m=0; m < TDim; m++)
                {
                    matrix_index_i = a*TDim + i;
                    matrix_index_j = r*TDim + m;

                    double value = 0.0;

                    value += dn_dx(a,i)*TauTwoDerivatives(r,m);
                    value += dn_dx(r,m)*TauTwoDerivatives(a,i);
                    
                    resultant_matrix(matrix_index_i, matrix_index_j) = 0.5*value*divU*volume;
                }
            }
        }
    }
    
    WriteMatrix(11, resultant_matrix, pElement, rCurrentProcessInfo);
}

template<unsigned int TDim>
void Matrix_12(
    VMSAdjointElement<TDim>* pElement,
    array_1d< double, TDim >& Velocity,
    double Density,
    double DynamicViscosity, 
    ProcessInfo& rCurrentProcessInfo,
    Element::MatrixType& rTauOneDerivatives
)
{
    constexpr unsigned int TNumNodes = TDim + 1;
    constexpr unsigned int TCoordLocalSize = TDim * TNumNodes;

    Element::MatrixType resultant_matrix;

    resultant_matrix.resize(TCoordLocalSize, TCoordLocalSize);
    resultant_matrix.clear();

    BoundedMatrix<double, TNumNodes, TDim> dn_dx;
    array_1d< double, TNumNodes > N;
    double volume;
    GeometryUtils::CalculateGeometryData(pElement->GetGeometry(),dn_dx,N,volume);

    unsigned int matrix_index_i, matrix_index_j;

    BoundedVector<double, TNumNodes> nodal_pressure_vector;
    for (unsigned int iNode = 0; iNode < TNumNodes; iNode++)
         nodal_pressure_vector[iNode] = pElement->GetGeometry()[iNode].FastGetSolutionStepValue(PRESSURE);

    for (unsigned int a=0; a < TNumNodes; a++)
    {
        for (unsigned int r=0; r < TNumNodes; r++)
        {
            for (unsigned int i=0; i < TDim; i++)
            {
                for (unsigned int m=0; m < TDim; m++)
                {
                    matrix_index_i = a*TDim + i;
                    matrix_index_j = r*TDim + m;

                    double u_dot_delta_Na = 0.0;
                    for (unsigned int k = 0; k < TDim; k++)
                        u_dot_delta_Na += Velocity[k]*dn_dx(a,k);
                    
                    double u_dot_delta_Nr = 0.0;
                    for (unsigned int k = 0; k < TDim; k++)
                        u_dot_delta_Nr += Velocity[k]*dn_dx(r,k);
                    
                    double daba_i_P = 0.0;
                    for (unsigned int k=0; k < TNumNodes; k++)
                        daba_i_P += dn_dx(k,i)*nodal_pressure_vector[k];
                    double daba_m_P = 0.0;
                    for (unsigned int k=0; k < TNumNodes; k++)
                        daba_m_P += dn_dx(k,m)*nodal_pressure_vector[k];                    

                    double value = Density*u_dot_delta_Na*rTauOneDerivatives(r,m)*daba_i_P;
                    value += Density*u_dot_delta_Nr*rTauOneDerivatives(a,i)*daba_m_P;
                    
                    resultant_matrix(matrix_index_i, matrix_index_j) = 0.5*value*volume;
                }
            }
        }
    }
    
    WriteMatrix(12, resultant_matrix, pElement, rCurrentProcessInfo);
}

template<unsigned int TDim>
void Matrix_13(
    VMSAdjointElement<TDim>* pElement,
    array_1d< double, TDim >& Velocity,
    double Density,
    double DynamicViscosity, 
    ProcessInfo& rCurrentProcessInfo,
    double TauOne
)
{
    constexpr unsigned int TNumNodes = TDim + 1;
    constexpr unsigned int TCoordLocalSize = TDim * TNumNodes;

    Element::MatrixType resultant_matrix;

    resultant_matrix.resize(TCoordLocalSize, TCoordLocalSize);
    resultant_matrix.clear();

    BoundedMatrix<double, TNumNodes, TDim> dn_dx;
    array_1d< double, TNumNodes > N;
    double volume;
    GeometryUtils::CalculateGeometryData(pElement->GetGeometry(),dn_dx,N,volume);

    unsigned int matrix_index_i, matrix_index_j;

    BoundedVector<double, TNumNodes> nodal_pressure_vector;
    for (unsigned int iNode = 0; iNode < TNumNodes; iNode++)
         nodal_pressure_vector[iNode] = pElement->GetGeometry()[iNode].FastGetSolutionStepValue(PRESSURE);

    for (unsigned int a=0; a < TNumNodes; a++)
    {
        for (unsigned int r=0; r < TNumNodes; r++)
        {
            for (unsigned int i=0; i < TDim; i++)
            {
                for (unsigned int m=0; m < TDim; m++)
                {
                    matrix_index_i = a*TDim + i;
                    matrix_index_j = r*TDim + m;

                    double daba_i_P = 0.0;
                    for (unsigned int k=0; k < TNumNodes; k++)
                        daba_i_P += dn_dx(k,i)*nodal_pressure_vector[k];
                    double daba_m_P = 0.0;
                    for (unsigned int k=0; k < TNumNodes; k++)
                        daba_m_P += dn_dx(k,m)*nodal_pressure_vector[k];                    

                    double value = 0.0;
                    value += Density*N[r]*dn_dx(a,m)*TauOne*daba_i_P;
                    value += Density*N[a]*dn_dx(r,i)*TauOne*daba_m_P;
                    
                    resultant_matrix(matrix_index_i, matrix_index_j) = 0.5*value*volume;
                }
            }
        }
    }
    
    WriteMatrix(13, resultant_matrix, pElement, rCurrentProcessInfo);
}

template<unsigned int TDim>
void Matrix_14(
    VMSAdjointElement<TDim>* pElement,
    array_1d< double, TDim >& Velocity,
    double Density,
    double DynamicViscosity, 
    ProcessInfo& rCurrentProcessInfo,
    double TauOne
)
{
    constexpr unsigned int TNumNodes = TDim + 1;

    Element::MatrixType resultant_matrix;

    resultant_matrix.resize(TNumNodes, TNumNodes);
    resultant_matrix.clear();

    BoundedMatrix<double, TNumNodes, TDim> dn_dx;
    array_1d< double, TNumNodes > N;
    double volume;
    GeometryUtils::CalculateGeometryData(pElement->GetGeometry(),dn_dx,N,volume);

    unsigned int matrix_index_i, matrix_index_j;

    for (unsigned int a=0; a < TNumNodes; a++)
    {
        for (unsigned int r=0; r < TNumNodes; r++)
        {

            matrix_index_i = a;
            matrix_index_j = r;

            double value = 0.0;

            for (unsigned int k=0; k < TDim; k++)
                value += dn_dx(a,k)*dn_dx(r,k);
            
            resultant_matrix(matrix_index_i, matrix_index_j) = value*TauOne*volume;
        }
    }
    
    WriteMatrix(14, resultant_matrix, pElement, rCurrentProcessInfo);
}

template<unsigned int TDim>
void NumericalDiffusionMatrix(
    VMSAdjointElement<TDim>* pElement,
    array_1d< double, TDim >& Velocity,
    ProcessInfo& rCurrentProcessInfo
)
{
    constexpr unsigned int TNumNodes = TDim + 1;
    constexpr unsigned int TBlockSize = TDim + 1;
    constexpr unsigned int TCoordLocalSize = TBlockSize * TNumNodes;

    Element::MatrixType resultant_matrix;

    resultant_matrix.resize(TCoordLocalSize, TCoordLocalSize);
    resultant_matrix.clear();

    BoundedMatrix<double, TNumNodes, TDim> dn_dx;
    array_1d< double, TNumNodes > N;
    double volume;
    GeometryUtils::CalculateGeometryData(pElement->GetGeometry(),dn_dx,N,volume);

    for (unsigned int a = 0; a < TNumNodes; ++a)
    {
        for (unsigned int b = 0; b < TNumNodes; ++b)
        {
            // (dN_a/dx_k dN_b/dx_k)
            double value = 0.0;
            for (unsigned int k = 0; k < TDim; k++)
                value += dn_dx(a,k) * dn_dx(b,k);
            
            value *= volume;

            for (unsigned int i=0; i < TBlockSize; i++)
                resultant_matrix(a*TBlockSize+i,b*TBlockSize+i) = value;
        }
    }

    WriteMatrix(18, resultant_matrix, pElement, rCurrentProcessInfo);
}

template<unsigned int TDim>
void VMSAdjointElement<TDim>::ProcessSymmetricMatrices(ProcessInfo& rCurrentProcessInfo)
{
    BoundedMatrix<double, TNumNodes, TDim> dn_dx;
    array_1d< double, TNumNodes > N;
    double volume;
    GeometryUtils::CalculateGeometryData(this->GetGeometry(),dn_dx,N,volume);
    
    // Density
    double Density;
    this->EvaluateInPoint(Density,DENSITY,N);

    // Dynamic viscosity
    double Viscosity;
    this->EvaluateInPoint(Viscosity,VISCOSITY,N);
    Viscosity *= Density;


    // u
    array_1d< double, TDim > Velocity;
    this->EvaluateInPoint(Velocity,VELOCITY,N);

    double TauOne, TauTwo;
    double VelNorm = norm_2(Velocity);
    double ElemSize = this->CalculateElementSize(volume);
    this->CalculateStabilizationParameters(TauOne,TauTwo,VelNorm,ElemSize,
            Density,Viscosity,rCurrentProcessInfo);

    array_1d< double, TCoordLocalSize > DetJDerivatives;
    this->CalculateDeterminantOfJacobianDerivatives(DetJDerivatives);

    Element::MatrixType TauOneDerivatives;
    TauOneDerivatives.resize(TNumNodes, TDim);
    TauOneDerivatives.clear();
    
    Element::MatrixType TauTwoDerivatives;
    TauTwoDerivatives.resize(TNumNodes, TDim);
    TauTwoDerivatives.clear();

    for (unsigned int i=0; i < TNumNodes; i++)
        for (unsigned int j=0; j < TDim; j++)
        {
            double DetJDeriv = DetJDerivatives[i*TDim+j];
            
            double TauOneDeriv, TauTwoDeriv;
            this->CalculateStabilizationParametersDerivative(
                    TauOneDeriv,TauTwoDeriv,TauOne,TauTwo,VelNorm,ElemSize,Density,
                    Viscosity,DetJDeriv);

            TauOneDerivatives(i,j) = TauOneDeriv;
            TauTwoDerivatives(i,j) = TauTwoDeriv;
        }

    Matrix_1<TDim>(this, Velocity, Density, Viscosity, rCurrentProcessInfo);
    Matrix_2<TDim>(this, Velocity, Density, Viscosity, rCurrentProcessInfo);
    Matrix_3<TDim>(this, Velocity, Density, Viscosity, rCurrentProcessInfo, TauOne);
    Matrix_4<TDim>(this, Velocity, Density, Viscosity, rCurrentProcessInfo, TauOneDerivatives);
    Matrix_5<TDim>(this, Velocity, Density, Viscosity, rCurrentProcessInfo, TauOne);
    Matrix_6<TDim>(this, Velocity, Density, Viscosity, rCurrentProcessInfo, TauOne);
    Matrix_7<TDim>(this, Velocity, Density, Viscosity, rCurrentProcessInfo, TauTwo);
    Matrix_8<TDim>(this, Velocity, Density, Viscosity, rCurrentProcessInfo);
    Matrix_9<TDim>(this, Velocity, Density, Viscosity, rCurrentProcessInfo);
    Matrix_10<TDim>(this, Velocity, Density, Viscosity, rCurrentProcessInfo);
    Matrix_11<TDim>(this, Velocity, Density, Viscosity, rCurrentProcessInfo, TauTwoDerivatives);
    Matrix_12<TDim>(this, Velocity, Density, Viscosity, rCurrentProcessInfo, TauOneDerivatives);
    Matrix_13<TDim>(this, Velocity, Density, Viscosity, rCurrentProcessInfo, TauOne);
    Matrix_14<TDim>(this, Velocity, Density, Viscosity, rCurrentProcessInfo, TauOne);
    NumericalDiffusionMatrix<TDim>(this, Velocity, rCurrentProcessInfo);

    Element::MatrixType resultant_matrix;

    this->CalculatePrimalGradientOfVMSSteadyTerm(resultant_matrix, rCurrentProcessInfo);

    resultant_matrix = 0.5*(resultant_matrix + boost::numeric::ublas::trans(resultant_matrix));
    
    WriteMatrix(17, resultant_matrix, this, rCurrentProcessInfo);
    
}

