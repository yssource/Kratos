//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Elisa Magliozzi
//

#include "custom_elements/compressible_navier_stokes.h"

namespace Kratos {

template<>
void CompressibleNavierStokes<3>::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    unsigned int Dim = 3;
    unsigned int Dimes = Dim+2;
    unsigned int NumNodes = 4;
    unsigned int DofSize  = NumNodes*(Dimes);

    if (rResult.size() != DofSize)
        rResult.resize(DofSize, false);

    for(unsigned int i=0; i<NumNodes; i++)
    {
        rResult[i*(Dimes)  ]  =  this->GetGeometry()[i].GetDof(DENSITY).EquationId();
        rResult[i*(Dimes)+1]  =  this->GetGeometry()[i].GetDof(MOMENT_X).EquationId();
        rResult[i*(Dimes)+2]  =  this->GetGeometry()[i].GetDof(MOMENT_Y).EquationId();
        rResult[i*(Dimes)+3]  =  this->GetGeometry()[i].GetDof(MOMENT_Z).EquationId(); 
        rResult[i*(Dimes)+4]  =  this->GetGeometry()[i].GetDof(TOTAL_ENERGY).EquationId();
    }

    KRATOS_CATCH("")
}


template<>
void CompressibleNavierStokes<2>::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    unsigned int Dim = 2;
    unsigned int Dimes = Dim+2;
    unsigned int NumNodes = 3;
    unsigned int DofSize  = NumNodes*(Dimes);

    if (rResult.size() != DofSize)
        rResult.resize(DofSize, false);

    for(unsigned int i=0; i<NumNodes; i++)
    {
       rResult[i*(Dimes)  ]  =  this->GetGeometry()[i].GetDof(DENSITY).EquationId();
        rResult[i*(Dimes)+1]  =  this->GetGeometry()[i].GetDof(MOMENT_X).EquationId();
        rResult[i*(Dimes)+2]  =  this->GetGeometry()[i].GetDof(MOMENT_Y).EquationId();
        rResult[i*(Dimes)+3]  =  this->GetGeometry()[i].GetDof(TOTAL_ENERGY).EquationId();
    }

    KRATOS_CATCH("")
}


template<>
void CompressibleNavierStokes<3>::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    unsigned int Dim = 3;
    unsigned int Dimes = Dim+2;
    unsigned int NumNodes = 4;
    unsigned int DofSize  = NumNodes*(Dimes);

    if (ElementalDofList.size() != DofSize)
        ElementalDofList.resize(DofSize);

    for(unsigned int i=0; i<NumNodes; i++)
    {
        ElementalDofList[i*(Dimes)  ]  =  this->GetGeometry()[i].pGetDof(DENSITY);
        ElementalDofList[i*(Dimes)+1]  =  this->GetGeometry()[i].pGetDof(MOMENT_X);
        ElementalDofList[i*(Dimes)+2]  =  this->GetGeometry()[i].pGetDof(MOMENT_Y);
        ElementalDofList[i*(Dimes)+3]  =  this->GetGeometry()[i].pGetDof(MOMENT_Z);
        ElementalDofList[i*(Dimes)+4]  =  this->GetGeometry()[i].pGetDof(TOTAL_ENERGY);
    }

    KRATOS_CATCH("");
}


template<>
void CompressibleNavierStokes<2>::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    unsigned int Dim = 2;
    unsigned int Dimes = Dim+2;
    unsigned int NumNodes = 3;
    unsigned int DofSize  = NumNodes*(Dimes);

    if (ElementalDofList.size() != DofSize)
        ElementalDofList.resize(DofSize);

    for(unsigned int i=0; i<NumNodes; i++)
    {
        ElementalDofList[i*(Dimes)  ]  =  this->GetGeometry()[i].pGetDof(DENSITY);
        ElementalDofList[i*(Dimes)+1]  =  this->GetGeometry()[i].pGetDof(MOMENT_X);
        ElementalDofList[i*(Dimes)+2]  =  this->GetGeometry()[i].pGetDof(MOMENT_Y);
        ElementalDofList[i*(Dimes)+3]  =  this->GetGeometry()[i].pGetDof(TOTAL_ENERGY);
    }

    KRATOS_CATCH("");
}


template<>
void CompressibleNavierStokes<3>::ComputeGaussPointLHSContribution(bounded_matrix<double,20,20>& lhs, const ElementDataStruct& data)
{
    const int nnodes = 4;
    const int dim = 3;
    const int dimes = dim+2;
    
    const double h = data.h;                                // Characteristic element size
    
    const double& bdf0 = data.bdf0;
    const double& bdf1 = data.bdf1;
    const double& bdf2 = data.bdf2;
    
    const bounded_matrix<double,nnodes,dimes>& U = data.U;
    const bounded_matrix<double,nnodes,dimes>& Un = data.Un;
    const bounded_matrix<double,nnodes,dimes>& Unn = data.Unn;
    const bounded_matrix<double,nnodes,dim>& f_ext = data.f_ext;
//     const bounded_matrix<double,nnodes,dim>& vmesh = data.vmesh;
    const array_1d<double,nnodes>& r = data.r;
    
    const double mu = data.mu;
    const double nu = data.nu;
    const double lambda = data.lambda;
    const double cv = data.cv;
    const double y = data.y;
    
    // Get shape function values
    const array_1d<double,nnodes>& N = data.N;
    const bounded_matrix<double,nnodes,dim>& DN = data.DN_DX;
    

    
    // Stabilization parameters
    const double c1 = 4.0;
    const double c2 = 2.0;
   
    //substitute_lhs_3D

}


template<>
void CompressibleNavierStokes<2>::ComputeGaussPointLHSContribution(bounded_matrix<double,12,12>& lhs, const ElementDataStruct& data)
{
    const int nnodes = 3;
    const int dim = 2;
    const int dimes = dim+2;
  
    const double h = data.h;                                // Characteristic element size
    
    const double& bdf0 = data.bdf0;
    const double& bdf1 = data.bdf1;
    const double& bdf2 = data.bdf2;
    
    const bounded_matrix<double,nnodes,dimes>& U = data.U;
    const bounded_matrix<double,nnodes,dimes>& Un = data.Un;
    const bounded_matrix<double,nnodes,dimes>& Unn = data.Unn;
    const bounded_matrix<double,nnodes,dim>& f_ext = data.f_ext;
    //     const bounded_matrix<double,nnodes,dim>& vmesh = data.vmesh;
    const array_1d<double,nnodes>& r = data.r;

    const double mu = data.mu;
    const double nu = data.nu;
    const double lambda = data.lambda;
    const double cv = data.cv;
    const double y = data.y;
    
    // Get shape function values
    const array_1d<double,nnodes>& N = data.N;
    const bounded_matrix<double,nnodes,dim>& DN = data.DN_DX;

    // Stabilization parameters
    const double c1 = 4.0;
    const double c2 = 2.0;

    //substitute_lhs_2D

}


template<>
void CompressibleNavierStokes<3>::ComputeGaussPointRHSContribution(array_1d<double,20>& rhs, const ElementDataStruct& data)
{
    const int nnodes = 4;
    const int dim = 3;
    const int dimes = dim+2;

    const double h = data.h;                                // Characteristic element size

    const double& bdf0 = data.bdf0;
    const double& bdf1 = data.bdf1;
    const double& bdf2 = data.bdf2;
    
    const bounded_matrix<double,nnodes,dimes>& U = data.U;
    const bounded_matrix<double,nnodes,dimes>& Un = data.Un;
    const bounded_matrix<double,nnodes,dimes>& Unn = data.Unn;
    const bounded_matrix<double,nnodes,dim>& f_ext = data.f_ext;
    const array_1d<double,nnodes>& r = data.r;

    const double mu = data.mu;
    const double nu = data.nu;
    const double lambda = data.lambda;
    const double cv = data.cv;
    const double y = data.y;
    
    // Get shape function values
    const array_1d<double,nnodes>& N = data.N;
    const bounded_matrix<double,nnodes,dim>& DN = data.DN_DX;

    // Auxiliary variables used in the calculation of the RHS
    const array_1d<double,dimes> U_gauss = prod(trans(U), N);
    const array_1d<double,dim> f_gauss = prod(trans(f_ext), N);
    const array_1d<double,dim> grad_U = prod(trans(DN), U);
    const double& r_gauss = inner_prod(data.N, data.r);
    
    const array_1d<double,dimes> accel_gauss = bdf0*U_gauss+bdf1*prod(trans(Un), N)+bdf2*prod(trans(Unn), N);
   
    // Stabilization parameters
    const double stab_c1 = 4.0;
    const double stab_c2 = 2.0;
    
    //substitute_rhs_3D
}


template<>
void CompressibleNavierStokes<2>::ComputeGaussPointRHSContribution(array_1d<double,12>& rhs, const ElementDataStruct& data)
{
    const int nnodes = 3;
    const int dim = 2;
    const int dimes = dim+2;
    
    const double h = data.h;                                // Characteristic element size
    
    const double& bdf0 = data.bdf0;
    const double& bdf1 = data.bdf1;
    const double& bdf2 = data.bdf2;

    const bounded_matrix<double,nnodes,dimes>& U = data.U;
    const bounded_matrix<double,nnodes,dimes>& Un = data.Un;
    const bounded_matrix<double,nnodes,dimes>& Unn = data.Unn;
    const bounded_matrix<double,nnodes,dim>& f_ext = data.f_ext;
    const array_1d<double,nnodes>& r = data.r;
        
    const double mu = data.mu;
    const double nu = data.nu;
    const double lambda = data.lambda;
    const double cv = data.cv;
    const double y = data.y;

    // Get shape function values
    const array_1d<double,nnodes>& N = data.N;
    const bounded_matrix<double,nnodes,dim>& DN = data.DN_DX;

    // Auxiliary variables used in the calculation of the RHS
    const array_1d<double,dimes> U_gauss = prod(trans(U), N);
    const array_1d<double,dim> f_gauss = prod(trans(f_ext), N);
    const array_1d<double,dim> grad_U = prod(trans(DN), U);
    const double& r_gauss = inner_prod(data.N, data.r);
    
    const array_1d<double,dimes> accel_gauss = bdf0*U_gauss+bdf1*prod(trans(Un), N)+bdf2*prod(trans(Unn), N);
   
    // Stabilization parameters
    const double c1 = 4.0;
    const double c2 = 2.0;

    //substitute_rhs_2D
}


template<>
double CompressibleNavierStokes<3>::SubscaleErrorEstimate(const ElementDataStruct& data)
{
    const int nnodes = 4;
    const int dim = 3;
    const int dimes = dim+2;

    const double h = data.h;                                // Characteristic element size

    const double& bdf0 = data.bdf0;
    const double& bdf1 = data.bdf1;
    const double& bdf2 = data.bdf2;

    const bounded_matrix<double,nnodes,dimes>& U = data.U;
    const bounded_matrix<double,nnodes,dimes>& Un = data.Un;
    const bounded_matrix<double,nnodes,dimes>& Unn = data.Unn;
    const bounded_matrix<double,nnodes,dim>& f_ext = data.f_ext;
    const array_1d<double,nnodes>& r = data.r;
          
    const double mu = data.mu;
    const double nu = data.nu;
    const double lambda = data.lambda;
    const double cv = data.cv;
    const double y = data.y;

    // Get shape function values
    const array_1d<double,nnodes>& N = data.N;
    const bounded_matrix<double,nnodes,dim>& DN = data.DN_DX;

    // Auxiliary variables used in the calculation of the error estimator
    array_1d<double,dimes> U_s_gauss; //WHAT IS THIS FOR?
    const array_1d<double,dimes> U_gauss = prod(trans(U), N);
    const array_1d<double,dim> f_gauss = prod(trans(f_ext), N);
    const array_1d<double,dimes> grad_U = prod(trans(DN), U);
    const double& r_gauss = inner_prod(data.N, data.r);
    
    const array_1d<double,dimes> accel_gauss = bdf0*U_gauss+bdf1*prod(trans(Un), N)+bdf2*prod(trans(Unn), N);
    
    // Stabilization parameters
    const double c1 = 4.0;
    const double c2 = 2.0;
   
    // Gauss point velocity subscale value computation
    //substitute_gausspt_subscale_3D

    const double U_gauss_norm = norm_2(U_gauss);
    const double U_s_gauss_norm = norm_2(U_s_gauss);

    return U_s_gauss_norm/U_gauss_norm;
}


template<>
double CompressibleNavierStokes<2>::SubscaleErrorEstimate(const ElementDataStruct& data)
{
    const int nnodes = 3;
    const int dim = 2;
    const int dimes = dim+2;

   const double h = data.h;                                // Characteristic element size

    const double& bdf0 = data.bdf0;
    const double& bdf1 = data.bdf1;
    const double& bdf2 = data.bdf2;

    const bounded_matrix<double,nnodes,dimes>& U = data.U;
    const bounded_matrix<double,nnodes,dimes>& Un = data.Un;
    const bounded_matrix<double,nnodes,dimes>& Unn = data.Unn;
    const bounded_matrix<double,nnodes,dim>& f_ext = data.f_ext;
    const array_1d<double,nnodes>& r = data.r;
          
    const double mu = data.mu;
    const double nu = data.nu;
    const double lambda = data.lambda;
    const double cv = data.cv;
    const double y = data.y;

    // Get shape function values
    const array_1d<double,nnodes>& N = data.N;
    const bounded_matrix<double,nnodes,dim>& DN = data.DN_DX;

    // Auxiliary variables used in the calculation of the error estimator
    array_1d<double,dimes> U_s_gauss; //WHAT IS THIS FOR?
    const array_1d<double,dimes> U_gauss = prod(trans(U), N);
    const array_1d<double,dim> f_gauss = prod(trans(f_ext), N);
    const array_1d<double,dimes> grad_U = prod(trans(DN), U);
    const double& r_gauss = inner_prod(data.N, data.r);
    
    const array_1d<double,dimes> accel_gauss = bdf0*U_gauss+bdf1*prod(trans(Un), N)+bdf2*prod(trans(Unn), N);

    // Gauss point velocity subscale value computation
    //substitute_gausspt_subscale_2D

    const double U_gauss_norm = norm_2(U_gauss);
    const double U_s_gauss_norm = norm_2(U_s_gauss);

    return U_s_gauss_norm/U_gauss_norm;
}

}
