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
    unsigned int BlockSize = Dim+2;
    unsigned int NumNodes = 4;
    unsigned int DofSize  = NumNodes*(BlockSize);

    if (rResult.size() != DofSize)
        rResult.resize(DofSize, false);

    for(unsigned int i=0; i<NumNodes; i++)
    {
        rResult[i*(BlockSize)  ]  =  this->GetGeometry()[i].GetDof(DENSITY).EquationId();
        rResult[i*(BlockSize)+1]  =  this->GetGeometry()[i].GetDof(MOMENT_X).EquationId();
        rResult[i*(BlockSize)+2]  =  this->GetGeometry()[i].GetDof(MOMENT_Y).EquationId();
        rResult[i*(BlockSize)+3]  =  this->GetGeometry()[i].GetDof(MOMENT_Z).EquationId(); 
        rResult[i*(BlockSize)+4]  =  this->GetGeometry()[i].GetDof(TOTAL_ENERGY).EquationId();
    }

    KRATOS_CATCH("")
}


template<>
void CompressibleNavierStokes<2>::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    unsigned int Dim = 2;
    unsigned int BlockSize = Dim+2;
    unsigned int NumNodes = 3;
    unsigned int DofSize  = NumNodes*(BlockSize);

    if (rResult.size() != DofSize)
        rResult.resize(DofSize, false);

    for(unsigned int i=0; i<NumNodes; i++)
    {
       rResult[i*(BlockSize)  ]  =  this->GetGeometry()[i].GetDof(DENSITY).EquationId();
        rResult[i*(BlockSize)+1]  =  this->GetGeometry()[i].GetDof(MOMENT_X).EquationId();
        rResult[i*(BlockSize)+2]  =  this->GetGeometry()[i].GetDof(MOMENT_Y).EquationId();
        rResult[i*(BlockSize)+3]  =  this->GetGeometry()[i].GetDof(TOTAL_ENERGY).EquationId();
    }

    KRATOS_CATCH("")
}


template<>
void CompressibleNavierStokes<3>::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    unsigned int Dim = 3;
    unsigned int BlockSize = Dim+2;
    unsigned int NumNodes = 4;
    unsigned int DofSize  = NumNodes*(BlockSize);

    if (ElementalDofList.size() != DofSize)
        ElementalDofList.resize(DofSize);

    for(unsigned int i=0; i<NumNodes; i++)
    {
        ElementalDofList[i*(BlockSize)  ]  =  this->GetGeometry()[i].pGetDof(DENSITY);
        ElementalDofList[i*(BlockSize)+1]  =  this->GetGeometry()[i].pGetDof(MOMENT_X);
        ElementalDofList[i*(BlockSize)+2]  =  this->GetGeometry()[i].pGetDof(MOMENT_Y);
        ElementalDofList[i*(BlockSize)+3]  =  this->GetGeometry()[i].pGetDof(MOMENT_Z);
        ElementalDofList[i*(BlockSize)+4]  =  this->GetGeometry()[i].pGetDof(TOTAL_ENERGY);
    }

    KRATOS_CATCH("");
}


template<>
void CompressibleNavierStokes<2>::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    unsigned int Dim = 2;
    unsigned int BlockSize = Dim+2;
    unsigned int NumNodes = 3;
    unsigned int DofSize  = NumNodes*(BlockSize);

    if (ElementalDofList.size() != DofSize)
        ElementalDofList.resize(DofSize);

    for(unsigned int i=0; i<NumNodes; i++)
    {
        ElementalDofList[i*(BlockSize)  ]  =  this->GetGeometry()[i].pGetDof(DENSITY);
        ElementalDofList[i*(BlockSize)+1]  =  this->GetGeometry()[i].pGetDof(MOMENT_X);
        ElementalDofList[i*(BlockSize)+2]  =  this->GetGeometry()[i].pGetDof(MOMENT_Y);
        ElementalDofList[i*(BlockSize)+3]  =  this->GetGeometry()[i].pGetDof(TOTAL_ENERGY);
    }

    KRATOS_CATCH("");
}


template<>
void CompressibleNavierStokes<3>::ComputeGaussPointLHSContribution(bounded_matrix<double,20,20>& lhs, const ElementDataStruct& data)
{
    const int nnodes = 4;
    const int dim = 3;
    const int BlockSize = dim+2;
    const double h = data.h; 
 
    const double& bdf0 = data.bdf0;
    const double& bdf1 = data.bdf1;
    const double& bdf2 = data.bdf2;
    
    const bounded_matrix<double,nnodes,BlockSize>& U = data.U;
    const bounded_matrix<double,nnodes,BlockSize>& Un = data.Un;
    const bounded_matrix<double,nnodes,BlockSize>& Unn = data.Unn;
    const bounded_matrix<double,nnodes,dim>& f_ext = data.f_ext;
    const array_1d<double,nnodes>& r = data.r;   
    const double mu = data.mu;
    const double nu = data.nu;
    const double lambda = data.lambda;
    const double cv = data.cv;
    const double y = data.y;
    const double cp = cv*y;
 
    // Get shape function values
    const array_1d<double,nnodes>& N = data.N;
    const bounded_matrix<double,nnodes,dim>& DN = data.DN_DX;
    
    // Stabilization parameters
    const double stab_c1 = 4.0;
    const double stab_c2 = 2.0;
    
    const array_1d<double,BlockSize> U_gauss= prod(trans(U),N);
    
    double tmp = U_gauss(dim+1)/U_gauss(0);
    for(unsigned int ll=0; ll<dim; ll++)
        tmp -=(U_gauss(ll+1)*U_gauss(ll+1))/(2*U_gauss(0)*U_gauss(0));
    double c = sqrt(y*(y-1)*tmp);

    double tau1inv = 0.0;
    for(unsigned int ll=0; ll<dim; ll++)
        tau1inv += (U_gauss(ll+1)/U_gauss(0))*(U_gauss(ll+1)/U_gauss(0));
    tau1inv = (sqrt(tau1inv)+c)*stab_c2/h;
    double tau2inv = stab_c1*nu/(h*h)+tau1inv;
    double tau3inv = stab_c1*lambda/(U_gauss(0)*cp*h*h)+tau1inv;
        
    const double tau1 = 1/tau1inv;
    const double tau2 = 1/tau2inv;
    const double tau3 = 1/tau3inv;
    //substitute_lhs_3D

}


template<>
void CompressibleNavierStokes<2>::ComputeGaussPointLHSContribution(bounded_matrix<double,12,12>& lhs, const ElementDataStruct& data)
{
    const int nnodes = 3;
    const int dim = 2;
    const int BlockSize = dim+2;
    const double h = data.h; 

    const double& bdf0 = data.bdf0;
    const double& bdf1 = data.bdf1;
    const double& bdf2 = data.bdf2;
    
    const bounded_matrix<double,nnodes,BlockSize>& U = data.U;
    const bounded_matrix<double,nnodes,BlockSize>& Un = data.Un;
    const bounded_matrix<double,nnodes,BlockSize>& Unn = data.Unn;
    const bounded_matrix<double,nnodes,dim>& f_ext = data.f_ext;
    const array_1d<double,nnodes>& r = data.r;
    const double mu = data.mu;
    const double nu = data.nu;
    const double lambda = data.lambda;
    const double cv = data.cv;
    const double y = data.y;
    const double cp = cv*y;
 
    // Get shape function values
    const array_1d<double,nnodes>& N = data.N;
    const bounded_matrix<double,nnodes,dim>& DN = data.DN_DX;
    
    // Stabilization parameters
    const double stab_c1 = 4.0;
    const double stab_c2 = 2.0;
    
    const array_1d<double,BlockSize> U_gauss= prod(trans(U),N);
    
    double tmp = U_gauss(dim+1)/U_gauss(0);
    for(unsigned int ll=0; ll<dim; ll++)
        tmp -=(U_gauss(ll+1)*U_gauss(ll+1))/(2*U_gauss(0)*U_gauss(0));
    double c = sqrt(y*(y-1)*tmp);

    double tau1inv = 0.0;
    for(unsigned int ll=0; ll<dim; ll++)
        tau1inv += (U_gauss(ll+1)/U_gauss(0))*(U_gauss(ll+1)/U_gauss(0));
    tau1inv = (sqrt(tau1inv)+c)*stab_c2/h;
    double tau2inv = stab_c1*nu/(h*h)+tau1inv;
    double tau3inv = stab_c1*lambda/(U_gauss(0)*cp*h*h)+tau1inv;
        
    const double tau1 = 1/tau1inv;
    const double tau2 = 1/tau2inv;
    const double tau3 = 1/tau3inv;

    const double clhs0 =             -pow(N[0], 2)*bdf0;
const double clhs1 =             N[0]*bdf0;
const double clhs2 =             -N[1]*clhs1;
const double clhs3 =             -N[2]*clhs1;
const double clhs4 =             -pow(N[1], 2)*bdf0;
const double clhs5 =             -N[1]*N[2]*bdf0;
const double clhs6 =             -pow(N[2], 2)*bdf0;
            lhs(0,0)=clhs0;
            lhs(0,1)=0;
            lhs(0,2)=0;
            lhs(0,3)=0;
            lhs(0,4)=clhs2;
            lhs(0,5)=0;
            lhs(0,6)=0;
            lhs(0,7)=0;
            lhs(0,8)=clhs3;
            lhs(0,9)=0;
            lhs(0,10)=0;
            lhs(0,11)=0;
            lhs(1,0)=0;
            lhs(1,1)=clhs0;
            lhs(1,2)=0;
            lhs(1,3)=0;
            lhs(1,4)=0;
            lhs(1,5)=clhs2;
            lhs(1,6)=0;
            lhs(1,7)=0;
            lhs(1,8)=0;
            lhs(1,9)=clhs3;
            lhs(1,10)=0;
            lhs(1,11)=0;
            lhs(2,0)=0;
            lhs(2,1)=0;
            lhs(2,2)=clhs0;
            lhs(2,3)=0;
            lhs(2,4)=0;
            lhs(2,5)=0;
            lhs(2,6)=clhs2;
            lhs(2,7)=0;
            lhs(2,8)=0;
            lhs(2,9)=0;
            lhs(2,10)=clhs3;
            lhs(2,11)=0;
            lhs(3,0)=0;
            lhs(3,1)=0;
            lhs(3,2)=0;
            lhs(3,3)=clhs0;
            lhs(3,4)=0;
            lhs(3,5)=0;
            lhs(3,6)=0;
            lhs(3,7)=clhs2;
            lhs(3,8)=0;
            lhs(3,9)=0;
            lhs(3,10)=0;
            lhs(3,11)=clhs3;
            lhs(4,0)=clhs2;
            lhs(4,1)=0;
            lhs(4,2)=0;
            lhs(4,3)=0;
            lhs(4,4)=clhs4;
            lhs(4,5)=0;
            lhs(4,6)=0;
            lhs(4,7)=0;
            lhs(4,8)=clhs5;
            lhs(4,9)=0;
            lhs(4,10)=0;
            lhs(4,11)=0;
            lhs(5,0)=0;
            lhs(5,1)=clhs2;
            lhs(5,2)=0;
            lhs(5,3)=0;
            lhs(5,4)=0;
            lhs(5,5)=clhs4;
            lhs(5,6)=0;
            lhs(5,7)=0;
            lhs(5,8)=0;
            lhs(5,9)=clhs5;
            lhs(5,10)=0;
            lhs(5,11)=0;
            lhs(6,0)=0;
            lhs(6,1)=0;
            lhs(6,2)=clhs2;
            lhs(6,3)=0;
            lhs(6,4)=0;
            lhs(6,5)=0;
            lhs(6,6)=clhs4;
            lhs(6,7)=0;
            lhs(6,8)=0;
            lhs(6,9)=0;
            lhs(6,10)=clhs5;
            lhs(6,11)=0;
            lhs(7,0)=0;
            lhs(7,1)=0;
            lhs(7,2)=0;
            lhs(7,3)=clhs2;
            lhs(7,4)=0;
            lhs(7,5)=0;
            lhs(7,6)=0;
            lhs(7,7)=clhs4;
            lhs(7,8)=0;
            lhs(7,9)=0;
            lhs(7,10)=0;
            lhs(7,11)=clhs5;
            lhs(8,0)=clhs3;
            lhs(8,1)=0;
            lhs(8,2)=0;
            lhs(8,3)=0;
            lhs(8,4)=clhs5;
            lhs(8,5)=0;
            lhs(8,6)=0;
            lhs(8,7)=0;
            lhs(8,8)=clhs6;
            lhs(8,9)=0;
            lhs(8,10)=0;
            lhs(8,11)=0;
            lhs(9,0)=0;
            lhs(9,1)=clhs3;
            lhs(9,2)=0;
            lhs(9,3)=0;
            lhs(9,4)=0;
            lhs(9,5)=clhs5;
            lhs(9,6)=0;
            lhs(9,7)=0;
            lhs(9,8)=0;
            lhs(9,9)=clhs6;
            lhs(9,10)=0;
            lhs(9,11)=0;
            lhs(10,0)=0;
            lhs(10,1)=0;
            lhs(10,2)=clhs3;
            lhs(10,3)=0;
            lhs(10,4)=0;
            lhs(10,5)=0;
            lhs(10,6)=clhs5;
            lhs(10,7)=0;
            lhs(10,8)=0;
            lhs(10,9)=0;
            lhs(10,10)=clhs6;
            lhs(10,11)=0;
            lhs(11,0)=0;
            lhs(11,1)=0;
            lhs(11,2)=0;
            lhs(11,3)=clhs3;
            lhs(11,4)=0;
            lhs(11,5)=0;
            lhs(11,6)=0;
            lhs(11,7)=clhs5;
            lhs(11,8)=0;
            lhs(11,9)=0;
            lhs(11,10)=0;
            lhs(11,11)=clhs6;


}


template<>
void CompressibleNavierStokes<3>::ComputeGaussPointRHSContribution(array_1d<double,20>& rhs, const ElementDataStruct& data)
{
    const int nnodes = 4;
    const int dim = 3;
    const int BlockSize = dim+2;
    const double h = data.h; 

    const double& bdf0 = data.bdf0;
    const double& bdf1 = data.bdf1;
    const double& bdf2 = data.bdf2;
    
    const bounded_matrix<double,nnodes,BlockSize>& U = data.U;
    const bounded_matrix<double,nnodes,BlockSize>& Un = data.Un;
    const bounded_matrix<double,nnodes,BlockSize>& Unn = data.Unn;
    const bounded_matrix<double,nnodes,dim>& f_ext = data.f_ext;
    const array_1d<double,nnodes>& r = data.r;
    const double mu = data.mu;
    const double nu = data.nu;
    const double lambda = data.lambda;
    const double cv = data.cv;
    const double y = data.y;
    const double cp = cv*y;
    
    
    // Get shape function values
    const array_1d<double,nnodes>& N = data.N;
    const bounded_matrix<double,nnodes,dim>& DN = data.DN_DX;

    // Auxiliary variables used in the calculation of the RHS
    const array_1d<double,BlockSize> U_gauss = prod(trans(U), N);
    const array_1d<double,dim> f_gauss = prod(trans(f_ext), N);
    const bounded_matrix<double,dim,BlockSize> grad_U = prod(trans(DN), U);
    const array_1d<double,BlockSize> accel_gauss = bdf0*U_gauss+bdf1*prod(trans(Un), N)+bdf2*prod(trans(Unn), N);
    
    // Stabilization parameters
    const double stab_c1 = 4.0;
    const double stab_c2 = 2.0;
    double tmp = U_gauss(dim+1)/U_gauss(0);
    for(unsigned int ll=0; ll<dim; ll++)
        tmp -=(U_gauss(ll+1)*U_gauss(ll+1))/(2*U_gauss(0)*U_gauss(0));
    double c = sqrt(y*(y-1)*tmp);

    double tau1inv = 0.0;
    for(unsigned int ll=0; ll<dim; ll++)
        tau1inv += (U_gauss(ll+1)/U_gauss(0))*(U_gauss(ll+1)/U_gauss(0));
    tau1inv = (sqrt(tau1inv)+c)*stab_c2/h;
    double tau2inv = stab_c1*nu/(h*h)+tau1inv;
    double tau3inv = stab_c1*lambda/(U_gauss(0)*cp*h*h)+tau1inv;
        
    const double tau1 = 1/tau1inv;
    const double tau2 = 1/tau2inv;
    const double tau3 = 1/tau3inv;
    //substitute_rhs_3D
}


template<>
void CompressibleNavierStokes<2>::ComputeGaussPointRHSContribution(array_1d<double,12>& rhs, const ElementDataStruct& data)
{
    const int nnodes = 3;
    const int dim = 2;
    const int BlockSize = dim+2;
    const double h = data.h;

    const double& bdf0 = data.bdf0;
    const double& bdf1 = data.bdf1;
    const double& bdf2 = data.bdf2;

    const bounded_matrix<double,nnodes,BlockSize>& U = data.U;
    const bounded_matrix<double,nnodes,BlockSize>& Un = data.Un;
    const bounded_matrix<double,nnodes,BlockSize>& Unn = data.Unn;
    const bounded_matrix<double,nnodes,dim>& f_ext = data.f_ext;
    const array_1d<double,nnodes>& r = data.r;
    const double mu = data.mu;
    const double nu = data.nu;
    const double lambda = data.lambda;
    const double cv = data.cv;
    const double y = data.y;
    const double cp = cv*y;
    

    // Get shape function values
    const array_1d<double,nnodes>& N = data.N;
    const bounded_matrix<double,nnodes,dim>& DN = data.DN_DX;

    // Auxiliary variables used in the calculation of the RHS
    const array_1d<double,BlockSize> U_gauss = prod(trans(U), N);
    const array_1d<double,dim> f_gauss = prod(trans(f_ext), N);
    const bounded_matrix<double,dim,BlockSize> grad_U = prod(trans(DN), U);
    const array_1d<double,BlockSize> accel_gauss = bdf0*U_gauss+bdf1*prod(trans(Un), N)+bdf2*prod(trans(Unn), N);
    
    // Stabilization parameters
    const double stab_c1 = 4.0;
    const double stab_c2 = 2.0;
    double tmp = U_gauss(dim+1)/U_gauss(0);
    for(unsigned int ll=0; ll<dim; ll++)
        tmp -=(U_gauss(ll+1)*U_gauss(ll+1))/(2*U_gauss(0)*U_gauss(0));
    double c = sqrt(y*(y-1)*tmp);

    double tau1inv = 0.0;
    for(unsigned int ll=0; ll<dim; ll++)
        tau1inv += (U_gauss(ll+1)/U_gauss(0))*(U_gauss(ll+1)/U_gauss(0));
    tau1inv = (sqrt(tau1inv)+c)*stab_c2/h;
    double tau2inv = stab_c1*nu/(h*h)+tau1inv;
    double tau3inv = stab_c1*lambda/(U_gauss(0)*cp*h*h)+tau1inv;
        
    const double tau1 = 1/tau1inv;
    const double tau2 = 1/tau2inv;
    const double tau3 = 1/tau3inv;

    const double crhs0 =             N[0]*(U(0,0)*bdf0 + Un(0,0)*bdf1 + Unn(0,0)*bdf2) + N[1]*(U(1,0)*bdf0 + Un(1,0)*bdf1 + Unn(1,0)*bdf2) + N[2]*(U(2,0)*bdf0 + Un(2,0)*bdf1 + Unn(2,0)*bdf2);
const double crhs1 =             N[0]*(U(0,1)*bdf0 + Un(0,1)*bdf1 + Unn(0,1)*bdf2) + N[1]*(U(1,1)*bdf0 + Un(1,1)*bdf1 + Unn(1,1)*bdf2) + N[2]*(U(2,1)*bdf0 + Un(2,1)*bdf1 + Unn(2,1)*bdf2);
const double crhs2 =             N[0]*(U(0,2)*bdf0 + Un(0,2)*bdf1 + Unn(0,2)*bdf2) + N[1]*(U(1,2)*bdf0 + Un(1,2)*bdf1 + Unn(1,2)*bdf2) + N[2]*(U(2,2)*bdf0 + Un(2,2)*bdf1 + Unn(2,2)*bdf2);
const double crhs3 =             N[0]*(U(0,3)*bdf0 + Un(0,3)*bdf1 + Unn(0,3)*bdf2) + N[1]*(U(1,3)*bdf0 + Un(1,3)*bdf1 + Unn(1,3)*bdf2) + N[2]*(U(2,3)*bdf0 + Un(2,3)*bdf1 + Unn(2,3)*bdf2);
            rhs[0]=N[0]*crhs0;
            rhs[1]=N[0]*crhs1;
            rhs[2]=N[0]*crhs2;
            rhs[3]=N[0]*crhs3;
            rhs[4]=N[1]*crhs0;
            rhs[5]=N[1]*crhs1;
            rhs[6]=N[1]*crhs2;
            rhs[7]=N[1]*crhs3;
            rhs[8]=N[2]*crhs0;
            rhs[9]=N[2]*crhs1;
            rhs[10]=N[2]*crhs2;
            rhs[11]=N[2]*crhs3;

}

/*
template<>
double CompressibleNavierStokes<3>::SubscaleErrorEstimate(const ElementDataStruct& data)
{
    const int nnodes = 4;
    const int dim = 3;
    const int BlockSize = dim+2;

    const double h = data.h;                                // Characteristic element size

    const double& bdf0 = data.bdf0;
    const double& bdf1 = data.bdf1;
    const double& bdf2 = data.bdf2;

    const bounded_matrix<double,nnodes,BlockSize>& U = data.U;
    const bounded_matrix<double,nnodes,BlockSize>& Un = data.Un;
    const bounded_matrix<double,nnodes,BlockSize>& Unn = data.Unn;
    const bounded_matrix<double,nnodes,dim>& f_ext = data.f_ext;
    const array_1d<double,nnodes>& r = data.r;
          
    const double mu = data.mu;
    const double nu = data.nu;
    const double lambda = data.lambda;
    const double cv = data.cv;
    const double cp = data.cp;
    const double y = data.y;

    // Get shape function values
    const array_1d<double,nnodes>& N = data.N;
    const bounded_matrix<double,nnodes,dim>& DN = data.DN_DX;

    // Auxiliary variables used in the calculation of the error estimator
    array_1d<double,BlockSize> U_s_gauss; //WHAT IS THIS FOR?
    const array_1d<double,BlockSize> U_gauss = prod(trans(U), N);
    const array_1d<double,dim> f_gauss = prod(trans(f_ext), N);
    const bounded_matrix<double,dim,BlockSize> grad_U = prod(trans(DN), U);
    const double& r_gauss = inner_prod(data.N, data.r);
    
    const array_1d<double,BlockSize> accel_gauss = bdf0*U_gauss+bdf1*prod(trans(Un), N)+bdf2*prod(trans(Unn), N);
    
    // Stabilization parameters
    const double stab_c1 = 4.0;
    const double stab_c2 = 2.0;
   
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
    const int BlockSize = dim+2;

   const double h = data.h;                                // Characteristic element size

    const double& bdf0 = data.bdf0;
    const double& bdf1 = data.bdf1;
    const double& bdf2 = data.bdf2;

    const bounded_matrix<double,nnodes,BlockSize>& U = data.U;
    const bounded_matrix<double,nnodes,BlockSize>& Un = data.Un;
    const bounded_matrix<double,nnodes,BlockSize>& Unn = data.Unn;
    const bounded_matrix<double,nnodes,dim>& f_ext = data.f_ext;
    const array_1d<double,nnodes>& r = data.r;
          
    const double mu = data.mu;
    const double nu = data.nu;
    const double lambda = data.lambda;
    const double cv = data.cv;
    const double cp = data.cp;
    const double y = data.y;

    // Get shape function values
    const array_1d<double,nnodes>& N = data.N;
    const bounded_matrix<double,nnodes,dim>& DN = data.DN_DX;

    // Auxiliary variables used in the calculation of the error estimator
    array_1d<double,BlockSize> U_s_gauss; //WHAT IS THIS FOR?
    const array_1d<double,BlockSize> U_gauss = prod(trans(U), N);
    const array_1d<double,dim> f_gauss = prod(trans(f_ext), N);
    const bounded_matrix<double,dim,BlockSize> grad_U = prod(trans(DN), U);
    const double& r_gauss = inner_prod(data.N, data.r);
    
    const array_1d<double,BlockSize> accel_gauss = bdf0*U_gauss+bdf1*prod(trans(Un), N)+bdf2*prod(trans(Unn), N);

    // Gauss point velocity subscale value computation
    //substitute_gausspt_subscale_2D

    const double U_gauss_norm = norm_2(U_gauss);
    const double U_s_gauss_norm = norm_2(U_s_gauss);

    return U_s_gauss_norm/U_gauss_norm;
}
*/

}
