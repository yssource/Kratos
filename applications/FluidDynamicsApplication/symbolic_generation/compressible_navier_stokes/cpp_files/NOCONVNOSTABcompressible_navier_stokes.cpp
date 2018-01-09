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
        rResult[i*(BlockSize)+1]  =  this->GetGeometry()[i].GetDof(MOMENTUM_X).EquationId();
        rResult[i*(BlockSize)+2]  =  this->GetGeometry()[i].GetDof(MOMENTUM_Y).EquationId();
        rResult[i*(BlockSize)+3]  =  this->GetGeometry()[i].GetDof(MOMENTUM_Z).EquationId(); 
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
        rResult[i*(BlockSize)+1]  =  this->GetGeometry()[i].GetDof(MOMENTUM_X).EquationId();
        rResult[i*(BlockSize)+2]  =  this->GetGeometry()[i].GetDof(MOMENTUM_Y).EquationId();
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
        ElementalDofList[i*(BlockSize)+1]  =  this->GetGeometry()[i].pGetDof(MOMENTUM_X);
        ElementalDofList[i*(BlockSize)+2]  =  this->GetGeometry()[i].pGetDof(MOMENTUM_Y);
        ElementalDofList[i*(BlockSize)+3]  =  this->GetGeometry()[i].pGetDof(MOMENTUM_Z);
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
        ElementalDofList[i*(BlockSize)+1]  =  this->GetGeometry()[i].pGetDof(MOMENTUM_X);
        ElementalDofList[i*(BlockSize)+2]  =  this->GetGeometry()[i].pGetDof(MOMENTUM_Y);
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

    const double clhs0 =             pow(N[0], 2);
const double clhs1 =             bdf0*clhs0;
const double clhs2 =             -clhs1;
const double clhs3 =             N[0]*bdf0;
const double clhs4 =             N[1]*clhs3;
const double clhs5 =             -clhs4;
const double clhs6 =             N[2]*clhs3;
const double clhs7 =             -clhs6;
const double clhs8 =             N[0]*f_ext(0,0) + N[1]*f_ext(1,0) + N[2]*f_ext(2,0);
const double clhs9 =             clhs0*clhs8;
const double clhs10 =             N[0]*U(0,0) + N[1]*U(1,0) + N[2]*U(2,0);
const double clhs11 =             pow(clhs10, -2);
const double clhs12 =             (1.0L/3.0L)*DN(0,0)*clhs11*mu;
const double clhs13 =             4*DN(0,0);
const double clhs14 =             N[0]*U(0,1) + N[1]*U(1,1) + N[2]*U(2,1);
const double clhs15 =             3*DN(0,1);
const double clhs16 =             N[0]*U(0,2) + N[1]*U(1,2) + N[2]*U(2,2);
const double clhs17 =             DN(0,0)*U(0,1) + DN(1,0)*U(1,1) + DN(2,0)*U(2,1);
const double clhs18 =             N[0]*clhs17;
const double clhs19 =             DN(0,1)*U(0,2) + DN(1,1)*U(1,2) + DN(2,1)*U(2,2);
const double clhs20 =             N[0]*clhs19;
const double clhs21 =             3*clhs20;
const double clhs22 =             DN(0,0)*U(0,0) + DN(1,0)*U(1,0) + DN(2,0)*U(2,0);
const double clhs23 =             1.0/clhs10;
const double clhs24 =             N[0]*clhs14*clhs22*clhs23;
const double clhs25 =             6*N[0]*U(0,2) + 6*N[1]*U(1,2) + 6*N[2]*U(2,2);
const double clhs26 =             DN(0,1)*U(0,0) + DN(1,1)*U(1,0) + DN(2,1)*U(2,0);
const double clhs27 =             N[0]*clhs23;
const double clhs28 =             clhs26*clhs27;
const double clhs29 =             clhs13*clhs14 + clhs15*clhs16 + 4*clhs18 + clhs21 - 8*clhs24 - clhs25*clhs28;
const double clhs30 =             (1.0L/3.0L)*DN(0,1)*clhs11*mu;
const double clhs31 =             2*DN(0,0);
const double clhs32 =             clhs14*clhs15;
const double clhs33 =             DN(0,0)*U(0,2) + DN(1,0)*U(1,2) + DN(2,0)*U(2,2);
const double clhs34 =             N[0]*clhs33;
const double clhs35 =             2*clhs34;
const double clhs36 =             DN(0,1)*U(0,1) + DN(1,1)*U(1,1) + DN(2,1)*U(2,1);
const double clhs37 =             N[0]*clhs36;
const double clhs38 =             3*clhs37;
const double clhs39 =             clhs22*clhs27;
const double clhs40 =             4*clhs39;
const double clhs41 =             6*N[0]*U(0,1) + 6*N[1]*U(1,1) + 6*N[2]*U(2,1);
const double clhs42 =             -clhs16*clhs31 + clhs16*clhs40 - clhs28*clhs41 + clhs29 + clhs32 - clhs35 + clhs38;
const double clhs43 =             DN(0,0)*clhs23*mu;
const double clhs44 =             DN(0,0) - clhs39;
const double clhs45 =             clhs43*clhs44;
const double clhs46 =             (1.0L/3.0L)*DN(0,1)*clhs23*mu;
const double clhs47 =             3*clhs28;
const double clhs48 =             clhs13 + clhs15 - clhs40 - clhs47;
const double clhs49 =             (1.0L/3.0L)*clhs23*mu;
const double clhs50 =             3*DN(0,0);
const double clhs51 =             DN(0,1) - clhs28;
const double clhs52 =             2*clhs39;
const double clhs53 =             -clhs15 + clhs31 + clhs47 - clhs52;
const double clhs54 =             N[0]*N[1];
const double clhs55 =             clhs54*clhs8;
const double clhs56 =             4*DN(1,0);
const double clhs57 =             3*DN(1,1);
const double clhs58 =             N[1]*clhs17;
const double clhs59 =             N[1]*clhs19;
const double clhs60 =             3*clhs59;
const double clhs61 =             8*N[0]*U(0,1) + 8*N[1]*U(1,1) + 8*N[2]*U(2,1);
const double clhs62 =             N[1]*clhs23;
const double clhs63 =             clhs22*clhs62;
const double clhs64 =             clhs26*clhs62;
const double clhs65 =             clhs14*clhs56 + clhs16*clhs57 - clhs25*clhs64 + 4*clhs58 + clhs60 - clhs61*clhs63;
const double clhs66 =             2*DN(1,0);
const double clhs67 =             clhs14*clhs57;
const double clhs68 =             N[1]*clhs33;
const double clhs69 =             2*clhs68;
const double clhs70 =             N[1]*clhs36;
const double clhs71 =             3*clhs70;
const double clhs72 =             4*clhs63;
const double clhs73 =             -clhs16*clhs66 + clhs16*clhs72 - clhs41*clhs64 + clhs65 + clhs67 - clhs69 + clhs71;
const double clhs74 =             DN(1,0) - clhs63;
const double clhs75 =             clhs43*clhs74;
const double clhs76 =             3*clhs64;
const double clhs77 =             clhs56 + clhs57 - clhs72 - clhs76;
const double clhs78 =             DN(1,1) - clhs64;
const double clhs79 =             2*clhs63;
const double clhs80 =             -clhs57 + clhs66 + clhs76 - clhs79;
const double clhs81 =             N[0]*N[2];
const double clhs82 =             clhs8*clhs81;
const double clhs83 =             4*DN(2,0);
const double clhs84 =             3*DN(2,1);
const double clhs85 =             N[2]*clhs17;
const double clhs86 =             N[2]*clhs19;
const double clhs87 =             3*clhs86;
const double clhs88 =             N[2]*clhs23;
const double clhs89 =             clhs22*clhs88;
const double clhs90 =             clhs26*clhs88;
const double clhs91 =             clhs14*clhs83 + clhs16*clhs84 - clhs25*clhs90 - clhs61*clhs89 + 4*clhs85 + clhs87;
const double clhs92 =             2*DN(2,0);
const double clhs93 =             clhs14*clhs84;
const double clhs94 =             N[2]*clhs33;
const double clhs95 =             2*clhs94;
const double clhs96 =             N[2]*clhs36;
const double clhs97 =             3*clhs96;
const double clhs98 =             4*clhs89;
const double clhs99 =             -clhs16*clhs92 + clhs16*clhs98 - clhs41*clhs90 + clhs91 + clhs93 - clhs95 + clhs97;
const double clhs100 =             DN(2,0) - clhs89;
const double clhs101 =             clhs100*clhs43;
const double clhs102 =             3*clhs90;
const double clhs103 =             -clhs102 + clhs83 + clhs84 - clhs98;
const double clhs104 =             DN(2,1) - clhs90;
const double clhs105 =             2*clhs89;
const double clhs106 =             clhs102 - clhs105 - clhs84 + clhs92;
const double clhs107 =             N[0]*f_ext(0,1) + N[1]*f_ext(1,1) + N[2]*f_ext(2,1);
const double clhs108 =             clhs0*clhs107;
const double clhs109 =             clhs16*clhs50;
const double clhs110 =             2*DN(0,1);
const double clhs111 =             3*clhs34;
const double clhs112 =             2*clhs37;
const double clhs113 =             4*clhs28;
const double clhs114 =             clhs109 - clhs110*clhs14 + clhs111 - clhs112 + clhs113*clhs14 - clhs25*clhs39;
const double clhs115 =             4*DN(0,1);
const double clhs116 =             3*clhs18;
const double clhs117 =             8*N[0]*U(0,2) + 8*N[1]*U(1,2) + 8*N[2]*U(2,2);
const double clhs118 =             clhs114 + clhs115*clhs16 + clhs116 - clhs117*clhs28 + clhs14*clhs50 + 4*clhs20 - 6*clhs24;
const double clhs119 =             3*clhs39;
const double clhs120 =             2*clhs28;
const double clhs121 =             clhs110 + clhs119 - clhs120 - clhs50;
const double clhs122 =             -clhs113 + clhs115 - clhs119 + clhs50;
const double clhs123 =             clhs107*clhs54;
const double clhs124 =             3*DN(1,0);
const double clhs125 =             clhs124*clhs16;
const double clhs126 =             2*DN(1,1);
const double clhs127 =             3*clhs68;
const double clhs128 =             2*clhs70;
const double clhs129 =             4*clhs64;
const double clhs130 =             clhs125 - clhs126*clhs14 + clhs127 - clhs128 + clhs129*clhs14 - clhs25*clhs63;
const double clhs131 =             4*DN(1,1);
const double clhs132 =             3*clhs58;
const double clhs133 =             -clhs117*clhs64 + clhs124*clhs14 + clhs130 + clhs131*clhs16 + clhs132 - clhs41*clhs63 + 4*clhs59;
const double clhs134 =             3*clhs63;
const double clhs135 =             2*clhs64;
const double clhs136 =             -clhs124 + clhs126 + clhs134 - clhs135;
const double clhs137 =             clhs124 - clhs129 + clhs131 - clhs134;
const double clhs138 =             clhs107*clhs81;
const double clhs139 =             3*DN(2,0);
const double clhs140 =             clhs139*clhs16;
const double clhs141 =             2*DN(2,1);
const double clhs142 =             3*clhs94;
const double clhs143 =             2*clhs96;
const double clhs144 =             4*clhs90;
const double clhs145 =             -clhs14*clhs141 + clhs14*clhs144 + clhs140 + clhs142 - clhs143 - clhs25*clhs89;
const double clhs146 =             4*DN(2,1);
const double clhs147 =             3*clhs85;
const double clhs148 =             -clhs117*clhs90 + clhs139*clhs14 + clhs145 + clhs146*clhs16 + clhs147 - clhs41*clhs89 + 4*clhs86;
const double clhs149 =             3*clhs89;
const double clhs150 =             2*clhs90;
const double clhs151 =             -clhs139 + clhs141 + clhs149 - clhs150;
const double clhs152 =             clhs139 - clhs144 + clhs146 - clhs149;
const double clhs153 =             N[0]*r[0] + N[1]*r[1] + N[2]*r[2];
const double clhs154 =             (1.0L/3.0L)*DN(0,0)*clhs11;
const double clhs155 =             1.0/cv;
const double clhs156 =             3*N[0]*clhs155*lambda;
const double clhs157 =             DN(0,0)*U(0,3) + DN(1,0)*U(1,3) + DN(2,0)*U(2,3);
const double clhs158 =             clhs14*clhs16*clhs23*mu;
const double clhs159 =             4*clhs16*clhs23*mu;
const double clhs160 =             6*clhs14*clhs23*mu;
const double clhs161 =             clhs155*lambda;
const double clhs162 =             -clhs161 + mu;
const double clhs163 =             clhs162*clhs23*clhs25;
const double clhs164 =             4*mu;
const double clhs165 =             3*clhs161;
const double clhs166 =             clhs164 - clhs165;
const double clhs167 =             2*clhs14*clhs166*clhs23;
const double clhs168 =             3*N[0]*clhs11*clhs14*clhs16*mu;
const double clhs169 =             clhs165*(N[0]*U(0,3) + N[1]*U(1,3) + N[2]*U(2,3));
const double clhs170 =             pow(clhs14, 2)*clhs23;
const double clhs171 =             3*mu;
const double clhs172 =             pow(clhs16, 2)*clhs23;
const double clhs173 =             -clhs165*clhs170;
const double clhs174 =             -clhs165*clhs172;
const double clhs175 =             clhs164*clhs170 + clhs171*clhs172 + clhs173 + clhs174;
const double clhs176 =             clhs169 + clhs175;
const double clhs177 =             clhs161*(2*N[0]*U(0,3) + 2*N[1]*U(1,3) + 2*N[2]*U(2,3));
const double clhs178 =             clhs175 + clhs177;
const double clhs179 =             -DN(0,0)*clhs176 - DN(0,1)*clhs158 + clhs119*clhs178 - clhs156*clhs157 + clhs159*clhs37 - clhs160*clhs20 - clhs163*clhs34 - clhs167*clhs18 + clhs168*clhs26;
const double clhs180 =             (1.0L/3.0L)*DN(0,1)*clhs11;
const double clhs181 =             DN(0,1)*U(0,3) + DN(1,1)*U(1,3) + DN(2,1)*U(2,3);
const double clhs182 =             6*clhs16*clhs23*mu;
const double clhs183 =             4*clhs14*clhs23*mu;
const double clhs184 =             clhs162*clhs23*clhs41;
const double clhs185 =             2*clhs16*clhs166*clhs23;
const double clhs186 =             clhs164*clhs172 + clhs170*clhs171 + clhs173 + clhs174;
const double clhs187 =             clhs169 + clhs186;
const double clhs188 =             clhs177 + clhs186;
const double clhs189 =             -DN(0,0)*clhs158 - DN(0,1)*clhs187 - clhs156*clhs181 + clhs168*clhs22 + clhs179 - clhs18*clhs182 + clhs183*clhs34 - clhs184*clhs37 - clhs185*clhs20 + clhs188*clhs47;
const double clhs190 =             clhs16*mu;
const double clhs191 =             clhs14*clhs166;
const double clhs192 =             DN(0,0)*clhs191 - clhs110*clhs190 + clhs166*clhs18 - clhs190*clhs28 - clhs191*clhs52 + clhs21*mu;
const double clhs193 =             -N[0]*clhs162*clhs23*clhs26*clhs41 + clhs109*mu + clhs162*clhs32 + clhs162*clhs38 - clhs190*clhs39 + clhs192 - clhs35*mu;
const double clhs194 =             clhs14*mu;
const double clhs195 =             -N[0]*clhs162*clhs22*clhs23*clhs25 + clhs109*clhs162 + clhs111*clhs162 - clhs112*mu - clhs194*clhs28 + clhs32*mu;
const double clhs196 =             clhs16*clhs166;
const double clhs197 =             DN(0,1)*clhs196 + clhs116*mu - clhs120*clhs196 + clhs166*clhs20 - clhs194*clhs31 - clhs194*clhs39 + clhs195;
const double clhs198 =             DN(0,0)*clhs155*clhs23*lambda;
const double clhs199 =             DN(0,1)*clhs155*clhs23*lambda;
const double clhs200 =             clhs44 + clhs51;
const double clhs201 =             clhs153*clhs54;
const double clhs202 =             3*N[1]*clhs155*lambda;
const double clhs203 =             3*N[1]*clhs11*clhs14*clhs16*mu;
const double clhs204 =             -DN(1,0)*clhs176 - DN(1,1)*clhs158 + clhs134*clhs178 - clhs157*clhs202 + clhs159*clhs70 - clhs160*clhs59 - clhs163*clhs68 - clhs167*clhs58 + clhs203*clhs26;
const double clhs205 =             -DN(1,0)*clhs158 - DN(1,1)*clhs187 - clhs181*clhs202 - clhs182*clhs58 + clhs183*clhs68 - clhs184*clhs70 - clhs185*clhs59 + clhs188*clhs76 + clhs203*clhs22 + clhs204;
const double clhs206 =             DN(1,0)*clhs191 - clhs126*clhs190 + clhs166*clhs58 - clhs190*clhs64 - clhs191*clhs79 + clhs60*mu;
const double clhs207 =             -N[1]*clhs162*clhs23*clhs26*clhs41 + clhs125*mu + clhs162*clhs67 + clhs162*clhs71 - clhs190*clhs63 + clhs206 - clhs69*mu;
const double clhs208 =             -N[1]*clhs162*clhs22*clhs23*clhs25 + clhs125*clhs162 + clhs127*clhs162 - clhs128*mu - clhs194*clhs64 + clhs67*mu;
const double clhs209 =             DN(1,1)*clhs196 + clhs132*mu - clhs135*clhs196 + clhs166*clhs59 - clhs194*clhs63 - clhs194*clhs66 + clhs208;
const double clhs210 =             clhs74 + clhs78;
const double clhs211 =             clhs153*clhs81;
const double clhs212 =             3*N[2]*clhs155*lambda;
const double clhs213 =             3*N[2]*clhs11*clhs14*clhs16*mu;
const double clhs214 =             -DN(2,0)*clhs176 - DN(2,1)*clhs158 + clhs149*clhs178 - clhs157*clhs212 + clhs159*clhs96 - clhs160*clhs86 - clhs163*clhs94 - clhs167*clhs85 + clhs213*clhs26;
const double clhs215 =             -DN(2,0)*clhs158 - DN(2,1)*clhs187 + clhs102*clhs188 - clhs181*clhs212 - clhs182*clhs85 + clhs183*clhs94 - clhs184*clhs96 - clhs185*clhs86 + clhs213*clhs22 + clhs214;
const double clhs216 =             DN(2,0)*clhs191 - clhs105*clhs191 - clhs141*clhs190 + clhs166*clhs85 - clhs190*clhs90 + clhs87*mu;
const double clhs217 =             -N[2]*clhs162*clhs23*clhs26*clhs41 + clhs140*mu + clhs162*clhs93 + clhs162*clhs97 - clhs190*clhs89 + clhs216 - clhs95*mu;
const double clhs218 =             -N[2]*clhs162*clhs22*clhs23*clhs25 + clhs140*clhs162 + clhs142*clhs162 - clhs143*mu - clhs194*clhs90 + clhs93*mu;
const double clhs219 =             DN(2,1)*clhs196 + clhs147*mu - clhs150*clhs196 + clhs166*clhs86 - clhs194*clhs89 - clhs194*clhs92 + clhs218;
const double clhs220 =             clhs100 + clhs104;
const double clhs221 =             pow(N[1], 2);
const double clhs222 =             bdf0*clhs221;
const double clhs223 =             -clhs222;
const double clhs224 =             N[1]*N[2];
const double clhs225 =             bdf0*clhs224;
const double clhs226 =             -clhs225;
const double clhs227 =             (1.0L/3.0L)*DN(1,0)*clhs11*mu;
const double clhs228 =             (1.0L/3.0L)*DN(1,1)*clhs11*mu;
const double clhs229 =             DN(1,0)*clhs23*mu;
const double clhs230 =             clhs229*clhs44;
const double clhs231 =             (1.0L/3.0L)*DN(1,1)*clhs23*mu;
const double clhs232 =             clhs221*clhs8;
const double clhs233 =             clhs229*clhs74;
const double clhs234 =             clhs224*clhs8;
const double clhs235 =             clhs100*clhs229;
const double clhs236 =             clhs107*clhs221;
const double clhs237 =             clhs107*clhs224;
const double clhs238 =             (1.0L/3.0L)*DN(1,0)*clhs11;
const double clhs239 =             (1.0L/3.0L)*DN(1,1)*clhs11;
const double clhs240 =             DN(1,0)*clhs155*clhs23*lambda;
const double clhs241 =             DN(1,1)*clhs155*clhs23*lambda;
const double clhs242 =             clhs153*clhs224;
const double clhs243 =             pow(N[2], 2);
const double clhs244 =             bdf0*clhs243;
const double clhs245 =             -clhs244;
const double clhs246 =             (1.0L/3.0L)*DN(2,0)*clhs11*mu;
const double clhs247 =             (1.0L/3.0L)*DN(2,1)*clhs11*mu;
const double clhs248 =             DN(2,0)*clhs23*mu;
const double clhs249 =             clhs248*clhs44;
const double clhs250 =             (1.0L/3.0L)*DN(2,1)*clhs23*mu;
const double clhs251 =             clhs248*clhs74;
const double clhs252 =             clhs243*clhs8;
const double clhs253 =             clhs100*clhs248;
const double clhs254 =             clhs107*clhs243;
const double clhs255 =             (1.0L/3.0L)*DN(2,0)*clhs11;
const double clhs256 =             (1.0L/3.0L)*DN(2,1)*clhs11;
const double clhs257 =             DN(2,0)*clhs155*clhs23*lambda;
const double clhs258 =             DN(2,1)*clhs155*clhs23*lambda;
            lhs(0,0)=clhs2;
            lhs(0,1)=0;
            lhs(0,2)=0;
            lhs(0,3)=0;
            lhs(0,4)=clhs5;
            lhs(0,5)=0;
            lhs(0,6)=0;
            lhs(0,7)=0;
            lhs(0,8)=clhs7;
            lhs(0,9)=0;
            lhs(0,10)=0;
            lhs(0,11)=0;
            lhs(1,0)=clhs12*clhs29 + clhs30*clhs42 + clhs9;
            lhs(1,1)=clhs2 - 4.0L/3.0L*clhs45 - clhs46*clhs48;
            lhs(1,2)=clhs49*(DN(0,1)*clhs53 - clhs50*clhs51);
            lhs(1,3)=0;
            lhs(1,4)=clhs12*clhs65 + clhs30*clhs73 + clhs55;
            lhs(1,5)=-clhs46*clhs77 + clhs5 - 4.0L/3.0L*clhs75;
            lhs(1,6)=clhs49*(DN(0,1)*clhs80 - clhs50*clhs78);
            lhs(1,7)=0;
            lhs(1,8)=clhs12*clhs91 + clhs30*clhs99 + clhs82;
            lhs(1,9)=-4.0L/3.0L*clhs101 - clhs103*clhs46 + clhs7;
            lhs(1,10)=clhs49*(DN(0,1)*clhs106 - clhs104*clhs50);
            lhs(1,11)=0;
            lhs(2,0)=clhs108 + clhs114*clhs12 + clhs118*clhs30;
            lhs(2,1)=clhs49*(DN(0,1)*clhs121 + clhs31*clhs51);
            lhs(2,2)=-clhs122*clhs46 + clhs2 - clhs45;
            lhs(2,3)=0;
            lhs(2,4)=clhs12*clhs130 + clhs123 + clhs133*clhs30;
            lhs(2,5)=clhs49*(DN(0,1)*clhs136 + clhs31*clhs78);
            lhs(2,6)=-clhs137*clhs46 + clhs5 - clhs75;
            lhs(2,7)=0;
            lhs(2,8)=clhs12*clhs145 + clhs138 + clhs148*clhs30;
            lhs(2,9)=clhs49*(DN(0,1)*clhs151 + clhs104*clhs31);
            lhs(2,10)=-clhs101 - clhs152*clhs46 + clhs7;
            lhs(2,11)=0;
            lhs(3,0)=clhs0*clhs153 - clhs154*clhs179 - clhs180*clhs189;
            lhs(3,1)=-clhs154*clhs192 - clhs180*clhs193 + clhs9;
            lhs(3,2)=clhs108 - clhs154*clhs195 - clhs180*clhs197;
            lhs(3,3)=-clhs1 - clhs198*clhs44 - clhs199*clhs200;
            lhs(3,4)=-clhs154*clhs204 - clhs180*clhs205 + clhs201;
            lhs(3,5)=-clhs154*clhs206 - clhs180*clhs207 + clhs55;
            lhs(3,6)=clhs123 - clhs154*clhs208 - clhs180*clhs209;
            lhs(3,7)=-clhs198*clhs74 - clhs199*clhs210 - clhs4;
            lhs(3,8)=-clhs154*clhs214 - clhs180*clhs215 + clhs211;
            lhs(3,9)=-clhs154*clhs216 - clhs180*clhs217 + clhs82;
            lhs(3,10)=clhs138 - clhs154*clhs218 - clhs180*clhs219;
            lhs(3,11)=-clhs100*clhs198 - clhs199*clhs220 - clhs6;
            lhs(4,0)=clhs5;
            lhs(4,1)=0;
            lhs(4,2)=0;
            lhs(4,3)=0;
            lhs(4,4)=clhs223;
            lhs(4,5)=0;
            lhs(4,6)=0;
            lhs(4,7)=0;
            lhs(4,8)=clhs226;
            lhs(4,9)=0;
            lhs(4,10)=0;
            lhs(4,11)=0;
            lhs(5,0)=clhs227*clhs29 + clhs228*clhs42 + clhs55;
            lhs(5,1)=-4.0L/3.0L*clhs230 - clhs231*clhs48 + clhs5;
            lhs(5,2)=clhs49*(DN(1,1)*clhs53 - clhs124*clhs51);
            lhs(5,3)=0;
            lhs(5,4)=clhs227*clhs65 + clhs228*clhs73 + clhs232;
            lhs(5,5)=clhs223 - clhs231*clhs77 - 4.0L/3.0L*clhs233;
            lhs(5,6)=clhs49*(DN(1,1)*clhs80 - clhs124*clhs78);
            lhs(5,7)=0;
            lhs(5,8)=clhs227*clhs91 + clhs228*clhs99 + clhs234;
            lhs(5,9)=-clhs103*clhs231 + clhs226 - 4.0L/3.0L*clhs235;
            lhs(5,10)=clhs49*(DN(1,1)*clhs106 - clhs104*clhs124);
            lhs(5,11)=0;
            lhs(6,0)=clhs114*clhs227 + clhs118*clhs228 + clhs123;
            lhs(6,1)=clhs49*(DN(1,1)*clhs121 + clhs51*clhs66);
            lhs(6,2)=-clhs122*clhs231 - clhs230 + clhs5;
            lhs(6,3)=0;
            lhs(6,4)=clhs130*clhs227 + clhs133*clhs228 + clhs236;
            lhs(6,5)=clhs49*(DN(1,1)*clhs136 + clhs66*clhs78);
            lhs(6,6)=-clhs137*clhs231 + clhs223 - clhs233;
            lhs(6,7)=0;
            lhs(6,8)=clhs145*clhs227 + clhs148*clhs228 + clhs237;
            lhs(6,9)=clhs49*(DN(1,1)*clhs151 + clhs104*clhs66);
            lhs(6,10)=-clhs152*clhs231 + clhs226 - clhs235;
            lhs(6,11)=0;
            lhs(7,0)=-clhs179*clhs238 - clhs189*clhs239 + clhs201;
            lhs(7,1)=-clhs192*clhs238 - clhs193*clhs239 + clhs55;
            lhs(7,2)=clhs123 - clhs195*clhs238 - clhs197*clhs239;
            lhs(7,3)=-clhs200*clhs241 - clhs240*clhs44 - clhs4;
            lhs(7,4)=clhs153*clhs221 - clhs204*clhs238 - clhs205*clhs239;
            lhs(7,5)=-clhs206*clhs238 - clhs207*clhs239 + clhs232;
            lhs(7,6)=-clhs208*clhs238 - clhs209*clhs239 + clhs236;
            lhs(7,7)=-clhs210*clhs241 - clhs222 - clhs240*clhs74;
            lhs(7,8)=-clhs214*clhs238 - clhs215*clhs239 + clhs242;
            lhs(7,9)=-clhs216*clhs238 - clhs217*clhs239 + clhs234;
            lhs(7,10)=-clhs218*clhs238 - clhs219*clhs239 + clhs237;
            lhs(7,11)=-clhs100*clhs240 - clhs220*clhs241 - clhs225;
            lhs(8,0)=clhs7;
            lhs(8,1)=0;
            lhs(8,2)=0;
            lhs(8,3)=0;
            lhs(8,4)=clhs226;
            lhs(8,5)=0;
            lhs(8,6)=0;
            lhs(8,7)=0;
            lhs(8,8)=clhs245;
            lhs(8,9)=0;
            lhs(8,10)=0;
            lhs(8,11)=0;
            lhs(9,0)=clhs246*clhs29 + clhs247*clhs42 + clhs82;
            lhs(9,1)=-4.0L/3.0L*clhs249 - clhs250*clhs48 + clhs7;
            lhs(9,2)=clhs49*(DN(2,1)*clhs53 - clhs139*clhs51);
            lhs(9,3)=0;
            lhs(9,4)=clhs234 + clhs246*clhs65 + clhs247*clhs73;
            lhs(9,5)=clhs226 - clhs250*clhs77 - 4.0L/3.0L*clhs251;
            lhs(9,6)=clhs49*(DN(2,1)*clhs80 - clhs139*clhs78);
            lhs(9,7)=0;
            lhs(9,8)=clhs246*clhs91 + clhs247*clhs99 + clhs252;
            lhs(9,9)=-clhs103*clhs250 + clhs245 - 4.0L/3.0L*clhs253;
            lhs(9,10)=clhs49*(DN(2,1)*clhs106 - clhs104*clhs139);
            lhs(9,11)=0;
            lhs(10,0)=clhs114*clhs246 + clhs118*clhs247 + clhs138;
            lhs(10,1)=clhs49*(DN(2,1)*clhs121 + clhs51*clhs92);
            lhs(10,2)=-clhs122*clhs250 - clhs249 + clhs7;
            lhs(10,3)=0;
            lhs(10,4)=clhs130*clhs246 + clhs133*clhs247 + clhs237;
            lhs(10,5)=clhs49*(DN(2,1)*clhs136 + clhs78*clhs92);
            lhs(10,6)=-clhs137*clhs250 + clhs226 - clhs251;
            lhs(10,7)=0;
            lhs(10,8)=clhs145*clhs246 + clhs148*clhs247 + clhs254;
            lhs(10,9)=clhs49*(DN(2,1)*clhs151 + clhs104*clhs92);
            lhs(10,10)=-clhs152*clhs250 + clhs245 - clhs253;
            lhs(10,11)=0;
            lhs(11,0)=-clhs179*clhs255 - clhs189*clhs256 + clhs211;
            lhs(11,1)=-clhs192*clhs255 - clhs193*clhs256 + clhs82;
            lhs(11,2)=clhs138 - clhs195*clhs255 - clhs197*clhs256;
            lhs(11,3)=-clhs200*clhs258 - clhs257*clhs44 - clhs6;
            lhs(11,4)=-clhs204*clhs255 - clhs205*clhs256 + clhs242;
            lhs(11,5)=-clhs206*clhs255 - clhs207*clhs256 + clhs234;
            lhs(11,6)=-clhs208*clhs255 - clhs209*clhs256 + clhs237;
            lhs(11,7)=-clhs210*clhs258 - clhs225 - clhs257*clhs74;
            lhs(11,8)=clhs153*clhs243 - clhs214*clhs255 - clhs215*clhs256;
            lhs(11,9)=-clhs216*clhs255 - clhs217*clhs256 + clhs252;
            lhs(11,10)=-clhs218*clhs255 - clhs219*clhs256 + clhs254;
            lhs(11,11)=-clhs100*clhs257 - clhs220*clhs258 - clhs244;


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
const double crhs1 =             N[0]*f_ext(0,0) + N[1]*f_ext(1,0) + N[2]*f_ext(2,0);
const double crhs2 =             N[0]*U(0,0) + N[1]*U(1,0) + N[2]*U(2,0);
const double crhs3 =             N[0]*crhs2;
const double crhs4 =             N[0]*(U(0,1)*bdf0 + Un(0,1)*bdf1 + Unn(0,1)*bdf2) + N[1]*(U(1,1)*bdf0 + Un(1,1)*bdf1 + Unn(1,1)*bdf2) + N[2]*(U(2,1)*bdf0 + Un(2,1)*bdf1 + Unn(2,1)*bdf2);
const double crhs5 =             1.0/crhs2;
const double crhs6 =             (1.0L/3.0L)*DN(0,0)*crhs5*mu;
const double crhs7 =             DN(0,0)*U(0,1);
const double crhs8 =             4*crhs7;
const double crhs9 =             DN(0,1)*U(0,2);
const double crhs10 =             3*crhs9;
const double crhs11 =             DN(1,0)*U(1,1);
const double crhs12 =             4*crhs11;
const double crhs13 =             DN(1,1)*U(1,2);
const double crhs14 =             3*crhs13;
const double crhs15 =             DN(2,0)*U(2,1);
const double crhs16 =             4*crhs15;
const double crhs17 =             DN(2,1)*U(2,2);
const double crhs18 =             3*crhs17;
const double crhs19 =             DN(0,0)*U(0,0) + DN(1,0)*U(1,0) + DN(2,0)*U(2,0);
const double crhs20 =             N[0]*U(0,1) + N[1]*U(1,1) + N[2]*U(2,1);
const double crhs21 =             crhs19*crhs20*crhs5;
const double crhs22 =             4*crhs21;
const double crhs23 =             DN(0,1)*U(0,0) + DN(1,1)*U(1,0) + DN(2,1)*U(2,0);
const double crhs24 =             N[0]*U(0,2) + N[1]*U(1,2) + N[2]*U(2,2);
const double crhs25 =             crhs23*crhs24*crhs5;
const double crhs26 =             3*crhs25;
const double crhs27 =             -crhs10 - crhs12 - crhs14 - crhs16 - crhs18 + crhs22 + crhs26 - crhs8;
const double crhs28 =             (1.0L/3.0L)*DN(0,1)*crhs5*mu;
const double crhs29 =             DN(0,0)*U(0,2);
const double crhs30 =             DN(0,1)*U(0,1);
const double crhs31 =             DN(1,0)*U(1,2);
const double crhs32 =             DN(1,1)*U(1,1);
const double crhs33 =             DN(2,0)*U(2,2);
const double crhs34 =             DN(2,1)*U(2,1);
const double crhs35 =             crhs19*crhs24*crhs5;
const double crhs36 =             crhs20*crhs23*crhs5;
const double crhs37 =             crhs10 + crhs12 + crhs14 + crhs16 + crhs18 - crhs22 - crhs26 - 2*crhs29 + 3*crhs30 - 2*crhs31 + 3*crhs32 - 2*crhs33 + 3*crhs34 + 2*crhs35 - 3*crhs36 + crhs8;
const double crhs38 =             N[0]*f_ext(0,1) + N[1]*f_ext(1,1) + N[2]*f_ext(2,1);
const double crhs39 =             N[0]*(U(0,2)*bdf0 + Un(0,2)*bdf1 + Unn(0,2)*bdf2) + N[1]*(U(1,2)*bdf0 + Un(1,2)*bdf1 + Unn(1,2)*bdf2) + N[2]*(U(2,2)*bdf0 + Un(2,2)*bdf1 + Unn(2,2)*bdf2);
const double crhs40 =             3*crhs29;
const double crhs41 =             2*crhs30;
const double crhs42 =             3*crhs31;
const double crhs43 =             2*crhs32;
const double crhs44 =             3*crhs33;
const double crhs45 =             2*crhs34;
const double crhs46 =             3*crhs35;
const double crhs47 =             2*crhs36;
const double crhs48 =             -crhs40 + crhs41 - crhs42 + crhs43 - crhs44 + crhs45 + crhs46 - crhs47;
const double crhs49 =             3*crhs11 + 4*crhs13 + 3*crhs15 + 4*crhs17 - 3*crhs21 - 4*crhs25 + crhs40 - crhs41 + crhs42 - crhs43 + crhs44 - crhs45 - crhs46 + crhs47 + 3*crhs7 + 4*crhs9;
const double crhs50 =             N[0]*(U(0,3)*bdf0 + Un(0,3)*bdf1 + Unn(0,3)*bdf2) + N[1]*(U(1,3)*bdf0 + Un(1,3)*bdf1 + Unn(1,3)*bdf2) + N[2]*(U(2,3)*bdf0 + Un(2,3)*bdf1 + Unn(2,3)*bdf2);
const double crhs51 =             crhs1*crhs20 + crhs2*(N[0]*r[0] + N[1]*r[1] + N[2]*r[2]) + crhs24*crhs38;
const double crhs52 =             3*(lambda/cv);
const double crhs53 =             crhs30 + crhs32 + crhs34;
const double crhs54 =             crhs13 + crhs17 + crhs9;
const double crhs55 =             crhs20*crhs5;
const double crhs56 =             3*mu;
const double crhs57 =             -crhs52;
const double crhs58 =             crhs56 + crhs57;
const double crhs59 =             crhs29 + crhs31 + crhs33;
const double crhs60 =             4*mu;
const double crhs61 =             crhs57 + crhs60;
const double crhs62 =             crhs11 + crhs15 + crhs7;
const double crhs63 =             crhs20*crhs24*mu/pow(crhs2, 2);
const double crhs64 =             pow(crhs20, 2)*crhs5;
const double crhs65 =             pow(crhs24, 2)*crhs5;
const double crhs66 =             -crhs52*crhs64 - crhs52*crhs65 + crhs52*(N[0]*U(0,3) + N[1]*U(1,3) + N[2]*U(2,3));
const double crhs67 =             -crhs19*crhs5*(crhs56*crhs65 + crhs60*crhs64 + crhs66) - crhs23*crhs63 - 2*crhs24*crhs5*crhs53*mu + crhs24*crhs5*crhs58*crhs59 + crhs52*(DN(0,0)*U(0,3) + DN(1,0)*U(1,3) + DN(2,0)*U(2,3)) + 3*crhs54*crhs55*mu + crhs55*crhs61*crhs62;
const double crhs68 =             (1.0L/3.0L)*crhs5*crhs67;
const double crhs69 =             crhs24*crhs5;
const double crhs70 =             (1.0L/3.0L)*crhs5*(-crhs19*crhs63 + crhs20*crhs5*crhs53*crhs58 - crhs23*crhs5*(crhs56*crhs64 + crhs60*crhs65 + crhs66) + crhs52*(DN(0,1)*U(0,3) + DN(1,1)*U(1,3) + DN(2,1)*U(2,3)) + crhs54*crhs61*crhs69 - 2*crhs55*crhs59*mu + 3*crhs62*crhs69*mu + crhs67);
const double crhs71 =             N[1]*crhs2;
const double crhs72 =             (1.0L/3.0L)*DN(1,0)*crhs5*mu;
const double crhs73 =             (1.0L/3.0L)*DN(1,1)*crhs5*mu;
const double crhs74 =             N[2]*crhs2;
const double crhs75 =             (1.0L/3.0L)*DN(2,0)*crhs5*mu;
const double crhs76 =             (1.0L/3.0L)*DN(2,1)*crhs5*mu;
            rhs[0]=N[0]*crhs0;
            rhs[1]=N[0]*crhs4 - crhs1*crhs3 - crhs27*crhs6 + crhs28*crhs37;
            rhs[2]=N[0]*crhs39 + crhs28*crhs49 - crhs3*crhs38 - crhs48*crhs6;
            rhs[3]=DN(0,0)*crhs68 + DN(0,1)*crhs70 + N[0]*crhs50 - N[0]*crhs51;
            rhs[4]=N[1]*crhs0;
            rhs[5]=N[1]*crhs4 - crhs1*crhs71 - crhs27*crhs72 + crhs37*crhs73;
            rhs[6]=N[1]*crhs39 - crhs38*crhs71 - crhs48*crhs72 + crhs49*crhs73;
            rhs[7]=DN(1,0)*crhs68 + DN(1,1)*crhs70 + N[1]*crhs50 - N[1]*crhs51;
            rhs[8]=N[2]*crhs0;
            rhs[9]=N[2]*crhs4 - crhs1*crhs74 - crhs27*crhs75 + crhs37*crhs76;
            rhs[10]=N[2]*crhs39 - crhs38*crhs74 - crhs48*crhs75 + crhs49*crhs76;
            rhs[11]=DN(2,0)*crhs68 + DN(2,1)*crhs70 + N[2]*crhs50 - N[2]*crhs51;

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
