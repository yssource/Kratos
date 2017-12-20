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
const double clhs2 =             DN(0,0)*N[0];
const double clhs3 =             DN(0,1)*N[0];
const double clhs4 =             N[0]*bdf0;
const double clhs5 =             N[1]*clhs4;
const double clhs6 =             -clhs5;
const double clhs7 =             DN(1,0)*N[0];
const double clhs8 =             DN(1,1)*N[0];
const double clhs9 =             N[2]*clhs4;
const double clhs10 =             -clhs9;
const double clhs11 =             DN(2,0)*N[0];
const double clhs12 =             DN(2,1)*N[0];
const double clhs13 =             N[0]*f_ext(0,0) + N[1]*f_ext(1,0) + N[2]*f_ext(2,0);
const double clhs14 =             clhs0*clhs13;
const double clhs15 =             N[0]*U(0,0) + N[1]*U(1,0) + N[2]*U(2,0);
const double clhs16 =             pow(clhs15, -2);
const double clhs17 =             (1.0L/3.0L)*DN(0,0)*clhs16*mu;
const double clhs18 =             4*DN(0,0);
const double clhs19 =             N[0]*U(0,1) + N[1]*U(1,1) + N[2]*U(2,1);
const double clhs20 =             3*DN(0,1);
const double clhs21 =             N[0]*U(0,2) + N[1]*U(1,2) + N[2]*U(2,2);
const double clhs22 =             DN(0,0)*U(0,1) + DN(1,0)*U(1,1) + DN(2,0)*U(2,1);
const double clhs23 =             N[0]*clhs22;
const double clhs24 =             DN(0,1)*U(0,2) + DN(1,1)*U(1,2) + DN(2,1)*U(2,2);
const double clhs25 =             N[0]*clhs24;
const double clhs26 =             3*clhs25;
const double clhs27 =             DN(0,0)*U(0,0) + DN(1,0)*U(1,0) + DN(2,0)*U(2,0);
const double clhs28 =             1.0/clhs15;
const double clhs29 =             N[0]*clhs28;
const double clhs30 =             clhs27*clhs29;
const double clhs31 =             clhs19*clhs30;
const double clhs32 =             DN(0,1)*U(0,0) + DN(1,1)*U(1,0) + DN(2,1)*U(2,0);
const double clhs33 =             clhs29*clhs32;
const double clhs34 =             clhs21*clhs33;
const double clhs35 =             clhs18*clhs19 + clhs20*clhs21 + 4*clhs23 + clhs26 - 8*clhs31 - 6*clhs34;
const double clhs36 =             (1.0L/3.0L)*DN(0,1)*clhs16*mu;
const double clhs37 =             2*DN(0,0);
const double clhs38 =             clhs21*clhs37;
const double clhs39 =             clhs19*clhs20;
const double clhs40 =             DN(0,0)*U(0,2) + DN(1,0)*U(1,2) + DN(2,0)*U(2,2);
const double clhs41 =             N[0]*clhs40;
const double clhs42 =             2*clhs41;
const double clhs43 =             DN(0,1)*U(0,1) + DN(1,1)*U(1,1) + DN(2,1)*U(2,1);
const double clhs44 =             N[0]*clhs43;
const double clhs45 =             3*clhs44;
const double clhs46 =             4*clhs30;
const double clhs47 =             clhs21*clhs46;
const double clhs48 =             clhs19*clhs33;
const double clhs49 =             6*clhs48;
const double clhs50 =             clhs35 - clhs38 + clhs39 - clhs42 + clhs45 + clhs47 - clhs49;
const double clhs51 =             (1.0L/2.0L)*N[0]*clhs16;
const double clhs52 =             2*DN(0,1);
const double clhs53 =             clhs19*clhs52;
const double clhs54 =             clhs21*clhs53;
const double clhs55 =             2*clhs44;
const double clhs56 =             clhs21*clhs55;
const double clhs57 =             2*N[0]*U(0,1) + 2*N[1]*U(1,1) + 2*N[2]*U(2,1);
const double clhs58 =             y - 3;
const double clhs59 =             clhs23*clhs58;
const double clhs60 =             y - 1;
const double clhs61 =             clhs21*clhs60;
const double clhs62 =             clhs42*clhs61;
const double clhs63 =             4*clhs33;
const double clhs64 =             clhs19*clhs63;
const double clhs65 =             pow(clhs19, 2);
const double clhs66 =             clhs60*clhs65;
const double clhs67 =             pow(clhs21, 2);
const double clhs68 =             clhs60*clhs67;
const double clhs69 =             clhs66 + clhs68;
const double clhs70 =             -2*clhs65 + clhs69;
const double clhs71 =             2*clhs30;
const double clhs72 =             DN(0,0)*clhs70 + clhs21*clhs64 - clhs25*clhs57 - clhs54 - clhs56 + clhs57*clhs59 + clhs62 - clhs70*clhs71;
const double clhs73 =             DN(0,0)*clhs28*mu;
const double clhs74 =             DN(0,0) - clhs30;
const double clhs75 =             clhs73*clhs74;
const double clhs76 =             (1.0L/3.0L)*DN(0,1)*clhs28*mu;
const double clhs77 =             3*clhs33;
const double clhs78 =             clhs18 + clhs20 - clhs46 - clhs77;
const double clhs79 =             DN(0,1)*clhs21;
const double clhs80 =             clhs25 - clhs34 + clhs79;
const double clhs81 =             DN(0,0)*clhs19;
const double clhs82 =             clhs31*clhs58 - clhs58*clhs81 - clhs59 + clhs80;
const double clhs83 =             DN(0,0)*mu;
const double clhs84 =             DN(0,1) - clhs33;
const double clhs85 =             clhs83*clhs84;
const double clhs86 =             (1.0L/3.0L)*DN(0,1)*mu;
const double clhs87 =             -clhs20 + clhs37 - clhs71 + clhs77;
const double clhs88 =             DN(0,1)*clhs19;
const double clhs89 =             DN(0,0)*clhs21;
const double clhs90 =             clhs41*clhs60;
const double clhs91 =             -clhs30*clhs61 - clhs44 + clhs48 + clhs60*clhs89 - clhs88 + clhs90;
const double clhs92 =             N[0]*N[1];
const double clhs93 =             clhs13*clhs92;
const double clhs94 =             4*DN(1,0);
const double clhs95 =             3*DN(1,1);
const double clhs96 =             N[1]*clhs22;
const double clhs97 =             N[1]*clhs24;
const double clhs98 =             3*clhs97;
const double clhs99 =             N[1]*clhs28;
const double clhs100 =             clhs27*clhs99;
const double clhs101 =             clhs100*clhs19;
const double clhs102 =             clhs32*clhs99;
const double clhs103 =             clhs102*clhs21;
const double clhs104 =             -8*clhs101 - 6*clhs103 + clhs19*clhs94 + clhs21*clhs95 + 4*clhs96 + clhs98;
const double clhs105 =             2*DN(1,0);
const double clhs106 =             clhs105*clhs21;
const double clhs107 =             clhs19*clhs95;
const double clhs108 =             N[1]*clhs40;
const double clhs109 =             2*clhs108;
const double clhs110 =             N[1]*clhs43;
const double clhs111 =             3*clhs110;
const double clhs112 =             4*clhs100;
const double clhs113 =             clhs112*clhs21;
const double clhs114 =             clhs102*clhs19;
const double clhs115 =             6*clhs114;
const double clhs116 =             clhs104 - clhs106 + clhs107 - clhs109 + clhs111 + clhs113 - clhs115;
const double clhs117 =             2*DN(1,1);
const double clhs118 =             clhs117*clhs19;
const double clhs119 =             clhs118*clhs21;
const double clhs120 =             2*clhs110;
const double clhs121 =             clhs120*clhs21;
const double clhs122 =             clhs58*clhs96;
const double clhs123 =             clhs109*clhs61;
const double clhs124 =             4*clhs102;
const double clhs125 =             clhs124*clhs19;
const double clhs126 =             2*clhs100;
const double clhs127 =             DN(1,0)*clhs70 - clhs119 - clhs121 + clhs122*clhs57 + clhs123 + clhs125*clhs21 - clhs126*clhs70 - clhs57*clhs97;
const double clhs128 =             DN(1,0) - clhs100;
const double clhs129 =             clhs128*clhs73;
const double clhs130 =             3*clhs102;
const double clhs131 =             -clhs112 - clhs130 + clhs94 + clhs95;
const double clhs132 =             DN(1,1)*clhs21;
const double clhs133 =             -clhs103 + clhs132 + clhs97;
const double clhs134 =             DN(1,0)*clhs19;
const double clhs135 =             clhs101*clhs58 - clhs122 + clhs133 - clhs134*clhs58;
const double clhs136 =             DN(1,1) - clhs102;
const double clhs137 =             clhs136*clhs83;
const double clhs138 =             clhs105 - clhs126 + clhs130 - clhs95;
const double clhs139 =             DN(1,1)*clhs19;
const double clhs140 =             DN(1,0)*clhs21;
const double clhs141 =             clhs108*clhs60;
const double clhs142 =             -clhs100*clhs61 - clhs110 + clhs114 - clhs139 + clhs140*clhs60 + clhs141;
const double clhs143 =             N[0]*N[2];
const double clhs144 =             clhs13*clhs143;
const double clhs145 =             4*DN(2,0);
const double clhs146 =             3*DN(2,1);
const double clhs147 =             N[2]*clhs22;
const double clhs148 =             N[2]*clhs24;
const double clhs149 =             3*clhs148;
const double clhs150 =             N[2]*clhs28;
const double clhs151 =             clhs150*clhs27;
const double clhs152 =             clhs151*clhs19;
const double clhs153 =             clhs150*clhs32;
const double clhs154 =             clhs153*clhs21;
const double clhs155 =             clhs145*clhs19 + clhs146*clhs21 + 4*clhs147 + clhs149 - 8*clhs152 - 6*clhs154;
const double clhs156 =             2*DN(2,0);
const double clhs157 =             clhs156*clhs21;
const double clhs158 =             clhs146*clhs19;
const double clhs159 =             N[2]*clhs40;
const double clhs160 =             2*clhs159;
const double clhs161 =             N[2]*clhs43;
const double clhs162 =             3*clhs161;
const double clhs163 =             4*clhs151;
const double clhs164 =             clhs163*clhs21;
const double clhs165 =             clhs153*clhs19;
const double clhs166 =             6*clhs165;
const double clhs167 =             clhs155 - clhs157 + clhs158 - clhs160 + clhs162 + clhs164 - clhs166;
const double clhs168 =             2*DN(2,1);
const double clhs169 =             clhs168*clhs19;
const double clhs170 =             clhs169*clhs21;
const double clhs171 =             2*clhs161;
const double clhs172 =             clhs171*clhs21;
const double clhs173 =             clhs147*clhs58;
const double clhs174 =             clhs160*clhs61;
const double clhs175 =             4*clhs153;
const double clhs176 =             clhs175*clhs19;
const double clhs177 =             2*clhs151;
const double clhs178 =             DN(2,0)*clhs70 - clhs148*clhs57 - clhs170 - clhs172 + clhs173*clhs57 + clhs174 + clhs176*clhs21 - clhs177*clhs70;
const double clhs179 =             DN(2,0) - clhs151;
const double clhs180 =             clhs179*clhs73;
const double clhs181 =             3*clhs153;
const double clhs182 =             clhs145 + clhs146 - clhs163 - clhs181;
const double clhs183 =             DN(2,1)*clhs21;
const double clhs184 =             clhs148 - clhs154 + clhs183;
const double clhs185 =             DN(2,0)*clhs19;
const double clhs186 =             clhs152*clhs58 - clhs173 + clhs184 - clhs185*clhs58;
const double clhs187 =             DN(2,1) - clhs153;
const double clhs188 =             clhs187*clhs83;
const double clhs189 =             -clhs146 + clhs156 - clhs177 + clhs181;
const double clhs190 =             DN(2,1)*clhs19;
const double clhs191 =             DN(2,0)*clhs21;
const double clhs192 =             clhs159*clhs60;
const double clhs193 =             -clhs151*clhs61 - clhs161 + clhs165 - clhs190 + clhs191*clhs60 + clhs192;
const double clhs194 =             N[0]*f_ext(0,1) + N[1]*f_ext(1,1) + N[2]*f_ext(2,1);
const double clhs195 =             clhs0*clhs194;
const double clhs196 =             3*DN(0,0);
const double clhs197 =             clhs196*clhs21;
const double clhs198 =             3*clhs41;
const double clhs199 =             clhs21*clhs30;
const double clhs200 =             6*clhs199;
const double clhs201 =             clhs197 + clhs198 - clhs200 - clhs53 - clhs55 + clhs64;
const double clhs202 =             4*DN(0,1);
const double clhs203 =             3*clhs23;
const double clhs204 =             clhs19*clhs196 + clhs201 + clhs202*clhs21 + clhs203 + 4*clhs25 - 6*clhs31 - 8*clhs34;
const double clhs205 =             clhs19*clhs38;
const double clhs206 =             2*N[0]*U(0,2) + 2*N[1]*U(1,2) + 2*N[2]*U(2,2);
const double clhs207 =             clhs19*clhs42;
const double clhs208 =             clhs25*clhs58;
const double clhs209 =             clhs19*clhs60;
const double clhs210 =             clhs209*clhs55;
const double clhs211 =             -2*clhs67 + clhs69;
const double clhs212 =             2*clhs33;
const double clhs213 =             DN(0,1)*clhs211 + clhs19*clhs47 - clhs205 + clhs206*clhs208 - clhs206*clhs23 - clhs207 + clhs210 - clhs211*clhs212;
const double clhs214 =             3*clhs30;
const double clhs215 =             clhs196 - clhs214;
const double clhs216 =             clhs212 + clhs215 - clhs52;
const double clhs217 =             clhs44*clhs60;
const double clhs218 =             -clhs199 + clhs209*clhs33 - clhs217 + clhs41 - clhs60*clhs88 + clhs89;
const double clhs219 =             clhs202 + clhs215 - clhs63;
const double clhs220 =             clhs23 - clhs31 + clhs81;
const double clhs221 =             -clhs208 + clhs220 + clhs34*clhs58 - clhs58*clhs79;
const double clhs222 =             clhs194*clhs92;
const double clhs223 =             3*DN(1,0);
const double clhs224 =             clhs21*clhs223;
const double clhs225 =             3*clhs108;
const double clhs226 =             clhs100*clhs21;
const double clhs227 =             6*clhs226;
const double clhs228 =             -clhs118 - clhs120 + clhs125 + clhs224 + clhs225 - clhs227;
const double clhs229 =             4*DN(1,1);
const double clhs230 =             3*clhs96;
const double clhs231 =             -6*clhs101 - 8*clhs103 + clhs19*clhs223 + clhs21*clhs229 + clhs228 + clhs230 + 4*clhs97;
const double clhs232 =             clhs106*clhs19;
const double clhs233 =             clhs109*clhs19;
const double clhs234 =             clhs58*clhs97;
const double clhs235 =             clhs120*clhs209;
const double clhs236 =             2*clhs102;
const double clhs237 =             DN(1,1)*clhs211 + clhs113*clhs19 + clhs206*clhs234 - clhs206*clhs96 - clhs211*clhs236 - clhs232 - clhs233 + clhs235;
const double clhs238 =             3*clhs100;
const double clhs239 =             clhs223 - clhs238;
const double clhs240 =             -clhs117 + clhs236 + clhs239;
const double clhs241 =             clhs110*clhs60;
const double clhs242 =             clhs102*clhs209 + clhs108 - clhs139*clhs60 + clhs140 - clhs226 - clhs241;
const double clhs243 =             -clhs124 + clhs229 + clhs239;
const double clhs244 =             -clhs101 + clhs134 + clhs96;
const double clhs245 =             clhs103*clhs58 - clhs132*clhs58 - clhs234 + clhs244;
const double clhs246 =             clhs143*clhs194;
const double clhs247 =             3*DN(2,0);
const double clhs248 =             clhs21*clhs247;
const double clhs249 =             3*clhs159;
const double clhs250 =             clhs151*clhs21;
const double clhs251 =             6*clhs250;
const double clhs252 =             -clhs169 - clhs171 + clhs176 + clhs248 + clhs249 - clhs251;
const double clhs253 =             4*DN(2,1);
const double clhs254 =             3*clhs147;
const double clhs255 =             4*clhs148 - 6*clhs152 - 8*clhs154 + clhs19*clhs247 + clhs21*clhs253 + clhs252 + clhs254;
const double clhs256 =             clhs157*clhs19;
const double clhs257 =             clhs160*clhs19;
const double clhs258 =             clhs148*clhs58;
const double clhs259 =             clhs171*clhs209;
const double clhs260 =             2*clhs153;
const double clhs261 =             DN(2,1)*clhs211 - clhs147*clhs206 + clhs164*clhs19 + clhs206*clhs258 - clhs211*clhs260 - clhs256 - clhs257 + clhs259;
const double clhs262 =             3*clhs151;
const double clhs263 =             clhs247 - clhs262;
const double clhs264 =             -clhs168 + clhs260 + clhs263;
const double clhs265 =             clhs161*clhs60;
const double clhs266 =             clhs153*clhs209 + clhs159 - clhs190*clhs60 + clhs191 - clhs250 - clhs265;
const double clhs267 =             -clhs175 + clhs253 + clhs263;
const double clhs268 =             clhs147 - clhs152 + clhs185;
const double clhs269 =             clhs154*clhs58 - clhs183*clhs58 - clhs258 + clhs268;
const double clhs270 =             N[0]*r[0] + N[1]*r[1] + N[2]*r[2];
const double clhs271 =             (1.0L/3.0L)*DN(0,0)*clhs16;
const double clhs272 =             1.0/cv;
const double clhs273 =             3*N[0]*clhs272*lambda;
const double clhs274 =             DN(0,0)*U(0,3) + DN(1,0)*U(1,3) + DN(2,0)*U(2,3);
const double clhs275 =             clhs21*clhs28*mu;
const double clhs276 =             4*clhs21*clhs28*mu;
const double clhs277 =             6*clhs19*clhs28*mu;
const double clhs278 =             clhs272*lambda;
const double clhs279 =             -clhs278 + mu;
const double clhs280 =             6*clhs21*clhs279*clhs28;
const double clhs281 =             4*mu;
const double clhs282 =             3*clhs278;
const double clhs283 =             clhs281 - clhs282;
const double clhs284 =             clhs23*clhs283;
const double clhs285 =             clhs28*clhs57;
const double clhs286 =             3*N[0]*clhs16*clhs19*clhs21*mu;
const double clhs287 =             N[0]*U(0,3);
const double clhs288 =             N[1]*U(1,3);
const double clhs289 =             N[2]*U(2,3);
const double clhs290 =             clhs282*(clhs287 + clhs288 + clhs289);
const double clhs291 =             clhs28*clhs65;
const double clhs292 =             3*mu;
const double clhs293 =             clhs28*clhs67;
const double clhs294 =             -clhs282*clhs291;
const double clhs295 =             -clhs282*clhs293;
const double clhs296 =             clhs281*clhs291 + clhs292*clhs293 + clhs294 + clhs295;
const double clhs297 =             clhs290 + clhs296;
const double clhs298 =             2*clhs287;
const double clhs299 =             2*clhs288;
const double clhs300 =             2*clhs289;
const double clhs301 =             clhs298 + clhs299 + clhs300;
const double clhs302 =             clhs278*clhs301;
const double clhs303 =             clhs296 + clhs302;
const double clhs304 =             -DN(0,0)*clhs297 + clhs214*clhs303 - clhs25*clhs277 - clhs273*clhs274 - clhs275*clhs88 + clhs276*clhs44 - clhs280*clhs41 - clhs284*clhs285 + clhs286*clhs32;
const double clhs305 =             2*N[0]*y;
const double clhs306 =             clhs274*clhs305;
const double clhs307 =             DN(0,1)*U(0,3) + DN(1,1)*U(1,3) + DN(2,1)*U(2,3);
const double clhs308 =             clhs305*clhs307;
const double clhs309 =             4*N[0]*U(0,1) + 4*N[1]*U(1,1) + 4*N[2]*U(2,1);
const double clhs310 =             clhs21*clhs28*clhs309;
const double clhs311 =             clhs28*clhs68;
const double clhs312 =             -clhs298;
const double clhs313 =             -clhs299;
const double clhs314 =             -clhs300;
const double clhs315 =             clhs301*clhs60;
const double clhs316 =             -clhs315;
const double clhs317 =             clhs28*clhs66;
const double clhs318 =             clhs28*clhs60;
const double clhs319 =             clhs318*(clhs65 + clhs67);
const double clhs320 =             clhs312 + clhs313 + clhs314 + clhs316 + clhs317 + clhs319;
const double clhs321 =             clhs311 + clhs320;
const double clhs322 =             clhs311 + clhs312 + clhs313 + clhs314 + clhs316;
const double clhs323 =             5*clhs317 + clhs319 + clhs322;
const double clhs324 =             5*clhs311 + clhs320;
const double clhs325 =             clhs317 + 2*clhs319 + clhs322;
const double clhs326 =             clhs325*clhs57;
const double clhs327 =             clhs206*clhs325;
const double clhs328 =             -clhs19*clhs306 - clhs21*clhs308 + clhs217*clhs310 + clhs23*clhs323 + clhs25*clhs324 - clhs30*clhs326 + clhs310*clhs90 + clhs321*clhs79 + clhs321*clhs81 - clhs327*clhs33;
const double clhs329 =             (1.0L/3.0L)*DN(0,1)*clhs16;
const double clhs330 =             6*clhs21*clhs28*mu;
const double clhs331 =             4*clhs19*clhs28*mu;
const double clhs332 =             6*clhs19*clhs279*clhs28;
const double clhs333 =             clhs25*clhs283;
const double clhs334 =             clhs206*clhs28;
const double clhs335 =             clhs281*clhs293 + clhs291*clhs292 + clhs294 + clhs295;
const double clhs336 =             clhs290 + clhs335;
const double clhs337 =             clhs302 + clhs335;
const double clhs338 =             -DN(0,1)*clhs336 - clhs23*clhs330 + clhs27*clhs286 - clhs273*clhs307 - clhs275*clhs81 + clhs304 + clhs331*clhs41 - clhs332*clhs44 - clhs333*clhs334 + clhs337*clhs77;
const double clhs339 =             clhs21*mu;
const double clhs340 =             clhs283*clhs57;
const double clhs341 =             clhs26*mu + clhs283*clhs81 + clhs284 - clhs30*clhs340 - clhs33*clhs339 - clhs339*clhs52;
const double clhs342 =             clhs197*mu + clhs279*clhs39 + clhs279*clhs45 - clhs279*clhs49 - clhs30*clhs339 + clhs341 - clhs42*mu;
const double clhs343 =             (1.0L/2.0L)*N[0]*clhs28;
const double clhs344 =             6*clhs19*clhs28*clhs60;
const double clhs345 =             N[0]*clhs16*clhs21*clhs309*clhs60;
const double clhs346 =             clhs298 + clhs299 + clhs300 + clhs315;
const double clhs347 =             -clhs311 - 3*clhs317 + clhs346;
const double clhs348 =             -DN(0,0)*clhs347 + N[0]*clhs24*clhs318*clhs57 + clhs23*clhs344 + clhs28*clhs62 - clhs30*clhs323 - clhs306 + clhs318*clhs54 + clhs318*clhs56 - clhs32*clhs345;
const double clhs349 =             clhs19*mu;
const double clhs350 =             clhs197*clhs279 + clhs198*clhs279 - clhs200*clhs279 - clhs33*clhs349 + clhs39*mu - clhs55*mu;
const double clhs351 =             clhs206*clhs283;
const double clhs352 =             clhs203*mu + clhs283*clhs79 - clhs30*clhs349 - clhs33*clhs351 + clhs333 - clhs349*clhs37 + clhs350;
const double clhs353 =             6*clhs21*clhs28*clhs60;
const double clhs354 =             -3*clhs311 - clhs317 + clhs346;
const double clhs355 =             -DN(0,1)*clhs354 + N[0]*clhs206*clhs22*clhs318 + clhs205*clhs318 + clhs207*clhs318 + clhs210*clhs28 + clhs25*clhs353 - clhs27*clhs345 - clhs308 - clhs324*clhs33;
const double clhs356 =             DN(0,0)*clhs272*clhs28*lambda;
const double clhs357 =             DN(0,1)*clhs272*clhs28*lambda;
const double clhs358 =             clhs74 + clhs84;
const double clhs359 =             N[0]*clhs28*y;
const double clhs360 =             clhs220 + clhs80;
const double clhs361 =             clhs270*clhs92;
const double clhs362 =             3*N[1]*clhs272*lambda;
const double clhs363 =             clhs283*clhs96;
const double clhs364 =             3*N[1]*clhs16*clhs19*clhs21*mu;
const double clhs365 =             -DN(1,0)*clhs297 - clhs108*clhs280 + clhs110*clhs276 - clhs139*clhs275 + clhs238*clhs303 - clhs274*clhs362 - clhs277*clhs97 - clhs285*clhs363 + clhs32*clhs364;
const double clhs366 =             2*N[1]*y;
const double clhs367 =             clhs274*clhs366;
const double clhs368 =             clhs307*clhs366;
const double clhs369 =             -clhs100*clhs326 - clhs102*clhs327 + clhs132*clhs321 + clhs134*clhs321 + clhs141*clhs310 - clhs19*clhs367 - clhs21*clhs368 + clhs241*clhs310 + clhs323*clhs96 + clhs324*clhs97;
const double clhs370 =             clhs283*clhs97;
const double clhs371 =             -DN(1,1)*clhs336 + clhs108*clhs331 - clhs110*clhs332 + clhs130*clhs337 - clhs134*clhs275 + clhs27*clhs364 - clhs307*clhs362 - clhs330*clhs96 - clhs334*clhs370 + clhs365;
const double clhs372 =             -clhs100*clhs340 - clhs102*clhs339 - clhs117*clhs339 + clhs134*clhs283 + clhs363 + clhs98*mu;
const double clhs373 =             -clhs100*clhs339 + clhs107*clhs279 - clhs109*mu + clhs111*clhs279 - clhs115*clhs279 + clhs224*mu + clhs372;
const double clhs374 =             N[1]*clhs16*clhs21*clhs309*clhs60;
const double clhs375 =             -DN(1,0)*clhs347 + N[1]*clhs24*clhs318*clhs57 - clhs100*clhs323 + clhs119*clhs318 + clhs121*clhs318 + clhs123*clhs28 - clhs32*clhs374 + clhs344*clhs96 - clhs367;
const double clhs376 =             -clhs102*clhs349 + clhs107*mu - clhs120*mu + clhs224*clhs279 + clhs225*clhs279 - clhs227*clhs279;
const double clhs377 =             -clhs100*clhs349 - clhs102*clhs351 - clhs105*clhs349 + clhs132*clhs283 + clhs230*mu + clhs370 + clhs376;
const double clhs378 =             -DN(1,1)*clhs354 + N[1]*clhs206*clhs22*clhs318 - clhs102*clhs324 + clhs232*clhs318 + clhs233*clhs318 + clhs235*clhs28 - clhs27*clhs374 + clhs353*clhs97 - clhs368;
const double clhs379 =             clhs128 + clhs136;
const double clhs380 =             clhs133 + clhs244;
const double clhs381 =             clhs143*clhs270;
const double clhs382 =             3*N[2]*clhs272*lambda;
const double clhs383 =             clhs147*clhs283;
const double clhs384 =             3*N[2]*clhs16*clhs19*clhs21*mu;
const double clhs385 =             -DN(2,0)*clhs297 - clhs148*clhs277 - clhs159*clhs280 + clhs161*clhs276 - clhs190*clhs275 + clhs262*clhs303 - clhs274*clhs382 - clhs285*clhs383 + clhs32*clhs384;
const double clhs386 =             2*N[2]*y;
const double clhs387 =             clhs274*clhs386;
const double clhs388 =             clhs307*clhs386;
const double clhs389 =             clhs147*clhs323 + clhs148*clhs324 - clhs151*clhs326 - clhs153*clhs327 + clhs183*clhs321 + clhs185*clhs321 - clhs19*clhs387 + clhs192*clhs310 - clhs21*clhs388 + clhs265*clhs310;
const double clhs390 =             clhs148*clhs283;
const double clhs391 =             -DN(2,1)*clhs336 - clhs147*clhs330 + clhs159*clhs331 - clhs161*clhs332 + clhs181*clhs337 - clhs185*clhs275 + clhs27*clhs384 - clhs307*clhs382 - clhs334*clhs390 + clhs385;
const double clhs392 =             clhs149*mu - clhs151*clhs340 - clhs153*clhs339 - clhs168*clhs339 + clhs185*clhs283 + clhs383;
const double clhs393 =             -clhs151*clhs339 + clhs158*clhs279 - clhs160*mu + clhs162*clhs279 - clhs166*clhs279 + clhs248*mu + clhs392;
const double clhs394 =             N[2]*clhs16*clhs21*clhs309*clhs60;
const double clhs395 =             -DN(2,0)*clhs347 + N[2]*clhs24*clhs318*clhs57 + clhs147*clhs344 - clhs151*clhs323 + clhs170*clhs318 + clhs172*clhs318 + clhs174*clhs28 - clhs32*clhs394 - clhs387;
const double clhs396 =             -clhs153*clhs349 + clhs158*mu - clhs171*mu + clhs248*clhs279 + clhs249*clhs279 - clhs251*clhs279;
const double clhs397 =             -clhs151*clhs349 - clhs153*clhs351 - clhs156*clhs349 + clhs183*clhs283 + clhs254*mu + clhs390 + clhs396;
const double clhs398 =             -DN(2,1)*clhs354 + N[2]*clhs206*clhs22*clhs318 + clhs148*clhs353 - clhs153*clhs324 + clhs256*clhs318 + clhs257*clhs318 + clhs259*clhs28 - clhs27*clhs394 - clhs388;
const double clhs399 =             clhs179 + clhs187;
const double clhs400 =             clhs184 + clhs268;
const double clhs401 =             DN(0,0)*N[1];
const double clhs402 =             DN(0,1)*N[1];
const double clhs403 =             pow(N[1], 2);
const double clhs404 =             bdf0*clhs403;
const double clhs405 =             DN(1,0)*N[1];
const double clhs406 =             DN(1,1)*N[1];
const double clhs407 =             N[1]*N[2];
const double clhs408 =             bdf0*clhs407;
const double clhs409 =             -clhs408;
const double clhs410 =             DN(2,0)*N[1];
const double clhs411 =             DN(2,1)*N[1];
const double clhs412 =             (1.0L/3.0L)*DN(1,0)*clhs16*mu;
const double clhs413 =             (1.0L/3.0L)*DN(1,1)*clhs16*mu;
const double clhs414 =             (1.0L/2.0L)*N[1]*clhs16;
const double clhs415 =             DN(1,0)*clhs28*mu;
const double clhs416 =             clhs415*clhs74;
const double clhs417 =             (1.0L/3.0L)*DN(1,1)*clhs28*mu;
const double clhs418 =             DN(1,0)*mu;
const double clhs419 =             clhs418*clhs84;
const double clhs420 =             (1.0L/3.0L)*DN(1,1)*mu;
const double clhs421 =             clhs13*clhs403;
const double clhs422 =             clhs128*clhs415;
const double clhs423 =             clhs136*clhs418;
const double clhs424 =             clhs13*clhs407;
const double clhs425 =             clhs179*clhs415;
const double clhs426 =             clhs187*clhs418;
const double clhs427 =             clhs194*clhs403;
const double clhs428 =             clhs194*clhs407;
const double clhs429 =             (1.0L/3.0L)*DN(1,0)*clhs16;
const double clhs430 =             (1.0L/3.0L)*DN(1,1)*clhs16;
const double clhs431 =             (1.0L/2.0L)*N[1]*clhs28;
const double clhs432 =             DN(1,0)*clhs272*clhs28*lambda;
const double clhs433 =             DN(1,1)*clhs272*clhs28*lambda;
const double clhs434 =             N[1]*clhs28*y;
const double clhs435 =             clhs270*clhs407;
const double clhs436 =             DN(0,0)*N[2];
const double clhs437 =             DN(0,1)*N[2];
const double clhs438 =             DN(1,0)*N[2];
const double clhs439 =             DN(1,1)*N[2];
const double clhs440 =             pow(N[2], 2);
const double clhs441 =             bdf0*clhs440;
const double clhs442 =             DN(2,0)*N[2];
const double clhs443 =             DN(2,1)*N[2];
const double clhs444 =             (1.0L/3.0L)*DN(2,0)*clhs16*mu;
const double clhs445 =             (1.0L/3.0L)*DN(2,1)*clhs16*mu;
const double clhs446 =             (1.0L/2.0L)*N[2]*clhs16;
const double clhs447 =             DN(2,0)*clhs28*mu;
const double clhs448 =             clhs447*clhs74;
const double clhs449 =             (1.0L/3.0L)*DN(2,1)*clhs28*mu;
const double clhs450 =             DN(2,0)*mu;
const double clhs451 =             clhs450*clhs84;
const double clhs452 =             (1.0L/3.0L)*DN(2,1)*mu;
const double clhs453 =             clhs128*clhs447;
const double clhs454 =             clhs136*clhs450;
const double clhs455 =             clhs13*clhs440;
const double clhs456 =             clhs179*clhs447;
const double clhs457 =             clhs187*clhs450;
const double clhs458 =             clhs194*clhs440;
const double clhs459 =             (1.0L/3.0L)*DN(2,0)*clhs16;
const double clhs460 =             (1.0L/3.0L)*DN(2,1)*clhs16;
const double clhs461 =             (1.0L/2.0L)*N[2]*clhs28;
const double clhs462 =             DN(2,0)*clhs272*clhs28*lambda;
const double clhs463 =             DN(2,1)*clhs272*clhs28*lambda;
const double clhs464 =             N[2]*clhs28*y;
            lhs(0,0)=-clhs1;
            lhs(0,1)=-clhs2;
            lhs(0,2)=-clhs3;
            lhs(0,3)=0;
            lhs(0,4)=clhs6;
            lhs(0,5)=-clhs7;
            lhs(0,6)=-clhs8;
            lhs(0,7)=0;
            lhs(0,8)=clhs10;
            lhs(0,9)=-clhs11;
            lhs(0,10)=-clhs12;
            lhs(0,11)=0;
            lhs(1,0)=clhs14 + clhs17*clhs35 + clhs36*clhs50 - clhs51*clhs72;
            lhs(1,1)=-clhs1 - clhs29*clhs82 - 4.0L/3.0L*clhs75 - clhs76*clhs78;
            lhs(1,2)=clhs28*(N[0]*clhs91 - clhs85 + clhs86*clhs87);
            lhs(1,3)=-clhs2*clhs60;
            lhs(1,4)=clhs104*clhs17 + clhs116*clhs36 - clhs127*clhs51 + clhs93;
            lhs(1,5)=-4.0L/3.0L*clhs129 - clhs131*clhs76 - clhs135*clhs29 - clhs5;
            lhs(1,6)=clhs28*(N[0]*clhs142 - clhs137 + clhs138*clhs86);
            lhs(1,7)=-clhs60*clhs7;
            lhs(1,8)=clhs144 + clhs155*clhs17 + clhs167*clhs36 - clhs178*clhs51;
            lhs(1,9)=-4.0L/3.0L*clhs180 - clhs182*clhs76 - clhs186*clhs29 - clhs9;
            lhs(1,10)=clhs28*(N[0]*clhs193 - clhs188 + clhs189*clhs86);
            lhs(1,11)=-clhs11*clhs60;
            lhs(2,0)=clhs17*clhs201 + clhs195 + clhs204*clhs36 - clhs213*clhs51;
            lhs(2,1)=clhs28*(-N[0]*clhs218 - clhs216*clhs86 + (2.0L/3.0L)*clhs85);
            lhs(2,2)=-clhs1 - clhs219*clhs76 - clhs221*clhs29 - clhs75;
            lhs(2,3)=-clhs3*clhs60;
            lhs(2,4)=clhs17*clhs228 + clhs222 + clhs231*clhs36 - clhs237*clhs51;
            lhs(2,5)=clhs28*(-N[0]*clhs242 + (2.0L/3.0L)*clhs137 - clhs240*clhs86);
            lhs(2,6)=-clhs129 - clhs243*clhs76 - clhs245*clhs29 - clhs5;
            lhs(2,7)=-clhs60*clhs8;
            lhs(2,8)=clhs17*clhs252 + clhs246 + clhs255*clhs36 - clhs261*clhs51;
            lhs(2,9)=clhs28*(-N[0]*clhs266 + (2.0L/3.0L)*clhs188 - clhs264*clhs86);
            lhs(2,10)=-clhs180 - clhs267*clhs76 - clhs269*clhs29 - clhs9;
            lhs(2,11)=-clhs12*clhs60;
            lhs(3,0)=clhs0*clhs270 - clhs271*clhs304 - clhs328*clhs51 - clhs329*clhs338;
            lhs(3,1)=clhs14 - clhs271*clhs341 - clhs329*clhs342 + clhs343*clhs348;
            lhs(3,2)=clhs195 - clhs271*clhs350 - clhs329*clhs352 + clhs343*clhs355;
            lhs(3,3)=-clhs1 - clhs356*clhs74 - clhs357*clhs358 - clhs359*clhs360;
            lhs(3,4)=-clhs271*clhs365 - clhs329*clhs371 + clhs361 - clhs369*clhs51;
            lhs(3,5)=-clhs271*clhs372 - clhs329*clhs373 + clhs343*clhs375 + clhs93;
            lhs(3,6)=clhs222 - clhs271*clhs376 - clhs329*clhs377 + clhs343*clhs378;
            lhs(3,7)=-clhs128*clhs356 - clhs357*clhs379 - clhs359*clhs380 - clhs5;
            lhs(3,8)=-clhs271*clhs385 - clhs329*clhs391 + clhs381 - clhs389*clhs51;
            lhs(3,9)=clhs144 - clhs271*clhs392 - clhs329*clhs393 + clhs343*clhs395;
            lhs(3,10)=clhs246 - clhs271*clhs396 - clhs329*clhs397 + clhs343*clhs398;
            lhs(3,11)=-clhs179*clhs356 - clhs357*clhs399 - clhs359*clhs400 - clhs9;
            lhs(4,0)=clhs6;
            lhs(4,1)=-clhs401;
            lhs(4,2)=-clhs402;
            lhs(4,3)=0;
            lhs(4,4)=-clhs404;
            lhs(4,5)=-clhs405;
            lhs(4,6)=-clhs406;
            lhs(4,7)=0;
            lhs(4,8)=clhs409;
            lhs(4,9)=-clhs410;
            lhs(4,10)=-clhs411;
            lhs(4,11)=0;
            lhs(5,0)=clhs35*clhs412 + clhs413*clhs50 - clhs414*clhs72 + clhs93;
            lhs(5,1)=-4.0L/3.0L*clhs416 - clhs417*clhs78 - clhs5 - clhs82*clhs99;
            lhs(5,2)=clhs28*(N[1]*clhs91 - clhs419 + clhs420*clhs87);
            lhs(5,3)=-clhs401*clhs60;
            lhs(5,4)=clhs104*clhs412 + clhs116*clhs413 - clhs127*clhs414 + clhs421;
            lhs(5,5)=-clhs131*clhs417 - clhs135*clhs99 - clhs404 - 4.0L/3.0L*clhs422;
            lhs(5,6)=clhs28*(N[1]*clhs142 + clhs138*clhs420 - clhs423);
            lhs(5,7)=-clhs405*clhs60;
            lhs(5,8)=clhs155*clhs412 + clhs167*clhs413 - clhs178*clhs414 + clhs424;
            lhs(5,9)=-clhs182*clhs417 - clhs186*clhs99 - clhs408 - 4.0L/3.0L*clhs425;
            lhs(5,10)=clhs28*(N[1]*clhs193 + clhs189*clhs420 - clhs426);
            lhs(5,11)=-clhs410*clhs60;
            lhs(6,0)=clhs201*clhs412 + clhs204*clhs413 - clhs213*clhs414 + clhs222;
            lhs(6,1)=clhs28*(-N[1]*clhs218 - clhs216*clhs420 + (2.0L/3.0L)*clhs419);
            lhs(6,2)=-clhs219*clhs417 - clhs221*clhs99 - clhs416 - clhs5;
            lhs(6,3)=-clhs402*clhs60;
            lhs(6,4)=clhs228*clhs412 + clhs231*clhs413 - clhs237*clhs414 + clhs427;
            lhs(6,5)=clhs28*(-N[1]*clhs242 - clhs240*clhs420 + (2.0L/3.0L)*clhs423);
            lhs(6,6)=-clhs243*clhs417 - clhs245*clhs99 - clhs404 - clhs422;
            lhs(6,7)=-clhs406*clhs60;
            lhs(6,8)=clhs252*clhs412 + clhs255*clhs413 - clhs261*clhs414 + clhs428;
            lhs(6,9)=clhs28*(-N[1]*clhs266 - clhs264*clhs420 + (2.0L/3.0L)*clhs426);
            lhs(6,10)=-clhs267*clhs417 - clhs269*clhs99 - clhs408 - clhs425;
            lhs(6,11)=-clhs411*clhs60;
            lhs(7,0)=-clhs304*clhs429 - clhs328*clhs414 - clhs338*clhs430 + clhs361;
            lhs(7,1)=-clhs341*clhs429 - clhs342*clhs430 + clhs348*clhs431 + clhs93;
            lhs(7,2)=clhs222 - clhs350*clhs429 - clhs352*clhs430 + clhs355*clhs431;
            lhs(7,3)=-clhs358*clhs433 - clhs360*clhs434 - clhs432*clhs74 - clhs5;
            lhs(7,4)=clhs270*clhs403 - clhs365*clhs429 - clhs369*clhs414 - clhs371*clhs430;
            lhs(7,5)=-clhs372*clhs429 - clhs373*clhs430 + clhs375*clhs431 + clhs421;
            lhs(7,6)=-clhs376*clhs429 - clhs377*clhs430 + clhs378*clhs431 + clhs427;
            lhs(7,7)=-clhs128*clhs432 - clhs379*clhs433 - clhs380*clhs434 - clhs404;
            lhs(7,8)=-clhs385*clhs429 - clhs389*clhs414 - clhs391*clhs430 + clhs435;
            lhs(7,9)=-clhs392*clhs429 - clhs393*clhs430 + clhs395*clhs431 + clhs424;
            lhs(7,10)=-clhs396*clhs429 - clhs397*clhs430 + clhs398*clhs431 + clhs428;
            lhs(7,11)=-clhs179*clhs432 - clhs399*clhs433 - clhs400*clhs434 - clhs408;
            lhs(8,0)=clhs10;
            lhs(8,1)=-clhs436;
            lhs(8,2)=-clhs437;
            lhs(8,3)=0;
            lhs(8,4)=clhs409;
            lhs(8,5)=-clhs438;
            lhs(8,6)=-clhs439;
            lhs(8,7)=0;
            lhs(8,8)=-clhs441;
            lhs(8,9)=-clhs442;
            lhs(8,10)=-clhs443;
            lhs(8,11)=0;
            lhs(9,0)=clhs144 + clhs35*clhs444 + clhs445*clhs50 - clhs446*clhs72;
            lhs(9,1)=-clhs150*clhs82 - 4.0L/3.0L*clhs448 - clhs449*clhs78 - clhs9;
            lhs(9,2)=clhs28*(N[2]*clhs91 - clhs451 + clhs452*clhs87);
            lhs(9,3)=-clhs436*clhs60;
            lhs(9,4)=clhs104*clhs444 + clhs116*clhs445 - clhs127*clhs446 + clhs424;
            lhs(9,5)=-clhs131*clhs449 - clhs135*clhs150 - clhs408 - 4.0L/3.0L*clhs453;
            lhs(9,6)=clhs28*(N[2]*clhs142 + clhs138*clhs452 - clhs454);
            lhs(9,7)=-clhs438*clhs60;
            lhs(9,8)=clhs155*clhs444 + clhs167*clhs445 - clhs178*clhs446 + clhs455;
            lhs(9,9)=-clhs150*clhs186 - clhs182*clhs449 - clhs441 - 4.0L/3.0L*clhs456;
            lhs(9,10)=clhs28*(N[2]*clhs193 + clhs189*clhs452 - clhs457);
            lhs(9,11)=-clhs442*clhs60;
            lhs(10,0)=clhs201*clhs444 + clhs204*clhs445 - clhs213*clhs446 + clhs246;
            lhs(10,1)=clhs28*(-N[2]*clhs218 - clhs216*clhs452 + (2.0L/3.0L)*clhs451);
            lhs(10,2)=-clhs150*clhs221 - clhs219*clhs449 - clhs448 - clhs9;
            lhs(10,3)=-clhs437*clhs60;
            lhs(10,4)=clhs228*clhs444 + clhs231*clhs445 - clhs237*clhs446 + clhs428;
            lhs(10,5)=clhs28*(-N[2]*clhs242 - clhs240*clhs452 + (2.0L/3.0L)*clhs454);
            lhs(10,6)=-clhs150*clhs245 - clhs243*clhs449 - clhs408 - clhs453;
            lhs(10,7)=-clhs439*clhs60;
            lhs(10,8)=clhs252*clhs444 + clhs255*clhs445 - clhs261*clhs446 + clhs458;
            lhs(10,9)=clhs28*(-N[2]*clhs266 - clhs264*clhs452 + (2.0L/3.0L)*clhs457);
            lhs(10,10)=-clhs150*clhs269 - clhs267*clhs449 - clhs441 - clhs456;
            lhs(10,11)=-clhs443*clhs60;
            lhs(11,0)=-clhs304*clhs459 - clhs328*clhs446 - clhs338*clhs460 + clhs381;
            lhs(11,1)=clhs144 - clhs341*clhs459 - clhs342*clhs460 + clhs348*clhs461;
            lhs(11,2)=clhs246 - clhs350*clhs459 - clhs352*clhs460 + clhs355*clhs461;
            lhs(11,3)=-clhs358*clhs463 - clhs360*clhs464 - clhs462*clhs74 - clhs9;
            lhs(11,4)=-clhs365*clhs459 - clhs369*clhs446 - clhs371*clhs460 + clhs435;
            lhs(11,5)=-clhs372*clhs459 - clhs373*clhs460 + clhs375*clhs461 + clhs424;
            lhs(11,6)=-clhs376*clhs459 - clhs377*clhs460 + clhs378*clhs461 + clhs428;
            lhs(11,7)=-clhs128*clhs462 - clhs379*clhs463 - clhs380*clhs464 - clhs408;
            lhs(11,8)=clhs270*clhs440 - clhs385*clhs459 - clhs389*clhs446 - clhs391*clhs460;
            lhs(11,9)=-clhs392*clhs459 - clhs393*clhs460 + clhs395*clhs461 + clhs455;
            lhs(11,10)=-clhs396*clhs459 - clhs397*clhs460 + clhs398*clhs461 + clhs458;
            lhs(11,11)=-clhs179*clhs462 - clhs399*clhs463 - clhs400*clhs464 - clhs441;


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

    const double crhs0 =             DN(0,1)*U(0,2);
const double crhs1 =             DN(1,1)*U(1,2);
const double crhs2 =             DN(2,1)*U(2,2);
const double crhs3 =             crhs0 + crhs1 + crhs2;
const double crhs4 =             DN(0,0)*U(0,1);
const double crhs5 =             DN(1,0)*U(1,1);
const double crhs6 =             DN(2,0)*U(2,1);
const double crhs7 =             N[0]*(U(0,0)*bdf0 + Un(0,0)*bdf1 + Unn(0,0)*bdf2) + N[1]*(U(1,0)*bdf0 + Un(1,0)*bdf1 + Unn(1,0)*bdf2) + N[2]*(U(2,0)*bdf0 + Un(2,0)*bdf1 + Unn(2,0)*bdf2) + crhs3 + crhs4 + crhs5 + crhs6;
const double crhs8 =             N[0]*f_ext(0,0) + N[1]*f_ext(1,0) + N[2]*f_ext(2,0);
const double crhs9 =             N[0]*U(0,0) + N[1]*U(1,0) + N[2]*U(2,0);
const double crhs10 =             N[0]*crhs9;
const double crhs11 =             N[0]*(U(0,1)*bdf0 + Un(0,1)*bdf1 + Unn(0,1)*bdf2) + N[1]*(U(1,1)*bdf0 + Un(1,1)*bdf1 + Unn(1,1)*bdf2) + N[2]*(U(2,1)*bdf0 + Un(2,1)*bdf1 + Unn(2,1)*bdf2);
const double crhs12 =             1.0/crhs9;
const double crhs13 =             (1.0L/3.0L)*DN(0,0)*crhs12*mu;
const double crhs14 =             DN(0,0)*U(0,0) + DN(1,0)*U(1,0) + DN(2,0)*U(2,0);
const double crhs15 =             N[0]*U(0,1) + N[1]*U(1,1) + N[2]*U(2,1);
const double crhs16 =             crhs12*crhs14*crhs15;
const double crhs17 =             N[0]*U(0,2) + N[1]*U(1,2) + N[2]*U(2,2);
const double crhs18 =             DN(0,1)*U(0,0) + DN(1,1)*U(1,0) + DN(2,1)*U(2,0);
const double crhs19 =             3*crhs12*crhs18;
const double crhs20 =             -3*crhs0 - 3*crhs1 + 4*crhs16 + crhs17*crhs19 - 3*crhs2 - 4*crhs4 - 4*crhs5 - 4*crhs6;
const double crhs21 =             (1.0L/3.0L)*DN(0,1)*crhs12*mu;
const double crhs22 =             DN(0,0)*U(0,2);
const double crhs23 =             DN(0,1)*U(0,1);
const double crhs24 =             DN(1,0)*U(1,2);
const double crhs25 =             DN(1,1)*U(1,1);
const double crhs26 =             DN(2,0)*U(2,2);
const double crhs27 =             DN(2,1)*U(2,1);
const double crhs28 =             crhs12*crhs14*crhs17;
const double crhs29 =             crhs15*crhs19 + crhs20 + 2*crhs22 - 3*crhs23 + 2*crhs24 - 3*crhs25 + 2*crhs26 - 3*crhs27 - 2*crhs28;
const double crhs30 =             (1.0L/2.0L)*N[0];
const double crhs31 =             2*y;
const double crhs32 =             crhs31 - 2;
const double crhs33 =             DN(0,0)*U(0,3) + DN(1,0)*U(1,3) + DN(2,0)*U(2,3);
const double crhs34 =             crhs23 + crhs25 + crhs27;
const double crhs35 =             2*crhs12*crhs17;
const double crhs36 =             crhs34*crhs35;
const double crhs37 =             2*crhs12*crhs15;
const double crhs38 =             y - 3;
const double crhs39 =             crhs4 + crhs5 + crhs6;
const double crhs40 =             crhs22 + crhs24 + crhs26;
const double crhs41 =             crhs12*crhs17*crhs40;
const double crhs42 =             pow(crhs9, -2);
const double crhs43 =             2*crhs15*crhs17*crhs42;
const double crhs44 =             pow(crhs15, 2);
const double crhs45 =             y - 1;
const double crhs46 =             crhs44*crhs45;
const double crhs47 =             pow(crhs17, 2);
const double crhs48 =             crhs45*crhs47;
const double crhs49 =             crhs46 + crhs48;
const double crhs50 =             crhs14*crhs42*(-2*crhs44 + crhs49) - crhs18*crhs43 + crhs3*crhs37 + crhs32*crhs33 - crhs32*crhs41 + crhs36 - crhs37*crhs38*crhs39;
const double crhs51 =             N[0]*f_ext(0,1) + N[1]*f_ext(1,1) + N[2]*f_ext(2,1);
const double crhs52 =             N[0]*(U(0,2)*bdf0 + Un(0,2)*bdf1 + Unn(0,2)*bdf2) + N[1]*(U(1,2)*bdf0 + Un(1,2)*bdf1 + Unn(1,2)*bdf2) + N[2]*(U(2,2)*bdf0 + Un(2,2)*bdf1 + Unn(2,2)*bdf2);
const double crhs53 =             -crhs18*crhs37 - 3*crhs22 + 2*crhs23 - 3*crhs24 + 2*crhs25 - 3*crhs26 + 2*crhs27 + 3*crhs28;
const double crhs54 =             crhs12*crhs17*crhs18;
const double crhs55 =             -4*crhs0 - 4*crhs1 + 3*crhs16 - 4*crhs2 - 3*crhs4 - 3*crhs5 + crhs53 + 4*crhs54 - 3*crhs6;
const double crhs56 =             DN(0,1)*U(0,3) + DN(1,1)*U(1,3) + DN(2,1)*U(2,3);
const double crhs57 =             crhs37*crhs40;
const double crhs58 =             crhs12*crhs15*crhs34;
const double crhs59 =             -crhs14*crhs43 + crhs18*crhs42*(-2*crhs47 + crhs49) - crhs3*crhs35*crhs38 + crhs32*crhs56 - crhs32*crhs58 + crhs35*crhs39 + crhs57;
const double crhs60 =             N[0]*(U(0,3)*bdf0 + Un(0,3)*bdf1 + Unn(0,3)*bdf2) + N[1]*(U(1,3)*bdf0 + Un(1,3)*bdf1 + Unn(1,3)*bdf2) + N[2]*(U(2,3)*bdf0 + Un(2,3)*bdf1 + Unn(2,3)*bdf2);
const double crhs61 =             crhs15*crhs8 + crhs17*crhs51 + crhs9*(N[0]*r[0] + N[1]*r[1] + N[2]*r[2]);
const double crhs62 =             3*(lambda/cv);
const double crhs63 =             crhs12*crhs15;
const double crhs64 =             3*mu;
const double crhs65 =             -crhs62;
const double crhs66 =             crhs64 + crhs65;
const double crhs67 =             4*mu;
const double crhs68 =             crhs65 + crhs67;
const double crhs69 =             crhs15*crhs17*crhs42*mu;
const double crhs70 =             crhs12*crhs44;
const double crhs71 =             crhs12*crhs47;
const double crhs72 =             N[0]*U(0,3);
const double crhs73 =             N[1]*U(1,3);
const double crhs74 =             N[2]*U(2,3);
const double crhs75 =             crhs72 + crhs73 + crhs74;
const double crhs76 =             -crhs62*crhs70 - crhs62*crhs71 + crhs62*crhs75;
const double crhs77 =             -crhs12*crhs14*(crhs64*crhs71 + crhs67*crhs70 + crhs76) - crhs18*crhs69 + 3*crhs3*crhs63*mu + crhs33*crhs62 - crhs36*mu + crhs39*crhs63*crhs68 + crhs41*crhs66;
const double crhs78 =             (1.0L/3.0L)*crhs12*crhs77;
const double crhs79 =             crhs12*crhs17;
const double crhs80 =             -2*crhs72;
const double crhs81 =             -2*crhs73;
const double crhs82 =             -2*crhs74;
const double crhs83 =             -crhs32*crhs75;
const double crhs84 =             crhs12*crhs48;
const double crhs85 =             crhs12*crhs46;
const double crhs86 =             crhs80 + crhs81 + crhs82 + crhs83 + crhs85;
const double crhs87 =             crhs12*crhs45*(crhs44 + crhs47) + crhs84 + crhs86;
const double crhs88 =             crhs12*(-crhs12*crhs15*crhs17*crhs32*crhs40 + crhs15*crhs31*crhs33 - crhs15*crhs32*crhs34*crhs79 + crhs16*crhs87 + crhs17*crhs31*crhs56 - crhs3*(3*crhs84 + crhs86) - crhs39*(crhs80 + crhs81 + crhs82 + crhs83 + crhs84 + 3*crhs85) + crhs54*crhs87);
const double crhs89 =             (1.0L/3.0L)*crhs12*(-crhs12*crhs18*(crhs64*crhs70 + crhs67*crhs71 + crhs76) - crhs14*crhs69 + crhs3*crhs68*crhs79 + 3*crhs39*crhs79*mu + crhs56*crhs62 - crhs57*mu + crhs58*crhs66 + crhs77);
const double crhs90 =             N[1]*crhs9;
const double crhs91 =             (1.0L/3.0L)*DN(1,0)*crhs12*mu;
const double crhs92 =             (1.0L/3.0L)*DN(1,1)*crhs12*mu;
const double crhs93 =             (1.0L/2.0L)*N[1];
const double crhs94 =             N[2]*crhs9;
const double crhs95 =             (1.0L/3.0L)*DN(2,0)*crhs12*mu;
const double crhs96 =             (1.0L/3.0L)*DN(2,1)*crhs12*mu;
const double crhs97 =             (1.0L/2.0L)*N[2];
            rhs[0]=N[0]*crhs7;
            rhs[1]=N[0]*crhs11 - crhs10*crhs8 - crhs13*crhs20 - crhs21*crhs29 + crhs30*crhs50;
            rhs[2]=N[0]*crhs52 - crhs10*crhs51 - crhs13*crhs53 - crhs21*crhs55 + crhs30*crhs59;
            rhs[3]=DN(0,0)*crhs78 + DN(0,1)*crhs89 + N[0]*crhs60 - N[0]*crhs61 + crhs30*crhs88;
            rhs[4]=N[1]*crhs7;
            rhs[5]=N[1]*crhs11 - crhs20*crhs91 - crhs29*crhs92 + crhs50*crhs93 - crhs8*crhs90;
            rhs[6]=N[1]*crhs52 - crhs51*crhs90 - crhs53*crhs91 - crhs55*crhs92 + crhs59*crhs93;
            rhs[7]=DN(1,0)*crhs78 + DN(1,1)*crhs89 + N[1]*crhs60 - N[1]*crhs61 + crhs88*crhs93;
            rhs[8]=N[2]*crhs7;
            rhs[9]=N[2]*crhs11 - crhs20*crhs95 - crhs29*crhs96 + crhs50*crhs97 - crhs8*crhs94;
            rhs[10]=N[2]*crhs52 - crhs51*crhs94 - crhs53*crhs95 - crhs55*crhs96 + crhs59*crhs97;
            rhs[11]=DN(2,0)*crhs78 + DN(2,1)*crhs89 + N[2]*crhs60 - N[2]*crhs61 + crhs88*crhs97;

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
