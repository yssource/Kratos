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
const double clhs1 =             -bdf0*clhs0;
const double clhs2 =             DN(0,0)*N[0];
const double clhs3 =             DN(0,1)*N[0];
const double clhs4 =             N[0]*bdf0;
const double clhs5 =             -N[1]*clhs4;
const double clhs6 =             DN(1,0)*N[0];
const double clhs7 =             DN(1,1)*N[0];
const double clhs8 =             -N[2]*clhs4;
const double clhs9 =             DN(2,0)*N[0];
const double clhs10 =             DN(2,1)*N[0];
const double clhs11 =             N[0]*f_ext(0,0) + N[1]*f_ext(1,0) + N[2]*f_ext(2,0);
const double clhs12 =             clhs0*clhs11;
const double clhs13 =             N[0]*U(0,0) + N[1]*U(1,0) + N[2]*U(2,0);
const double clhs14 =             pow(clhs13, -2);
const double clhs15 =             (2.0L/3.0L)*DN(0,0)*clhs14*mu;
const double clhs16 =             2*DN(0,0);
const double clhs17 =             N[0]*U(0,1) + N[1]*U(1,1) + N[2]*U(2,1);
const double clhs18 =             clhs16*clhs17;
const double clhs19 =             N[0]*U(0,2) + N[1]*U(1,2) + N[2]*U(2,2);
const double clhs20 =             DN(0,1)*clhs19;
const double clhs21 =             DN(0,0)*U(0,1) + DN(1,0)*U(1,1) + DN(2,0)*U(2,1);
const double clhs22 =             N[0]*clhs21;
const double clhs23 =             2*clhs22;
const double clhs24 =             DN(0,1)*U(0,2) + DN(1,1)*U(1,2) + DN(2,1)*U(2,2);
const double clhs25 =             N[0]*clhs24;
const double clhs26 =             DN(0,0)*U(0,0) + DN(1,0)*U(1,0) + DN(2,0)*U(2,0);
const double clhs27 =             1.0/clhs13;
const double clhs28 =             N[0]*clhs27;
const double clhs29 =             clhs26*clhs28;
const double clhs30 =             4*clhs29;
const double clhs31 =             clhs17*clhs30;
const double clhs32 =             2*N[0]*U(0,2) + 2*N[1]*U(1,2) + 2*N[2]*U(2,2);
const double clhs33 =             DN(0,1)*U(0,0) + DN(1,1)*U(1,0) + DN(2,1)*U(2,0);
const double clhs34 =             clhs28*clhs33;
const double clhs35 =             clhs18 - clhs20 + clhs23 - clhs25 - clhs31 + clhs32*clhs34;
const double clhs36 =             (1.0L/3.0L)*DN(0,1)*clhs14*mu;
const double clhs37 =             4*DN(0,0);
const double clhs38 =             2*DN(0,1);
const double clhs39 =             clhs19*clhs38;
const double clhs40 =             2*clhs25;
const double clhs41 =             clhs17*clhs29;
const double clhs42 =             4*clhs34;
const double clhs43 =             clhs19*clhs42;
const double clhs44 =             3*DN(0,0);
const double clhs45 =             clhs19*clhs44;
const double clhs46 =             3*DN(0,1);
const double clhs47 =             clhs17*clhs46;
const double clhs48 =             DN(0,0)*U(0,2) + DN(1,0)*U(1,2) + DN(2,0)*U(2,2);
const double clhs49 =             N[0]*clhs48;
const double clhs50 =             3*clhs49;
const double clhs51 =             DN(0,1)*U(0,1) + DN(1,1)*U(1,1) + DN(2,1)*U(2,1);
const double clhs52 =             N[0]*clhs51;
const double clhs53 =             3*clhs52;
const double clhs54 =             clhs19*clhs29;
const double clhs55 =             6*clhs54;
const double clhs56 =             clhs17*clhs34;
const double clhs57 =             6*clhs56;
const double clhs58 =             clhs45 + clhs47 + clhs50 + clhs53 - clhs55 - clhs57;
const double clhs59 =             clhs17*clhs37 + 4*clhs22 - clhs39 - clhs40 - 8*clhs41 + clhs43 + clhs58;
const double clhs60 =             (1.0L/2.0L)*N[0]*clhs14;
const double clhs61 =             clhs17*clhs39;
const double clhs62 =             clhs17*clhs40;
const double clhs63 =             y - 3;
const double clhs64 =             clhs17*clhs63;
const double clhs65 =             y - 1;
const double clhs66 =             clhs49*clhs65;
const double clhs67 =             pow(clhs17, 2);
const double clhs68 =             clhs65*clhs67;
const double clhs69 =             pow(clhs19, 2);
const double clhs70 =             clhs65*clhs69;
const double clhs71 =             clhs68 + clhs70;
const double clhs72 =             -2*clhs67 + clhs71;
const double clhs73 =             2*clhs29;
const double clhs74 =             DN(0,0)*clhs72 + clhs17*clhs43 + clhs23*clhs64 - clhs32*clhs52 + clhs32*clhs66 - clhs61 - clhs62 - clhs72*clhs73;
const double clhs75 =             DN(0,0)*clhs27*mu;
const double clhs76 =             DN(0,0) - clhs29;
const double clhs77 =             clhs75*clhs76;
const double clhs78 =             (1.0L/3.0L)*DN(0,1)*clhs27*mu;
const double clhs79 =             3*clhs34;
const double clhs80 =             -clhs46 + clhs79;
const double clhs81 =             clhs30 - clhs37 + clhs80;
const double clhs82 =             clhs19*clhs34;
const double clhs83 =             clhs20 + clhs25 - clhs82;
const double clhs84 =             DN(0,0)*clhs17;
const double clhs85 =             -clhs22*clhs63 + clhs29*clhs64 - clhs63*clhs84 + clhs83;
const double clhs86 =             (1.0L/3.0L)*clhs27;
const double clhs87 =             mu*(DN(0,1) - clhs34);
const double clhs88 =             DN(0,1)*mu;
const double clhs89 =             3*clhs29;
const double clhs90 =             2*clhs34;
const double clhs91 =             -clhs38 + clhs44 - clhs89 + clhs90;
const double clhs92 =             3*N[0];
const double clhs93 =             DN(0,1)*clhs17;
const double clhs94 =             DN(0,0)*clhs19;
const double clhs95 =             -clhs52 - clhs54*clhs65 + clhs56 + clhs65*clhs94 + clhs66 - clhs93;
const double clhs96 =             N[0]*N[1];
const double clhs97 =             clhs11*clhs96;
const double clhs98 =             2*DN(1,0);
const double clhs99 =             clhs17*clhs98;
const double clhs100 =             DN(1,1)*clhs19;
const double clhs101 =             N[1]*clhs21;
const double clhs102 =             2*clhs101;
const double clhs103 =             N[1]*clhs24;
const double clhs104 =             N[1]*clhs27;
const double clhs105 =             clhs104*clhs26;
const double clhs106 =             4*clhs105;
const double clhs107 =             clhs106*clhs17;
const double clhs108 =             clhs104*clhs33;
const double clhs109 =             -clhs100 + clhs102 - clhs103 - clhs107 + clhs108*clhs32 + clhs99;
const double clhs110 =             4*DN(1,0);
const double clhs111 =             2*DN(1,1);
const double clhs112 =             clhs111*clhs19;
const double clhs113 =             2*clhs103;
const double clhs114 =             clhs105*clhs17;
const double clhs115 =             4*clhs108;
const double clhs116 =             clhs115*clhs19;
const double clhs117 =             3*DN(1,0);
const double clhs118 =             clhs117*clhs19;
const double clhs119 =             3*DN(1,1);
const double clhs120 =             clhs119*clhs17;
const double clhs121 =             N[1]*clhs48;
const double clhs122 =             3*clhs121;
const double clhs123 =             N[1]*clhs51;
const double clhs124 =             3*clhs123;
const double clhs125 =             clhs105*clhs19;
const double clhs126 =             6*clhs125;
const double clhs127 =             clhs108*clhs17;
const double clhs128 =             6*clhs127;
const double clhs129 =             clhs118 + clhs120 + clhs122 + clhs124 - clhs126 - clhs128;
const double clhs130 =             4*clhs101 + clhs110*clhs17 - clhs112 - clhs113 - 8*clhs114 + clhs116 + clhs129;
const double clhs131 =             clhs112*clhs17;
const double clhs132 =             clhs113*clhs17;
const double clhs133 =             clhs121*clhs65;
const double clhs134 =             2*clhs105;
const double clhs135 =             DN(1,0)*clhs72 + clhs102*clhs64 + clhs116*clhs17 - clhs123*clhs32 - clhs131 - clhs132 + clhs133*clhs32 - clhs134*clhs72;
const double clhs136 =             DN(1,0) - clhs105;
const double clhs137 =             clhs136*clhs75;
const double clhs138 =             3*clhs108;
const double clhs139 =             -clhs119 + clhs138;
const double clhs140 =             clhs106 - clhs110 + clhs139;
const double clhs141 =             clhs108*clhs19;
const double clhs142 =             clhs100 + clhs103 - clhs141;
const double clhs143 =             DN(1,0)*clhs17;
const double clhs144 =             -clhs101*clhs63 + clhs105*clhs64 + clhs142 - clhs143*clhs63;
const double clhs145 =             mu*(DN(1,1) - clhs108);
const double clhs146 =             3*clhs105;
const double clhs147 =             2*clhs108;
const double clhs148 =             -clhs111 + clhs117 - clhs146 + clhs147;
const double clhs149 =             DN(1,1)*clhs17;
const double clhs150 =             DN(1,0)*clhs19;
const double clhs151 =             -clhs123 - clhs125*clhs65 + clhs127 + clhs133 - clhs149 + clhs150*clhs65;
const double clhs152 =             N[0]*N[2];
const double clhs153 =             clhs11*clhs152;
const double clhs154 =             2*DN(2,0);
const double clhs155 =             clhs154*clhs17;
const double clhs156 =             DN(2,1)*clhs19;
const double clhs157 =             N[2]*clhs21;
const double clhs158 =             2*clhs157;
const double clhs159 =             N[2]*clhs24;
const double clhs160 =             N[2]*clhs27;
const double clhs161 =             clhs160*clhs26;
const double clhs162 =             4*clhs161;
const double clhs163 =             clhs162*clhs17;
const double clhs164 =             clhs160*clhs33;
const double clhs165 =             clhs155 - clhs156 + clhs158 - clhs159 - clhs163 + clhs164*clhs32;
const double clhs166 =             4*DN(2,0);
const double clhs167 =             2*DN(2,1);
const double clhs168 =             clhs167*clhs19;
const double clhs169 =             2*clhs159;
const double clhs170 =             clhs161*clhs17;
const double clhs171 =             4*clhs164;
const double clhs172 =             clhs171*clhs19;
const double clhs173 =             3*DN(2,0);
const double clhs174 =             clhs173*clhs19;
const double clhs175 =             3*DN(2,1);
const double clhs176 =             clhs17*clhs175;
const double clhs177 =             N[2]*clhs48;
const double clhs178 =             3*clhs177;
const double clhs179 =             N[2]*clhs51;
const double clhs180 =             3*clhs179;
const double clhs181 =             clhs161*clhs19;
const double clhs182 =             6*clhs181;
const double clhs183 =             clhs164*clhs17;
const double clhs184 =             6*clhs183;
const double clhs185 =             clhs174 + clhs176 + clhs178 + clhs180 - clhs182 - clhs184;
const double clhs186 =             4*clhs157 + clhs166*clhs17 - clhs168 - clhs169 - 8*clhs170 + clhs172 + clhs185;
const double clhs187 =             clhs168*clhs17;
const double clhs188 =             clhs169*clhs17;
const double clhs189 =             clhs177*clhs65;
const double clhs190 =             2*clhs161;
const double clhs191 =             DN(2,0)*clhs72 + clhs158*clhs64 + clhs17*clhs172 - clhs179*clhs32 - clhs187 - clhs188 + clhs189*clhs32 - clhs190*clhs72;
const double clhs192 =             DN(2,0) - clhs161;
const double clhs193 =             clhs192*clhs75;
const double clhs194 =             3*clhs164;
const double clhs195 =             -clhs175 + clhs194;
const double clhs196 =             clhs162 - clhs166 + clhs195;
const double clhs197 =             clhs164*clhs19;
const double clhs198 =             clhs156 + clhs159 - clhs197;
const double clhs199 =             DN(2,0)*clhs17;
const double clhs200 =             -clhs157*clhs63 + clhs161*clhs64 + clhs198 - clhs199*clhs63;
const double clhs201 =             mu*(DN(2,1) - clhs164);
const double clhs202 =             3*clhs161;
const double clhs203 =             2*clhs164;
const double clhs204 =             -clhs167 + clhs173 - clhs202 + clhs203;
const double clhs205 =             DN(2,1)*clhs17;
const double clhs206 =             DN(2,0)*clhs19;
const double clhs207 =             -clhs179 - clhs181*clhs65 + clhs183 + clhs189 - clhs205 + clhs206*clhs65;
const double clhs208 =             N[0]*f_ext(0,1) + N[1]*f_ext(1,1) + N[2]*f_ext(2,1);
const double clhs209 =             clhs0*clhs208;
const double clhs210 =             DN(0,0)*clhs14*mu;
const double clhs211 =             clhs49 + clhs94;
const double clhs212 =             2*N[0]*U(0,1) + 2*N[1]*U(1,1) + 2*N[2]*U(2,1);
const double clhs213 =             clhs211 - clhs212*clhs34 - clhs29*clhs32 + clhs52 + clhs93;
const double clhs214 =             4*DN(0,1);
const double clhs215 =             -clhs18 + clhs19*clhs214 - clhs23 + 4*clhs25 + clhs31 + clhs58 - 8*clhs82;
const double clhs216 =             clhs18*clhs19;
const double clhs217 =             clhs19*clhs23;
const double clhs218 =             clhs19*clhs63;
const double clhs219 =             clhs52*clhs65;
const double clhs220 =             -2*clhs69 + clhs71;
const double clhs221 =             DN(0,1)*clhs220 + clhs19*clhs31 + clhs212*clhs219 - clhs212*clhs49 - clhs216 - clhs217 + clhs218*clhs40 - clhs220*clhs90;
const double clhs222 =             clhs16 - clhs73 + clhs80;
const double clhs223 =             clhs211 - clhs219 - clhs54 + clhs56*clhs65 - clhs65*clhs93;
const double clhs224 =             -clhs214 + clhs42 - clhs44 + clhs89;
const double clhs225 =             clhs22 - clhs41 + clhs84;
const double clhs226 =             -clhs20*clhs63 + clhs218*clhs34 + clhs225 - clhs25*clhs63;
const double clhs227 =             clhs208*clhs96;
const double clhs228 =             clhs121 + clhs150;
const double clhs229 =             -clhs105*clhs32 - clhs108*clhs212 + clhs123 + clhs149 + clhs228;
const double clhs230 =             4*DN(1,1);
const double clhs231 =             -clhs102 + 4*clhs103 + clhs107 + clhs129 - 8*clhs141 + clhs19*clhs230 - clhs99;
const double clhs232 =             clhs19*clhs99;
const double clhs233 =             clhs102*clhs19;
const double clhs234 =             clhs123*clhs65;
const double clhs235 =             DN(1,1)*clhs220 + clhs107*clhs19 + clhs113*clhs218 - clhs121*clhs212 - clhs147*clhs220 + clhs212*clhs234 - clhs232 - clhs233;
const double clhs236 =             -clhs134 + clhs139 + clhs98;
const double clhs237 =             -clhs125 + clhs127*clhs65 - clhs149*clhs65 + clhs228 - clhs234;
const double clhs238 =             clhs115 - clhs117 + clhs146 - clhs230;
const double clhs239 =             clhs101 - clhs114 + clhs143;
const double clhs240 =             -clhs100*clhs63 - clhs103*clhs63 + clhs108*clhs218 + clhs239;
const double clhs241 =             clhs152*clhs208;
const double clhs242 =             clhs177 + clhs206;
const double clhs243 =             -clhs161*clhs32 - clhs164*clhs212 + clhs179 + clhs205 + clhs242;
const double clhs244 =             4*DN(2,1);
const double clhs245 =             -clhs155 - clhs158 + 4*clhs159 + clhs163 + clhs185 + clhs19*clhs244 - 8*clhs197;
const double clhs246 =             clhs155*clhs19;
const double clhs247 =             clhs158*clhs19;
const double clhs248 =             clhs179*clhs65;
const double clhs249 =             DN(2,1)*clhs220 + clhs163*clhs19 + clhs169*clhs218 - clhs177*clhs212 - clhs203*clhs220 + clhs212*clhs248 - clhs246 - clhs247;
const double clhs250 =             clhs154 - clhs190 + clhs195;
const double clhs251 =             -clhs181 + clhs183*clhs65 - clhs205*clhs65 + clhs242 - clhs248;
const double clhs252 =             clhs171 - clhs173 + clhs202 - clhs244;
const double clhs253 =             clhs157 - clhs170 + clhs199;
const double clhs254 =             -clhs156*clhs63 - clhs159*clhs63 + clhs164*clhs218 + clhs253;
const double clhs255 =             N[0]*r[0] + N[1]*r[1] + N[2]*r[2];
const double clhs256 =             (1.0L/3.0L)*DN(0,0)*clhs14;
const double clhs257 =             1.0/cv;
const double clhs258 =             3*N[0]*clhs257*lambda;
const double clhs259 =             DN(0,0)*U(0,3) + DN(1,0)*U(1,3) + DN(2,0)*U(2,3);
const double clhs260 =             clhs19*clhs27*mu;
const double clhs261 =             6*clhs19*clhs27*mu;
const double clhs262 =             4*clhs17*clhs27*mu;
const double clhs263 =             clhs257*lambda;
const double clhs264 =             -clhs263 + mu;
const double clhs265 =             6*clhs19*clhs264*clhs27;
const double clhs266 =             4*mu;
const double clhs267 =             3*clhs263;
const double clhs268 =             clhs266 - clhs267;
const double clhs269 =             clhs17*clhs268*clhs27;
const double clhs270 =             3*N[0]*clhs14*clhs17*clhs19*mu;
const double clhs271 =             N[0]*U(0,3);
const double clhs272 =             N[1]*U(1,3);
const double clhs273 =             N[2]*U(2,3);
const double clhs274 =             clhs267*(clhs271 + clhs272 + clhs273);
const double clhs275 =             clhs27*clhs67;
const double clhs276 =             3*mu;
const double clhs277 =             clhs27*clhs69;
const double clhs278 =             -clhs267*clhs275;
const double clhs279 =             -clhs267*clhs277;
const double clhs280 =             clhs266*clhs275 + clhs276*clhs277 + clhs278 + clhs279;
const double clhs281 =             clhs274 + clhs280;
const double clhs282 =             2*clhs271;
const double clhs283 =             2*clhs272;
const double clhs284 =             2*clhs273;
const double clhs285 =             clhs282 + clhs283 + clhs284;
const double clhs286 =             clhs263*clhs285;
const double clhs287 =             clhs280 + clhs286;
const double clhs288 =             DN(0,0)*clhs281 + clhs23*clhs269 - clhs25*clhs262 + clhs258*clhs259 + clhs260*clhs93 + clhs261*clhs52 + clhs265*clhs49 - clhs270*clhs33 - clhs287*clhs89;
const double clhs289 =             2*N[0]*y;
const double clhs290 =             clhs259*clhs289;
const double clhs291 =             DN(0,1)*U(0,3) + DN(1,1)*U(1,3) + DN(2,1)*U(2,3);
const double clhs292 =             clhs289*clhs291;
const double clhs293 =             4*clhs17*clhs19*clhs27;
const double clhs294 =             clhs27*clhs70;
const double clhs295 =             -clhs282;
const double clhs296 =             -clhs283;
const double clhs297 =             -clhs284;
const double clhs298 =             clhs285*clhs65;
const double clhs299 =             -clhs298;
const double clhs300 =             clhs27*clhs68;
const double clhs301 =             clhs27*clhs65;
const double clhs302 =             clhs301*(clhs67 + clhs69);
const double clhs303 =             clhs295 + clhs296 + clhs297 + clhs299 + clhs300 + clhs302;
const double clhs304 =             clhs294 + clhs303;
const double clhs305 =             clhs294 + clhs295 + clhs296 + clhs297 + clhs299;
const double clhs306 =             5*clhs300 + clhs302 + clhs305;
const double clhs307 =             5*clhs294 + clhs303;
const double clhs308 =             clhs282 + clhs283 + clhs284 - clhs294 + clhs298 - clhs300 - 2*clhs302;
const double clhs309 =             clhs212*clhs308;
const double clhs310 =             N[0]*clhs27*clhs32*clhs33;
const double clhs311 =             -clhs17*clhs290 - clhs19*clhs292 + clhs20*clhs304 + clhs219*clhs293 + clhs22*clhs306 + clhs25*clhs307 + clhs29*clhs309 + clhs293*clhs66 + clhs304*clhs84 + clhs308*clhs310;
const double clhs312 =             (1.0L/3.0L)*DN(0,1)*clhs14;
const double clhs313 =             4*clhs19*clhs27*mu;
const double clhs314 =             6*clhs17*clhs27*mu;
const double clhs315 =             6*clhs17*clhs264*clhs27;
const double clhs316 =             clhs19*clhs268*clhs27;
const double clhs317 =             clhs266*clhs277 + clhs275*clhs276 + clhs278 + clhs279;
const double clhs318 =             clhs274 + clhs317;
const double clhs319 =             clhs286 + clhs317;
const double clhs320 =             DN(0,1)*clhs318 - clhs22*clhs313 + clhs258*clhs291 - clhs26*clhs270 + clhs260*clhs84 + clhs288 + clhs314*clhs49 + clhs315*clhs52 + clhs316*clhs40 - clhs319*clhs79;
const double clhs321 =             clhs19*mu;
const double clhs322 =             clhs212*clhs268;
const double clhs323 =             -clhs22*clhs268 - clhs268*clhs84 + clhs29*clhs322 + clhs321*clhs34 - clhs321*clhs46 + clhs40*mu;
const double clhs324 =             clhs16*clhs321 - clhs264*clhs47 - clhs264*clhs53 + clhs264*clhs57 + clhs29*clhs321 + clhs323 - clhs50*mu;
const double clhs325 =             (1.0L/2.0L)*N[0]*clhs27;
const double clhs326 =             6*clhs17*clhs27*clhs65;
const double clhs327 =             clhs27*clhs32;
const double clhs328 =             4*N[0]*clhs14*clhs17*clhs19*clhs65;
const double clhs329 =             3*clhs300 + clhs305;
const double clhs330 =             -DN(0,0)*clhs329 - N[0]*clhs27*clhs32*clhs48*clhs65 - clhs219*clhs327 - clhs22*clhs326 + clhs29*clhs306 + clhs290 - clhs301*clhs61 - clhs301*clhs62 + clhs328*clhs33;
const double clhs331 =             clhs17*mu;
const double clhs332 =             -clhs264*clhs45 - clhs264*clhs50 + clhs264*clhs55 + clhs331*clhs34 + clhs331*clhs38 - clhs53*mu;
const double clhs333 =             -clhs20*clhs268 + clhs23*mu - clhs25*clhs268 + clhs268*clhs310 + clhs29*clhs331 - clhs331*clhs44 + clhs332;
const double clhs334 =             clhs212*clhs27;
const double clhs335 =             6*clhs19*clhs27*clhs65;
const double clhs336 =             3*clhs294 + clhs295 + clhs296 + clhs297 + clhs299 + clhs300;
const double clhs337 =             -DN(0,1)*clhs336 - N[0]*clhs212*clhs27*clhs51*clhs65 - clhs216*clhs301 - clhs217*clhs301 - clhs25*clhs335 + clhs26*clhs328 + clhs292 + clhs307*clhs34 - clhs334*clhs66;
const double clhs338 =             DN(0,0)*clhs257*clhs27*lambda;
const double clhs339 =             DN(0,1)*clhs257*clhs27*lambda;
const double clhs340 =             -DN(0,0) - DN(0,1) + clhs29 + clhs34;
const double clhs341 =             N[0]*clhs27*y;
const double clhs342 =             clhs225 + clhs83;
const double clhs343 =             clhs255*clhs96;
const double clhs344 =             3*N[1]*clhs257*lambda;
const double clhs345 =             3*N[1]*clhs14*clhs17*clhs19*mu;
const double clhs346 =             DN(1,0)*clhs281 + clhs102*clhs269 - clhs103*clhs262 + clhs121*clhs265 + clhs123*clhs261 - clhs146*clhs287 + clhs149*clhs260 + clhs259*clhs344 - clhs33*clhs345;
const double clhs347 =             2*N[1]*y;
const double clhs348 =             clhs259*clhs347;
const double clhs349 =             clhs291*clhs347;
const double clhs350 =             N[1]*clhs27*clhs32*clhs33;
const double clhs351 =             clhs100*clhs304 + clhs101*clhs306 + clhs103*clhs307 + clhs105*clhs309 + clhs133*clhs293 + clhs143*clhs304 - clhs17*clhs348 - clhs19*clhs349 + clhs234*clhs293 + clhs308*clhs350;
const double clhs352 =             DN(1,1)*clhs318 - clhs101*clhs313 + clhs113*clhs316 + clhs121*clhs314 + clhs123*clhs315 - clhs138*clhs319 + clhs143*clhs260 - clhs26*clhs345 + clhs291*clhs344 + clhs346;
const double clhs353 =             -clhs101*clhs268 + clhs105*clhs322 + clhs108*clhs321 + clhs113*mu - clhs119*clhs321 - clhs143*clhs268;
const double clhs354 =             clhs105*clhs321 - clhs120*clhs264 - clhs122*mu - clhs124*clhs264 + clhs128*clhs264 + clhs321*clhs98 + clhs353;
const double clhs355 =             4*N[1]*clhs14*clhs17*clhs19*clhs65;
const double clhs356 =             -DN(1,0)*clhs329 - N[1]*clhs27*clhs32*clhs48*clhs65 - clhs101*clhs326 + clhs105*clhs306 - clhs131*clhs301 - clhs132*clhs301 - clhs234*clhs327 + clhs33*clhs355 + clhs348;
const double clhs357 =             clhs108*clhs331 + clhs111*clhs331 - clhs118*clhs264 - clhs122*clhs264 - clhs124*mu + clhs126*clhs264;
const double clhs358 =             -clhs100*clhs268 + clhs102*mu - clhs103*clhs268 + clhs105*clhs331 - clhs117*clhs331 + clhs268*clhs350 + clhs357;
const double clhs359 =             -DN(1,1)*clhs336 - N[1]*clhs212*clhs27*clhs51*clhs65 - clhs103*clhs335 + clhs108*clhs307 - clhs133*clhs334 - clhs232*clhs301 - clhs233*clhs301 + clhs26*clhs355 + clhs349;
const double clhs360 =             -DN(1,0) - DN(1,1) + clhs105 + clhs108;
const double clhs361 =             clhs142 + clhs239;
const double clhs362 =             clhs152*clhs255;
const double clhs363 =             3*N[2]*clhs257*lambda;
const double clhs364 =             3*N[2]*clhs14*clhs17*clhs19*mu;
const double clhs365 =             DN(2,0)*clhs281 + clhs158*clhs269 - clhs159*clhs262 + clhs177*clhs265 + clhs179*clhs261 - clhs202*clhs287 + clhs205*clhs260 + clhs259*clhs363 - clhs33*clhs364;
const double clhs366 =             2*N[2]*y;
const double clhs367 =             clhs259*clhs366;
const double clhs368 =             clhs291*clhs366;
const double clhs369 =             N[2]*clhs27*clhs32*clhs33;
const double clhs370 =             clhs156*clhs304 + clhs157*clhs306 + clhs159*clhs307 + clhs161*clhs309 - clhs17*clhs367 + clhs189*clhs293 - clhs19*clhs368 + clhs199*clhs304 + clhs248*clhs293 + clhs308*clhs369;
const double clhs371 =             DN(2,1)*clhs318 - clhs157*clhs313 + clhs169*clhs316 + clhs177*clhs314 + clhs179*clhs315 - clhs194*clhs319 + clhs199*clhs260 - clhs26*clhs364 + clhs291*clhs363 + clhs365;
const double clhs372 =             -clhs157*clhs268 + clhs161*clhs322 + clhs164*clhs321 + clhs169*mu - clhs175*clhs321 - clhs199*clhs268;
const double clhs373 =             clhs154*clhs321 + clhs161*clhs321 - clhs176*clhs264 - clhs178*mu - clhs180*clhs264 + clhs184*clhs264 + clhs372;
const double clhs374 =             4*N[2]*clhs14*clhs17*clhs19*clhs65;
const double clhs375 =             -DN(2,0)*clhs329 - N[2]*clhs27*clhs32*clhs48*clhs65 - clhs157*clhs326 + clhs161*clhs306 - clhs187*clhs301 - clhs188*clhs301 - clhs248*clhs327 + clhs33*clhs374 + clhs367;
const double clhs376 =             clhs164*clhs331 + clhs167*clhs331 - clhs174*clhs264 - clhs178*clhs264 - clhs180*mu + clhs182*clhs264;
const double clhs377 =             -clhs156*clhs268 + clhs158*mu - clhs159*clhs268 + clhs161*clhs331 - clhs173*clhs331 + clhs268*clhs369 + clhs376;
const double clhs378 =             -DN(2,1)*clhs336 - N[2]*clhs212*clhs27*clhs51*clhs65 - clhs159*clhs335 + clhs164*clhs307 - clhs189*clhs334 - clhs246*clhs301 - clhs247*clhs301 + clhs26*clhs374 + clhs368;
const double clhs379 =             -DN(2,0) - DN(2,1) + clhs161 + clhs164;
const double clhs380 =             clhs198 + clhs253;
const double clhs381 =             DN(0,0)*N[1];
const double clhs382 =             DN(0,1)*N[1];
const double clhs383 =             pow(N[1], 2);
const double clhs384 =             -bdf0*clhs383;
const double clhs385 =             DN(1,0)*N[1];
const double clhs386 =             DN(1,1)*N[1];
const double clhs387 =             N[1]*N[2];
const double clhs388 =             -bdf0*clhs387;
const double clhs389 =             DN(2,0)*N[1];
const double clhs390 =             DN(2,1)*N[1];
const double clhs391 =             (2.0L/3.0L)*DN(1,0)*clhs14*mu;
const double clhs392 =             (1.0L/3.0L)*DN(1,1)*clhs14*mu;
const double clhs393 =             (1.0L/2.0L)*N[1]*clhs14;
const double clhs394 =             DN(1,0)*clhs27*mu;
const double clhs395 =             clhs394*clhs76;
const double clhs396 =             (1.0L/3.0L)*DN(1,1)*clhs27*mu;
const double clhs397 =             DN(1,1)*mu;
const double clhs398 =             3*N[1];
const double clhs399 =             clhs11*clhs383;
const double clhs400 =             clhs136*clhs394;
const double clhs401 =             clhs11*clhs387;
const double clhs402 =             clhs192*clhs394;
const double clhs403 =             DN(1,0)*clhs14*mu;
const double clhs404 =             clhs208*clhs383;
const double clhs405 =             clhs208*clhs387;
const double clhs406 =             (1.0L/3.0L)*DN(1,0)*clhs14;
const double clhs407 =             (1.0L/3.0L)*DN(1,1)*clhs14;
const double clhs408 =             (1.0L/2.0L)*N[1]*clhs27;
const double clhs409 =             DN(1,0)*clhs257*clhs27*lambda;
const double clhs410 =             DN(1,1)*clhs257*clhs27*lambda;
const double clhs411 =             N[1]*clhs27*y;
const double clhs412 =             clhs255*clhs387;
const double clhs413 =             DN(0,0)*N[2];
const double clhs414 =             DN(0,1)*N[2];
const double clhs415 =             DN(1,0)*N[2];
const double clhs416 =             DN(1,1)*N[2];
const double clhs417 =             pow(N[2], 2);
const double clhs418 =             -bdf0*clhs417;
const double clhs419 =             DN(2,0)*N[2];
const double clhs420 =             DN(2,1)*N[2];
const double clhs421 =             (2.0L/3.0L)*DN(2,0)*clhs14*mu;
const double clhs422 =             (1.0L/3.0L)*DN(2,1)*clhs14*mu;
const double clhs423 =             (1.0L/2.0L)*N[2]*clhs14;
const double clhs424 =             DN(2,0)*clhs27*mu;
const double clhs425 =             clhs424*clhs76;
const double clhs426 =             (1.0L/3.0L)*DN(2,1)*clhs27*mu;
const double clhs427 =             DN(2,1)*mu;
const double clhs428 =             3*N[2];
const double clhs429 =             clhs136*clhs424;
const double clhs430 =             clhs11*clhs417;
const double clhs431 =             clhs192*clhs424;
const double clhs432 =             DN(2,0)*clhs14*mu;
const double clhs433 =             clhs208*clhs417;
const double clhs434 =             (1.0L/3.0L)*DN(2,0)*clhs14;
const double clhs435 =             (1.0L/3.0L)*DN(2,1)*clhs14;
const double clhs436 =             (1.0L/2.0L)*N[2]*clhs27;
const double clhs437 =             DN(2,0)*clhs257*clhs27*lambda;
const double clhs438 =             DN(2,1)*clhs257*clhs27*lambda;
const double clhs439 =             N[2]*clhs27*y;
            lhs(0,0)=clhs1;
            lhs(0,1)=-clhs2;
            lhs(0,2)=-clhs3;
            lhs(0,3)=0;
            lhs(0,4)=clhs5;
            lhs(0,5)=-clhs6;
            lhs(0,6)=-clhs7;
            lhs(0,7)=0;
            lhs(0,8)=clhs8;
            lhs(0,9)=-clhs9;
            lhs(0,10)=-clhs10;
            lhs(0,11)=0;
            lhs(1,0)=clhs12 - clhs15*clhs35 - clhs36*clhs59 - clhs60*clhs74;
            lhs(1,1)=clhs1 - clhs28*clhs85 + (4.0L/3.0L)*clhs77 - clhs78*clhs81;
            lhs(1,2)=clhs86*(-clhs16*clhs87 + clhs88*clhs91 + clhs92*clhs95);
            lhs(1,3)=-clhs2*clhs65;
            lhs(1,4)=-clhs109*clhs15 - clhs130*clhs36 - clhs135*clhs60 + clhs97;
            lhs(1,5)=(4.0L/3.0L)*clhs137 - clhs140*clhs78 - clhs144*clhs28 + clhs5;
            lhs(1,6)=clhs86*(-clhs145*clhs16 + clhs148*clhs88 + clhs151*clhs92);
            lhs(1,7)=-clhs6*clhs65;
            lhs(1,8)=-clhs15*clhs165 + clhs153 - clhs186*clhs36 - clhs191*clhs60;
            lhs(1,9)=(4.0L/3.0L)*clhs193 - clhs196*clhs78 - clhs200*clhs28 + clhs8;
            lhs(1,10)=clhs86*(-clhs16*clhs201 + clhs204*clhs88 + clhs207*clhs92);
            lhs(1,11)=-clhs65*clhs9;
            lhs(2,0)=clhs209 - clhs210*clhs213 - clhs215*clhs36 - clhs221*clhs60;
            lhs(2,1)=clhs86*(-clhs222*clhs88 - clhs223*clhs92 + clhs44*clhs87);
            lhs(2,2)=clhs1 - clhs224*clhs78 - clhs226*clhs28 + clhs77;
            lhs(2,3)=-clhs3*clhs65;
            lhs(2,4)=-clhs210*clhs229 + clhs227 - clhs231*clhs36 - clhs235*clhs60;
            lhs(2,5)=clhs86*(clhs145*clhs44 - clhs236*clhs88 - clhs237*clhs92);
            lhs(2,6)=clhs137 - clhs238*clhs78 - clhs240*clhs28 + clhs5;
            lhs(2,7)=-clhs65*clhs7;
            lhs(2,8)=-clhs210*clhs243 + clhs241 - clhs245*clhs36 - clhs249*clhs60;
            lhs(2,9)=clhs86*(clhs201*clhs44 - clhs250*clhs88 - clhs251*clhs92);
            lhs(2,10)=clhs193 - clhs252*clhs78 - clhs254*clhs28 + clhs8;
            lhs(2,11)=-clhs10*clhs65;
            lhs(3,0)=clhs0*clhs255 - clhs256*clhs288 - clhs311*clhs60 - clhs312*clhs320;
            lhs(3,1)=clhs12 - clhs256*clhs323 - clhs312*clhs324 - clhs325*clhs330;
            lhs(3,2)=clhs209 - clhs256*clhs332 - clhs312*clhs333 - clhs325*clhs337;
            lhs(3,3)=clhs1 + clhs338*clhs76 - clhs339*clhs340 - clhs341*clhs342;
            lhs(3,4)=-clhs256*clhs346 - clhs312*clhs352 + clhs343 - clhs351*clhs60;
            lhs(3,5)=-clhs256*clhs353 - clhs312*clhs354 - clhs325*clhs356 + clhs97;
            lhs(3,6)=clhs227 - clhs256*clhs357 - clhs312*clhs358 - clhs325*clhs359;
            lhs(3,7)=clhs136*clhs338 - clhs339*clhs360 - clhs341*clhs361 + clhs5;
            lhs(3,8)=-clhs256*clhs365 - clhs312*clhs371 + clhs362 - clhs370*clhs60;
            lhs(3,9)=clhs153 - clhs256*clhs372 - clhs312*clhs373 - clhs325*clhs375;
            lhs(3,10)=clhs241 - clhs256*clhs376 - clhs312*clhs377 - clhs325*clhs378;
            lhs(3,11)=clhs192*clhs338 - clhs339*clhs379 - clhs341*clhs380 + clhs8;
            lhs(4,0)=clhs5;
            lhs(4,1)=-clhs381;
            lhs(4,2)=-clhs382;
            lhs(4,3)=0;
            lhs(4,4)=clhs384;
            lhs(4,5)=-clhs385;
            lhs(4,6)=-clhs386;
            lhs(4,7)=0;
            lhs(4,8)=clhs388;
            lhs(4,9)=-clhs389;
            lhs(4,10)=-clhs390;
            lhs(4,11)=0;
            lhs(5,0)=-clhs35*clhs391 - clhs392*clhs59 - clhs393*clhs74 + clhs97;
            lhs(5,1)=-clhs104*clhs85 + (4.0L/3.0L)*clhs395 - clhs396*clhs81 + clhs5;
            lhs(5,2)=clhs86*(clhs397*clhs91 + clhs398*clhs95 - clhs87*clhs98);
            lhs(5,3)=-clhs381*clhs65;
            lhs(5,4)=-clhs109*clhs391 - clhs130*clhs392 - clhs135*clhs393 + clhs399;
            lhs(5,5)=-clhs104*clhs144 - clhs140*clhs396 + clhs384 + (4.0L/3.0L)*clhs400;
            lhs(5,6)=clhs86*(-clhs145*clhs98 + clhs148*clhs397 + clhs151*clhs398);
            lhs(5,7)=-clhs385*clhs65;
            lhs(5,8)=-clhs165*clhs391 - clhs186*clhs392 - clhs191*clhs393 + clhs401;
            lhs(5,9)=-clhs104*clhs200 - clhs196*clhs396 + clhs388 + (4.0L/3.0L)*clhs402;
            lhs(5,10)=clhs86*(-clhs201*clhs98 + clhs204*clhs397 + clhs207*clhs398);
            lhs(5,11)=-clhs389*clhs65;
            lhs(6,0)=-clhs213*clhs403 - clhs215*clhs392 - clhs221*clhs393 + clhs227;
            lhs(6,1)=clhs86*(clhs117*clhs87 - clhs222*clhs397 - clhs223*clhs398);
            lhs(6,2)=-clhs104*clhs226 - clhs224*clhs396 + clhs395 + clhs5;
            lhs(6,3)=-clhs382*clhs65;
            lhs(6,4)=-clhs229*clhs403 - clhs231*clhs392 - clhs235*clhs393 + clhs404;
            lhs(6,5)=clhs86*(clhs117*clhs145 - clhs236*clhs397 - clhs237*clhs398);
            lhs(6,6)=-clhs104*clhs240 - clhs238*clhs396 + clhs384 + clhs400;
            lhs(6,7)=-clhs386*clhs65;
            lhs(6,8)=-clhs243*clhs403 - clhs245*clhs392 - clhs249*clhs393 + clhs405;
            lhs(6,9)=clhs86*(clhs117*clhs201 - clhs250*clhs397 - clhs251*clhs398);
            lhs(6,10)=-clhs104*clhs254 - clhs252*clhs396 + clhs388 + clhs402;
            lhs(6,11)=-clhs390*clhs65;
            lhs(7,0)=-clhs288*clhs406 - clhs311*clhs393 - clhs320*clhs407 + clhs343;
            lhs(7,1)=-clhs323*clhs406 - clhs324*clhs407 - clhs330*clhs408 + clhs97;
            lhs(7,2)=clhs227 - clhs332*clhs406 - clhs333*clhs407 - clhs337*clhs408;
            lhs(7,3)=-clhs340*clhs410 - clhs342*clhs411 + clhs409*clhs76 + clhs5;
            lhs(7,4)=clhs255*clhs383 - clhs346*clhs406 - clhs351*clhs393 - clhs352*clhs407;
            lhs(7,5)=-clhs353*clhs406 - clhs354*clhs407 - clhs356*clhs408 + clhs399;
            lhs(7,6)=-clhs357*clhs406 - clhs358*clhs407 - clhs359*clhs408 + clhs404;
            lhs(7,7)=clhs136*clhs409 - clhs360*clhs410 - clhs361*clhs411 + clhs384;
            lhs(7,8)=-clhs365*clhs406 - clhs370*clhs393 - clhs371*clhs407 + clhs412;
            lhs(7,9)=-clhs372*clhs406 - clhs373*clhs407 - clhs375*clhs408 + clhs401;
            lhs(7,10)=-clhs376*clhs406 - clhs377*clhs407 - clhs378*clhs408 + clhs405;
            lhs(7,11)=clhs192*clhs409 - clhs379*clhs410 - clhs380*clhs411 + clhs388;
            lhs(8,0)=clhs8;
            lhs(8,1)=-clhs413;
            lhs(8,2)=-clhs414;
            lhs(8,3)=0;
            lhs(8,4)=clhs388;
            lhs(8,5)=-clhs415;
            lhs(8,6)=-clhs416;
            lhs(8,7)=0;
            lhs(8,8)=clhs418;
            lhs(8,9)=-clhs419;
            lhs(8,10)=-clhs420;
            lhs(8,11)=0;
            lhs(9,0)=clhs153 - clhs35*clhs421 - clhs422*clhs59 - clhs423*clhs74;
            lhs(9,1)=-clhs160*clhs85 + (4.0L/3.0L)*clhs425 - clhs426*clhs81 + clhs8;
            lhs(9,2)=clhs86*(-clhs154*clhs87 + clhs427*clhs91 + clhs428*clhs95);
            lhs(9,3)=-clhs413*clhs65;
            lhs(9,4)=-clhs109*clhs421 - clhs130*clhs422 - clhs135*clhs423 + clhs401;
            lhs(9,5)=-clhs140*clhs426 - clhs144*clhs160 + clhs388 + (4.0L/3.0L)*clhs429;
            lhs(9,6)=clhs86*(-clhs145*clhs154 + clhs148*clhs427 + clhs151*clhs428);
            lhs(9,7)=-clhs415*clhs65;
            lhs(9,8)=-clhs165*clhs421 - clhs186*clhs422 - clhs191*clhs423 + clhs430;
            lhs(9,9)=-clhs160*clhs200 - clhs196*clhs426 + clhs418 + (4.0L/3.0L)*clhs431;
            lhs(9,10)=clhs86*(-clhs154*clhs201 + clhs204*clhs427 + clhs207*clhs428);
            lhs(9,11)=-clhs419*clhs65;
            lhs(10,0)=-clhs213*clhs432 - clhs215*clhs422 - clhs221*clhs423 + clhs241;
            lhs(10,1)=clhs86*(clhs173*clhs87 - clhs222*clhs427 - clhs223*clhs428);
            lhs(10,2)=-clhs160*clhs226 - clhs224*clhs426 + clhs425 + clhs8;
            lhs(10,3)=-clhs414*clhs65;
            lhs(10,4)=-clhs229*clhs432 - clhs231*clhs422 - clhs235*clhs423 + clhs405;
            lhs(10,5)=clhs86*(clhs145*clhs173 - clhs236*clhs427 - clhs237*clhs428);
            lhs(10,6)=-clhs160*clhs240 - clhs238*clhs426 + clhs388 + clhs429;
            lhs(10,7)=-clhs416*clhs65;
            lhs(10,8)=-clhs243*clhs432 - clhs245*clhs422 - clhs249*clhs423 + clhs433;
            lhs(10,9)=clhs86*(clhs173*clhs201 - clhs250*clhs427 - clhs251*clhs428);
            lhs(10,10)=-clhs160*clhs254 - clhs252*clhs426 + clhs418 + clhs431;
            lhs(10,11)=-clhs420*clhs65;
            lhs(11,0)=-clhs288*clhs434 - clhs311*clhs423 - clhs320*clhs435 + clhs362;
            lhs(11,1)=clhs153 - clhs323*clhs434 - clhs324*clhs435 - clhs330*clhs436;
            lhs(11,2)=clhs241 - clhs332*clhs434 - clhs333*clhs435 - clhs337*clhs436;
            lhs(11,3)=-clhs340*clhs438 - clhs342*clhs439 + clhs437*clhs76 + clhs8;
            lhs(11,4)=-clhs346*clhs434 - clhs351*clhs423 - clhs352*clhs435 + clhs412;
            lhs(11,5)=-clhs353*clhs434 - clhs354*clhs435 - clhs356*clhs436 + clhs401;
            lhs(11,6)=-clhs357*clhs434 - clhs358*clhs435 - clhs359*clhs436 + clhs405;
            lhs(11,7)=clhs136*clhs437 - clhs360*clhs438 - clhs361*clhs439 + clhs388;
            lhs(11,8)=clhs255*clhs417 - clhs365*clhs434 - clhs370*clhs423 - clhs371*clhs435;
            lhs(11,9)=-clhs372*clhs434 - clhs373*clhs435 - clhs375*clhs436 + clhs430;
            lhs(11,10)=-clhs376*clhs434 - clhs377*clhs435 - clhs378*clhs436 + clhs433;
            lhs(11,11)=clhs192*clhs437 - clhs379*clhs438 - clhs380*clhs439 + clhs418;


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
const double crhs13 =             DN(0,0)*crhs12*mu;
const double crhs14 =             DN(0,0)*U(0,0) + DN(1,0)*U(1,0) + DN(2,0)*U(2,0);
const double crhs15 =             N[0]*U(0,1) + N[1]*U(1,1) + N[2]*U(2,1);
const double crhs16 =             DN(0,1)*U(0,0) + DN(1,1)*U(1,0) + DN(2,1)*U(2,0);
const double crhs17 =             N[0]*U(0,2) + N[1]*U(1,2) + N[2]*U(2,2);
const double crhs18 =             -4.0L/3.0L*DN(0,0)*U(0,1) + (2.0L/3.0L)*DN(0,1)*U(0,2) - 4.0L/3.0L*DN(1,0)*U(1,1) + (2.0L/3.0L)*DN(1,1)*U(1,2) - 4.0L/3.0L*DN(2,0)*U(2,1) + (2.0L/3.0L)*DN(2,1)*U(2,2) + (4.0L/3.0L)*crhs12*crhs14*crhs15 - 2.0L/3.0L*crhs12*crhs16*crhs17;
const double crhs19 =             (1.0L/3.0L)*DN(0,1)*crhs12*mu;
const double crhs20 =             DN(0,0)*U(0,2);
const double crhs21 =             3*crhs20;
const double crhs22 =             DN(0,1)*U(0,1);
const double crhs23 =             3*crhs22;
const double crhs24 =             DN(1,0)*U(1,2);
const double crhs25 =             3*crhs24;
const double crhs26 =             DN(1,1)*U(1,1);
const double crhs27 =             3*crhs26;
const double crhs28 =             DN(2,0)*U(2,2);
const double crhs29 =             3*crhs28;
const double crhs30 =             DN(2,1)*U(2,1);
const double crhs31 =             3*crhs30;
const double crhs32 =             crhs12*crhs14*crhs15;
const double crhs33 =             crhs12*crhs17;
const double crhs34 =             crhs14*crhs33;
const double crhs35 =             3*crhs34;
const double crhs36 =             crhs12*crhs15;
const double crhs37 =             crhs16*crhs36;
const double crhs38 =             3*crhs37;
const double crhs39 =             crhs16*crhs33;
const double crhs40 =             2*crhs0 + 2*crhs1 + 2*crhs2 - crhs21 - crhs23 - crhs25 - crhs27 - crhs29 - crhs31 + 4*crhs32 + crhs35 + crhs38 - 2*crhs39 - 4*crhs4 - 4*crhs5 - 4*crhs6;
const double crhs41 =             (1.0L/2.0L)*N[0];
const double crhs42 =             2*y;
const double crhs43 =             crhs42 - 2;
const double crhs44 =             DN(0,0)*U(0,3) + DN(1,0)*U(1,3) + DN(2,0)*U(2,3);
const double crhs45 =             crhs22 + crhs26 + crhs30;
const double crhs46 =             2*crhs12*crhs17;
const double crhs47 =             2*crhs12*crhs15;
const double crhs48 =             crhs3*crhs47;
const double crhs49 =             y - 3;
const double crhs50 =             crhs4 + crhs5 + crhs6;
const double crhs51 =             crhs20 + crhs24 + crhs28;
const double crhs52 =             crhs12*crhs17*crhs51;
const double crhs53 =             pow(crhs9, -2);
const double crhs54 =             2*crhs15*crhs17*crhs53;
const double crhs55 =             pow(crhs15, 2);
const double crhs56 =             y - 1;
const double crhs57 =             crhs55*crhs56;
const double crhs58 =             pow(crhs17, 2);
const double crhs59 =             crhs56*crhs58;
const double crhs60 =             crhs57 + crhs59;
const double crhs61 =             crhs14*crhs53*(-2*crhs55 + crhs60) - crhs16*crhs54 + crhs43*crhs44 - crhs43*crhs52 + crhs45*crhs46 - crhs47*crhs49*crhs50 + crhs48;
const double crhs62 =             N[0]*f_ext(0,1) + N[1]*f_ext(1,1) + N[2]*f_ext(2,1);
const double crhs63 =             N[0]*(U(0,2)*bdf0 + Un(0,2)*bdf1 + Unn(0,2)*bdf2) + N[1]*(U(1,2)*bdf0 + Un(1,2)*bdf1 + Unn(1,2)*bdf2) + N[2]*(U(2,2)*bdf0 + Un(2,2)*bdf1 + Unn(2,2)*bdf2);
const double crhs64 =             -crhs20 - crhs22 - crhs24 - crhs26 - crhs28 - crhs30 + crhs34 + crhs37;
const double crhs65 =             4*crhs0 + 4*crhs1 + 4*crhs2 + crhs21 + crhs23 + crhs25 + crhs27 + crhs29 + crhs31 + 2*crhs32 - crhs35 - crhs38 - 4*crhs39 - 2*crhs4 - 2*crhs5 - 2*crhs6;
const double crhs66 =             DN(0,1)*U(0,3) + DN(1,1)*U(1,3) + DN(2,1)*U(2,3);
const double crhs67 =             crhs46*crhs50;
const double crhs68 =             crhs12*crhs15*crhs45;
const double crhs69 =             -crhs14*crhs54 + crhs16*crhs53*(-2*crhs58 + crhs60) - crhs3*crhs46*crhs49 + crhs43*crhs66 - crhs43*crhs68 + crhs47*crhs51 + crhs67;
const double crhs70 =             N[0]*(U(0,3)*bdf0 + Un(0,3)*bdf1 + Unn(0,3)*bdf2) + N[1]*(U(1,3)*bdf0 + Un(1,3)*bdf1 + Unn(1,3)*bdf2) + N[2]*(U(2,3)*bdf0 + Un(2,3)*bdf1 + Unn(2,3)*bdf2);
const double crhs71 =             crhs15*crhs8 + crhs17*crhs62 + crhs9*(N[0]*r[0] + N[1]*r[1] + N[2]*r[2]);
const double crhs72 =             3*(lambda/cv);
const double crhs73 =             3*mu;
const double crhs74 =             crhs12*crhs17*crhs45;
const double crhs75 =             -crhs72;
const double crhs76 =             crhs73 + crhs75;
const double crhs77 =             4*mu;
const double crhs78 =             crhs75 + crhs77;
const double crhs79 =             crhs15*crhs17*crhs53*mu;
const double crhs80 =             crhs12*crhs55;
const double crhs81 =             crhs12*crhs58;
const double crhs82 =             N[0]*U(0,3);
const double crhs83 =             N[1]*U(1,3);
const double crhs84 =             N[2]*U(2,3);
const double crhs85 =             crhs82 + crhs83 + crhs84;
const double crhs86 =             -crhs72*crhs80 - crhs72*crhs81 + crhs72*crhs85;
const double crhs87 =             -crhs12*crhs14*(crhs73*crhs81 + crhs77*crhs80 + crhs86) - crhs16*crhs79 + crhs36*crhs50*crhs78 + crhs44*crhs72 - crhs48*mu + crhs52*crhs76 + crhs73*crhs74;
const double crhs88 =             (1.0L/3.0L)*crhs12*crhs87;
const double crhs89 =             -2*crhs82;
const double crhs90 =             -2*crhs83;
const double crhs91 =             -2*crhs84;
const double crhs92 =             -crhs43*crhs85;
const double crhs93 =             crhs12*crhs59;
const double crhs94 =             crhs12*crhs57;
const double crhs95 =             crhs89 + crhs90 + crhs91 + crhs92 + crhs94;
const double crhs96 =             crhs12*crhs56*(crhs55 + crhs58) + crhs93 + crhs95;
const double crhs97 =             crhs12*(-crhs12*crhs15*crhs17*crhs43*crhs51 + crhs15*crhs42*crhs44 - crhs15*crhs43*crhs74 + crhs17*crhs42*crhs66 - crhs3*(3*crhs93 + crhs95) + crhs32*crhs96 + crhs39*crhs96 - crhs50*(crhs89 + crhs90 + crhs91 + crhs92 + crhs93 + 3*crhs94));
const double crhs98 =             (1.0L/3.0L)*crhs12*(-crhs12*crhs16*(crhs73*crhs80 + crhs77*crhs81 + crhs86) - crhs14*crhs79 + crhs3*crhs33*crhs78 + 3*crhs36*crhs51*mu + crhs66*crhs72 - crhs67*mu + crhs68*crhs76 + crhs87);
const double crhs99 =             N[1]*crhs9;
const double crhs100 =             DN(1,0)*crhs12*mu;
const double crhs101 =             (1.0L/3.0L)*DN(1,1)*crhs12*mu;
const double crhs102 =             (1.0L/2.0L)*N[1];
const double crhs103 =             N[2]*crhs9;
const double crhs104 =             DN(2,0)*crhs12*mu;
const double crhs105 =             (1.0L/3.0L)*DN(2,1)*crhs12*mu;
const double crhs106 =             (1.0L/2.0L)*N[2];
            rhs[0]=N[0]*crhs7;
            rhs[1]=N[0]*crhs11 - crhs10*crhs8 + crhs13*crhs18 + crhs19*crhs40 + crhs41*crhs61;
            rhs[2]=N[0]*crhs63 - crhs10*crhs62 + crhs13*crhs64 - crhs19*crhs65 + crhs41*crhs69;
            rhs[3]=-DN(0,0)*crhs88 - DN(0,1)*crhs98 + N[0]*crhs70 - N[0]*crhs71 + crhs41*crhs97;
            rhs[4]=N[1]*crhs7;
            rhs[5]=N[1]*crhs11 + crhs100*crhs18 + crhs101*crhs40 + crhs102*crhs61 - crhs8*crhs99;
            rhs[6]=N[1]*crhs63 + crhs100*crhs64 - crhs101*crhs65 + crhs102*crhs69 - crhs62*crhs99;
            rhs[7]=-DN(1,0)*crhs88 - DN(1,1)*crhs98 + N[1]*crhs70 - N[1]*crhs71 + crhs102*crhs97;
            rhs[8]=N[2]*crhs7;
            rhs[9]=N[2]*crhs11 - crhs103*crhs8 + crhs104*crhs18 + crhs105*crhs40 + crhs106*crhs61;
            rhs[10]=N[2]*crhs63 - crhs103*crhs62 + crhs104*crhs64 - crhs105*crhs65 + crhs106*crhs69;
            rhs[11]=-DN(2,0)*crhs88 - DN(2,1)*crhs98 + N[2]*crhs70 - N[2]*crhs71 + crhs106*crhs97;

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
