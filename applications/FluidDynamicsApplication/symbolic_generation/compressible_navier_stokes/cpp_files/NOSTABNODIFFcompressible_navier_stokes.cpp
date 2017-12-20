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

    const double clhs0 =             DN(0,0)*N[0];
const double clhs1 =             DN(0,1)*N[0];
const double clhs2 =             N[0]*bdf0;
const double clhs3 =             -N[1]*clhs2;
const double clhs4 =             DN(1,0)*N[0];
const double clhs5 =             DN(1,1)*N[0];
const double clhs6 =             -N[2]*clhs2;
const double clhs7 =             DN(2,0)*N[0];
const double clhs8 =             DN(2,1)*N[0];
const double clhs9 =             N[0]*f_ext(0,0) + N[1]*f_ext(1,0) + N[2]*f_ext(2,0);
const double clhs10 =             N[0]*clhs9;
const double clhs11 =             N[0]*U(0,0) + N[1]*U(1,0) + N[2]*U(2,0);
const double clhs12 =             pow(clhs11, -2);
const double clhs13 =             (1.0L/2.0L)*clhs12;
const double clhs14 =             N[0]*U(0,1) + N[1]*U(1,1) + N[2]*U(2,1);
const double clhs15 =             DN(0,1)*clhs14;
const double clhs16 =             2*N[0]*U(0,2) + 2*N[1]*U(1,2) + 2*N[2]*U(2,2);
const double clhs17 =             DN(0,1)*U(0,1) + DN(1,1)*U(1,1) + DN(2,1)*U(2,1);
const double clhs18 =             N[0]*clhs17;
const double clhs19 =             DN(0,1)*U(0,2) + DN(1,1)*U(1,2) + DN(2,1)*U(2,2);
const double clhs20 =             N[0]*clhs19;
const double clhs21 =             2*N[0]*U(0,1) + 2*N[1]*U(1,1) + 2*N[2]*U(2,1);
const double clhs22 =             y - 3;
const double clhs23 =             DN(0,0)*U(0,1) + DN(1,0)*U(1,1) + DN(2,0)*U(2,1);
const double clhs24 =             N[0]*clhs23;
const double clhs25 =             clhs22*clhs24;
const double clhs26 =             DN(0,0)*U(0,2) + DN(1,0)*U(1,2) + DN(2,0)*U(2,2);
const double clhs27 =             N[0]*clhs26;
const double clhs28 =             2*y - 2;
const double clhs29 =             N[0]*U(0,2) + N[1]*U(1,2) + N[2]*U(2,2);
const double clhs30 =             clhs28*clhs29;
const double clhs31 =             4*N[0]*U(0,2) + 4*N[1]*U(1,2) + 4*N[2]*U(2,2);
const double clhs32 =             DN(0,1)*U(0,0) + DN(1,1)*U(1,0) + DN(2,1)*U(2,0);
const double clhs33 =             1.0/clhs11;
const double clhs34 =             N[0]*clhs32*clhs33;
const double clhs35 =             clhs14*clhs34;
const double clhs36 =             pow(clhs14, 2);
const double clhs37 =             y - 1;
const double clhs38 =             clhs36*clhs37;
const double clhs39 =             pow(clhs29, 2);
const double clhs40 =             clhs37*clhs39;
const double clhs41 =             clhs38 + clhs40;
const double clhs42 =             -2*clhs36 + clhs41;
const double clhs43 =             DN(0,0)*U(0,0) + DN(1,0)*U(1,0) + DN(2,0)*U(2,0);
const double clhs44 =             N[0]*clhs33*clhs43;
const double clhs45 =             2*clhs36*clhs37;
const double clhs46 =             2*clhs37*clhs39;
const double clhs47 =             -4*clhs36 + clhs45 + clhs46;
const double clhs48 =             clhs10 - clhs13*(DN(0,0)*clhs42 - clhs15*clhs16 - clhs16*clhs18 - clhs20*clhs21 + clhs21*clhs25 + clhs27*clhs30 + clhs31*clhs35 - clhs44*clhs47);
const double clhs49 =             DN(0,1)*clhs29;
const double clhs50 =             clhs29*clhs34;
const double clhs51 =             clhs20 + clhs49 - clhs50;
const double clhs52 =             DN(0,0)*clhs14;
const double clhs53 =             clhs14*clhs44;
const double clhs54 =             clhs2 + clhs33*(-clhs22*clhs52 + clhs22*clhs53 - clhs25 + clhs51);
const double clhs55 =             N[0]*clhs33;
const double clhs56 =             DN(0,0)*clhs29;
const double clhs57 =             clhs27*clhs37;
const double clhs58 =             clhs29*clhs44;
const double clhs59 =             -clhs15 - clhs18 + clhs35 + clhs37*clhs56 - clhs37*clhs58 + clhs57;
const double clhs60 =             N[1]*clhs9;
const double clhs61 =             DN(1,1)*clhs14;
const double clhs62 =             N[1]*clhs17;
const double clhs63 =             N[1]*clhs19;
const double clhs64 =             N[1]*clhs23;
const double clhs65 =             clhs22*clhs64;
const double clhs66 =             N[1]*clhs26;
const double clhs67 =             N[1]*clhs32*clhs33;
const double clhs68 =             clhs14*clhs67;
const double clhs69 =             N[1]*clhs33*clhs43;
const double clhs70 =             -clhs13*(DN(1,0)*clhs42 - clhs16*clhs61 - clhs16*clhs62 - clhs21*clhs63 + clhs21*clhs65 + clhs30*clhs66 + clhs31*clhs68 - clhs47*clhs69) + clhs60;
const double clhs71 =             N[1]*bdf0;
const double clhs72 =             DN(1,1)*clhs29;
const double clhs73 =             clhs29*clhs67;
const double clhs74 =             clhs63 + clhs72 - clhs73;
const double clhs75 =             DN(1,0)*clhs14;
const double clhs76 =             clhs14*clhs69;
const double clhs77 =             clhs33*(-clhs22*clhs75 + clhs22*clhs76 - clhs65 + clhs74) + clhs71;
const double clhs78 =             DN(1,0)*clhs29;
const double clhs79 =             clhs37*clhs66;
const double clhs80 =             clhs29*clhs69;
const double clhs81 =             clhs37*clhs78 - clhs37*clhs80 - clhs61 - clhs62 + clhs68 + clhs79;
const double clhs82 =             N[2]*clhs9;
const double clhs83 =             DN(2,1)*clhs14;
const double clhs84 =             N[2]*clhs17;
const double clhs85 =             N[2]*clhs19;
const double clhs86 =             N[2]*clhs23;
const double clhs87 =             clhs22*clhs86;
const double clhs88 =             N[2]*clhs26;
const double clhs89 =             N[2]*clhs32*clhs33;
const double clhs90 =             clhs14*clhs89;
const double clhs91 =             N[2]*clhs33*clhs43;
const double clhs92 =             -clhs13*(DN(2,0)*clhs42 - clhs16*clhs83 - clhs16*clhs84 - clhs21*clhs85 + clhs21*clhs87 + clhs30*clhs88 + clhs31*clhs90 - clhs47*clhs91) + clhs82;
const double clhs93 =             N[2]*bdf0;
const double clhs94 =             DN(2,1)*clhs29;
const double clhs95 =             clhs29*clhs89;
const double clhs96 =             clhs85 + clhs94 - clhs95;
const double clhs97 =             DN(2,0)*clhs14;
const double clhs98 =             clhs14*clhs91;
const double clhs99 =             clhs33*(-clhs22*clhs97 + clhs22*clhs98 - clhs87 + clhs96) + clhs93;
const double clhs100 =             DN(2,0)*clhs29;
const double clhs101 =             clhs37*clhs88;
const double clhs102 =             clhs29*clhs91;
const double clhs103 =             clhs100*clhs37 + clhs101 - clhs102*clhs37 - clhs83 - clhs84 + clhs90;
const double clhs104 =             N[0]*f_ext(0,1) + N[1]*f_ext(1,1) + N[2]*f_ext(2,1);
const double clhs105 =             N[0]*clhs104;
const double clhs106 =             clhs20*clhs22;
const double clhs107 =             clhs14*clhs28;
const double clhs108 =             -2*clhs39 + clhs41;
const double clhs109 =             -4*clhs39 + clhs45 + clhs46;
const double clhs110 =             clhs105 + clhs13*(-DN(0,1)*clhs108 - clhs106*clhs16 - clhs107*clhs18 + clhs109*clhs34 + clhs16*clhs24 + clhs16*clhs52 + clhs21*clhs27 - clhs31*clhs53);
const double clhs111 =             clhs18*clhs37;
const double clhs112 =             -clhs111 - clhs15*clhs37 + clhs27 + clhs35*clhs37 + clhs56 - clhs58;
const double clhs113 =             clhs24 + clhs52 - clhs53;
const double clhs114 =             clhs2 + clhs33*(-clhs106 + clhs113 - clhs22*clhs49 + clhs22*clhs50);
const double clhs115 =             N[1]*clhs104;
const double clhs116 =             clhs22*clhs63;
const double clhs117 =             clhs115 + clhs13*(-DN(1,1)*clhs108 - clhs107*clhs62 + clhs109*clhs67 - clhs116*clhs16 + clhs16*clhs64 + clhs16*clhs75 + clhs21*clhs66 - clhs31*clhs76);
const double clhs118 =             clhs37*clhs62;
const double clhs119 =             -clhs118 - clhs37*clhs61 + clhs37*clhs68 + clhs66 + clhs78 - clhs80;
const double clhs120 =             clhs64 + clhs75 - clhs76;
const double clhs121 =             clhs33*(-clhs116 + clhs120 - clhs22*clhs72 + clhs22*clhs73) + clhs71;
const double clhs122 =             N[2]*clhs104;
const double clhs123 =             clhs22*clhs85;
const double clhs124 =             clhs122 + clhs13*(-DN(2,1)*clhs108 - clhs107*clhs84 + clhs109*clhs89 - clhs123*clhs16 + clhs16*clhs86 + clhs16*clhs97 + clhs21*clhs88 - clhs31*clhs98);
const double clhs125 =             clhs37*clhs84;
const double clhs126 =             clhs100 - clhs102 - clhs125 - clhs37*clhs83 + clhs37*clhs90 + clhs88;
const double clhs127 =             clhs86 + clhs97 - clhs98;
const double clhs128 =             clhs33*(-clhs123 + clhs127 - clhs22*clhs94 + clhs22*clhs95) + clhs93;
const double clhs129 =             N[0]*r[0] + N[1]*r[1] + N[2]*r[2];
const double clhs130 =             2*N[0]*y;
const double clhs131 =             DN(0,0)*U(0,3) + DN(1,0)*U(1,3) + DN(2,0)*U(2,3);
const double clhs132 =             clhs130*clhs131;
const double clhs133 =             DN(0,1)*U(0,3) + DN(1,1)*U(1,3) + DN(2,1)*U(2,3);
const double clhs134 =             clhs130*clhs133;
const double clhs135 =             clhs14*clhs31*clhs33;
const double clhs136 =             clhs33*clhs40;
const double clhs137 =             N[0]*U(0,3);
const double clhs138 =             -2*clhs137;
const double clhs139 =             N[1]*U(1,3);
const double clhs140 =             -2*clhs139;
const double clhs141 =             N[2]*U(2,3);
const double clhs142 =             -2*clhs141;
const double clhs143 =             -clhs28*(clhs137 + clhs139 + clhs141);
const double clhs144 =             clhs33*clhs38;
const double clhs145 =             clhs33*(clhs36 + clhs39);
const double clhs146 =             clhs145*clhs37;
const double clhs147 =             clhs138 + clhs140 + clhs142 + clhs143 + clhs144 + clhs146;
const double clhs148 =             clhs136 + clhs147;
const double clhs149 =             clhs136 + clhs138 + clhs140 + clhs142 + clhs143;
const double clhs150 =             5*clhs144 + clhs146 + clhs149;
const double clhs151 =             5*clhs136 + clhs147;
const double clhs152 =             clhs138 + clhs140 + clhs142 + clhs143 + clhs144;
const double clhs153 =             clhs136 + clhs145*clhs28 + clhs152;
const double clhs154 =             clhs153*clhs21;
const double clhs155 =             clhs153*clhs16;
const double clhs156 =             N[0]*clhs129 - clhs13*(clhs111*clhs135 - clhs132*clhs14 - clhs134*clhs29 + clhs135*clhs57 + clhs148*clhs49 + clhs148*clhs52 + clhs150*clhs24 + clhs151*clhs20 - clhs154*clhs44 - clhs155*clhs34);
const double clhs157 =             (1.0L/2.0L)*clhs33;
const double clhs158 =             clhs28*clhs29*clhs33;
const double clhs159 =             6*clhs14*clhs33*clhs37;
const double clhs160 =             clhs14*clhs28*clhs33;
const double clhs161 =             N[0]*clhs12*clhs14*clhs31*clhs37;
const double clhs162 =             3*clhs144 + clhs149;
const double clhs163 =             clhs10 + clhs157*(DN(0,0)*clhs162 + N[0]*clhs26*clhs28*clhs29*clhs33 - clhs132 + clhs15*clhs158 - clhs150*clhs44 + clhs158*clhs18 + clhs159*clhs24 + clhs160*clhs20 - clhs161*clhs32);
const double clhs164 =             6*clhs29*clhs33*clhs37;
const double clhs165 =             3*clhs136 + clhs152;
const double clhs166 =             clhs105 + clhs157*(DN(0,1)*clhs165 + N[0]*clhs14*clhs17*clhs28*clhs33 - clhs134 - clhs151*clhs34 + clhs158*clhs24 + clhs158*clhs52 + clhs160*clhs27 - clhs161*clhs43 + clhs164*clhs20);
const double clhs167 =             clhs33*y;
const double clhs168 =             clhs167*(clhs113 + clhs51) + clhs2;
const double clhs169 =             2*N[1]*y;
const double clhs170 =             clhs131*clhs169;
const double clhs171 =             clhs133*clhs169;
const double clhs172 =             N[1]*clhs129 - clhs13*(clhs118*clhs135 + clhs135*clhs79 - clhs14*clhs170 + clhs148*clhs72 + clhs148*clhs75 + clhs150*clhs64 + clhs151*clhs63 - clhs154*clhs69 - clhs155*clhs67 - clhs171*clhs29);
const double clhs173 =             N[1]*clhs12*clhs14*clhs31*clhs37;
const double clhs174 =             clhs157*(DN(1,0)*clhs162 + N[1]*clhs26*clhs28*clhs29*clhs33 - clhs150*clhs69 + clhs158*clhs61 + clhs158*clhs62 + clhs159*clhs64 + clhs160*clhs63 - clhs170 - clhs173*clhs32) + clhs60;
const double clhs175 =             clhs115 + clhs157*(DN(1,1)*clhs165 + N[1]*clhs14*clhs17*clhs28*clhs33 - clhs151*clhs67 + clhs158*clhs64 + clhs158*clhs75 + clhs160*clhs66 + clhs164*clhs63 - clhs171 - clhs173*clhs43);
const double clhs176 =             clhs167*(clhs120 + clhs74) + clhs71;
const double clhs177 =             2*N[2]*y;
const double clhs178 =             clhs131*clhs177;
const double clhs179 =             clhs133*clhs177;
const double clhs180 =             N[2]*clhs129 - clhs13*(clhs101*clhs135 + clhs125*clhs135 - clhs14*clhs178 + clhs148*clhs94 + clhs148*clhs97 + clhs150*clhs86 + clhs151*clhs85 - clhs154*clhs91 - clhs155*clhs89 - clhs179*clhs29);
const double clhs181 =             N[2]*clhs12*clhs14*clhs31*clhs37;
const double clhs182 =             clhs157*(DN(2,0)*clhs162 + N[2]*clhs26*clhs28*clhs29*clhs33 - clhs150*clhs91 + clhs158*clhs83 + clhs158*clhs84 + clhs159*clhs86 + clhs160*clhs85 - clhs178 - clhs181*clhs32) + clhs82;
const double clhs183 =             clhs122 + clhs157*(DN(2,1)*clhs165 + N[2]*clhs14*clhs17*clhs28*clhs33 - clhs151*clhs89 + clhs158*clhs86 + clhs158*clhs97 + clhs160*clhs88 + clhs164*clhs85 - clhs179 - clhs181*clhs43);
const double clhs184 =             clhs167*(clhs127 + clhs96) + clhs93;
const double clhs185 =             DN(0,0)*N[1];
const double clhs186 =             DN(0,1)*N[1];
const double clhs187 =             DN(1,0)*N[1];
const double clhs188 =             DN(1,1)*N[1];
const double clhs189 =             -N[2]*clhs71;
const double clhs190 =             DN(2,0)*N[1];
const double clhs191 =             DN(2,1)*N[1];
const double clhs192 =             N[1]*clhs33;
const double clhs193 =             DN(0,0)*N[2];
const double clhs194 =             DN(0,1)*N[2];
const double clhs195 =             DN(1,0)*N[2];
const double clhs196 =             DN(1,1)*N[2];
const double clhs197 =             DN(2,0)*N[2];
const double clhs198 =             DN(2,1)*N[2];
const double clhs199 =             N[2]*clhs33;
            lhs(0,0)=-pow(N[0], 2)*bdf0;
            lhs(0,1)=-clhs0;
            lhs(0,2)=-clhs1;
            lhs(0,3)=0;
            lhs(0,4)=clhs3;
            lhs(0,5)=-clhs4;
            lhs(0,6)=-clhs5;
            lhs(0,7)=0;
            lhs(0,8)=clhs6;
            lhs(0,9)=-clhs7;
            lhs(0,10)=-clhs8;
            lhs(0,11)=0;
            lhs(1,0)=N[0]*clhs48;
            lhs(1,1)=-N[0]*clhs54;
            lhs(1,2)=clhs55*clhs59;
            lhs(1,3)=-clhs0*clhs37;
            lhs(1,4)=N[0]*clhs70;
            lhs(1,5)=-N[0]*clhs77;
            lhs(1,6)=clhs55*clhs81;
            lhs(1,7)=-clhs37*clhs4;
            lhs(1,8)=N[0]*clhs92;
            lhs(1,9)=-N[0]*clhs99;
            lhs(1,10)=clhs103*clhs55;
            lhs(1,11)=-clhs37*clhs7;
            lhs(2,0)=N[0]*clhs110;
            lhs(2,1)=-clhs112*clhs55;
            lhs(2,2)=-N[0]*clhs114;
            lhs(2,3)=-clhs1*clhs37;
            lhs(2,4)=N[0]*clhs117;
            lhs(2,5)=-clhs119*clhs55;
            lhs(2,6)=-N[0]*clhs121;
            lhs(2,7)=-clhs37*clhs5;
            lhs(2,8)=N[0]*clhs124;
            lhs(2,9)=-clhs126*clhs55;
            lhs(2,10)=-N[0]*clhs128;
            lhs(2,11)=-clhs37*clhs8;
            lhs(3,0)=N[0]*clhs156;
            lhs(3,1)=N[0]*clhs163;
            lhs(3,2)=N[0]*clhs166;
            lhs(3,3)=-N[0]*clhs168;
            lhs(3,4)=N[0]*clhs172;
            lhs(3,5)=N[0]*clhs174;
            lhs(3,6)=N[0]*clhs175;
            lhs(3,7)=-N[0]*clhs176;
            lhs(3,8)=N[0]*clhs180;
            lhs(3,9)=N[0]*clhs182;
            lhs(3,10)=N[0]*clhs183;
            lhs(3,11)=-N[0]*clhs184;
            lhs(4,0)=clhs3;
            lhs(4,1)=-clhs185;
            lhs(4,2)=-clhs186;
            lhs(4,3)=0;
            lhs(4,4)=-pow(N[1], 2)*bdf0;
            lhs(4,5)=-clhs187;
            lhs(4,6)=-clhs188;
            lhs(4,7)=0;
            lhs(4,8)=clhs189;
            lhs(4,9)=-clhs190;
            lhs(4,10)=-clhs191;
            lhs(4,11)=0;
            lhs(5,0)=N[1]*clhs48;
            lhs(5,1)=-N[1]*clhs54;
            lhs(5,2)=clhs192*clhs59;
            lhs(5,3)=-clhs185*clhs37;
            lhs(5,4)=N[1]*clhs70;
            lhs(5,5)=-N[1]*clhs77;
            lhs(5,6)=clhs192*clhs81;
            lhs(5,7)=-clhs187*clhs37;
            lhs(5,8)=N[1]*clhs92;
            lhs(5,9)=-N[1]*clhs99;
            lhs(5,10)=clhs103*clhs192;
            lhs(5,11)=-clhs190*clhs37;
            lhs(6,0)=N[1]*clhs110;
            lhs(6,1)=-clhs112*clhs192;
            lhs(6,2)=-N[1]*clhs114;
            lhs(6,3)=-clhs186*clhs37;
            lhs(6,4)=N[1]*clhs117;
            lhs(6,5)=-clhs119*clhs192;
            lhs(6,6)=-N[1]*clhs121;
            lhs(6,7)=-clhs188*clhs37;
            lhs(6,8)=N[1]*clhs124;
            lhs(6,9)=-clhs126*clhs192;
            lhs(6,10)=-N[1]*clhs128;
            lhs(6,11)=-clhs191*clhs37;
            lhs(7,0)=N[1]*clhs156;
            lhs(7,1)=N[1]*clhs163;
            lhs(7,2)=N[1]*clhs166;
            lhs(7,3)=-N[1]*clhs168;
            lhs(7,4)=N[1]*clhs172;
            lhs(7,5)=N[1]*clhs174;
            lhs(7,6)=N[1]*clhs175;
            lhs(7,7)=-N[1]*clhs176;
            lhs(7,8)=N[1]*clhs180;
            lhs(7,9)=N[1]*clhs182;
            lhs(7,10)=N[1]*clhs183;
            lhs(7,11)=-N[1]*clhs184;
            lhs(8,0)=clhs6;
            lhs(8,1)=-clhs193;
            lhs(8,2)=-clhs194;
            lhs(8,3)=0;
            lhs(8,4)=clhs189;
            lhs(8,5)=-clhs195;
            lhs(8,6)=-clhs196;
            lhs(8,7)=0;
            lhs(8,8)=-pow(N[2], 2)*bdf0;
            lhs(8,9)=-clhs197;
            lhs(8,10)=-clhs198;
            lhs(8,11)=0;
            lhs(9,0)=N[2]*clhs48;
            lhs(9,1)=-N[2]*clhs54;
            lhs(9,2)=clhs199*clhs59;
            lhs(9,3)=-clhs193*clhs37;
            lhs(9,4)=N[2]*clhs70;
            lhs(9,5)=-N[2]*clhs77;
            lhs(9,6)=clhs199*clhs81;
            lhs(9,7)=-clhs195*clhs37;
            lhs(9,8)=N[2]*clhs92;
            lhs(9,9)=-N[2]*clhs99;
            lhs(9,10)=clhs103*clhs199;
            lhs(9,11)=-clhs197*clhs37;
            lhs(10,0)=N[2]*clhs110;
            lhs(10,1)=-clhs112*clhs199;
            lhs(10,2)=-N[2]*clhs114;
            lhs(10,3)=-clhs194*clhs37;
            lhs(10,4)=N[2]*clhs117;
            lhs(10,5)=-clhs119*clhs199;
            lhs(10,6)=-N[2]*clhs121;
            lhs(10,7)=-clhs196*clhs37;
            lhs(10,8)=N[2]*clhs124;
            lhs(10,9)=-clhs126*clhs199;
            lhs(10,10)=-N[2]*clhs128;
            lhs(10,11)=-clhs198*clhs37;
            lhs(11,0)=N[2]*clhs156;
            lhs(11,1)=N[2]*clhs163;
            lhs(11,2)=N[2]*clhs166;
            lhs(11,3)=-N[2]*clhs168;
            lhs(11,4)=N[2]*clhs172;
            lhs(11,5)=N[2]*clhs174;
            lhs(11,6)=N[2]*clhs175;
            lhs(11,7)=-N[2]*clhs176;
            lhs(11,8)=N[2]*clhs180;
            lhs(11,9)=N[2]*clhs182;
            lhs(11,10)=N[2]*clhs183;
            lhs(11,11)=-N[2]*clhs184;


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

    const double crhs0 =             DN(0,1)*U(0,2) + DN(1,1)*U(1,2) + DN(2,1)*U(2,2);
const double crhs1 =             DN(0,0)*U(0,1);
const double crhs2 =             DN(1,0)*U(1,1);
const double crhs3 =             DN(2,0)*U(2,1);
const double crhs4 =             N[0]*(U(0,0)*bdf0 + Un(0,0)*bdf1 + Unn(0,0)*bdf2) + N[1]*(U(1,0)*bdf0 + Un(1,0)*bdf1 + Unn(1,0)*bdf2) + N[2]*(U(2,0)*bdf0 + Un(2,0)*bdf1 + Unn(2,0)*bdf2) + crhs0 + crhs1 + crhs2 + crhs3;
const double crhs5 =             y - 1;
const double crhs6 =             DN(0,0)*U(0,3) + DN(1,0)*U(1,3) + DN(2,0)*U(2,3);
const double crhs7 =             N[0]*U(0,0) + N[1]*U(1,0) + N[2]*U(2,0);
const double crhs8 =             N[0]*f_ext(0,0) + N[1]*f_ext(1,0) + N[2]*f_ext(2,0);
const double crhs9 =             DN(0,1)*U(0,1) + DN(1,1)*U(1,1) + DN(2,1)*U(2,1);
const double crhs10 =             N[0]*U(0,2) + N[1]*U(1,2) + N[2]*U(2,2);
const double crhs11 =             1.0/crhs7;
const double crhs12 =             crhs10*crhs11;
const double crhs13 =             crhs12*crhs9;
const double crhs14 =             N[0]*U(0,1) + N[1]*U(1,1) + N[2]*U(2,1);
const double crhs15 =             crhs11*crhs14;
const double crhs16 =             y - 3;
const double crhs17 =             crhs1 + crhs2 + crhs3;
const double crhs18 =             DN(0,0)*U(0,2) + DN(1,0)*U(1,2) + DN(2,0)*U(2,2);
const double crhs19 =             crhs10*crhs11*crhs18;
const double crhs20 =             DN(0,1)*U(0,0) + DN(1,1)*U(1,0) + DN(2,1)*U(2,0);
const double crhs21 =             pow(crhs7, -2);
const double crhs22 =             crhs10*crhs14*crhs21;
const double crhs23 =             (1.0L/2.0L)*crhs21;
const double crhs24 =             DN(0,0)*U(0,0) + DN(1,0)*U(1,0) + DN(2,0)*U(2,0);
const double crhs25 =             pow(crhs14, 2);
const double crhs26 =             crhs25*crhs5;
const double crhs27 =             pow(crhs10, 2);
const double crhs28 =             crhs27*crhs5;
const double crhs29 =             crhs26 + crhs28;
const double crhs30 =             N[0]*(U(0,1)*bdf0 + Un(0,1)*bdf1 + Unn(0,1)*bdf2) + N[1]*(U(1,1)*bdf0 + Un(1,1)*bdf1 + Unn(1,1)*bdf2) + N[2]*(U(2,1)*bdf0 + Un(2,1)*bdf1 + Unn(2,1)*bdf2) + crhs0*crhs15 + crhs13 - crhs15*crhs16*crhs17 - crhs19*crhs5 - crhs20*crhs22 + crhs23*crhs24*(-2*crhs25 + crhs29) + crhs5*crhs6 - crhs7*crhs8;
const double crhs31 =             DN(0,1)*U(0,3) + DN(1,1)*U(1,3) + DN(2,1)*U(2,3);
const double crhs32 =             N[0]*f_ext(0,1) + N[1]*f_ext(1,1) + N[2]*f_ext(2,1);
const double crhs33 =             N[0]*(U(0,2)*bdf0 + Un(0,2)*bdf1 + Unn(0,2)*bdf2) + N[1]*(U(1,2)*bdf0 + Un(1,2)*bdf1 + Unn(1,2)*bdf2) + N[2]*(U(2,2)*bdf0 + Un(2,2)*bdf1 + Unn(2,2)*bdf2) - crhs0*crhs12*crhs16 + crhs12*crhs17 + crhs15*crhs18 - crhs15*crhs5*crhs9 + crhs20*crhs23*(-2*crhs27 + crhs29) - crhs22*crhs24 + crhs31*crhs5 - crhs32*crhs7;
const double crhs34 =             2*y;
const double crhs35 =             crhs34 - 2;
const double crhs36 =             crhs14*crhs35;
const double crhs37 =             N[0]*U(0,3);
const double crhs38 =             -2*crhs37;
const double crhs39 =             N[1]*U(1,3);
const double crhs40 =             -2*crhs39;
const double crhs41 =             N[2]*U(2,3);
const double crhs42 =             -2*crhs41;
const double crhs43 =             -crhs35*(crhs37 + crhs39 + crhs41);
const double crhs44 =             crhs11*crhs28;
const double crhs45 =             crhs11*crhs26;
const double crhs46 =             crhs38 + crhs40 + crhs42 + crhs43 + crhs45;
const double crhs47 =             crhs11*crhs5*(crhs25 + crhs27) + crhs44 + crhs46;
const double crhs48 =             N[0]*(U(0,3)*bdf0 + Un(0,3)*bdf1 + Unn(0,3)*bdf2) + N[1]*(U(1,3)*bdf0 + Un(1,3)*bdf1 + Unn(1,3)*bdf2) + N[2]*(U(2,3)*bdf0 + Un(2,3)*bdf1 + Unn(2,3)*bdf2) - crhs10*crhs32 + (1.0L/2.0L)*crhs11*(-crhs0*(3*crhs44 + crhs46) + crhs10*crhs31*crhs34 + crhs12*crhs20*crhs47 - crhs13*crhs36 + crhs14*crhs34*crhs6 + crhs15*crhs24*crhs47 - crhs17*(crhs38 + crhs40 + crhs42 + crhs43 + crhs44 + 3*crhs45) - crhs19*crhs36) - crhs14*crhs8 - crhs7*(N[0]*r[0] + N[1]*r[1] + N[2]*r[2]);
            rhs[0]=N[0]*crhs4;
            rhs[1]=N[0]*crhs30;
            rhs[2]=N[0]*crhs33;
            rhs[3]=N[0]*crhs48;
            rhs[4]=N[1]*crhs4;
            rhs[5]=N[1]*crhs30;
            rhs[6]=N[1]*crhs33;
            rhs[7]=N[1]*crhs48;
            rhs[8]=N[2]*crhs4;
            rhs[9]=N[2]*crhs30;
            rhs[10]=N[2]*crhs33;
            rhs[11]=N[2]*crhs48;

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
