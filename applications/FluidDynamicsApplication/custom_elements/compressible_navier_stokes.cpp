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
    const double stab_c1 = 4.0;
    const double stab_c2 = 2.0;
   
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
    const double stab_c1 = 4.0;
    const double stab_c2 = 2.0;

    const double clhs0 =             pow(N[0], 2);
const double clhs1 =             1.0/stab_c2;
const double clhs2 =             N[0]*U(0,0);
const double clhs3 =             N[1]*U(1,0);
const double clhs4 =             N[2]*U(2,0);
const double clhs5 =             clhs2 + clhs3 + clhs4;
const double clhs6 =             1.0/clhs5;
const double clhs7 =             N[0]*U(0,1);
const double clhs8 =             N[1]*U(1,1);
const double clhs9 =             N[2]*U(2,1);
const double clhs10 =             clhs7 + clhs8 + clhs9;
const double clhs11 =             N[0]*U(0,2);
const double clhs12 =             N[1]*U(1,2);
const double clhs13 =             N[2]*U(2,2);
const double clhs14 =             clhs11 + clhs12 + clhs13;
const double clhs15 =             fabs(clhs6*(clhs10 + clhs14));
const double clhs16 =             y - 1;
const double clhs17 =             N[0]*U(0,3);
const double clhs18 =             2*clhs17;
const double clhs19 =             -clhs18;
const double clhs20 =             N[1]*U(1,3);
const double clhs21 =             2*clhs20;
const double clhs22 =             -clhs21;
const double clhs23 =             N[2]*U(2,3);
const double clhs24 =             2*clhs23;
const double clhs25 =             -clhs24;
const double clhs26 =             pow(clhs10, 2);
const double clhs27 =             clhs26*clhs6;
const double clhs28 =             pow(clhs14, 2);
const double clhs29 =             clhs28*clhs6;
const double clhs30 =             clhs27 + clhs29;
const double clhs31 =             2*clhs15 - clhs16*clhs6*y*(clhs19 + clhs22 + clhs25 + clhs30);
const double clhs32 =             1.0/clhs31;
const double clhs33 =             DN(0,0)*U(0,1);
const double clhs34 =             DN(0,1)*U(0,2);
const double clhs35 =             DN(1,0)*U(1,1);
const double clhs36 =             DN(1,1)*U(1,2);
const double clhs37 =             DN(2,0)*U(2,1);
const double clhs38 =             DN(2,1)*U(2,2);
const double clhs39 =             3*N[0];
const double clhs40 =             3*N[1];
const double clhs41 =             3*N[2];
const double clhs42 =             DN(0,0)*U(0,0) + DN(1,0)*U(1,0) + DN(2,0)*U(2,0);
const double clhs43 =             pow(clhs42, 2);
const double clhs44 =             pow(clhs5, -2);
const double clhs45 =             7*clhs44*mu;
const double clhs46 =             DN(0,1)*U(0,0) + DN(1,1)*U(1,0) + DN(2,1)*U(2,0);
const double clhs47 =             pow(clhs46, 2);
const double clhs48 =             1.0/cv;
const double clhs49 =             3*clhs44*clhs48*lambda;
const double clhs50 =             clhs33 + clhs35 + clhs37;
const double clhs51 =             4*clhs44*clhs50*mu;
const double clhs52 =             DN(0,0)*U(0,2) + DN(1,0)*U(1,2) + DN(2,0)*U(2,2);
const double clhs53 =             3*clhs44*clhs52;
const double clhs54 =             2*clhs42*clhs44*mu;
const double clhs55 =             DN(0,1)*U(0,1) + DN(1,1)*U(1,1) + DN(2,1)*U(2,1);
const double clhs56 =             5*clhs42*clhs44*mu;
const double clhs57 =             clhs34 + clhs36 + clhs38;
const double clhs58 =             clhs44*clhs46*clhs50*mu;
const double clhs59 =             clhs44*clhs52*mu;
const double clhs60 =             clhs46*mu;
const double clhs61 =             3*clhs44*clhs55;
const double clhs62 =             4*clhs44*clhs57*mu;
const double clhs63 =             8*N[0]*U(0,1) + 8*N[1]*U(1,1) + 8*N[2]*U(2,1);
const double clhs64 =             pow(clhs5, -3);
const double clhs65 =             clhs43*clhs64*mu;
const double clhs66 =             6*N[0]*U(0,2) + 6*N[1]*U(1,2) + 6*N[2]*U(2,2);
const double clhs67 =             6*N[0]*U(0,1) + 6*N[1]*U(1,1) + 6*N[2]*U(2,1);
const double clhs68 =             clhs47*clhs64*mu;
const double clhs69 =             8*N[0]*U(0,2) + 8*N[1]*U(1,2) + 8*N[2]*U(2,2);
const double clhs70 =             DN(0,0)*U(0,3) + DN(1,0)*U(1,3) + DN(2,0)*U(2,3);
const double clhs71 =             clhs42*clhs44;
const double clhs72 =             DN(0,1)*U(0,3) + DN(1,1)*U(1,3) + DN(2,1)*U(2,3);
const double clhs73 =             clhs44*clhs46;
const double clhs74 =             clhs48*lambda;
const double clhs75 =             -clhs74 + mu;
const double clhs76 =             3*clhs44*clhs52*clhs75;
const double clhs77 =             3*clhs44*clhs55*clhs75;
const double clhs78 =             4*mu;
const double clhs79 =             3*clhs74;
const double clhs80 =             clhs78 - clhs79;
const double clhs81 =             clhs44*clhs50*clhs80;
const double clhs82 =             clhs44*clhs57*clhs80;
const double clhs83 =             clhs14*clhs55;
const double clhs84 =             clhs42*clhs64*mu;
const double clhs85 =             clhs10*clhs57;
const double clhs86 =             clhs84*clhs85;
const double clhs87 =             clhs14*clhs50;
const double clhs88 =             clhs46*clhs64*mu;
const double clhs89 =             clhs10*clhs52;
const double clhs90 =             clhs88*clhs89;
const double clhs91 =             clhs74 - mu;
const double clhs92 =             clhs64*clhs91;
const double clhs93 =             clhs64*clhs67*clhs91;
const double clhs94 =             clhs10*clhs43;
const double clhs95 =             -clhs78 + clhs79;
const double clhs96 =             2*clhs64*clhs95;
const double clhs97 =             clhs14*clhs47;
const double clhs98 =             clhs42*clhs52*clhs64*clhs91;
const double clhs99 =             clhs46*clhs55;
const double clhs100 =             pow(clhs5, -4);
const double clhs101 =             clhs100*clhs42*mu;
const double clhs102 =             clhs28*clhs6*mu;
const double clhs103 =             clhs18 + clhs21 + clhs24;
const double clhs104 =             clhs103*clhs74 - clhs27*clhs79 - clhs29*clhs79;
const double clhs105 =             3*clhs102 + clhs104 + clhs27*clhs78;
const double clhs106 =             clhs26*clhs6*mu;
const double clhs107 =             clhs104 + 3*clhs106 + clhs29*clhs78;
const double clhs108 =             clhs107*clhs64;
const double clhs109 =             -clhs10*clhs101*clhs46*clhs66 - clhs10*clhs42*clhs50*clhs96 - 3*clhs105*clhs43*clhs64 - 3*clhs108*clhs47 - clhs14*clhs46*clhs57*clhs96 + 3*clhs33 + 3*clhs34 + 3*clhs35 + 3*clhs36 + 3*clhs37 + 3*clhs38 + clhs39*(U(0,0)*bdf0 + Un(0,0)*bdf1 + Unn(0,0)*bdf2) + clhs40*(U(1,0)*bdf0 + Un(1,0)*bdf1 + Unn(1,0)*bdf2) + clhs41*(U(2,0)*bdf0 + Un(2,0)*bdf1 + Unn(2,0)*bdf2) + clhs42*clhs51 + clhs42*clhs53*mu - clhs42*clhs76 - clhs42*clhs81 + clhs43*clhs45 + clhs43*clhs49 - clhs43*clhs66*clhs92 + clhs45*clhs47 + clhs46*clhs54 + clhs46*clhs62 - clhs46*clhs77 - clhs46*clhs82 + clhs47*clhs49 - clhs47*clhs93 + 3*clhs48*clhs70*clhs71*lambda + 3*clhs48*clhs72*clhs73*lambda + clhs55*clhs56 - clhs56*clhs57 - 5*clhs58 + clhs59*(5*DN(0,1)*U(0,0) + 5*DN(1,1)*U(1,0) + 5*DN(2,1)*U(2,0)) + clhs60*clhs61 - clhs63*clhs65 - clhs65*clhs66 - clhs66*clhs98 - clhs67*clhs68 - clhs68*clhs69 + clhs83*clhs84 + clhs86 + clhs87*clhs88 + clhs90 - clhs93*clhs99 - clhs94*clhs96 - clhs96*clhs97;
const double clhs110 =             2*DN(0,0)*mu;
const double clhs111 =             pow(DN(0,0), 2);
const double clhs112 =             7*mu;
const double clhs113 =             pow(DN(0,1), 2);
const double clhs114 =             6*DN(0,0)*N[0];
const double clhs115 =             clhs10*clhs114;
const double clhs116 =             6*DN(0,1)*N[0];
const double clhs117 =             clhs116*clhs14;
const double clhs118 =             3*clhs0;
const double clhs119 =             clhs118*clhs50;
const double clhs120 =             clhs118*clhs52;
const double clhs121 =             clhs118*clhs55;
const double clhs122 =             clhs118*clhs57;
const double clhs123 =             3*clhs0*y;
const double clhs124 =             y - 3;
const double clhs125 =             6*DN(0,0)*clhs14*clhs16;
const double clhs126 =             6*DN(0,1)*clhs10*clhs16;
const double clhs127 =             14*DN(0,0)*mu;
const double clhs128 =             N[0]*clhs42*clhs6;
const double clhs129 =             8*DN(0,0)*mu;
const double clhs130 =             N[0]*clhs50*clhs6;
const double clhs131 =             6*DN(0,0)*N[0]*clhs52*clhs6;
const double clhs132 =             N[0]*clhs46*clhs6;
const double clhs133 =             N[0]*clhs55;
const double clhs134 =             10*DN(0,0)*clhs6*mu;
const double clhs135 =             N[0]*clhs57;
const double clhs136 =             2*DN(0,1)*mu;
const double clhs137 =             10*N[0]*clhs52;
const double clhs138 =             14*DN(0,1)*mu;
const double clhs139 =             6*DN(0,1)*N[0]*clhs55*clhs6;
const double clhs140 =             8*DN(0,1)*mu;
const double clhs141 =             clhs111*mu;
const double clhs142 =             6*clhs14*clhs6;
const double clhs143 =             clhs113*mu;
const double clhs144 =             6*clhs10*clhs6;
const double clhs145 =             6*DN(0,0)*clhs48*lambda;
const double clhs146 =             N[0]*clhs6;
const double clhs147 =             6*DN(0,1)*clhs48*lambda;
const double clhs148 =             12*DN(0,0)*clhs10*clhs14*clhs146;
const double clhs149 =             2*clhs6*clhs95;
const double clhs150 =             clhs0*clhs42;
const double clhs151 =             clhs144*clhs150;
const double clhs152 =             clhs142*clhs50;
const double clhs153 =             clhs144*clhs52;
const double clhs154 =             clhs0*clhs46;
const double clhs155 =             clhs142*clhs154;
const double clhs156 =             clhs142*clhs55;
const double clhs157 =             clhs144*clhs57;
const double clhs158 =             DN(0,1)*mu;
const double clhs159 =             6*DN(0,0);
const double clhs160 =             clhs10*clhs14*clhs44;
const double clhs161 =             clhs159*clhs160;
const double clhs162 =             24*N[0]*U(0,1) + 24*N[1]*U(1,1) + 24*N[2]*U(2,1);
const double clhs163 =             DN(0,0)*N[0]*clhs42*clhs44*mu;
const double clhs164 =             18*N[0]*U(0,2) + 18*N[1]*U(1,2) + 18*N[2]*U(2,2);
const double clhs165 =             30*DN(0,0)*N[0]*clhs44*clhs46*mu;
const double clhs166 =             DN(0,0)*clhs44*mu;
const double clhs167 =             3*N[0]*clhs14*clhs55;
const double clhs168 =             3*N[0]*clhs10*clhs57;
const double clhs169 =             30*DN(0,1)*N[0]*clhs42*clhs44*mu;
const double clhs170 =             DN(0,1)*clhs44*mu;
const double clhs171 =             3*N[0]*clhs14*clhs50;
const double clhs172 =             3*N[0]*clhs10*clhs52;
const double clhs173 =             18*N[0]*U(0,1) + 18*N[1]*U(1,1) + 18*N[2]*U(2,1);
const double clhs174 =             DN(0,1)*N[0]*clhs44*clhs46*mu;
const double clhs175 =             24*N[0]*U(0,2) + 24*N[1]*U(1,2) + 24*N[2]*U(2,2);
const double clhs176 =             clhs0*y;
const double clhs177 =             12*N[0]*clhs6;
const double clhs178 =             DN(0,1)*clhs10*clhs14*clhs16;
const double clhs179 =             6*y - 6;
const double clhs180 =             clhs14*clhs179*clhs6;
const double clhs181 =             18*clhs0*clhs6;
const double clhs182 =             clhs0*clhs179*clhs6;
const double clhs183 =             clhs10*clhs179*clhs6;
const double clhs184 =             -y + 3;
const double clhs185 =             clhs144*clhs184*clhs50;
const double clhs186 =             clhs142*clhs184*clhs57;
const double clhs187 =             18*N[0];
const double clhs188 =             clhs0*clhs10*clhs164;
const double clhs189 =             12*N[0]*U(0,1) + 12*N[1]*U(1,1) + 12*N[2]*U(2,1);
const double clhs190 =             clhs0*clhs14*clhs16*clhs173;
const double clhs191 =             clhs44*clhs52;
const double clhs192 =             clhs44*clhs55;
const double clhs193 =             2*clhs26;
const double clhs194 =             clhs16*clhs26;
const double clhs195 =             clhs16*clhs28;
const double clhs196 =             -clhs194 - clhs195;
const double clhs197 =             clhs193 + clhs196;
const double clhs198 =             clhs197*clhs6;
const double clhs199 =             2*clhs28;
const double clhs200 =             clhs196 + clhs199;
const double clhs201 =             clhs200*clhs6;
const double clhs202 =             9*clhs0;
const double clhs203 =             clhs194 + clhs195;
const double clhs204 =             -clhs193 + clhs203;
const double clhs205 =             clhs204*clhs42*clhs44;
const double clhs206 =             -clhs199 + clhs203;
const double clhs207 =             clhs206*clhs44*clhs46;
const double clhs208 =             clhs103*clhs16;
const double clhs209 =             -clhs208;
const double clhs210 =             clhs195*clhs6;
const double clhs211 =             clhs19 + clhs209 + clhs210 + clhs22 + clhs25;
const double clhs212 =             clhs194*clhs6;
const double clhs213 =             clhs16*clhs6*(clhs26 + clhs28);
const double clhs214 =             clhs211 + 5*clhs212 + clhs213;
const double clhs215 =             DN(0,0)*clhs214;
const double clhs216 =             clhs19 + clhs209 + clhs212 + clhs213 + clhs22 + clhs25;
const double clhs217 =             5*clhs210 + clhs216;
const double clhs218 =             DN(0,1)*clhs217;
const double clhs219 =             3*clhs6;
const double clhs220 =             -clhs212;
const double clhs221 =             -clhs210;
const double clhs222 =             -2*clhs213;
const double clhs223 =             clhs18 + clhs208 + clhs21 + clhs221 + clhs222 + clhs24;
const double clhs224 =             clhs220 + clhs223;
const double clhs225 =             clhs224*clhs6;
const double clhs226 =             -7*clhs212 + clhs223;
const double clhs227 =             clhs226*clhs6;
const double clhs228 =             clhs18 + clhs208 + clhs21 + clhs220 + clhs24;
const double clhs229 =             -7*clhs210 + clhs222 + clhs228;
const double clhs230 =             clhs229*clhs6;
const double clhs231 =             8*clhs106;
const double clhs232 =             6*clhs102;
const double clhs233 =             clhs79*(clhs17 + clhs20 + clhs23);
const double clhs234 =             6*clhs48*lambda;
const double clhs235 =             clhs234*clhs27;
const double clhs236 =             clhs234*clhs29;
const double clhs237 =             -clhs233 + clhs235 + clhs236;
const double clhs238 =             -clhs231 - clhs232 + clhs237;
const double clhs239 =             6*clhs106;
const double clhs240 =             8*clhs102;
const double clhs241 =             clhs237 - clhs239 - clhs240;
const double clhs242 =             clhs211 + clhs212 + 3*clhs213;
const double clhs243 =             9*clhs0*clhs242;
const double clhs244 =             -12*DN(0,0)*N[0]*clhs10*clhs14*clhs46*clhs64*mu + 6*DN(0,0)*N[0]*clhs10*clhs44*clhs50*clhs80 + 18*DN(0,0)*N[0]*clhs14*clhs44*clhs52*clhs75 + 2*DN(0,0)*N[0]*clhs50*clhs6*clhs95 + DN(0,0)*clhs14*clhs187*clhs42*clhs44*clhs75 + 6*DN(0,0)*clhs146*clhs48*clhs70*lambda + 18*DN(0,1)*N[0]*clhs10*clhs44*clhs55*clhs75 - DN(0,1)*N[0]*clhs14*clhs189*clhs42*clhs64*mu + 6*DN(0,1)*N[0]*clhs14*clhs44*clhs57*clhs80 + 2*DN(0,1)*N[0]*clhs57*clhs6*clhs95 + 12*DN(0,1)*clhs10*clhs14*clhs146 + DN(0,1)*clhs10*clhs187*clhs44*clhs46*clhs75 - DN(0,1)*clhs110 + 10*DN(0,1)*clhs130*mu - DN(0,1)*clhs137*clhs6*mu + 6*DN(0,1)*clhs146*clhs48*clhs72*lambda + N[0]*clhs125 + N[0]*clhs126 + N[0]*clhs140*clhs57*clhs6 + clhs0*clhs152 + clhs0*clhs153 + clhs0*clhs156 + clhs0*clhs157 + clhs0*clhs185 + clhs0*clhs186 + 3*clhs0*clhs226*clhs42*clhs6 + 3*clhs0*clhs229*clhs46*clhs6 + clhs10*clhs111*clhs149 - clhs10*clhs116 + 8*clhs10*clhs141*clhs6 + clhs10*clhs16*clhs181*clhs50 + clhs10*clhs165 - clhs10*clhs169 + 6*clhs10*clhs176*clhs6*clhs70 + clhs10*clhs243*clhs42*clhs44 + clhs105*clhs111*clhs219 + clhs107*clhs113*clhs219 + clhs110*clhs132 - clhs111*clhs112 + clhs111*clhs142*clhs91 - clhs111*clhs79 - clhs112*clhs113 + clhs113*clhs14*clhs149 + clhs113*clhs144*clhs91 - clhs113*clhs79 - clhs114*clhs14 + clhs114*clhs198 + clhs114*clhs238*clhs42*clhs44 + clhs115*clhs124 + clhs115*clhs225 + clhs115*clhs42*clhs44*clhs80 - clhs115*y - clhs115 + clhs116*clhs201 + clhs116*clhs241*clhs44*clhs46 + clhs117*clhs124 + clhs117*clhs225 + clhs117*clhs44*clhs46*clhs80 - clhs117*y - clhs117 + clhs119*clhs124 + clhs119*clhs227 - clhs119*y - clhs119 + clhs120*clhs16 - clhs120 + clhs121*clhs16 - clhs121 + clhs122*clhs124 + clhs122*clhs230 - clhs122*y - clhs122 - clhs123*clhs70 - clhs123*clhs72 + clhs127*clhs128 + clhs128*clhs136 + clhs128*clhs145 + clhs129*clhs130 + clhs131*clhs91 + clhs131*mu + clhs132*clhs138 + clhs132*clhs147 - clhs133*clhs134 + clhs134*clhs135 + clhs139*clhs91 + clhs139*mu + 8*clhs14*clhs143*clhs6 + clhs14*clhs16*clhs181*clhs57 - clhs14*clhs165 + clhs14*clhs169 + 6*clhs14*clhs176*clhs6*clhs72 + clhs14*clhs243*clhs44*clhs46 + clhs141*clhs142 + clhs142*clhs150 + clhs143*clhs144 + clhs144*clhs154 + clhs148*clhs16 + clhs148 - clhs150*clhs180 + clhs151*clhs184 + clhs151*y + clhs151 - clhs154*clhs183 + clhs155*clhs184 + clhs155*y + clhs155 + clhs158*clhs161 - clhs162*clhs163 - clhs163*clhs164 + clhs166*clhs167 + clhs166*clhs168 + clhs170*clhs171 + clhs170*clhs172 - clhs173*clhs174 - clhs174*clhs175 + clhs177*clhs178 + clhs182*clhs83 + clhs182*clhs85 + clhs182*clhs87 + clhs182*clhs89 - clhs188*clhs71 - clhs188*clhs73 - clhs190*clhs191 - clhs190*clhs192 - clhs190*clhs71 - clhs190*clhs73 + clhs202*clhs205 + clhs202*clhs207 + clhs215*clhs39 + clhs218*clhs39;
const double clhs245 =             N[0]*clhs44;
const double clhs246 =             clhs16*y*(-clhs17 - clhs20 - clhs23 + clhs30);
const double clhs247 =             // Not supported in C:
// re
re(clhs2);
const double clhs248 =             // Not supported in C:
// re
re(clhs3);
const double clhs249 =             // Not supported in C:
// re
re(clhs4);
const double clhs250 =             clhs247 + clhs248 + clhs249;
const double clhs251 =             // Not supported in C:
// im
im(clhs2);
const double clhs252 =             // Not supported in C:
// im
im(clhs3);
const double clhs253 =             // Not supported in C:
// im
im(clhs4);
const double clhs254 =             clhs251 + clhs252 + clhs253;
const double clhs255 =             pow(clhs250, 2) + pow(clhs254, 2);
const double clhs256 =             1/(clhs15*pow(clhs255, 2));
const double clhs257 =             // Not supported in C:
// re
// re
// re
// re
// re
// re
re(clhs11) + re(clhs12) + re(clhs13) + re(clhs7) + re(clhs8) + re(clhs9);
const double clhs258 =             clhs250*clhs257;
const double clhs259 =             // Not supported in C:
// im
// im
// im
// im
// im
// im
im(clhs11) + im(clhs12) + im(clhs13) + im(clhs7) + im(clhs8) + im(clhs9);
const double clhs260 =             clhs254*clhs259;
const double clhs261 =             clhs258 + clhs260;
const double clhs262 =             // Not supported in C:
// Derivative
Derivative(clhs247, U(0,0));
const double clhs263 =             // Not supported in C:
// Derivative
Derivative(clhs251, U(0,0));
const double clhs264 =             1.0/clhs255;
const double clhs265 =             2*clhs264*(clhs250*clhs262 + clhs254*clhs263);
const double clhs266 =             clhs250*clhs259;
const double clhs267 =             clhs254*clhs257;
const double clhs268 =             clhs266 - clhs267;
const double clhs269 =             clhs245*clhs246 + clhs256*(clhs261*(clhs257*clhs262 - clhs258*clhs265 + clhs259*clhs263 - clhs260*clhs265) + clhs268*(-clhs257*clhs263 + clhs259*clhs262 - clhs265*clhs266 + clhs265*clhs267));
const double clhs270 =             pow(clhs31, -2);
const double clhs271 =             14*DN(0,0)*clhs44*mu;
const double clhs272 =             2*DN(0,0)*clhs44*mu;
const double clhs273 =             14*DN(0,1)*clhs44*mu;
const double clhs274 =             16*DN(0,0)*clhs10*clhs64*mu;
const double clhs275 =             12*DN(0,0)*clhs14*clhs64*mu;
const double clhs276 =             12*DN(0,1)*clhs10*clhs64*mu;
const double clhs277 =             16*DN(0,1)*clhs14*clhs64*mu;
const double clhs278 =             6*DN(0,0)*clhs10*clhs100;
const double clhs279 =             6*DN(0,1)*clhs10;
const double clhs280 =             clhs100*clhs14*clhs42*mu;
const double clhs281 =             clhs42*clhs64;
const double clhs282 =             6*DN(0,0)*clhs105;
const double clhs283 =             6*DN(0,1)*clhs107*clhs64;
const double clhs284 =             -clhs136*clhs71 + clhs14*clhs278*clhs46*mu - clhs145*clhs71 - clhs147*clhs73 - clhs271*clhs42 - clhs272*clhs46 - clhs273*clhs46 + clhs274*clhs42 + clhs275*clhs42 + clhs276*clhs46 + clhs277*clhs46 + clhs279*clhs280 + clhs281*clhs282 + clhs283*clhs46;
const double clhs285 =             6*DN(0,1);
const double clhs286 =             6*DN(0,0)*clhs6;
const double clhs287 =             clhs10*clhs286;
const double clhs288 =             6*DN(0,1)*clhs6;
const double clhs289 =             clhs14*clhs288;
const double clhs290 =             N[0]*clhs50;
const double clhs291 =             6*clhs6;
const double clhs292 =             clhs290*clhs291;
const double clhs293 =             N[0]*clhs52;
const double clhs294 =             clhs135*clhs291;
const double clhs295 =             clhs44*clhs50;
const double clhs296 =             DN(0,0)*clhs44*clhs55*mu;
const double clhs297 =             DN(0,0)*clhs44*clhs57*mu;
const double clhs298 =             DN(0,1)*clhs44*clhs50*mu;
const double clhs299 =             DN(0,1)*clhs44*clhs52*mu;
const double clhs300 =             clhs44*clhs57;
const double clhs301 =             6*N[0]*clhs6*y;
const double clhs302 =             clhs179*clhs6;
const double clhs303 =             clhs44*clhs55*clhs75;
const double clhs304 =             DN(0,0)*clhs81;
const double clhs305 =             DN(0,1)*clhs82;
const double clhs306 =             clhs160*clhs285;
const double clhs307 =             6*N[0]*clhs42*clhs44;
const double clhs308 =             clhs10*clhs307;
const double clhs309 =             6*clhs44;
const double clhs310 =             clhs14*clhs290;
const double clhs311 =             clhs309*clhs310;
const double clhs312 =             clhs10*clhs293;
const double clhs313 =             clhs309*clhs312;
const double clhs314 =             6*N[0]*clhs44*clhs46;
const double clhs315 =             clhs14*clhs314;
const double clhs316 =             6*clhs14*clhs44;
const double clhs317 =             clhs133*clhs316;
const double clhs318 =             6*clhs10*clhs44;
const double clhs319 =             clhs135*clhs318;
const double clhs320 =             20*DN(0,1)*U(0,0) + 20*DN(1,1)*U(1,0) + 20*DN(2,1)*U(2,0);
const double clhs321 =             DN(0,0)*clhs10*clhs64*mu;
const double clhs322 =             DN(0,0)*clhs14*clhs64*mu;
const double clhs323 =             2*DN(0,0)*clhs64*mu;
const double clhs324 =             20*DN(0,0)*U(0,0) + 20*DN(1,0)*U(1,0) + 20*DN(2,0)*U(2,0);
const double clhs325 =             DN(0,1)*clhs10*clhs64*mu;
const double clhs326 =             DN(0,1)*clhs14*clhs64*mu;
const double clhs327 =             2*DN(0,1)*clhs64*mu;
const double clhs328 =             6*N[0]*y;
const double clhs329 =             6*clhs10*clhs124*clhs44;
const double clhs330 =             clhs290*clhs329;
const double clhs331 =             6*clhs124*clhs14*clhs44;
const double clhs332 =             clhs135*clhs331;
const double clhs333 =             clhs14*clhs16;
const double clhs334 =             18*clhs44;
const double clhs335 =             clhs179*clhs44;
const double clhs336 =             clhs10*clhs16;
const double clhs337 =             N[0]*clhs14*clhs55;
const double clhs338 =             clhs10*clhs179*clhs44;
const double clhs339 =             DN(0,0)*clhs14*clhs91;
const double clhs340 =             12*clhs42*clhs64;
const double clhs341 =             12*clhs52*clhs64;
const double clhs342 =             12*DN(0,1)*clhs10*clhs64*clhs91;
const double clhs343 =             4*DN(0,0)*clhs10*clhs64*clhs95;
const double clhs344 =             4*DN(0,1)*clhs14*clhs64*clhs95;
const double clhs345 =             N[0]*clhs14*clhs189*clhs64;
const double clhs346 =             clhs345*clhs42;
const double clhs347 =             clhs345*clhs46;
const double clhs348 =             clhs14*clhs16*clhs64;
const double clhs349 =             12*clhs14*clhs16;
const double clhs350 =             N[0]*clhs14*clhs16*clhs55;
const double clhs351 =             3*clhs44;
const double clhs352 =             DN(0,0)*clhs204;
const double clhs353 =             clhs351*clhs352;
const double clhs354 =             DN(0,1)*clhs206;
const double clhs355 =             clhs351*clhs354;
const double clhs356 =             6*N[0]*clhs42*clhs64;
const double clhs357 =             6*N[0]*clhs64;
const double clhs358 =             DN(0,0)*clhs10;
const double clhs359 =             3*clhs358;
const double clhs360 =             clhs44*(clhs210 + clhs216);
const double clhs361 =             DN(0,1)*clhs14;
const double clhs362 =             3*clhs361;
const double clhs363 =             3*N[0]*clhs214;
const double clhs364 =             3*clhs290;
const double clhs365 =             3*N[0]*clhs217;
const double clhs366 =             3*clhs135;
const double clhs367 =             6*DN(0,0)*clhs191*clhs75 + DN(0,0)*clhs219*(clhs18 + clhs208 + clhs21 - 3*clhs212 + clhs221 + clhs24) + DN(0,1)*clhs219*(-3*clhs210 + clhs228) - N[0]*clhs10*clhs16*clhs334*clhs50 + N[0]*clhs10*clhs349*clhs52*clhs64 - N[0]*clhs14*clhs16*clhs334*clhs57 + N[0]*clhs189*clhs348*clhs42 + N[0]*clhs189*clhs348*clhs46 - clhs10*clhs125*clhs44 + clhs10*clhs224*clhs356 + clhs10*clhs288 - clhs10*clhs314 - clhs10*clhs328*clhs44*clhs70 + clhs124*clhs308 + clhs124*clhs315 - clhs125*clhs6 - clhs126*clhs14*clhs44 - clhs126*clhs6 - clhs129*clhs295 + clhs133*clhs291 - clhs133*clhs302 - clhs135*clhs338 + clhs14*clhs224*clhs357*clhs46 + clhs14*clhs286 - clhs14*clhs307 - clhs14*clhs328*clhs44*clhs72 - clhs140*clhs300 - clhs145*clhs44*clhs70 - clhs147*clhs44*clhs72 + clhs159*clhs16 - clhs159*clhs59 + clhs159 + clhs16*clhs285 - clhs161 + clhs184*clhs287 + clhs184*clhs289 + clhs184*clhs292 + clhs184*clhs294 + clhs189*clhs350*clhs64 + clhs197*clhs356 + clhs200*clhs357*clhs46 + clhs214*clhs364*clhs44 + clhs217*clhs366*clhs44 + clhs284 + clhs285*clhs303 - clhs285*clhs44*clhs55*mu + clhs285 + clhs287*y + clhs287 + clhs289*y + clhs289 + clhs291*clhs293 + clhs292*y + clhs292 - clhs293*clhs302 + clhs294*y + clhs294 + 10*clhs296 - 10*clhs297 - 10*clhs298 + 10*clhs299 + clhs301*clhs70 + clhs301*clhs72 + 2*clhs304 + 2*clhs305 - clhs306 + clhs307*clhs333 - clhs308*y - clhs308 - clhs310*clhs335 - clhs311 - clhs312*clhs335 - clhs313 + clhs314*clhs336 - clhs315*y - clhs315 - clhs317 - clhs319 - clhs320*clhs321 + clhs320*clhs322 - clhs323*clhs83 - clhs323*clhs85 + clhs324*clhs325 - clhs324*clhs326 - clhs327*clhs87 - clhs327*clhs89 + clhs330 + clhs332 - clhs335*clhs337 + clhs339*clhs340 + clhs339*clhs341 + clhs342*clhs46 + clhs342*clhs55 + clhs343*clhs42 + clhs343*clhs50 + clhs344*clhs46 + clhs344*clhs57 + clhs346 + clhs347 + clhs353 + clhs355 + clhs359*clhs360 + clhs360*clhs362 + clhs363*clhs71 + clhs365*clhs73;
const double clhs368 =             (2.0L/9.0L)*clhs1*clhs109*clhs270*clhs367*h;
const double clhs369 =             3*bdf0;
const double clhs370 =             14*N[0]*clhs64*mu;
const double clhs371 =             3*clhs44*clhs48*clhs70*lambda;
const double clhs372 =             3*clhs44*clhs48*clhs72*lambda;
const double clhs373 =             6*N[0]*clhs48*clhs64*lambda;
const double clhs374 =             clhs322*clhs55;
const double clhs375 =             clhs321*clhs57;
const double clhs376 =             clhs326*clhs50;
const double clhs377 =             clhs325*clhs52;
const double clhs378 =             N[0]*clhs42*clhs50*clhs64*mu;
const double clhs379 =             N[0]*clhs42*clhs52*clhs64*mu;
const double clhs380 =             4*clhs42*clhs46*clhs64*mu;
const double clhs381 =             10*clhs42*clhs64*mu;
const double clhs382 =             10*clhs46*clhs64*mu;
const double clhs383 =             N[0]*clhs46*clhs55*clhs64*mu;
const double clhs384 =             N[0]*clhs46*clhs57*clhs64*mu;
const double clhs385 =             N[0]*clhs100*clhs43*mu;
const double clhs386 =             N[0]*clhs100*clhs47*mu;
const double clhs387 =             clhs42*clhs70;
const double clhs388 =             clhs46*clhs72;
const double clhs389 =             12*DN(0,0)*clhs14*clhs64*clhs75;
const double clhs390 =             clhs52*clhs64;
const double clhs391 =             12*DN(0,1)*clhs10*clhs64*clhs75;
const double clhs392 =             clhs55*clhs64*clhs75;
const double clhs393 =             6*clhs42*clhs64*clhs75;
const double clhs394 =             6*clhs46*clhs64*clhs75;
const double clhs395 =             4*DN(0,0)*clhs10*clhs64*clhs80;
const double clhs396 =             -6*clhs48*lambda + 8*mu;
const double clhs397 =             clhs396*clhs50*clhs64;
const double clhs398 =             4*DN(0,1)*clhs14*clhs64*clhs80;
const double clhs399 =             clhs396*clhs57*clhs64;
const double clhs400 =             N[0]*clhs100*clhs75;
const double clhs401 =             clhs164*clhs43;
const double clhs402 =             clhs173*clhs47;
const double clhs403 =             clhs396*clhs42*clhs64;
const double clhs404 =             clhs396*clhs46*clhs64;
const double clhs405 =             clhs100*clhs46*mu;
const double clhs406 =             6*N[0]*clhs100*clhs80;
const double clhs407 =             N[0]*clhs100*clhs42*clhs52*clhs75;
const double clhs408 =             N[0]*clhs100*clhs46*clhs55*clhs75;
const double clhs409 =             N[0]*clhs100*clhs42*clhs50*clhs80;
const double clhs410 =             N[0]*clhs100*clhs46*clhs57*clhs80;
const double clhs411 =             pow(clhs5, -5);
const double clhs412 =             clhs14*clhs162*clhs411*clhs42*clhs46*mu;
const double clhs413 =             6*N[0]*clhs100;
const double clhs414 =             clhs233 - clhs235 - clhs236;
const double clhs415 =             clhs231 + clhs232 + clhs414;
const double clhs416 =             clhs415*clhs43;
const double clhs417 =             clhs239 + clhs240 + clhs414;
const double clhs418 =             clhs417*clhs47;
const double clhs419 =             -6*DN(0,0)*clhs14*clhs390*clhs75 - DN(0,0)*clhs371 - DN(0,0)*clhs51 - DN(0,0)*clhs53*mu + DN(0,0)*clhs76 - DN(0,1)*clhs372 - DN(0,1)*clhs62 + DN(0,1)*clhs77 - N[0]*clhs369 + N[0]*clhs380 - N[0]*clhs412 + clhs101*clhs167 + clhs101*clhs168 + clhs133*clhs381 - clhs133*clhs394 - clhs135*clhs381 - clhs135*clhs404 + clhs137*clhs88 - clhs158*clhs61 - clhs162*clhs385 - clhs164*clhs385 + clhs164*clhs407 + clhs171*clhs405 + clhs172*clhs405 - clhs173*clhs386 + clhs173*clhs408 - clhs175*clhs386 - clhs279*clhs392 + clhs284 - clhs290*clhs382 - clhs290*clhs403 - clhs293*clhs393 - 5*clhs296 + 5*clhs297 + 5*clhs298 - 5*clhs299 + clhs304 + clhs305 - clhs358*clhs397 - clhs361*clhs399 + clhs370*clhs43 + clhs370*clhs47 + clhs373*clhs387 + clhs373*clhs388 + clhs373*clhs43 + clhs373*clhs47 - clhs374 - clhs375 - clhs376 - clhs377 + 8*clhs378 + 6*clhs379 + 6*clhs383 + 8*clhs384 - clhs389*clhs42 - clhs391*clhs46 - clhs395*clhs42 - clhs398*clhs46 + clhs400*clhs401 + clhs400*clhs402 + clhs406*clhs94 + clhs406*clhs97 + clhs409*clhs67 + clhs410*clhs66 - clhs413*clhs416 - clhs413*clhs418;
const double clhs420 =             (1.0L/9.0L)*clhs1*clhs32*clhs367*h;
const double clhs421 =             N[0]*bdf0;
const double clhs422 =             N[1]*clhs44;
const double clhs423 =             // Not supported in C:
// Derivative
Derivative(clhs248, U(1,0));
const double clhs424 =             // Not supported in C:
// Derivative
Derivative(clhs252, U(1,0));
const double clhs425 =             2*clhs264*(clhs250*clhs423 + clhs254*clhs424);
const double clhs426 =             clhs246*clhs422 + clhs256*(clhs261*(clhs257*clhs423 - clhs258*clhs425 + clhs259*clhs424 - clhs260*clhs425) + clhs268*(-clhs257*clhs424 + clhs259*clhs423 - clhs266*clhs425 + clhs267*clhs425));
const double clhs427 =             14*DN(1,0)*mu;
const double clhs428 =             DN(1,0)*mu;
const double clhs429 =             2*DN(1,0)*mu;
const double clhs430 =             DN(1,0)*clhs44*clhs55*mu;
const double clhs431 =             5*clhs44*clhs57;
const double clhs432 =             2*DN(1,1)*mu;
const double clhs433 =             5*DN(1,1);
const double clhs434 =             clhs44*clhs50*mu;
const double clhs435 =             14*DN(1,1)*mu;
const double clhs436 =             14*N[1]*clhs64*mu;
const double clhs437 =             6*DN(1,0)*clhs48*lambda;
const double clhs438 =             6*DN(1,1)*clhs48*lambda;
const double clhs439 =             6*N[1]*clhs48*clhs64*lambda;
const double clhs440 =             16*DN(1,0)*clhs10;
const double clhs441 =             12*DN(1,0)*clhs14;
const double clhs442 =             DN(1,0)*clhs14;
const double clhs443 =             clhs55*clhs64*mu;
const double clhs444 =             clhs442*clhs443;
const double clhs445 =             DN(1,0)*clhs10;
const double clhs446 =             clhs57*clhs64*mu;
const double clhs447 =             clhs445*clhs446;
const double clhs448 =             DN(1,1)*clhs14;
const double clhs449 =             clhs50*clhs64*mu;
const double clhs450 =             clhs448*clhs449;
const double clhs451 =             DN(1,1)*clhs10;
const double clhs452 =             clhs52*clhs64*mu;
const double clhs453 =             clhs451*clhs452;
const double clhs454 =             12*DN(1,1)*clhs10;
const double clhs455 =             16*DN(1,1)*clhs14;
const double clhs456 =             N[1]*clhs42*clhs50*clhs64*mu;
const double clhs457 =             N[1]*clhs42*clhs52*clhs64*mu;
const double clhs458 =             N[1]*clhs55;
const double clhs459 =             N[1]*clhs57;
const double clhs460 =             N[1]*clhs50;
const double clhs461 =             N[1]*clhs52;
const double clhs462 =             N[1]*clhs46*clhs55*clhs64*mu;
const double clhs463 =             N[1]*clhs46*clhs57*clhs64*mu;
const double clhs464 =             N[1]*clhs100*clhs43*mu;
const double clhs465 =             N[1]*clhs100*clhs47*mu;
const double clhs466 =             12*DN(1,0)*clhs14*clhs64*clhs75;
const double clhs467 =             6*clhs52*clhs64*clhs75;
const double clhs468 =             12*DN(1,1)*clhs10*clhs64*clhs75;
const double clhs469 =             6*clhs55*clhs64*clhs75;
const double clhs470 =             4*DN(1,0)*clhs10*clhs64*clhs80;
const double clhs471 =             4*DN(1,1)*clhs14*clhs64*clhs80;
const double clhs472 =             N[1]*clhs100*clhs75;
const double clhs473 =             6*DN(1,0)*clhs10*clhs14;
const double clhs474 =             DN(1,1)*clhs10*clhs66;
const double clhs475 =             3*clhs458;
const double clhs476 =             3*clhs459;
const double clhs477 =             clhs10*clhs100*clhs42*mu;
const double clhs478 =             3*clhs460;
const double clhs479 =             clhs100*clhs14*clhs46*mu;
const double clhs480 =             3*clhs461;
const double clhs481 =             clhs10*clhs100*clhs46*mu;
const double clhs482 =             6*N[1]*clhs100*clhs80;
const double clhs483 =             N[1]*clhs100*clhs42*clhs52*clhs75;
const double clhs484 =             clhs100*clhs173*clhs46*clhs75;
const double clhs485 =             N[1]*clhs100*clhs42*clhs50*clhs80;
const double clhs486 =             clhs100*clhs46*clhs66*clhs80;
const double clhs487 =             6*DN(1,0)*clhs105;
const double clhs488 =             6*DN(1,1)*clhs107*clhs64;
const double clhs489 =             6*N[1]*clhs100;
const double clhs490 =             -DN(1,0)*clhs371 - DN(1,0)*clhs51 + DN(1,0)*clhs76 + DN(1,0)*clhs81 - DN(1,1)*clhs372 - DN(1,1)*clhs61*mu - DN(1,1)*clhs62 + DN(1,1)*clhs77 + DN(1,1)*clhs82 - N[1]*clhs369 + N[1]*clhs380 - N[1]*clhs412 + clhs101*clhs474 - clhs162*clhs464 - clhs164*clhs464 + clhs164*clhs483 - clhs173*clhs465 - clhs175*clhs465 + clhs280*clhs475 + clhs281*clhs487 + clhs381*clhs458 - clhs381*clhs459 - clhs382*clhs460 + clhs382*clhs461 + clhs387*clhs439 + clhs388*clhs439 - clhs393*clhs461 - clhs394*clhs458 - clhs397*clhs445 - clhs399*clhs448 + clhs401*clhs472 + clhs402*clhs472 - clhs403*clhs460 - clhs404*clhs459 + clhs405*clhs473 - clhs416*clhs489 - clhs418*clhs489 - clhs42*clhs466 - clhs42*clhs470 - clhs427*clhs71 + clhs428*clhs431 - clhs428*clhs53 - clhs429*clhs73 + clhs43*clhs436 + clhs43*clhs439 - 5*clhs430 - clhs432*clhs71 + clhs433*clhs434 - clhs433*clhs59 - clhs435*clhs73 + clhs436*clhs47 - clhs437*clhs71 - clhs438*clhs73 + clhs439*clhs47 + clhs440*clhs84 + clhs441*clhs84 - clhs442*clhs467 - clhs444 - clhs447 - clhs450 - clhs451*clhs469 - clhs453 + clhs454*clhs88 + clhs455*clhs88 + 8*clhs456 + 6*clhs457 + clhs458*clhs484 + clhs459*clhs486 - clhs46*clhs468 - clhs46*clhs471 + clhs46*clhs488 + 6*clhs462 + 8*clhs463 + clhs476*clhs477 + clhs478*clhs479 + clhs480*clhs481 + clhs482*clhs94 + clhs482*clhs97 + clhs485*clhs67;
const double clhs491 =             (1.0L/9.0L)*clhs1*clhs109*clhs32*clhs44*h;
const double clhs492 =             6*DN(0,0)*N[1];
const double clhs493 =             clhs10*clhs492;
const double clhs494 =             6*DN(0,1)*N[1];
const double clhs495 =             clhs14*clhs494;
const double clhs496 =             6*N[0];
const double clhs497 =             clhs445*clhs496;
const double clhs498 =             clhs442*clhs496;
const double clhs499 =             clhs451*clhs496;
const double clhs500 =             clhs448*clhs496;
const double clhs501 =             6*N[1];
const double clhs502 =             clhs290*clhs501;
const double clhs503 =             clhs293*clhs501;
const double clhs504 =             clhs133*clhs501;
const double clhs505 =             clhs135*clhs501;
const double clhs506 =             6*N[0]*N[1]*y;
const double clhs507 =             16*DN(0,0)*clhs6*mu;
const double clhs508 =             12*DN(0,0)*clhs6*mu;
const double clhs509 =             20*DN(0,0)*clhs6*mu;
const double clhs510 =             28*DN(0,0)*clhs42*clhs6*mu;
const double clhs511 =             4*DN(0,0)*clhs46*clhs6*mu;
const double clhs512 =             20*DN(0,1)*clhs6*mu;
const double clhs513 =             12*DN(0,1)*clhs6*mu;
const double clhs514 =             16*DN(0,1)*clhs6*mu;
const double clhs515 =             4*DN(0,1)*clhs42*clhs6*mu;
const double clhs516 =             28*DN(0,1)*clhs46*clhs6*mu;
const double clhs517 =             12*DN(0,0)*N[1]*clhs48*clhs6*lambda;
const double clhs518 =             12*DN(0,1)*N[1]*clhs48*clhs6*lambda;
const double clhs519 =             12*DN(0,0)*clhs14*clhs6*clhs91;
const double clhs520 =             12*DN(0,0)*clhs6*clhs91;
const double clhs521 =             12*DN(0,1)*clhs10*clhs6*clhs91;
const double clhs522 =             12*DN(0,1)*clhs6*clhs91;
const double clhs523 =             4*DN(0,0)*clhs10*clhs6*clhs95;
const double clhs524 =             4*clhs460;
const double clhs525 =             DN(0,0)*clhs6*clhs95;
const double clhs526 =             12*DN(0,0)*clhs10*clhs14*clhs6;
const double clhs527 =             N[1]*clhs526;
const double clhs528 =             4*DN(0,1)*clhs14*clhs6*clhs95;
const double clhs529 =             4*clhs459;
const double clhs530 =             DN(0,1)*clhs6*clhs95;
const double clhs531 =             12*DN(0,1)*clhs10*clhs14*clhs6;
const double clhs532 =             clhs14*clhs445;
const double clhs533 =             clhs177*clhs532;
const double clhs534 =             12*N[0]*clhs14*clhs6;
const double clhs535 =             12*N[0]*N[1]*clhs42*clhs6;
const double clhs536 =             clhs10*clhs535;
const double clhs537 =             12*N[1]*clhs6;
const double clhs538 =             clhs310*clhs537;
const double clhs539 =             clhs312*clhs537;
const double clhs540 =             12*N[0]*N[1]*clhs46*clhs6;
const double clhs541 =             clhs14*clhs540;
const double clhs542 =             N[0]*clhs189*clhs57*clhs6;
const double clhs543 =             6*DN(0,0)*clhs14*clhs44*mu;
const double clhs544 =             48*N[0]*U(0,1) + 48*N[1]*U(1,1) + 48*N[2]*U(2,1);
const double clhs545 =             DN(0,0)*N[1]*clhs42*clhs44*mu;
const double clhs546 =             36*N[0]*U(0,2) + 36*N[1]*U(1,2) + 36*N[2]*U(2,2);
const double clhs547 =             60*DN(0,0)*N[1]*clhs44*clhs46*mu;
const double clhs548 =             6*DN(0,0)*clhs10*clhs44*mu;
const double clhs549 =             6*DN(0,1)*clhs44*mu;
const double clhs550 =             60*DN(0,1)*N[1]*clhs42*clhs44*mu;
const double clhs551 =             clhs14*clhs460;
const double clhs552 =             clhs10*clhs461;
const double clhs553 =             36*N[0]*U(0,1) + 36*N[1]*U(1,1) + 36*N[2]*U(2,1);
const double clhs554 =             DN(0,1)*N[1]*clhs44*clhs46*mu;
const double clhs555 =             48*N[0]*U(0,2) + 48*N[1]*U(1,2) + 48*N[2]*U(2,2);
const double clhs556 =             N[0]*N[1]*clhs6*y;
const double clhs557 =             clhs189*clhs70;
const double clhs558 =             12*N[0]*U(0,2) + 12*N[1]*U(1,2) + 12*N[2]*U(2,2);
const double clhs559 =             clhs558*clhs72;
const double clhs560 =             12*N[0]*clhs14*clhs16*clhs6;
const double clhs561 =             N[1]*clhs10*clhs16;
const double clhs562 =             36*N[0]*clhs50*clhs6;
const double clhs563 =             12*N[0]*clhs57*clhs6;
const double clhs564 =             36*N[0]*clhs14*clhs16*clhs57*clhs6;
const double clhs565 =             N[1]*clhs184;
const double clhs566 =             N[0]*clhs189*clhs50*clhs6;
const double clhs567 =             N[0]*clhs558*clhs57*clhs6;
const double clhs568 =             36*DN(0,0)*clhs14*clhs42*clhs44*clhs75;
const double clhs569 =             36*DN(0,0)*clhs14*clhs44*clhs75;
const double clhs570 =             36*DN(0,1)*clhs10*clhs44*clhs46*clhs75;
const double clhs571 =             36*DN(0,1)*clhs10*clhs44*clhs75;
const double clhs572 =             12*DN(0,0)*clhs10*clhs42*clhs44*clhs80;
const double clhs573 =             12*DN(0,0)*clhs10*clhs44*clhs80;
const double clhs574 =             12*DN(0,1)*clhs14*clhs44*clhs46*clhs80;
const double clhs575 =             12*DN(0,1)*clhs14*clhs44*clhs80;
const double clhs576 =             N[0]*N[1]*clhs10*clhs546;
const double clhs577 =             DN(0,0)*clhs14*clhs162*clhs46*clhs64*mu;
const double clhs578 =             DN(0,1)*clhs14*clhs162*clhs42*clhs64*mu;
const double clhs579 =             N[0]*N[1]*clhs14*clhs16*clhs553;
const double clhs580 =             36*N[0]*clhs10*clhs14*clhs16*clhs52;
const double clhs581 =             N[0]*clhs14*clhs16*clhs55*clhs553;
const double clhs582 =             6*N[0]*clhs197*clhs6;
const double clhs583 =             6*N[0]*clhs200*clhs6;
const double clhs584 =             18*N[0]*N[1];
const double clhs585 =             6*DN(0,0)*clhs105*clhs6;
const double clhs586 =             6*DN(0,1)*clhs107*clhs6;
const double clhs587 =             6*N[0]*clhs226*clhs42*clhs6;
const double clhs588 =             6*N[0]*clhs229*clhs46*clhs6;
const double clhs589 =             12*DN(0,0)*clhs238*clhs42*clhs44;
const double clhs590 =             12*DN(0,1)*clhs241*clhs44*clhs46;
const double clhs591 =             clhs173*clhs42*clhs44;
const double clhs592 =             N[0]*N[1]*clhs242;
const double clhs593 =             clhs164*clhs44*clhs46;
const double clhs594 =             -DN(1,0)*clhs127 - DN(1,0)*clhs136 - DN(1,0)*clhs145 + DN(1,0)*clhs363 + DN(1,0)*clhs519 + DN(1,0)*clhs523 + DN(1,0)*clhs582 + DN(1,0)*clhs585 - DN(1,1)*clhs110 - DN(1,1)*clhs138 - DN(1,1)*clhs147 + DN(1,1)*clhs365 + DN(1,1)*clhs521 + DN(1,1)*clhs528 + DN(1,1)*clhs583 + DN(1,1)*clhs586 + N[1]*clhs125 + N[1]*clhs126 + N[1]*clhs510 + N[1]*clhs511 + N[1]*clhs515 + N[1]*clhs516 + N[1]*clhs531 + N[1]*clhs542 + N[1]*clhs564 + N[1]*clhs568 + N[1]*clhs570 + N[1]*clhs572 + N[1]*clhs574 - N[1]*clhs577 - N[1]*clhs578 + N[1]*clhs587 + N[1]*clhs588 + N[1]*clhs589 + N[1]*clhs590 - clhs10*clhs494 + clhs10*clhs540 + clhs10*clhs547 - clhs10*clhs550 + clhs124*clhs493 + clhs124*clhs495 + clhs124*clhs497 + clhs124*clhs500 + clhs124*clhs502 + clhs124*clhs505 - clhs14*clhs492 + clhs14*clhs535 - clhs14*clhs547 + clhs14*clhs550 + clhs16*clhs498 + clhs16*clhs499 + clhs16*clhs503 + clhs16*clhs504 + clhs16*clhs527 + clhs16*clhs533 + clhs16*clhs538 + clhs16*clhs539 + clhs178*clhs537 + clhs184*clhs536 + clhs184*clhs541 + clhs198*clhs492 + clhs201*clhs494 + clhs205*clhs584 + clhs207*clhs584 + clhs215*clhs40 + clhs218*clhs40 + clhs225*clhs493 + clhs225*clhs495 + clhs225*clhs497 + clhs225*clhs500 + clhs227*clhs502 + clhs230*clhs505 - clhs333*clhs535 - clhs336*clhs540 + clhs337*clhs537 + clhs350*clhs537 + clhs42*clhs517 - clhs422*clhs580 - clhs422*clhs581 + clhs442*clhs508 - clhs442*clhs512 + clhs445*clhs507 + clhs445*clhs512 + clhs448*clhs509 + clhs448*clhs514 - clhs451*clhs509 + clhs451*clhs513 + clhs451*clhs534 + clhs451*clhs543 + clhs451*clhs560 - clhs458*clhs509 + clhs458*clhs513 + clhs458*clhs522 + clhs458*clhs543 + clhs458*clhs571 + clhs459*clhs509 + clhs459*clhs514 + clhs459*clhs548 + clhs459*clhs575 + clhs46*clhs518 + clhs460*clhs507 + clhs460*clhs512 + clhs460*clhs573 + clhs461*clhs508 - clhs461*clhs512 + clhs461*clhs520 + clhs461*clhs569 - clhs493*y - clhs493 - clhs495*y - clhs495 - clhs497*y - clhs497 - clhs498 - clhs499 - clhs500*y - clhs500 - clhs502*y - clhs502 - clhs503 - clhs504 - clhs505*y - clhs505 - clhs506*clhs70 - clhs506*clhs72 + clhs517*clhs70 + clhs518*clhs72 + clhs524*clhs525 + clhs527 + clhs529*clhs530 + clhs532*clhs549 + clhs533 + clhs536*y + clhs536 + clhs538 + clhs539 + clhs541*y + clhs541 - clhs544*clhs545 - clhs545*clhs546 + clhs549*clhs551 + clhs549*clhs552 - clhs553*clhs554 - clhs554*clhs555 + clhs556*clhs557 + clhs556*clhs559 + clhs561*clhs562 + clhs561*clhs563 + clhs565*clhs566 + clhs565*clhs567 - clhs576*clhs71 - clhs576*clhs73 - clhs579*clhs71 - clhs579*clhs73 + clhs591*clhs592 + clhs592*clhs593;
const double clhs595 =             N[2]*clhs44;
const double clhs596 =             // Not supported in C:
// Derivative
Derivative(clhs249, U(2,0));
const double clhs597 =             // Not supported in C:
// Derivative
Derivative(clhs253, U(2,0));
const double clhs598 =             2*clhs264*(clhs250*clhs596 + clhs254*clhs597);
const double clhs599 =             clhs246*clhs595 + clhs256*(clhs261*(clhs257*clhs596 - clhs258*clhs598 + clhs259*clhs597 - clhs260*clhs598) + clhs268*(-clhs257*clhs597 + clhs259*clhs596 - clhs266*clhs598 + clhs267*clhs598));
const double clhs600 =             14*DN(2,0);
const double clhs601 =             clhs42*clhs44*mu;
const double clhs602 =             DN(2,0)*mu;
const double clhs603 =             2*DN(2,0)*mu;
const double clhs604 =             DN(2,0)*clhs44*clhs55*mu;
const double clhs605 =             2*DN(2,1);
const double clhs606 =             5*DN(2,1);
const double clhs607 =             14*DN(2,1)*mu;
const double clhs608 =             14*N[2]*clhs64*mu;
const double clhs609 =             6*DN(2,0)*clhs48*lambda;
const double clhs610 =             6*DN(2,1)*clhs48*lambda;
const double clhs611 =             6*N[2]*clhs48*clhs64*lambda;
const double clhs612 =             16*DN(2,0)*clhs10;
const double clhs613 =             12*DN(2,0)*clhs14;
const double clhs614 =             DN(2,0)*clhs14;
const double clhs615 =             clhs443*clhs614;
const double clhs616 =             DN(2,0)*clhs10;
const double clhs617 =             clhs446*clhs616;
const double clhs618 =             DN(2,1)*clhs14;
const double clhs619 =             clhs449*clhs618;
const double clhs620 =             DN(2,1)*clhs10;
const double clhs621 =             clhs452*clhs620;
const double clhs622 =             12*DN(2,1)*clhs10;
const double clhs623 =             16*DN(2,1)*clhs14;
const double clhs624 =             N[2]*clhs42*clhs50*clhs64*mu;
const double clhs625 =             N[2]*clhs42*clhs52*clhs64*mu;
const double clhs626 =             N[2]*clhs55;
const double clhs627 =             N[2]*clhs57;
const double clhs628 =             N[2]*clhs50;
const double clhs629 =             N[2]*clhs52;
const double clhs630 =             N[2]*clhs46*clhs55*clhs64*mu;
const double clhs631 =             N[2]*clhs46*clhs57*clhs64*mu;
const double clhs632 =             N[2]*clhs100*clhs43*mu;
const double clhs633 =             N[2]*clhs100*clhs47*mu;
const double clhs634 =             12*DN(2,0)*clhs14*clhs64*clhs75;
const double clhs635 =             12*DN(2,1)*clhs10*clhs64*clhs75;
const double clhs636 =             4*DN(2,0)*clhs10*clhs64*clhs80;
const double clhs637 =             4*DN(2,1)*clhs14*clhs64*clhs80;
const double clhs638 =             N[2]*clhs100*clhs75;
const double clhs639 =             6*DN(2,0)*clhs10*clhs14;
const double clhs640 =             DN(2,1)*clhs10*clhs66;
const double clhs641 =             3*clhs626;
const double clhs642 =             3*clhs627;
const double clhs643 =             3*clhs628;
const double clhs644 =             3*clhs629;
const double clhs645 =             6*N[2]*clhs100*clhs80;
const double clhs646 =             N[2]*clhs100*clhs42*clhs52*clhs75;
const double clhs647 =             N[2]*clhs100*clhs42*clhs50*clhs80;
const double clhs648 =             6*DN(2,0)*clhs105;
const double clhs649 =             6*DN(2,1)*clhs107*clhs64;
const double clhs650 =             6*N[2]*clhs100;
const double clhs651 =             -DN(2,0)*clhs371 - DN(2,0)*clhs51 + DN(2,0)*clhs76 + DN(2,0)*clhs81 - DN(2,1)*clhs372 - DN(2,1)*clhs61*mu - DN(2,1)*clhs62 + DN(2,1)*clhs77 + DN(2,1)*clhs82 - N[2]*clhs369 + N[2]*clhs380 - N[2]*clhs412 + clhs101*clhs640 - clhs162*clhs632 - clhs164*clhs632 + clhs164*clhs646 - clhs173*clhs633 - clhs175*clhs633 + clhs280*clhs641 + clhs281*clhs648 + clhs381*clhs626 - clhs381*clhs627 - clhs382*clhs628 + clhs382*clhs629 + clhs387*clhs611 + clhs388*clhs611 - clhs393*clhs629 - clhs394*clhs626 - clhs397*clhs616 - clhs399*clhs618 + clhs401*clhs638 + clhs402*clhs638 - clhs403*clhs628 - clhs404*clhs627 + clhs405*clhs639 - clhs416*clhs650 - clhs418*clhs650 - clhs42*clhs634 - clhs42*clhs636 + clhs43*clhs608 + clhs43*clhs611 + clhs431*clhs602 + clhs434*clhs606 - clhs46*clhs635 - clhs46*clhs637 + clhs46*clhs649 - clhs467*clhs614 - clhs469*clhs620 + clhs47*clhs608 + clhs47*clhs611 + clhs477*clhs642 + clhs479*clhs643 + clhs481*clhs644 + clhs484*clhs626 + clhs486*clhs627 - clhs53*clhs602 - clhs59*clhs606 - clhs600*clhs601 - clhs601*clhs605 - clhs603*clhs73 - 5*clhs604 - clhs607*clhs73 - clhs609*clhs71 - clhs610*clhs73 + clhs612*clhs84 + clhs613*clhs84 - clhs615 - clhs617 - clhs619 - clhs621 + clhs622*clhs88 + clhs623*clhs88 + 8*clhs624 + 6*clhs625 + 6*clhs630 + 8*clhs631 + clhs645*clhs94 + clhs645*clhs97 + clhs647*clhs67;
const double clhs652 =             6*DN(0,0)*N[2];
const double clhs653 =             clhs10*clhs652;
const double clhs654 =             6*DN(0,1)*N[2];
const double clhs655 =             clhs14*clhs654;
const double clhs656 =             clhs496*clhs616;
const double clhs657 =             clhs496*clhs614;
const double clhs658 =             clhs496*clhs620;
const double clhs659 =             clhs496*clhs618;
const double clhs660 =             6*N[2];
const double clhs661 =             clhs290*clhs660;
const double clhs662 =             clhs293*clhs660;
const double clhs663 =             clhs133*clhs660;
const double clhs664 =             clhs135*clhs660;
const double clhs665 =             6*N[0]*N[2]*y;
const double clhs666 =             12*DN(0,0)*N[2]*clhs48*clhs6*lambda;
const double clhs667 =             12*DN(0,1)*N[2]*clhs48*clhs6*lambda;
const double clhs668 =             4*clhs628;
const double clhs669 =             N[2]*clhs526;
const double clhs670 =             4*clhs627;
const double clhs671 =             clhs14*clhs616;
const double clhs672 =             clhs177*clhs671;
const double clhs673 =             N[0]*N[2]*clhs42*clhs6;
const double clhs674 =             12*N[2]*clhs6;
const double clhs675 =             clhs310*clhs674;
const double clhs676 =             clhs312*clhs674;
const double clhs677 =             N[0]*N[2]*clhs46*clhs6;
const double clhs678 =             DN(0,0)*N[2]*clhs42*clhs44*mu;
const double clhs679 =             60*DN(0,0)*N[2]*clhs44*clhs46*mu;
const double clhs680 =             60*DN(0,1)*N[2]*clhs42*clhs44*mu;
const double clhs681 =             clhs14*clhs628;
const double clhs682 =             clhs10*clhs629;
const double clhs683 =             DN(0,1)*N[2]*clhs44*clhs46*mu;
const double clhs684 =             N[0]*N[2]*clhs189*clhs42*clhs6;
const double clhs685 =             N[0]*N[2]*clhs6*y;
const double clhs686 =             N[0]*N[2]*clhs46*clhs558*clhs6;
const double clhs687 =             N[2]*clhs10*clhs16;
const double clhs688 =             N[2]*clhs184;
const double clhs689 =             N[0]*N[2]*clhs10*clhs546;
const double clhs690 =             N[0]*N[2]*clhs14*clhs16*clhs553;
const double clhs691 =             18*N[0]*N[2];
const double clhs692 =             N[0]*N[2]*clhs242;
const double clhs693 =             -DN(2,0)*clhs127 - DN(2,0)*clhs136 - DN(2,0)*clhs145 + DN(2,0)*clhs363 + DN(2,0)*clhs519 + DN(2,0)*clhs523 + DN(2,0)*clhs582 + DN(2,0)*clhs585 - DN(2,1)*clhs110 - DN(2,1)*clhs138 - DN(2,1)*clhs147 + DN(2,1)*clhs365 + DN(2,1)*clhs521 + DN(2,1)*clhs528 + DN(2,1)*clhs583 + DN(2,1)*clhs586 + N[2]*clhs125 + N[2]*clhs126 + N[2]*clhs510 + N[2]*clhs511 + N[2]*clhs515 + N[2]*clhs516 + N[2]*clhs531 + N[2]*clhs542 + N[2]*clhs564 + N[2]*clhs568 + N[2]*clhs570 + N[2]*clhs572 + N[2]*clhs574 - N[2]*clhs577 - N[2]*clhs578 + N[2]*clhs587 + N[2]*clhs588 + N[2]*clhs589 + N[2]*clhs590 - 12*clhs10*clhs16*clhs677 - clhs10*clhs654 + clhs10*clhs679 - clhs10*clhs680 + clhs124*clhs653 + clhs124*clhs655 + clhs124*clhs656 + clhs124*clhs659 + clhs124*clhs661 + clhs124*clhs664 - clhs14*clhs652 - clhs14*clhs679 + clhs14*clhs680 + clhs16*clhs657 + clhs16*clhs658 + clhs16*clhs662 + clhs16*clhs663 + clhs16*clhs669 + clhs16*clhs672 + clhs16*clhs675 + clhs16*clhs676 + clhs178*clhs674 + clhs184*clhs684 + clhs184*clhs686 + clhs189*clhs673 + clhs189*clhs677 + clhs198*clhs652 + clhs201*clhs654 + clhs205*clhs691 + clhs207*clhs691 + clhs215*clhs41 + clhs218*clhs41 + clhs225*clhs653 + clhs225*clhs655 + clhs225*clhs656 + clhs225*clhs659 + clhs227*clhs661 + clhs230*clhs664 + clhs337*clhs674 - clhs349*clhs673 + clhs350*clhs674 + clhs42*clhs666 + clhs46*clhs667 + clhs507*clhs616 + clhs507*clhs628 + clhs508*clhs614 + clhs508*clhs629 + clhs509*clhs618 - clhs509*clhs620 - clhs509*clhs626 + clhs509*clhs627 - clhs512*clhs614 + clhs512*clhs616 + clhs512*clhs628 - clhs512*clhs629 + clhs513*clhs620 + clhs513*clhs626 + clhs514*clhs618 + clhs514*clhs627 + clhs520*clhs629 + clhs522*clhs626 + clhs525*clhs668 + clhs530*clhs670 + clhs534*clhs620 + clhs543*clhs620 + clhs543*clhs626 - clhs544*clhs678 - clhs546*clhs678 + clhs548*clhs627 + clhs549*clhs671 + clhs549*clhs681 + clhs549*clhs682 - clhs553*clhs683 - clhs555*clhs683 + clhs557*clhs685 + clhs558*clhs673 + clhs558*clhs677 + clhs559*clhs685 + clhs560*clhs620 + clhs562*clhs687 + clhs563*clhs687 + clhs566*clhs688 + clhs567*clhs688 + clhs569*clhs629 + clhs571*clhs626 + clhs573*clhs628 + clhs575*clhs627 - clhs580*clhs595 - clhs581*clhs595 + clhs591*clhs692 + clhs593*clhs692 - clhs653*y - clhs653 - clhs655*y - clhs655 - clhs656*y - clhs656 - clhs657 - clhs658 - clhs659*y - clhs659 - clhs661*y - clhs661 - clhs662 - clhs663 - clhs664*y - clhs664 - clhs665*clhs70 - clhs665*clhs72 + clhs666*clhs70 + clhs667*clhs72 + clhs669 + clhs672 + clhs675 + clhs676 + clhs684*y + clhs686*y - clhs689*clhs71 - clhs689*clhs73 - clhs690*clhs71 - clhs690*clhs73;
const double clhs694 =             N[0]*f_ext(0,0) + N[1]*f_ext(1,0) + N[2]*f_ext(2,0);
const double clhs695 =             (1.0L/3.0L)*DN(0,0)*clhs44*mu;
const double clhs696 =             4*clhs358;
const double clhs697 =             4*clhs290;
const double clhs698 =             6*N[0]*clhs46*clhs6;
const double clhs699 =             -clhs128*clhs63 - clhs14*clhs698 + clhs362 + clhs366 + clhs696 + clhs697;
const double clhs700 =             (1.0L/3.0L)*DN(0,1)*clhs44*mu;
const double clhs701 =             DN(0,0)*clhs14;
const double clhs702 =             DN(0,1)*clhs10;
const double clhs703 =             2*clhs293;
const double clhs704 =             4*N[0]*U(0,2) + 4*N[1]*U(1,2) + 4*N[2]*U(2,2);
const double clhs705 =             (1.0L/2.0L)*N[0]*clhs44;
const double clhs706 =             4*clhs135;
const double clhs707 =             2*y - 6;
const double clhs708 =             clhs10*clhs707;
const double clhs709 =             2*clhs16*clhs26;
const double clhs710 =             2*clhs16*clhs28;
const double clhs711 =             -4*clhs26 + clhs709 + clhs710;
const double clhs712 =             (4.0L/3.0L)*N[0]*clhs1*clhs109*clhs269*clhs270*h;
const double clhs713 =             (2.0L/3.0L)*N[0]*clhs1*clhs32*clhs694*h;
const double clhs714 =             -clhs434*(10*DN(0,1)*U(0,2) + 10*DN(1,1)*U(1,2) + 10*DN(2,1)*U(2,2)) + clhs59*(10*DN(0,1)*U(0,1) + 10*DN(1,1)*U(1,1) + 10*DN(2,1)*U(2,1));
const double clhs715 =             6*N[0]*U(0,0) + 6*N[1]*U(1,0) + 6*N[2]*U(2,0);
const double clhs716 =             pow(clhs50, 2);
const double clhs717 =             8*clhs44*mu;
const double clhs718 =             pow(clhs55, 2);
const double clhs719 =             6*clhs44*mu;
const double clhs720 =             6*clhs44*clhs75;
const double clhs721 =             clhs44*clhs50*clhs52*mu;
const double clhs722 =             14*clhs44*clhs55*mu;
const double clhs723 =             clhs44*clhs55*clhs57*mu;
const double clhs724 =             2*clhs44*clhs80;
const double clhs725 =             6*clhs42*clhs48*lambda;
const double clhs726 =             6*clhs48*clhs70*lambda;
const double clhs727 =             6*clhs46*clhs48*lambda;
const double clhs728 =             6*clhs48*clhs72*lambda;
const double clhs729 =             clhs10*clhs66;
const double clhs730 =             clhs189*clhs64*clhs91;
const double clhs731 =             clhs10*clhs50;
const double clhs732 =             16*clhs42*clhs64*mu;
const double clhs733 =             20*DN(0,1)*U(0,1) + 20*DN(1,1)*U(1,1) + 20*DN(2,1)*U(2,1);
const double clhs734 =             clhs42*clhs64*clhs733*mu;
const double clhs735 =             20*DN(0,0)*U(0,1) + 20*DN(1,0)*U(1,1) + 20*DN(2,0)*U(2,1);
const double clhs736 =             clhs46*clhs64*clhs735*mu;
const double clhs737 =             clhs50*clhs64;
const double clhs738 =             2*clhs57*clhs64*mu;
const double clhs739 =             4*N[0]*U(0,1) + 4*N[1]*U(1,1) + 4*N[2]*U(2,1);
const double clhs740 =             clhs14*clhs50*clhs91;
const double clhs741 =             clhs42*clhs739*clhs95;
const double clhs742 =             clhs64*clhs704*clhs95;
const double clhs743 =             clhs57*clhs64*clhs704*clhs95;
const double clhs744 =             6*clhs105*clhs42;
const double clhs745 =             clhs10*clhs734 - clhs10*clhs736 - clhs101*clhs14*clhs55*clhs67 - 6*clhs108*clhs46*clhs55 - clhs14*clhs405*clhs50*clhs67 - 12*clhs14*clhs50*clhs84 + 4*clhs14*clhs55*clhs737*mu - 16*clhs14*clhs55*clhs88 - clhs14*clhs734 + clhs14*clhs736 + clhs156 + clhs157 + clhs179*clhs70 - clhs180*clhs52 + clhs185 - clhs189*clhs46*clhs55*clhs64*mu - clhs192*clhs396*clhs57 + clhs192*clhs727 + clhs192*clhs728 + 3*clhs205 - 6*clhs295*clhs52*clhs75 + clhs295*clhs725 + clhs295*clhs726 - clhs340*clhs740 - clhs341*clhs740 + clhs434*(14*DN(0,0)*U(0,0) + 14*DN(1,0)*U(1,0) + 14*DN(2,0)*U(2,0)) + clhs46*clhs722 + clhs496*(U(0,1)*bdf0 + Un(0,1)*bdf1 + Unn(0,1)*bdf2) + clhs501*(U(1,1)*bdf0 + Un(1,1)*bdf1 + Unn(1,1)*bdf2) + clhs54*clhs55 + 2*clhs55*clhs64*clhs89*mu - clhs55*clhs743 + 2*clhs58 - clhs64*clhs716*clhs739*clhs95 + clhs660*(U(2,1)*bdf0 + Un(2,1)*bdf1 + Unn(2,1)*bdf2) - clhs694*clhs715 + clhs714 + clhs716*clhs717 - clhs716*clhs724 + clhs718*clhs719 - clhs718*clhs720 - clhs718*clhs730 + 6*clhs721 + 8*clhs723 - clhs729*clhs73 - clhs730*clhs99 - clhs731*clhs732 + clhs731*clhs738 - clhs737*clhs741 - clhs737*clhs744 - clhs742*clhs99;
const double clhs746 =             clhs31*stab_c2 + 2*nu*stab_c1/h;
const double clhs747 =             1.0/clhs746;
const double clhs748 =             (1.0L/9.0L)*clhs244*clhs44*clhs747*h;
const double clhs749 =             pow(clhs746, -2);
const double clhs750 =             (1.0L/9.0L)*clhs269*clhs367*clhs749*h*stab_c2;
const double clhs751 =             20*clhs57*clhs64*mu;
const double clhs752 =             20*clhs55*clhs64*mu;
const double clhs753 =             -clhs290*clhs751 + clhs293*clhs752;
const double clhs754 =             16*N[0]*clhs64*mu;
const double clhs755 =             12*N[0]*clhs64*mu;
const double clhs756 =             12*N[0]*clhs64*clhs75;
const double clhs757 =             4*clhs42*clhs64*mu;
const double clhs758 =             N[0]*clhs50*clhs52*clhs64*mu;
const double clhs759 =             4*clhs46*clhs64*mu;
const double clhs760 =             N[0]*clhs55*clhs57*clhs64*mu;
const double clhs761 =             4*N[0]*clhs64*clhs80;
const double clhs762 =             clhs14*clhs179*clhs44;
const double clhs763 =             12*clhs42*clhs48*clhs64*lambda;
const double clhs764 =             12*clhs48*clhs64*clhs70*lambda;
const double clhs765 =             12*clhs46*clhs48*clhs64*lambda;
const double clhs766 =             12*clhs48*clhs64*clhs72*lambda;
const double clhs767 =             12*clhs52*clhs64*clhs75;
const double clhs768 =             clhs553*clhs718;
const double clhs769 =             4*clhs57*clhs64*clhs80;
const double clhs770 =             clhs100*clhs42*clhs544*mu;
const double clhs771 =             36*clhs100*clhs42*mu;
const double clhs772 =             60*clhs10*clhs100*clhs42*mu;
const double clhs773 =             60*clhs10*clhs100*clhs46*mu;
const double clhs774 =             60*clhs100*clhs46*mu;
const double clhs775 =             12*clhs100*clhs55*mu;
const double clhs776 =             6*clhs10*clhs100*clhs57*mu;
const double clhs777 =             6*clhs100*clhs55*mu;
const double clhs778 =             clhs100*clhs46*clhs553*mu;
const double clhs779 =             N[0]*clhs100*clhs80;
const double clhs780 =             clhs189*clhs716;
const double clhs781 =             clhs42*clhs75;
const double clhs782 =             36*N[0]*clhs100*clhs14*clhs50;
const double clhs783 =             clhs52*clhs75;
const double clhs784 =             12*N[0]*clhs100*clhs14*clhs55*clhs80;
const double clhs785 =             clhs162*clhs411*clhs46*mu;
const double clhs786 =             clhs204*clhs42;
const double clhs787 =             12*clhs100*clhs415*clhs42;
const double clhs788 =             12*clhs100*clhs417*clhs46;
const double clhs789 =             (1.0L/18.0L)*clhs367*clhs747*h;
const double clhs790 =             N[0]*N[1];
const double clhs791 =             4*clhs445;
const double clhs792 =             8*clhs10*clhs42*clhs6;
const double clhs793 =             6*N[1]*clhs46*clhs6;
const double clhs794 =             -N[1]*clhs792 - clhs14*clhs793 + 3*clhs448 + clhs476 + clhs524 + clhs791;
const double clhs795 =             2*clhs461;
const double clhs796 =             clhs42*clhs6*clhs704;
const double clhs797 =             8*clhs10*clhs14*clhs46*clhs6;
const double clhs798 =             DN(1,0)*clhs204;
const double clhs799 =             clhs42*clhs6*clhs711;
const double clhs800 =             (4.0L/3.0L)*N[0]*clhs1*clhs109*clhs270*clhs426*h;
const double clhs801 =             (1.0L/9.0L)*clhs367*clhs426*clhs749*h*stab_c2;
const double clhs802 =             -clhs460*clhs751 + clhs461*clhs752;
const double clhs803 =             16*N[1]*clhs64*mu;
const double clhs804 =             12*N[1]*clhs64*mu;
const double clhs805 =             12*N[1]*clhs64*clhs75;
const double clhs806 =             20*clhs50*clhs64*mu;
const double clhs807 =             N[1]*clhs50*clhs52*clhs64*mu;
const double clhs808 =             N[1]*clhs55*clhs57*clhs64*mu;
const double clhs809 =             4*N[1]*clhs64*clhs80;
const double clhs810 =             N[1]*clhs14*clhs189*clhs64;
const double clhs811 =             clhs100*clhs55*mu;
const double clhs812 =             6*clhs100*clhs14*clhs50*mu;
const double clhs813 =             60*clhs100*clhs14*clhs42*mu;
const double clhs814 =             clhs100*clhs46*clhs555*mu;
const double clhs815 =             N[1]*clhs100*clhs80;
const double clhs816 =             36*N[1]*clhs100*clhs14*clhs50;
const double clhs817 =             clhs100*clhs46*clhs553*clhs75;
const double clhs818 =             clhs100*clhs46*clhs558*clhs80;
const double clhs819 =             clhs100*clhs558*clhs57*clhs80;
const double clhs820 =             clhs14*clhs162*clhs411*clhs42*mu;
const double clhs821 =             6*N[1]*clhs64;
const double clhs822 =             (1.0L/18.0L)*clhs44*clhs745*clhs747*h;
const double clhs823 =             N[0]*N[2];
const double clhs824 =             4*clhs616;
const double clhs825 =             6*N[2]*clhs46*clhs6;
const double clhs826 =             -N[2]*clhs792 - clhs14*clhs825 + 3*clhs618 + clhs642 + clhs668 + clhs824;
const double clhs827 =             2*clhs629;
const double clhs828 =             DN(2,0)*clhs204;
const double clhs829 =             (4.0L/3.0L)*N[0]*clhs1*clhs109*clhs270*clhs599*h;
const double clhs830 =             (1.0L/9.0L)*clhs367*clhs599*clhs749*h*stab_c2;
const double clhs831 =             -clhs628*clhs751 + clhs629*clhs752;
const double clhs832 =             16*N[2]*clhs64*mu;
const double clhs833 =             12*N[2]*clhs64*mu;
const double clhs834 =             12*N[2]*clhs64*clhs75;
const double clhs835 =             N[2]*clhs50*clhs52*clhs64*mu;
const double clhs836 =             N[2]*clhs55*clhs57*clhs64*mu;
const double clhs837 =             4*N[2]*clhs64*clhs80;
const double clhs838 =             N[2]*clhs189*clhs64;
const double clhs839 =             N[2]*clhs100*clhs80;
const double clhs840 =             36*N[2]*clhs100*clhs14*clhs50;
const double clhs841 =             6*N[2]*clhs64;
const double clhs842 =             N[0]*f_ext(0,1) + N[1]*f_ext(1,1) + N[2]*f_ext(2,1);
const double clhs843 =             2*clhs133;
const double clhs844 =             6*N[0]*clhs42*clhs6;
const double clhs845 =             clhs132*clhs739 - clhs14*clhs844 + 3*clhs293 + 3*clhs701 - 2*clhs702 - clhs843;
const double clhs846 =             clhs14*clhs707;
const double clhs847 =             -4*clhs28 + clhs709 + clhs710;
const double clhs848 =             (2.0L/3.0L)*N[0]*clhs1*clhs32*clhs842*h;
const double clhs849 =             pow(clhs52, 2);
const double clhs850 =             pow(clhs57, 2);
const double clhs851 =             14*clhs44*clhs52*mu;
const double clhs852 =             2*clhs44*clhs52*mu;
const double clhs853 =             6*DN(0,1)*U(0,2) + 6*DN(1,1)*U(1,2) + 6*DN(2,1)*U(2,2);
const double clhs854 =             2*DN(0,0)*U(0,2) + 2*DN(1,0)*U(1,2) + 2*DN(2,0)*U(2,2);
const double clhs855 =             20*DN(0,1)*U(0,2) + 20*DN(1,1)*U(1,2) + 20*DN(2,1)*U(2,2);
const double clhs856 =             20*DN(0,0)*U(0,2) + 20*DN(1,0)*U(1,2) + 20*DN(2,0)*U(2,2);
const double clhs857 =             clhs14*clhs46*clhs64*mu;
const double clhs858 =             12*clhs10*clhs57*clhs64*clhs91;
const double clhs859 =             clhs107*clhs64*clhs853;
const double clhs860 =             -clhs10*clhs101*clhs57*clhs66 + 4*clhs10*clhs390*clhs57*mu - clhs10*clhs405*clhs52*clhs66 - 12*clhs10*clhs57*clhs88 + clhs14*clhs55*clhs64*clhs854*mu - clhs14*clhs84*clhs855 + clhs152 + clhs153 + clhs179*clhs72 - clhs183*clhs55 + clhs186 + clhs191*clhs725 + clhs191*clhs726 + 3*clhs207 + clhs300*clhs727 + clhs300*clhs728 - clhs303*clhs853 - clhs390*clhs50*clhs739*clhs95 - clhs390*clhs741 - clhs390*clhs744 - clhs42*clhs52*clhs558*clhs64*mu + clhs42*clhs851 + 14*clhs44*clhs57*clhs60 - clhs46*clhs743 + clhs46*clhs852 - clhs46*clhs858 - clhs46*clhs859 + clhs496*(U(0,2)*bdf0 + Un(0,2)*bdf1 + Unn(0,2)*bdf2) + clhs501*(U(1,2)*bdf0 + Un(1,2)*bdf1 + Unn(1,2)*bdf2) + clhs54*clhs57 - clhs55*clhs858 - clhs558*clhs849*clhs92 - clhs558*clhs98 + clhs660*(U(2,2)*bdf0 + Un(2,2)*bdf1 + Unn(2,2)*bdf2) - clhs71*clhs729 + clhs714 - clhs715*clhs842 + clhs717*clhs850 + clhs719*clhs849 - clhs720*clhs849 + 8*clhs721 + 6*clhs723 - clhs724*clhs850 - clhs732*clhs89 + clhs738*clhs87 - clhs742*clhs850 - clhs81*clhs854 + clhs856*clhs857 - clhs857*(16*DN(0,1)*U(0,2) + 16*DN(1,1)*U(1,2) + 16*DN(2,1)*U(2,2)) + 20*clhs86 - 20*clhs90;
const double clhs861 =             12*clhs57*clhs64*clhs75;
const double clhs862 =             clhs52*clhs64*clhs80;
const double clhs863 =             clhs546*clhs849;
const double clhs864 =             clhs100*clhs14*mu;
const double clhs865 =             48*clhs100*clhs42*mu;
const double clhs866 =             clhs100*clhs42*clhs546*mu;
const double clhs867 =             clhs100*clhs853*mu;
const double clhs868 =             60*clhs14*clhs46*mu;
const double clhs869 =             N[0]*clhs100*clhs52;
const double clhs870 =             6*clhs14*clhs55*mu;
const double clhs871 =             12*clhs100*clhs57*mu;
const double clhs872 =             clhs558*clhs850;
const double clhs873 =             36*clhs10*clhs100*clhs57*clhs75;
const double clhs874 =             12*clhs100*clhs42*clhs80;
const double clhs875 =             clhs175*clhs411*clhs46*mu;
const double clhs876 =             clhs206*clhs46;
const double clhs877 =             12*clhs415*clhs42;
const double clhs878 =             2*clhs458;
const double clhs879 =             6*N[1]*clhs42*clhs6;
const double clhs880 =             clhs46*clhs6*clhs739;
const double clhs881 =             N[1]*clhs880 - clhs14*clhs879 + 3*clhs442 - 2*clhs451 + clhs480 - clhs878;
const double clhs882 =             8*clhs14*clhs46*clhs6;
const double clhs883 =             2*clhs10*clhs14*clhs42*clhs6;
const double clhs884 =             DN(1,1)*clhs206;
const double clhs885 =             clhs46*clhs6*clhs847;
const double clhs886 =             20*clhs52*clhs64*mu;
const double clhs887 =             clhs52*clhs64*clhs75;
const double clhs888 =             clhs100*clhs52*mu;
const double clhs889 =             N[1]*clhs100*clhs52;
const double clhs890 =             12*clhs10*clhs100*clhs52*clhs80;
const double clhs891 =             (1.0L/18.0L)*clhs44*clhs747*clhs860*h;
const double clhs892 =             2*clhs626;
const double clhs893 =             6*N[2]*clhs42*clhs6;
const double clhs894 =             N[2]*clhs880 - clhs14*clhs893 + 3*clhs614 - 2*clhs620 + clhs644 - clhs892;
const double clhs895 =             DN(2,1)*clhs206;
const double clhs896 =             N[2]*clhs100*clhs52;
            lhs(0,0)=-bdf0*clhs0 - 2.0L/9.0L*clhs1*clhs109*clhs244*clhs32*clhs44*h + clhs269*clhs368 + clhs419*clhs420;
            lhs(0,1)=-N[1]*clhs421 + clhs368*clhs426 + clhs420*clhs490 - clhs491*clhs594;
            lhs(0,2)=-N[2]*clhs421 + clhs368*clhs599 + clhs420*clhs651 - clhs491*clhs693;
            lhs(1,0)=clhs0*clhs694 + clhs419*clhs713 + clhs694*clhs712 + clhs695*clhs699 + clhs700*(-clhs10*clhs698 + clhs128*clhs704 + 3*clhs133 + clhs699 - 2*clhs701 + 3*clhs702 - clhs703) - clhs705*(N[0]*clhs14*clhs46*clhs6*clhs63 - clhs10*clhs706 - clhs128*clhs711 - clhs133*clhs704 + clhs290*clhs708 + clhs333*clhs703 + clhs352 - clhs702*clhs704) - clhs745*clhs748 + clhs745*clhs750 + clhs789*(-48*N[0]*clhs14*clhs405*clhs55 + clhs100*clhs14*clhs279*clhs50*mu - 60*clhs100*clhs337*clhs42*mu + clhs133*clhs757 + clhs133*clhs765 + clhs133*clhs766 - clhs133*clhs769 + clhs133*clhs772 - clhs133*clhs778 - clhs133*clhs788 - clhs136*clhs295 + clhs14*clhs278*clhs55*mu - clhs145*clhs295 - clhs147*clhs192 - clhs162*clhs337*clhs411*clhs42*mu + clhs189*clhs409 - clhs271*clhs50 - clhs272*clhs55 - clhs273*clhs55 + clhs274*clhs50 + clhs275*clhs50 + clhs276*clhs55 + clhs277*clhs55 + clhs282*clhs737 + clhs283*clhs55 + clhs290*clhs759 + clhs290*clhs763 + clhs290*clhs764 - clhs290*clhs767 - clhs290*clhs770 - clhs290*clhs773 + clhs290*clhs776 - clhs290*clhs787 - clhs293*clhs762 + clhs306 - clhs310*clhs771 + clhs310*clhs774 + clhs310*clhs775 - clhs310*clhs785 + clhs312*clhs777 + clhs317 + clhs319 - clhs321*clhs733 + clhs325*clhs735 - clhs330 - clhs347 - clhs353 + clhs357*clhs786 + 20*clhs374 - 20*clhs376 + 28*clhs378 + 28*clhs383 - clhs389*clhs50 - clhs391*clhs55 - clhs395*clhs50 - clhs398*clhs55 + clhs400*clhs768 + clhs408*clhs553 + clhs46*clhs784 + clhs496*clhs694 + clhs57*clhs784 + clhs716*clhs754 - clhs716*clhs761 + clhs718*clhs755 - clhs718*clhs756 + clhs753 + 12*clhs758 + 16*clhs760 + clhs779*clhs780 + clhs781*clhs782 + clhs782*clhs783);
            lhs(1,1)=clhs490*clhs713 - clhs594*clhs822 + clhs694*clhs790 + clhs694*clhs800 + clhs695*clhs794 + clhs700*(N[1]*clhs796 - clhs10*clhs793 - 2*clhs442 + 3*clhs451 + clhs475 + clhs794 - clhs795) - clhs705*(N[1]*clhs797 - N[1]*clhs799 - clhs10*clhs529 + clhs333*clhs795 - clhs451*clhs704 - clhs458*clhs704 + clhs460*clhs708 + clhs798) + clhs745*clhs801 + clhs789*(-DN(1,1)*clhs722 + clhs189*clhs485 - clhs192*clhs438 - clhs295*clhs427 - clhs295*clhs432 - clhs295*clhs437 + clhs316*clhs451 + clhs316*clhs458 + clhs318*clhs459 - clhs329*clhs460 - clhs351*clhs798 - clhs392*clhs454 - 2*clhs430 + clhs440*clhs449 + clhs441*clhs449 + clhs443*clhs454 + clhs443*clhs455 + 20*clhs444 - clhs445*clhs752 - 20*clhs450 + clhs451*clhs806 + clhs451*clhs812 + 28*clhs456 + clhs458*clhs757 + clhs458*clhs765 + clhs458*clhs766 - clhs458*clhs769 + clhs458*clhs772 - clhs458*clhs778 - clhs458*clhs788 - clhs458*clhs813 - clhs458*clhs814 + clhs458*clhs817 + clhs458*clhs818 + clhs458*clhs819 - clhs458*clhs820 - clhs46*clhs810 + clhs460*clhs759 + clhs460*clhs763 + clhs460*clhs764 - clhs460*clhs767 - clhs460*clhs770 - clhs460*clhs773 + clhs460*clhs776 - clhs460*clhs787 - clhs461*clhs762 + 28*clhs462 - clhs466*clhs50 - clhs470*clhs50 - clhs471*clhs55 + clhs472*clhs768 + clhs473*clhs811 + clhs487*clhs737 + clhs488*clhs55 + clhs501*clhs694 - clhs551*clhs771 + clhs551*clhs774 + clhs551*clhs775 - clhs551*clhs785 + clhs552*clhs777 + clhs716*clhs803 - clhs716*clhs809 + clhs718*clhs804 - clhs718*clhs805 + clhs780*clhs815 + clhs781*clhs816 + clhs783*clhs816 + clhs786*clhs821 + clhs802 + 12*clhs807 + 16*clhs808);
            lhs(1,2)=clhs651*clhs713 - clhs693*clhs822 + clhs694*clhs823 + clhs694*clhs829 + clhs695*clhs826 + clhs700*(N[2]*clhs796 - clhs10*clhs825 - 2*clhs614 + 3*clhs620 + clhs641 + clhs826 - clhs827) - clhs705*(N[2]*clhs797 - N[2]*clhs799 - clhs10*clhs670 + clhs333*clhs827 - clhs620*clhs704 - clhs626*clhs704 + clhs628*clhs708 + clhs828) + clhs745*clhs830 + clhs789*(-DN(2,1)*clhs722 - clhs14*clhs46*clhs838 + clhs189*clhs647 - clhs192*clhs610 - clhs295*clhs609 + clhs316*clhs620 + clhs316*clhs626 + clhs318*clhs627 - clhs329*clhs628 - clhs351*clhs828 - clhs392*clhs622 - clhs434*clhs600 - clhs434*clhs605 + clhs443*clhs622 + clhs443*clhs623 + clhs449*clhs612 + clhs449*clhs613 - clhs50*clhs634 - clhs50*clhs636 - clhs55*clhs637 + clhs55*clhs649 - 2*clhs604 + 20*clhs615 - clhs616*clhs752 - 20*clhs619 + clhs620*clhs806 + clhs620*clhs812 + 28*clhs624 + clhs626*clhs757 + clhs626*clhs765 + clhs626*clhs766 - clhs626*clhs769 + clhs626*clhs772 - clhs626*clhs778 - clhs626*clhs788 - clhs626*clhs813 - clhs626*clhs814 + clhs626*clhs817 + clhs626*clhs818 + clhs626*clhs819 - clhs626*clhs820 + clhs628*clhs759 + clhs628*clhs763 + clhs628*clhs764 - clhs628*clhs767 - clhs628*clhs770 - clhs628*clhs773 + clhs628*clhs776 - clhs628*clhs787 - clhs629*clhs762 + 28*clhs630 + clhs638*clhs768 + clhs639*clhs811 + clhs648*clhs737 + clhs660*clhs694 - clhs681*clhs771 + clhs681*clhs774 + clhs681*clhs775 - clhs681*clhs785 + clhs682*clhs777 + clhs716*clhs832 - clhs716*clhs837 + clhs718*clhs833 - clhs718*clhs834 + clhs780*clhs839 + clhs781*clhs840 + clhs783*clhs840 + clhs786*clhs841 + clhs831 + 12*clhs835 + 16*clhs836);
            lhs(2,0)=clhs0*clhs842 - clhs245*(2*clhs10*clhs128*clhs14 - clhs132*clhs847 + clhs135*clhs846 - clhs14*clhs358 - clhs310 - clhs312 + clhs336*clhs843 + clhs354) + clhs419*clhs848 + clhs695*clhs845 + clhs700*(-clhs10*clhs844 - clhs132*clhs69 + clhs359 + 4*clhs361 + clhs364 + clhs706 + clhs845) + clhs712*clhs842 - clhs748*clhs860 + clhs750*clhs860 + clhs789*(6*DN(0,0)*clhs10*clhs57*clhs864 + 6*DN(0,1)*clhs10*clhs52*clhs864 + 12*N[0]*clhs10*clhs100*clhs50*clhs52*clhs80 - clhs133*clhs338 - clhs133*clhs861 + clhs133*clhs873 + clhs135*clhs757 + clhs135*clhs765 + clhs135*clhs766 + clhs135*clhs772 - clhs135*clhs778 - clhs135*clhs788 - clhs135*clhs813 - clhs135*clhs814 + clhs135*clhs817 - clhs135*clhs820 - clhs136*clhs191 - clhs145*clhs191 - clhs147*clhs300 + clhs161 - clhs271*clhs52 - clhs272*clhs57 - clhs273*clhs57 + clhs274*clhs52 + clhs275*clhs52 + clhs276*clhs57 + clhs277*clhs57 + clhs282*clhs390 + clhs283*clhs57 + clhs293*clhs759 + clhs293*clhs763 + clhs293*clhs764 - clhs293*clhs866 + clhs310*clhs867 + clhs311 - clhs312*clhs774 - clhs312*clhs865 + clhs312*clhs871 + clhs312*clhs874 - clhs312*clhs875 + clhs313 + clhs322*clhs855 - clhs326*clhs856 - clhs332 - clhs346 - clhs355 + clhs357*clhs876 - 20*clhs375 + 20*clhs377 + 28*clhs379 + 28*clhs384 - clhs391*clhs57 - clhs398*clhs57 + clhs400*clhs863 + clhs407*clhs546 + clhs410*clhs558 + clhs496*clhs842 - clhs696*clhs862 - clhs697*clhs862 - clhs701*clhs767 + clhs753 + clhs754*clhs850 + clhs755*clhs849 - clhs756*clhs849 + 16*clhs758 + 12*clhs760 - clhs761*clhs850 + clhs779*clhs872 + clhs868*clhs869 + clhs869*clhs870 - clhs869*clhs877);
            lhs(2,1)=-clhs245*(N[1]*clhs883 - N[1]*clhs885 + clhs336*clhs878 + clhs459*clhs846 - clhs532 - clhs551 - clhs552 + clhs884) + clhs490*clhs848 - clhs594*clhs891 + clhs695*clhs881 + clhs700*(-N[1]*clhs882 - clhs10*clhs879 + 3*clhs445 + 4*clhs448 + clhs478 + clhs529 + clhs881) + clhs789*(-DN(1,0)*clhs851 - DN(1,1)*clhs852 + DN(1,1)*clhs859 - clhs191*clhs437 - clhs300*clhs429 - clhs300*clhs435 - clhs300*clhs438 + clhs309*clhs532 + clhs309*clhs551 + clhs309*clhs552 - clhs331*clhs459 - clhs338*clhs458 - clhs351*clhs884 + clhs390*clhs487 - clhs42*clhs810 + clhs440*clhs452 + clhs441*clhs452 - clhs441*clhs887 + clhs442*clhs751 + clhs446*clhs454 + clhs446*clhs455 - 20*clhs447 - clhs448*clhs886 + 20*clhs453 + 28*clhs457 - clhs458*clhs861 + clhs458*clhs873 + clhs459*clhs757 + clhs459*clhs765 + clhs459*clhs766 + clhs459*clhs772 - clhs459*clhs778 - clhs459*clhs788 - clhs459*clhs813 - clhs459*clhs814 + clhs459*clhs817 + clhs459*clhs818 - clhs459*clhs820 + clhs460*clhs890 + clhs461*clhs759 + clhs461*clhs763 + clhs461*clhs764 - clhs461*clhs866 + 28*clhs463 - clhs468*clhs57 - clhs471*clhs57 + clhs472*clhs863 + clhs474*clhs888 + clhs483*clhs546 + clhs501*clhs842 - clhs524*clhs862 + clhs532*clhs867 + clhs551*clhs867 - clhs552*clhs774 - clhs552*clhs865 + clhs552*clhs871 + clhs552*clhs874 - clhs552*clhs875 - clhs791*clhs862 + clhs802 + clhs803*clhs850 + clhs804*clhs849 - clhs805*clhs849 + 16*clhs807 + 12*clhs808 - clhs809*clhs850 + clhs815*clhs872 + clhs821*clhs876 + clhs868*clhs889 + clhs870*clhs889 - clhs877*clhs889) + clhs790*clhs842 + clhs800*clhs842 + clhs801*clhs860;
            lhs(2,2)=-clhs245*(N[2]*clhs883 - N[2]*clhs885 + clhs336*clhs892 + clhs627*clhs846 - clhs671 - clhs681 - clhs682 + clhs895) + clhs651*clhs848 - clhs693*clhs891 + clhs695*clhs894 + clhs700*(-N[2]*clhs882 - clhs10*clhs893 + 3*clhs616 + 4*clhs618 + clhs643 + clhs670 + clhs894) + clhs789*(-DN(2,0)*clhs851 - DN(2,1)*clhs852 + DN(2,1)*clhs859 - clhs14*clhs42*clhs838 - clhs191*clhs609 - clhs300*clhs603 - clhs300*clhs607 - clhs300*clhs610 + clhs309*clhs671 + clhs309*clhs681 + clhs309*clhs682 - clhs331*clhs627 - clhs338*clhs626 - clhs351*clhs895 + clhs390*clhs648 + clhs446*clhs622 + clhs446*clhs623 + clhs452*clhs612 + clhs452*clhs613 + clhs546*clhs646 - clhs57*clhs635 - clhs57*clhs637 - clhs613*clhs887 + clhs614*clhs751 - 20*clhs617 - clhs618*clhs886 + 20*clhs621 + 28*clhs625 - clhs626*clhs861 + clhs626*clhs873 + clhs627*clhs757 + clhs627*clhs765 + clhs627*clhs766 + clhs627*clhs772 - clhs627*clhs778 - clhs627*clhs788 - clhs627*clhs813 - clhs627*clhs814 + clhs627*clhs817 + clhs627*clhs818 - clhs627*clhs820 + clhs628*clhs890 + clhs629*clhs759 + clhs629*clhs763 + clhs629*clhs764 - clhs629*clhs866 + 28*clhs631 + clhs638*clhs863 + clhs640*clhs888 + clhs660*clhs842 - clhs668*clhs862 + clhs671*clhs867 + clhs681*clhs867 - clhs682*clhs774 - clhs682*clhs865 + clhs682*clhs871 + clhs682*clhs874 - clhs682*clhs875 - clhs824*clhs862 + clhs831 + clhs832*clhs850 + clhs833*clhs849 - clhs834*clhs849 + 16*clhs835 + 12*clhs836 - clhs837*clhs850 + clhs839*clhs872 + clhs841*clhs876 + clhs868*clhs896 + clhs870*clhs896 - clhs877*clhs896) + clhs823*clhs842 + clhs829*clhs842 + clhs830*clhs860;


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
    const bounded_matrix<double,dim,dimes> grad_U = prod(trans(DN), U);
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
    const bounded_matrix<double,dim,dimes> grad_U = prod(trans(DN), U);
    const double& r_gauss = inner_prod(data.N, data.r);
    
    const array_1d<double,dimes> accel_gauss = bdf0*U_gauss+bdf1*prod(trans(Un), N)+bdf2*prod(trans(Unn), N);
   
    // Stabilization parameters
    const double stab_c1 = 4.0;
    const double stab_c2 = 2.0;

    const double crhs0 =             DN(0,0)*U(0,1);
const double crhs1 =             DN(1,0)*U(1,1);
const double crhs2 =             DN(2,0)*U(2,1);
const double crhs3 =             crhs0 + crhs1 + crhs2;
const double crhs4 =             DN(0,1)*U(0,2);
const double crhs5 =             DN(1,1)*U(1,2);
const double crhs6 =             DN(2,1)*U(2,2);
const double crhs7 =             N[0]*(U(0,0)*bdf0 + Un(0,0)*bdf1 + Unn(0,0)*bdf2);
const double crhs8 =             N[1]*(U(1,0)*bdf0 + Un(1,0)*bdf1 + Unn(1,0)*bdf2);
const double crhs9 =             N[2]*(U(2,0)*bdf0 + Un(2,0)*bdf1 + Unn(2,0)*bdf2);
const double crhs10 =             1.0/stab_c2;
const double crhs11 =             N[0]*U(0,0) + N[1]*U(1,0) + N[2]*U(2,0);
const double crhs12 =             1.0/crhs11;
const double crhs13 =             N[0]*U(0,1) + N[1]*U(1,1) + N[2]*U(2,1);
const double crhs14 =             N[0]*U(0,2) + N[1]*U(1,2) + N[2]*U(2,2);
const double crhs15 =             y - 1;
const double crhs16 =             2*(N[0]*U(0,3));
const double crhs17 =             -crhs16;
const double crhs18 =             2*(N[1]*U(1,3));
const double crhs19 =             -crhs18;
const double crhs20 =             2*(N[2]*U(2,3));
const double crhs21 =             -crhs20;
const double crhs22 =             pow(crhs13, 2);
const double crhs23 =             crhs12*crhs22;
const double crhs24 =             pow(crhs14, 2);
const double crhs25 =             crhs12*crhs24;
const double crhs26 =             crhs12*crhs15*y*(crhs17 + crhs19 + crhs21 + crhs23 + crhs25) - 2*fabs(crhs12*(crhs13 + crhs14));
const double crhs27 =             1.0/crhs26;
const double crhs28 =             3*crhs0;
const double crhs29 =             3*crhs4;
const double crhs30 =             3*crhs1;
const double crhs31 =             3*crhs5;
const double crhs32 =             3*crhs2;
const double crhs33 =             3*crhs6;
const double crhs34 =             DN(0,0)*U(0,0) + DN(1,0)*U(1,0) + DN(2,0)*U(2,0);
const double crhs35 =             pow(crhs34, 2);
const double crhs36 =             pow(crhs11, -2);
const double crhs37 =             7*crhs36*mu;
const double crhs38 =             DN(0,1)*U(0,0) + DN(1,1)*U(1,0) + DN(2,1)*U(2,0);
const double crhs39 =             pow(crhs38, 2);
const double crhs40 =             1.0/cv;
const double crhs41 =             3*crhs36*crhs40*lambda;
const double crhs42 =             4*crhs36*mu;
const double crhs43 =             crhs3*crhs34;
const double crhs44 =             3*DN(0,0)*U(0,0) + 3*DN(1,0)*U(1,0) + 3*DN(2,0)*U(2,0);
const double crhs45 =             DN(0,0)*U(0,2);
const double crhs46 =             DN(1,0)*U(1,2);
const double crhs47 =             DN(2,0)*U(2,2);
const double crhs48 =             crhs45 + crhs46 + crhs47;
const double crhs49 =             crhs36*crhs48*mu;
const double crhs50 =             2*DN(0,1)*U(0,0) + 2*DN(1,1)*U(1,0) + 2*DN(2,1)*U(2,0);
const double crhs51 =             DN(0,1)*U(0,1);
const double crhs52 =             DN(1,1)*U(1,1);
const double crhs53 =             DN(2,1)*U(2,1);
const double crhs54 =             crhs51 + crhs52 + crhs53;
const double crhs55 =             5*crhs34*crhs36*mu;
const double crhs56 =             crhs4 + crhs5 + crhs6;
const double crhs57 =             5*DN(0,1)*U(0,0) + 5*DN(1,1)*U(1,0) + 5*DN(2,1)*U(2,0);
const double crhs58 =             crhs3*crhs36*mu;
const double crhs59 =             3*DN(0,1)*U(0,0) + 3*DN(1,1)*U(1,0) + 3*DN(2,1)*U(2,0);
const double crhs60 =             crhs36*crhs54*mu;
const double crhs61 =             crhs38*crhs56;
const double crhs62 =             pow(crhs11, -3);
const double crhs63 =             crhs35*crhs62*mu;
const double crhs64 =             6*N[0]*U(0,2) + 6*N[1]*U(1,2) + 6*N[2]*U(2,2);
const double crhs65 =             6*N[0]*U(0,1) + 6*N[1]*U(1,1) + 6*N[2]*U(2,1);
const double crhs66 =             crhs39*crhs62*mu;
const double crhs67 =             DN(0,0)*U(0,3) + DN(1,0)*U(1,3) + DN(2,0)*U(2,3);
const double crhs68 =             DN(0,1)*U(0,3) + DN(1,1)*U(1,3) + DN(2,1)*U(2,3);
const double crhs69 =             crhs40*lambda;
const double crhs70 =             -crhs69 + mu;
const double crhs71 =             crhs36*crhs48*crhs70;
const double crhs72 =             crhs62*crhs64*crhs70;
const double crhs73 =             4*mu;
const double crhs74 =             3*crhs69;
const double crhs75 =             crhs73 - crhs74;
const double crhs76 =             crhs36*crhs75;
const double crhs77 =             crhs14*crhs54;
const double crhs78 =             crhs34*crhs62*mu;
const double crhs79 =             crhs77*crhs78;
const double crhs80 =             crhs13*crhs56;
const double crhs81 =             crhs78*crhs80;
const double crhs82 =             crhs14*crhs3;
const double crhs83 =             crhs38*crhs62*mu;
const double crhs84 =             crhs82*crhs83;
const double crhs85 =             crhs13*crhs48;
const double crhs86 =             crhs83*crhs85;
const double crhs87 =             crhs34*crhs48;
const double crhs88 =             6*DN(0,1)*U(0,1) + 6*DN(1,1)*U(1,1) + 6*DN(2,1)*U(2,1);
const double crhs89 =             crhs38*crhs88;
const double crhs90 =             crhs62*crhs75;
const double crhs91 =             2*DN(0,0)*U(0,0) + 2*DN(1,0)*U(1,0) + 2*DN(2,0)*U(2,0);
const double crhs92 =             crhs38*crhs64;
const double crhs93 =             pow(crhs11, -4);
const double crhs94 =             3*mu;
const double crhs95 =             crhs16 + crhs18 + crhs20;
const double crhs96 =             -crhs23*crhs74 - crhs25*crhs74 + crhs69*crhs95;
const double crhs97 =             crhs23*crhs73 + crhs25*crhs94 + crhs96;
const double crhs98 =             crhs62*crhs97;
const double crhs99 =             crhs62*(crhs23*crhs94 + crhs25*crhs73 + crhs96);
const double crhs100 =             crhs13*crhs3*crhs90*crhs91 - crhs13*crhs34*crhs92*crhs93*mu + 2*crhs13*crhs35*crhs62*crhs75 + crhs13*crhs62*crhs70*crhs89 + 2*crhs14*crhs39*crhs62*crhs75 + crhs14*crhs50*crhs56*crhs90 + crhs28 + crhs29 + crhs30 + crhs31 + crhs32 + crhs33 + crhs34*crhs36*crhs50*mu + crhs34*crhs41*crhs67 + crhs35*crhs37 + crhs35*crhs41 + crhs35*crhs72 - 3*crhs35*crhs98 - crhs36*crhs54*crhs59*crhs70 + crhs37*crhs39 + crhs38*crhs41*crhs68 + crhs39*crhs41 + crhs39*crhs62*crhs65*crhs70 - 3*crhs39*crhs99 + crhs42*crhs43 + crhs42*crhs61 - crhs43*crhs76 + crhs44*crhs49 - crhs44*crhs71 + crhs49*crhs57 + crhs54*crhs55 - crhs55*crhs56 - crhs57*crhs58 + crhs59*crhs60 - crhs61*crhs76 - crhs63*crhs64 - crhs63*(8*N[0]*U(0,1) + 8*N[1]*U(1,1) + 8*N[2]*U(2,1)) - crhs65*crhs66 - crhs66*(8*N[0]*U(0,2) + 8*N[1]*U(1,2) + 8*N[2]*U(2,2)) + 3*crhs7 + crhs72*crhs87 + crhs79 + 3*crhs8 + crhs81 + crhs84 + crhs86 + 3*crhs9;
const double crhs101 =             6*DN(0,0);
const double crhs102 =             6*DN(0,1);
const double crhs103 =             crhs101*crhs15;
const double crhs104 =             crhs102*crhs15;
const double crhs105 =             6*DN(0,0)*crhs12;
const double crhs106 =             crhs105*crhs13;
const double crhs107 =             6*DN(0,1)*crhs12;
const double crhs108 =             crhs107*crhs14;
const double crhs109 =             6*N[0]*crhs12;
const double crhs110 =             crhs109*crhs3;
const double crhs111 =             crhs109*crhs48;
const double crhs112 =             crhs109*crhs54;
const double crhs113 =             crhs109*crhs56;
const double crhs114 =             14*DN(0,0)*U(0,0) + 14*DN(1,0)*U(1,0) + 14*DN(2,0)*U(2,0);
const double crhs115 =             DN(0,0)*crhs36*mu;
const double crhs116 =             8*DN(0,0)*U(0,1) + 8*DN(1,0)*U(1,1) + 8*DN(2,0)*U(2,1);
const double crhs117 =             10*DN(0,0)*crhs36*mu;
const double crhs118 =             DN(0,1)*crhs36*mu;
const double crhs119 =             10*DN(0,0)*U(0,1) + 10*DN(1,0)*U(1,1) + 10*DN(2,0)*U(2,1);
const double crhs120 =             10*crhs36*crhs48*mu;
const double crhs121 =             14*DN(0,1)*U(0,0) + 14*DN(1,1)*U(1,0) + 14*DN(2,1)*U(2,0);
const double crhs122 =             6*DN(0,1)*crhs54;
const double crhs123 =             8*DN(0,1)*U(0,2) + 8*DN(1,1)*U(1,2) + 8*DN(2,1)*U(2,2);
const double crhs124 =             6*N[0]*crhs12*y;
const double crhs125 =             y - 3;
const double crhs126 =             crhs12*crhs14;
const double crhs127 =             crhs12*crhs13;
const double crhs128 =             6*DN(0,0)*crhs36*crhs40*lambda;
const double crhs129 =             6*DN(0,1)*crhs36*crhs40*lambda;
const double crhs130 =             2*crhs3*crhs36*crhs75;
const double crhs131 =             crhs13*crhs14*crhs36;
const double crhs132 =             2*crhs36*crhs56*crhs75;
const double crhs133 =             6*N[0]*crhs34*crhs36;
const double crhs134 =             crhs13*crhs133;
const double crhs135 =             crhs133*crhs14;
const double crhs136 =             6*N[0]*crhs14*crhs36;
const double crhs137 =             crhs136*crhs3;
const double crhs138 =             6*N[0]*crhs13*crhs36;
const double crhs139 =             crhs138*crhs48;
const double crhs140 =             crhs138*crhs38;
const double crhs141 =             crhs136*crhs38;
const double crhs142 =             crhs138*crhs56;
const double crhs143 =             16*N[0]*U(0,1) + 16*N[1]*U(1,1) + 16*N[2]*U(2,1);
const double crhs144 =             DN(0,0)*crhs34*crhs62*mu;
const double crhs145 =             12*N[0]*U(0,2) + 12*N[1]*U(1,2) + 12*N[2]*U(2,2);
const double crhs146 =             20*DN(0,0)*crhs38*crhs62*mu;
const double crhs147 =             2*DN(0,0)*crhs62*mu;
const double crhs148 =             20*DN(0,1)*crhs34*crhs62*mu;
const double crhs149 =             2*DN(0,1)*crhs62*mu;
const double crhs150 =             12*N[0]*U(0,1) + 12*N[1]*U(1,1) + 12*N[2]*U(2,1);
const double crhs151 =             DN(0,1)*crhs38*crhs62*mu;
const double crhs152 =             16*N[0]*U(0,2) + 16*N[1]*U(1,2) + 16*N[2]*U(2,2);
const double crhs153 =             18*N[0]*crhs15*crhs36;
const double crhs154 =             crhs13*crhs3;
const double crhs155 =             crhs14*crhs56;
const double crhs156 =             DN(0,0)*crhs145*crhs62*crhs70;
const double crhs157 =             DN(0,1)*crhs150*crhs62*crhs70;
const double crhs158 =             crhs13*crhs34;
const double crhs159 =             4*DN(0,0)*crhs62*crhs75;
const double crhs160 =             4*DN(0,1)*crhs62*crhs75;
const double crhs161 =             N[0]*crhs145*crhs62;
const double crhs162 =             crhs13*crhs14*crhs38*crhs93*mu;
const double crhs163 =             crhs13*crhs14*crhs34*crhs93*mu;
const double crhs164 =             N[0]*crhs145*crhs15*crhs62;
const double crhs165 =             crhs13*crhs54;
const double crhs166 =             3*DN(0,0)*crhs36;
const double crhs167 =             crhs15*crhs22;
const double crhs168 =             crhs15*crhs24;
const double crhs169 =             crhs167 + crhs168;
const double crhs170 =             crhs169 - 2*crhs22;
const double crhs171 =             3*DN(0,1)*crhs36;
const double crhs172 =             crhs169 - 2*crhs24;
const double crhs173 =             6*N[0]*crhs62;
const double crhs174 =             3*crhs12;
const double crhs175 =             crhs12*crhs167;
const double crhs176 =             -crhs15*crhs95;
const double crhs177 =             crhs12*crhs168;
const double crhs178 =             crhs17 + crhs176 + crhs177 + crhs19 + crhs21;
const double crhs179 =             crhs17 + crhs175 + crhs176 + crhs19 + crhs21;
const double crhs180 =             crhs12*crhs15*(crhs22 + crhs24);
const double crhs181 =             crhs17 + crhs175 + crhs176 + crhs180 + crhs19 + crhs21;
const double crhs182 =             crhs177 + crhs181;
const double crhs183 =             N[0]*crhs36*(5*crhs175 + crhs178 + crhs180);
const double crhs184 =             N[0]*crhs36*(5*crhs177 + crhs181);
const double crhs185 =             N[0]*crhs62*(crhs177 + crhs179 + 2*crhs180);
const double crhs186 =             DN(0,0)*crhs130 - DN(0,0)*crhs174*(3*crhs175 + crhs178) + 6*DN(0,0)*crhs34*crhs98 + DN(0,1)*crhs120 + DN(0,1)*crhs132 - DN(0,1)*crhs174*(3*crhs177 + crhs179) + 6*DN(0,1)*crhs38*crhs99 + N[0]*crhs13*crhs145*crhs15*crhs34*crhs62 + N[0]*crhs13*crhs145*crhs15*crhs38*crhs62 - N[0]*crhs14*crhs15*crhs36*crhs88 - N[0]*crhs14*crhs36*crhs88 - crhs101*crhs131 + crhs101*crhs162 - crhs101*crhs49 + crhs101*crhs71 + crhs101 - crhs102*crhs131 + crhs102*crhs163 + crhs102 - crhs103*crhs126 - crhs103*crhs131 + crhs103 - crhs104*crhs127 - crhs104*crhs131 + crhs104 + crhs105*crhs14 - crhs106*crhs125 + crhs106*y + crhs106 + crhs107*crhs13 - crhs108*crhs125 + crhs108*y + crhs108 - crhs110*crhs125 + crhs110*y + crhs110 - crhs111*crhs15 + crhs111 - crhs112*crhs15 + crhs112 - crhs113*crhs125 + crhs113*y + crhs113 - crhs114*crhs115 - crhs115*crhs116 - crhs115*crhs50 + crhs117*crhs54 - crhs117*crhs56 - crhs118*crhs119 - crhs118*crhs121 - crhs118*crhs123 - crhs118*crhs91 + crhs122*crhs36*crhs70 - crhs122*crhs36*mu + crhs124*crhs67 + crhs124*crhs68 + crhs125*crhs134 + crhs125*crhs136*crhs56 + crhs125*crhs138*crhs3 + crhs125*crhs141 - crhs128*crhs34 - crhs128*crhs67 - crhs129*crhs38 - crhs129*crhs68 - crhs13*crhs146 + crhs13*crhs148 + crhs13*crhs161*crhs38 + crhs13*crhs166*crhs182 - crhs134*y - crhs134 + crhs135*crhs15 - crhs135 - crhs136*crhs68*y - crhs137*crhs15 - crhs137 - crhs138*crhs67*y - crhs139*crhs15 - crhs139 + crhs14*crhs146 - crhs14*crhs148 - crhs14*crhs160*crhs38 + crhs14*crhs171*crhs182 + crhs140*crhs15 - crhs140 - crhs141*y - crhs141 - crhs142*crhs15 - crhs142 + crhs143*crhs144 + crhs144*crhs145 - crhs147*crhs77 - crhs147*crhs80 - crhs149*crhs82 - crhs149*crhs85 + crhs150*crhs151 + crhs151*crhs152 - crhs153*crhs154 - crhs153*crhs155 - crhs154*crhs159 - crhs155*crhs160 - crhs156*crhs34 - crhs156*crhs48 - crhs157*crhs38 - crhs157*crhs54 - crhs158*crhs159 + crhs158*crhs161 + crhs164*crhs165 + crhs164*crhs85 + crhs166*crhs170 - crhs170*crhs173*crhs34 + crhs171*crhs172 - crhs172*crhs173*crhs38 + crhs183*crhs44 + crhs183*(crhs28 + crhs30 + crhs32) + crhs184*crhs59 + crhs184*(crhs29 + crhs31 + crhs33) - crhs185*crhs34*crhs65 - crhs185*crhs92;
const double crhs187 =             N[0]*f_ext(0,0) + N[1]*f_ext(1,0) + N[2]*f_ext(2,0);
const double crhs188 =             N[0]*crhs11;
const double crhs189 =             N[0]*(U(0,1)*bdf0 + Un(0,1)*bdf1 + Unn(0,1)*bdf2);
const double crhs190 =             N[1]*(U(1,1)*bdf0 + Un(1,1)*bdf1 + Unn(1,1)*bdf2);
const double crhs191 =             N[2]*(U(2,1)*bdf0 + Un(2,1)*bdf1 + Unn(2,1)*bdf2);
const double crhs192 =             (1.0L/3.0L)*DN(0,0)*crhs12*mu;
const double crhs193 =             crhs12*crhs13*crhs34;
const double crhs194 =             3*crhs12*crhs14;
const double crhs195 =             -4*crhs0 - 4*crhs1 + 4*crhs193 + crhs194*crhs38 - 4*crhs2 - crhs29 - crhs31 - crhs33;
const double crhs196 =             (1.0L/3.0L)*DN(0,1)*crhs12*mu;
const double crhs197 =             crhs15*crhs67;
const double crhs198 =             4*crhs12*crhs14;
const double crhs199 =             crhs12*crhs13*crhs56;
const double crhs200 =             2*DN(0,0)*U(0,1) + 2*DN(1,0)*U(1,1) + 2*DN(2,0)*U(2,1);
const double crhs201 =             crhs12*crhs125*crhs13;
const double crhs202 =             crhs12*crhs14*crhs15*crhs48;
const double crhs203 =             crhs13*crhs14*crhs36*crhs38;
const double crhs204 =             crhs170*crhs36;
const double crhs205 =             (2.0L/3.0L)*N[0]*crhs10*crhs100*crhs27*h;
const double crhs206 =             crhs36*crhs56*mu;
const double crhs207 =             -crhs119*crhs206 + crhs120*crhs54;
const double crhs208 =             6*N[0]*U(0,0) + 6*N[1]*U(1,0) + 6*N[2]*U(2,0);
const double crhs209 =             pow(crhs3, 2);
const double crhs210 =             8*crhs36*mu;
const double crhs211 =             pow(crhs54, 2);
const double crhs212 =             6*crhs36*mu;
const double crhs213 =             6*crhs36*crhs70;
const double crhs214 =             crhs36*crhs91*mu;
const double crhs215 =             6*DN(0,0)*U(0,1) + 6*DN(1,0)*U(1,1) + 6*DN(2,0)*U(2,1);
const double crhs216 =             2*crhs36*crhs75;
const double crhs217 =             6*crhs36*crhs40*lambda;
const double crhs218 =             crhs36*crhs40*crhs88*lambda;
const double crhs219 =             crhs150*crhs62*crhs70;
const double crhs220 =             crhs3*crhs34*crhs62*mu;
const double crhs221 =             20*crhs34*crhs62*mu;
const double crhs222 =             20*crhs38*crhs62*mu;
const double crhs223 =             crhs14*crhs54*crhs62;
const double crhs224 =             crhs62*mu;
const double crhs225 =             2*crhs51;
const double crhs226 =             2*crhs52;
const double crhs227 =             2*crhs53;
const double crhs228 =             crhs225 + crhs226 + crhs227;
const double crhs229 =             4*crhs13*crhs62*crhs75;
const double crhs230 =             crhs145*crhs62*crhs70;
const double crhs231 =             crhs54*crhs56;
const double crhs232 =             4*crhs14*crhs62*crhs75;
const double crhs233 =             (1.0L/18.0L)*crhs186*h/(-crhs26*stab_c2 + 2*nu*stab_c1/h);
const double crhs234 =             N[0]*f_ext(0,1) + N[1]*f_ext(1,1) + N[2]*f_ext(2,1);
const double crhs235 =             N[0]*(U(0,2)*bdf0 + Un(0,2)*bdf1 + Unn(0,2)*bdf2);
const double crhs236 =             N[1]*(U(1,2)*bdf0 + Un(1,2)*bdf1 + Unn(1,2)*bdf2);
const double crhs237 =             N[2]*(U(2,2)*bdf0 + Un(2,2)*bdf1 + Unn(2,2)*bdf2);
const double crhs238 =             -crhs127*crhs50 + crhs194*crhs34 + crhs225 + crhs226 + crhs227 - 3*crhs45 - 3*crhs46 - 3*crhs47;
const double crhs239 =             crhs15*crhs68;
const double crhs240 =             crhs126*crhs3;
const double crhs241 =             crhs127*crhs48;
const double crhs242 =             crhs12*crhs125*crhs14*crhs56;
const double crhs243 =             crhs12*crhs13*crhs15;
const double crhs244 =             crhs131*crhs34;
const double crhs245 =             crhs172*crhs36;
const double crhs246 =             pow(crhs48, 2);
const double crhs247 =             pow(crhs56, 2);
const double crhs248 =             6*crhs36*crhs40*crhs48*lambda;
const double crhs249 =             crhs34*crhs48*crhs62*mu;
const double crhs250 =             crhs48*crhs62;
const double crhs251 =             crhs38*crhs56*crhs62*mu;
const double crhs252 =             4*crhs13*crhs48*crhs62*crhs75;
            rhs[0]=N[0]*(crhs7 + crhs8 + crhs9) + N[0]*(crhs3 + 2*crhs4 + 2*crhs5 + 2*crhs6) - 1.0L/9.0L*crhs10*crhs100*crhs186*crhs27*h;
            rhs[1]=N[0]*(crhs189 + crhs190 + crhs191) + (1.0L/2.0L)*N[0]*(2*crhs197 + crhs198*crhs54 + 4*crhs199 - crhs200*crhs201 - 2*crhs202 - 4*crhs203 + crhs204*crhs34) - crhs187*crhs188 - crhs187*crhs205 - crhs192*crhs195 - crhs196*(-crhs126*crhs91 + crhs127*crhs59 + crhs195 + 2*crhs45 + 2*crhs46 + 2*crhs47 - 3*crhs51 - 3*crhs52 - 3*crhs53) + crhs233*(crhs114*crhs58 + crhs121*crhs60 + crhs123*crhs60 + crhs126*crhs88 + crhs13*crhs200*crhs224*crhs56 + crhs13*crhs224*crhs228*crhs48 - crhs132*crhs54 - 16*crhs14*crhs54*crhs83 - crhs143*crhs220 - crhs145*crhs220 - crhs150*crhs54*crhs83 - crhs154*crhs222 - crhs162*crhs215 - crhs163*crhs88 + crhs165*crhs221 - crhs187*crhs208 + 6*crhs189 + 6*crhs190 + 6*crhs191 + 6*crhs197 + 6*crhs199 - crhs201*crhs215 - 6*crhs202 - 6*crhs203 + crhs204*crhs44 + crhs207 + crhs209*crhs210 - crhs209*crhs216 + crhs209*crhs229 + crhs211*crhs212 - crhs211*crhs213 + crhs211*crhs219 + crhs214*crhs54 + crhs215*crhs36*crhs40*crhs67*lambda + crhs215*crhs49 - crhs215*crhs71 + crhs217*crhs43 + crhs218*crhs38 + crhs218*crhs68 + crhs219*crhs38*crhs54 + 4*crhs223*crhs3*mu + 4*crhs223*crhs38*crhs75 + crhs229*crhs43 + crhs230*crhs3*crhs48 + crhs230*crhs43 + crhs231*crhs232 - 6*crhs3*crhs34*crhs98 + crhs50*crhs58 - 20*crhs79 + 20*crhs84 - crhs89*crhs99);
            rhs[2]=N[0]*(crhs235 + crhs236 + crhs237) + N[0]*(-crhs228*crhs243 + 2*crhs239 + crhs240 + crhs241 - 2*crhs242 - crhs244 + crhs245*crhs38) - crhs188*crhs234 - crhs192*crhs238 - crhs196*(3*crhs193 + crhs198*crhs38 + crhs238 - crhs28 - crhs30 - crhs32 - 4*crhs4 - 4*crhs5 - 4*crhs6) - crhs205*crhs234 + crhs233*(crhs114*crhs49 + crhs116*crhs49 + crhs121*crhs206 + 4*crhs13*crhs250*crhs56*mu - crhs13*crhs34*crhs56*crhs64*crhs93*mu - crhs13*crhs48*crhs92*crhs93*mu - crhs130*crhs48 + crhs14*crhs222*crhs48 + 2*crhs14*crhs224*crhs3*crhs56 - crhs143*crhs249 - crhs145*crhs249 - crhs150*crhs251 - crhs152*crhs251 - crhs155*crhs221 + crhs207 - crhs208*crhs234 + crhs210*crhs247 + crhs212*crhs231 + crhs212*crhs246 - crhs213*crhs231 - crhs213*crhs246 + crhs214*crhs56 - crhs216*crhs247 + crhs217*crhs56*crhs68 + crhs217*crhs61 + crhs219*crhs231 + crhs219*crhs61 + 2*crhs223*crhs48*mu + crhs230*crhs246 + crhs230*crhs87 + crhs232*crhs247 + crhs232*crhs61 + 6*crhs235 + 6*crhs236 + 6*crhs237 + 6*crhs239 + 6*crhs240 + 6*crhs241 - 6*crhs242 - crhs243*crhs88 - 6*crhs244 + crhs245*crhs59 + crhs248*crhs34 + crhs248*crhs67 - 6*crhs250*crhs34*crhs97 + crhs252*crhs3 + crhs252*crhs34 - 6*crhs38*crhs56*crhs99 + crhs49*crhs50 + 20*crhs81 - 20*crhs86);

}

/*
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
    const bounded_matrix<double,dim,dimes> grad_U = prod(trans(DN), U);
    const double& r_gauss = inner_prod(data.N, data.r);
    
    const array_1d<double,dimes> accel_gauss = bdf0*U_gauss+bdf1*prod(trans(Un), N)+bdf2*prod(trans(Unn), N);
    
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
    const bounded_matrix<double,dim,dimes> grad_U = prod(trans(DN), U);
    const double& r_gauss = inner_prod(data.N, data.r);
    
    const array_1d<double,dimes> accel_gauss = bdf0*U_gauss+bdf1*prod(trans(Un), N)+bdf2*prod(trans(Unn), N);

    // Gauss point velocity subscale value computation
    //substitute_gausspt_subscale_2D

    const double U_gauss_norm = norm_2(U_gauss);
    const double U_s_gauss_norm = norm_2(U_s_gauss);

    return U_s_gauss_norm/U_gauss_norm;
}
*/

}
