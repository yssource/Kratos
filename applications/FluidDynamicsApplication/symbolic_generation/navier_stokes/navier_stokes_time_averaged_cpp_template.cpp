//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

#include "custom_elements/time_averaged_navier_stokes.h"

namespace Kratos {

template<>
void TimeAveragedNavierStokes<3>::ComputeGaussPointLHSContribution(
    BoundedMatrix<double,16,16>& lhs, 
    const ElementDataStruct& data){

    const int dim = 3;
    const int nnodes = 4;

    const double rho = inner_prod(data.N, data.rho);        // Density
    const double mu = inner_prod(data.N, data.mu);          // Dynamic viscosity
    const double h = data.h;                                // Characteristic element size
    const double c = data.c;                                // Wave velocity

    const double dt = data.dt;
    const double bdf0 = data.bdf0;
    const double dyn_tau = data.dyn_tau;

    const unsigned int  step = data.step;

    const BoundedMatrix<double,nnodes,dim>& y = data.y;
    const BoundedMatrix<double,nnodes,dim>& vmesh = data.vmesh;
    const BoundedMatrix<double,nnodes,dim>& vconv = y - vmesh;

    // Get constitutive matrix
    const Matrix& C = data.C;

    // Get shape function values
    const array_1d<double,nnodes>& N = data.N;
    const BoundedMatrix<double,nnodes,dim>& DN = data.DN_DX;

    // Stabilization parameters
    const double stab_c1 = 4.0;
    const double stab_c2 = 2.0;

    //substitute_lhs_3D
}


template<>
void TimeAveragedNavierStokes<2>::ComputeGaussPointLHSContribution(
    BoundedMatrix<double,9,9>& lhs,
    const ElementDataStruct& data){

    const int dim = 2;
    const int nnodes = 3;

    const double rho = inner_prod(data.N, data.rho);        // Density
    const double mu = inner_prod(data.N, data.mu);          // Dynamic viscosity
    const double h = data.h;                                // Characteristic element size
    const double c = data.c;                                // Wave velocity

    const double dt = data.dt;
    const double bdf0 = data.bdf0;
    const double dyn_tau = data.dyn_tau;

    const unsigned int  step = data.step;

    const BoundedMatrix<double,nnodes,dim>& y = data.y;
    const BoundedMatrix<double,nnodes,dim>& vmesh = data.vmesh;
    const BoundedMatrix<double,nnodes,dim>& vconv = y - vmesh;

    // Get constitutive matrix
    const Matrix& C = data.C;

    // Get shape function values
    const array_1d<double,nnodes>& N = data.N;
    const BoundedMatrix<double,nnodes,dim>& DN = data.DN_DX;

    // Stabilization parameters
    const double stab_c1 = 4.0;
    const double stab_c2 = 2.0;

    //substitute_lhs_2D
}


template<>
void TimeAveragedNavierStokes<3>::ComputeGaussPointRHSContribution(
    array_1d<double,16>& rhs, 
    const ElementDataStruct& data){

    const int dim = 3;
    const int nnodes = 4;
    const int strain_size = 6;

    const double rho = inner_prod(data.N, data.rho);        // Density
    const double mu = inner_prod(data.N, data.mu);          // Dynamic viscosity
    const double h = data.h;                                // Characteristic element size
    const double c = data.c;                                // Wave velocity

    const double dt = data.dt;
    const double bdf0 = data.bdf0;
    const double bdf1 = data.bdf1;
    const double bdf2 = data.bdf2;
    const double dyn_tau = data.dyn_tau;
    const unsigned int  step = data.step;

    const BoundedMatrix<double,nnodes,dim>& y = data.y;
    const BoundedMatrix<double,nnodes,dim>& yn = data.yn;
    const BoundedMatrix<double,nnodes,dim>& ynn = data.ynn;
    const BoundedMatrix<double,nnodes,dim>& ynnn = data.ynnn;
    const BoundedMatrix<double,nnodes,dim>& vmesh = data.vmesh;
    const BoundedMatrix<double,nnodes,dim>& vconv = y - vmesh;
    const BoundedMatrix<double,nnodes,dim>& f = data.f;
    const array_1d<double,nnodes>& x = data.x;
    const array_1d<double,nnodes>& xn = data.xn;
    const array_1d<double,nnodes>& xnn = data.xnn;
    const array_1d<double,nnodes>& xnnn = data.xnnn;
    const array_1d<double,strain_size>& stress = data.stress;

    // Get shape function values
    const array_1d<double,nnodes>& N = data.N;
    const BoundedMatrix<double,nnodes,dim>& DN = data.DN_DX;

    // Auxiliary variables used in the calculation of the RHS
    const array_1d<double,dim> f_gauss = prod(trans(f), N);
    const array_1d<double,dim> grad_p = prod(trans(DN), x);

    // Stabilization parameters
    const double stab_c1 = 4.0;
    const double stab_c2 = 2.0;

    //substitute_rhs_3D
}


template<>
void TimeAveragedNavierStokes<2>::ComputeGaussPointRHSContribution(
    array_1d<double,9>& rhs, 
    const ElementDataStruct& data){

    const int nnodes = 3;
    const int dim = 2;
    const int strain_size = 3;

    const double rho = inner_prod(data.N, data.rho);        // Density
    const double mu = inner_prod(data.N, data.mu);          // Dynamic viscosity
    const double h = data.h;                                // Characteristic element size
    const double c = data.c;                                // Wave velocity

    const double dt = data.dt;
    const double bdf0 = data.bdf0;
    const double bdf1 = data.bdf1;
    const double bdf2 = data.bdf2;
    const double dyn_tau = data.dyn_tau;
    const unsigned int  step = data.step;

    const BoundedMatrix<double,nnodes,dim>& y = data.y;
    const BoundedMatrix<double,nnodes,dim>& yn = data.yn;
    const BoundedMatrix<double,nnodes,dim>& ynn = data.ynn;
    const BoundedMatrix<double,nnodes,dim>& ynnn = data.ynnn;
    const BoundedMatrix<double,nnodes,dim>& vmesh = data.vmesh;
    const BoundedMatrix<double,nnodes,dim>& vconv = y - vmesh;
    const BoundedMatrix<double,nnodes,dim>& f = data.f;
    const array_1d<double,nnodes>& x = data.x;
    const array_1d<double,nnodes>& xn = data.xn;
    const array_1d<double,nnodes>& xnn = data.xnn;
    const array_1d<double,nnodes>& xnnn = data.xnnn;
    const array_1d<double,strain_size>& stress = data.stress;

    // Get shape function values
    const array_1d<double,nnodes>& N = data.N;
    const BoundedMatrix<double,nnodes,dim>& DN = data.DN_DX;

    // Auxiliary variables used in the calculation of the RHS
    const array_1d<double,dim> f_gauss = prod(trans(f), N);
    const array_1d<double,dim> grad_p = prod(trans(DN), x);

    // Stabilization parameters
    const double stab_c1 = 4.0;
    const double stab_c2 = 2.0;

    //substitute_rhs_2D
}

}
