// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Riccardo Rossi
//                   Ruben Zorrilla
//

// System includes


// External includes


// Project includes
#include "custom_conditions/fsi_line_load_condition_2d.h"

namespace Kratos
{
    //******************************* CONSTRUCTOR ****************************************
    //************************************************************************************

    FSILineLoadCondition2D::FSILineLoadCondition2D( IndexType NewId, GeometryType::Pointer pGeometry )
        : LineLoadCondition2D( NewId, pGeometry ) {}

    //************************************************************************************
    //************************************************************************************

    FSILineLoadCondition2D::FSILineLoadCondition2D( IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties )
        : LineLoadCondition2D( NewId, pGeometry, pProperties ) {}

    //********************************* CREATE *******************************************
    //************************************************************************************

    Condition::Pointer FSILineLoadCondition2D::Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties
        ) const {
        return boost::make_shared<FSILineLoadCondition2D>(NewId, pGeom, pProperties);
    }

    //************************************************************************************
    //************************************************************************************

    Condition::Pointer FSILineLoadCondition2D::Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties
        ) const {
        return boost::make_shared<FSILineLoadCondition2D>( NewId, GetGeometry().Create( ThisNodes ), pProperties );
    }

    //******************************* DESTRUCTOR *****************************************
    //************************************************************************************

    FSILineLoadCondition2D::~FSILineLoadCondition2D() {}

    //********************************* PUBLIC *******************************************
    //************************************************************************************

    void FSILineLoadCondition2D::CalculateMassMatrix(
        MatrixType& rMassMatrix,
        ProcessInfo& rCurrentProcessInfo) {

        // const double h = rCurrentProcessInfo[LENGHT]; //TO BE CHANGED
        // const double rho_fluid = rCurrentProcessInfo[DENSITY];  //TO BE CHANGED

        const double h = 2.5e3;            //TO BE CHANGED
        const double rho_fluid = 956.0;     //TO BE CHANGED

        // Get a reference to the condition geometry
        Geometry<Node<3>> &r_geometry = this->GetGeometry();
        const unsigned int n_nodes = r_geometry.size();
        const unsigned int domain_size = r_geometry.WorkingSpaceDimension();
        const unsigned int mat_size = n_nodes * domain_size;

        // Resize mass matrix
        if(rMassMatrix.size1() != mat_size || rMassMatrix.size2() != mat_size) {
            rMassMatrix.resize(mat_size,mat_size,false);
        }
        rMassMatrix.clear();

        // Get the condition shape functions
        Vector det_J;
        r_geometry.DeterminantOfJacobian(det_J, GeometryData::GI_GAUSS_2);
        const Matrix N_container = r_geometry.ShapeFunctionsValues(GeometryData::GI_GAUSS_2);
        const auto& integration_points = r_geometry.IntegrationPoints(GeometryData::GI_GAUSS_2);

        // Fill the mass matrix
        const unsigned int n_gauss = N_container.size1();
        for (unsigned int i_gauss = 0; i_gauss < n_gauss; ++i_gauss) {
            // Gauss pt. values
            const Vector aux_N = row(N_container, i_gauss);
            const double aux_weight = det_J(i_gauss)*integration_points[i_gauss].Weight();
            // Add current Gauss pt. contribution
            for(unsigned int i = 0; i < n_nodes; ++i) {
                for(unsigned int j = 0; j < n_nodes; ++j) {
                    for(unsigned int k = 0; k < domain_size; ++k) {
                        rMassMatrix(i*domain_size+k, j*domain_size+k) += aux_weight*aux_N(i)*aux_N(j);
                    }
                }
            }
        }

        // Multiphy by the fluid density and "attached layer" thickness
        rMassMatrix *= (h*rho_fluid);
    }

    void FSILineLoadCondition2D::CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
		VectorType& rRightHandSideVector,
		ProcessInfo& rCurrentProcessInfo) {

        // Get the condition constant parameters
        Geometry<Node<3>> &r_geometry = this->GetGeometry();
        const unsigned int n_nodes = r_geometry.size();
        const unsigned int domain_size = r_geometry.WorkingSpaceDimension();
        const unsigned int mat_size = n_nodes * domain_size;

        // Resize and initialize condition arrays
        if (rRightHandSideVector.size() != mat_size) {
            rRightHandSideVector.resize(mat_size, false);
        }

        if (rLeftHandSideMatrix.size1() != mat_size ||
            rLeftHandSideMatrix.size2() != mat_size) {
            rLeftHandSideMatrix.resize(mat_size, mat_size, false);
        }

        rLeftHandSideMatrix.clear();
        rRightHandSideVector.clear();

        // Call the base class to do the load integration
        LineLoadCondition2D::CalculateLocalSystem(
            rLeftHandSideMatrix,
            rRightHandSideVector,
            rCurrentProcessInfo);

        // Obtain the mass*acceleration at the old iteration and sum to the RHS,
        // so that the global effect of this condition is adding -M*(acc(it+1)-acc(it))
        Matrix mass_matrix;
        CalculateMassMatrix(mass_matrix,rCurrentProcessInfo);

        // Set the old acceleration vector from nodal data
        Vector acc_old_it(mat_size);
        for(unsigned int i_node = 0; i_node < n_nodes; ++i_node) {
            // Here assuming acceleration was stored in the non-historical
            // database at the end of the previous FSI non-linear iteration
            const auto& aux_acc_old_it = r_geometry[i_node].GetValue(ACCELERATION);
            for(unsigned int k = 0; k < domain_size; ++k) {
                acc_old_it(i_node*domain_size + k) = aux_acc_old_it(k);
            }
        }

        // Add the mass matrix times old acceleration product to the RHS
        // Note that the LHS contribution is added by the scheme
        noalias(rRightHandSideVector) += prod(mass_matrix, acc_old_it);
    }
} // Namespace Kratos
