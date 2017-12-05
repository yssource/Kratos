// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Riccardo Rossi
//                   Vicente Mataix
//

// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "custom_conditions/fsi_line_load_condition_2d.h"

namespace Kratos
{
    //******************************* CONSTRUCTOR ****************************************
    //************************************************************************************
    
    FSILineLoadCondition2D::FSILineLoadCondition2D( IndexType NewId, GeometryType::Pointer pGeometry )
        : BaseLoadCondition( NewId, pGeometry )
    {
        //DO NOT ADD DOFS HERE!!!
    }

    //************************************************************************************
    //************************************************************************************
    
    FSILineLoadCondition2D::FSILineLoadCondition2D( IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties )
        : BaseLoadCondition( NewId, pGeometry, pProperties )
    {
    }

    //********************************* CREATE *******************************************
    //************************************************************************************
    
    Condition::Pointer FSILineLoadCondition2D::Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties
        ) const
    {
        return boost::make_shared<FSILineLoadCondition2D>(NewId, pGeom, pProperties);
    }

    //************************************************************************************
    //************************************************************************************
    
    Condition::Pointer FSILineLoadCondition2D::Create( 
        IndexType NewId, 
        NodesArrayType const& ThisNodes,  
        PropertiesType::Pointer pProperties 
        ) const
    {
        return boost::make_shared<FSILineLoadCondition2D>( NewId, GetGeometry().Create( ThisNodes ), pProperties );
    }

    //******************************* DESTRUCTOR *****************************************
    //************************************************************************************
    
    FSILineLoadCondition2D::~FSILineLoadCondition2D()
    {
    }


    void FSILineLoadCondition2D::CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo)
    {
        const double h = rCurrentProcessInfo[LENGHT]; //TO BE CHANGED
        const double rho_fluid = rCurrentProcessInfo[DENSITY];  //TO BE CHANGED
        
        const double domain_size = GetGeometry().WorkingSpaceDimension();
        const double mat_size = domain_size * GetGeometry().size();
        if(rMassMatrix.size1() != mat_size && rMassMatrix.size2() != mat_size)
            rMassMatrix.resize(mat_size,mat_size,false);
        
        rMassMatrix.clear();
        
        Vector N;
        GetGeometry().ShapeFunctionValues(N);
        
        for(unsigned int i=0; i<GetGeometry().size(); ++i)
        {
            for(unsigned int k=0; k<domain_size; ++k)
                rMassMatrix(i*domain_size+k, i*domain_size+k) = N[i];
        }
    }
    
    void FSILineLoadCondition2D::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
				      VectorType& rRightHandSideVector,
				      ProcessInfo& rCurrentProcessInfo)
    {
        //call the base class to do the load integration
        LineLoadCondition2D::CalculateLocalSystem(rLeftHandSideMatrix,rRightHandSideVector,rCurrentProcessInfo); 
        
        //obtain the mass*acceleration at the old iteration and sum to the RHS, so that the global effect of this condition is adding -M*(acc(it+1)-acc(it))
        Matrix M;
        CalculateMassMatrix(M,rCurrentProcessInfo);
        
        Vector acc_old_it(mat_size);
        for(unsigned int i=0; i<GetGeometry().size(); ++i)
        {
            const auto& acc_old_it = GetGeometry()[i].GetValue(ACCELERATION); //here assuming acceleration was stored in the non-historical database at the end of the previous step
            for(unsigned int k=0; k<domain_size; ++k)
                acc_old_it[i*domain_size+k] = acc_old_it[k];
        }
        
        noalias(rRightHandSideVector) += prod(M, acc_old_it);
        
        
    }
} // Namespace Kratos


