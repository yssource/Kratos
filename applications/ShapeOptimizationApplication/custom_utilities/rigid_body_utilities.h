// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Willer Matthias, https://github.com/matthiaswiller
//
// ==============================================================================

#ifndef RIGID_BODY_UTILITIES_H
#define RIGID_BODY_UTILITIES_H

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------
#include <iostream>
#include <string>
#include <algorithm>

// ------------------------------------------------------------------------------
// External includes
// ------------------------------------------------------------------------------
#include <boost/python.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "includes/define.h"
#include "processes/process.h"
#include "includes/node.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "includes/kratos_flags.h"
#include "spaces/ublas_space.h"
#include "shape_optimization_application.h"
#include "utilities/svd_utils.h"

// ==============================================================================

namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
*/

class RigidBodyUtilities
{
public:
    ///@name Type Definitions
    ///@{

    // Type definitions for better reading later
    typedef array_1d<double,3> array_3d;
    typedef Node < 3 > NodeType;
    typedef Node < 3 > ::Pointer NodeTypePointer;
    typedef std::vector<NodeType::Pointer> NodeVector;
    typedef std::vector<NodeType::Pointer>::iterator NodeIterator;
    typedef std::vector<double>::iterator DoubleVectorIterator;
    typedef ModelPart::ConditionsContainerType ConditionsArrayType;

    // Type definitions for linear algebra including sparse systems
    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef SparseSpaceType::MatrixType SparseMatrixType;
    typedef SparseSpaceType::VectorType VectorType;

    /// Pointer definition of MapperVertexMorphing
    // KRATOS_CLASS_POINTER_DEFINITION(MapperVertexMorphing);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    RigidBodyUtilities( ModelPart& designSurface, Parameters optimizationSettings )
        : mrDesignSurface( designSurface ),
          mNumberOfDesignVariables(designSurface.Nodes().size())
    {
        CreateListOfRigidNodes();
    }

    /// Destructor.
    virtual ~RigidBodyUtilities()
    {
    }

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    // ==============================================================================
    void CreateListOfRigidNodes()
    {
        // for (auto & node_i : mrDesignSurface.Nodes())
        //     if(node_i.X() < -5.0)
        //         mNumberOfRigidNodes++;
        // mListOfRigidNodes.resize(mNumberOfRigidNodes)

        for (ModelPart::NodesContainerType::iterator node_it = mrDesignSurface.NodesBegin(); node_it != mrDesignSurface.NodesEnd(); ++node_it)
        {   
            NodeTypePointer pnode = *(node_it.base());

            if(pnode->X() < -5.0)
            {
                mNumberOfRigidNodes++;
                mListOfRigidNodes.push_back(pnode);
            }
        }
    }

    // --------------------------------------------------------------------------
    void FindRotationAndTranslationWithSVD()
    {
        // centroid A (member variable)
        Vector centroid_undeformed = ZeroVector(3);
        // centroid B (member variable)
        Vector centroid_deformed = ZeroVector(3);

        for(int node_index = 0 ; node_index<mListOfRigidNodes.size() ; node_index++)
        {
            // Get node information
            ModelPart::NodeType& node_i = *mListOfRigidNodes[node_index];

            array_3d& coord = node_i.Coordinates();
            array_3d& def = node_i.FastGetSolutionStepValue( SHAPE_UPDATE );

            centroid_undeformed(0) += coord[0];
            centroid_undeformed(1) += coord[1];
            centroid_undeformed(2) += coord[2];

            centroid_deformed(0) += coord[0] + def[0];
            centroid_deformed(1) += coord[1] + def[1];
            centroid_deformed(2) += coord[2] + def[2];
        }
        centroid_undeformed(0) = centroid_undeformed(0)/mNumberOfRigidNodes;
        centroid_undeformed(1) = centroid_undeformed(1)/mNumberOfRigidNodes;
        centroid_undeformed(2) = centroid_undeformed(2)/mNumberOfRigidNodes;
        centroid_deformed(0) = centroid_deformed(0)/mNumberOfRigidNodes;
        centroid_deformed(1) = centroid_deformed(1)/mNumberOfRigidNodes;
        centroid_deformed(2) = centroid_deformed(2)/mNumberOfRigidNodes;

        // H matrix (3x3)
        Matrix H = ZeroMatrix(3,3);
        for (int i=0;i<3;i++)
        {
            for(int j=0;j<3;j++)
            {
                for(int node_index = 0 ; node_index<mListOfRigidNodes.size() ; node_index++)
                {
                    // Get node information
                    ModelPart::NodeType& node_i = *mListOfRigidNodes[node_index];

                    array_3d& coord = node_i.Coordinates();
                    array_3d& def = node_i.FastGetSolutionStepValue( SHAPE_UPDATE );

                    H(i,j) += (coord[i] - centroid_undeformed(i))*(coord[j] + def[j] - centroid_deformed(j));
                }
            }
        }

        // U,V,S = svd(H)
        Matrix U = ZeroMatrix(3,3);
        Matrix V = ZeroMatrix(3,3);
        Matrix S = ZeroMatrix(3,3);

        Matrix R = prod(trans(H),H);

        SVDUtils<double>::JacobiSingularValueDecomposition(H, U, S, V);

        KRATOS_WATCH(U);

        // // R
        // Matrix R = prod(U,trans(V))

        double detR = MathUtils<double>::Det(R);

        // if (detR < 0) // special reflection case
        // {
        //     R(0,0) *= -1;
        //     R(0,1) *= -1;
        //     R(0,2) *= -1;
        // }

        // // t
        // Vector t = ZeroVector(3);
        // t = - R*centroid_undeformed + centroid_undeformed;

        // print(R);
        // print(t);
        KRATOS_WATCH(H);
    }

    // --------------------------------------------------------------------------
    void CorrectDesignUpdateWithRigidBodyConstraints()
    {
        FindRotationAndTranslationWithSVD();
    }

    // ==============================================================================

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "RigidBodyUtilities";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "RigidBodyUtilities";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
    }


    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{


    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    // // ==============================================================================
    // // Initialized by class constructor
    // // ==============================================================================
    ModelPart& mrDesignSurface;    
    const unsigned int mNumberOfDesignVariables;
    unsigned int mNumberOfRigidNodes = 0;
    NodeVector mListOfRigidNodes;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
//      RigidBodyUtilities& operator=(RigidBodyUtilities const& rOther);

    /// Copy constructor.
//      RigidBodyUtilities(RigidBodyUtilities const& rOther);


    ///@}

}; // Class RigidBodyUtilities

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // RIGID_BODY_UTILITIES_H
