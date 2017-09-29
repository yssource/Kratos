// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
//
// ==============================================================================

#ifndef CAD_RECONSTRUCTION_SOLVER_H
#define CAD_RECONSTRUCTION_SOLVER_H

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------

// ------------------------------------------------------------------------------
// External includes
// ------------------------------------------------------------------------------

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "includes/define.h"

// ==============================================================================

namespace Kratos
{
class CADReconstructionSolver
{
public:
    ///@name Type Definitions
    ///@{

    typedef UblasSpace<double, SparseMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, CompressedMatrix, Vector> CompressedSpaceType;
    typedef UblasSpace<double, Matrix, Vector> DenseSpaceType;
    typedef LinearSolver<CompressedSpaceType, DenseSpaceType > CompressedLinearSolverType;

    /// Pointer definition of CADReconstructionSolver
    KRATOS_CLASS_POINTER_DEFINITION(CADReconstructionSolver);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    CADReconstructionSolver( ReconstructionDataBase& reconstruction_data_base, 
                             ReconstructionConditionContainer& condition_container, 
                             CompressedLinearSolverType::Pointer linear_solver )
    : mrReconstructionDataBase( reconstruction_data_base ),
      mrReconstructionConditions( condition_container.GetReconstructionConditions() ),
      mrReconstructionConstraints( condition_container.GetConstraints() ),
      mpLinearSolver( linear_solver )     
    {      
    }

    /// Destructor.
    virtual ~CADReconstructionSolver()
    {
    }

    // --------------------------------------------------------------------------
    void InitializeEquationSystem()
    {
        IdentifyControlPointsRelevantForReconstruction();
        AssignEquationIdToControlPointsRelevantForReconstruction();
        InitializeSystemLHSAndRHS();        
    }

    // --------------------------------------------------------------------------
    void InitializeSystemLHSAndRHS()
    {
		mLHS.resize(3*mNumberOfRelevantControlPoints,3*mNumberOfRelevantControlPoints);
        mRHS.resize(3*mNumberOfRelevantControlPoints);
		mLHS.clear();
        mRHS.clear();
        
		mLHSWithoutConstraints.resize(3*mNumberOfRelevantControlPoints,3*mNumberOfRelevantControlPoints);
        mRHSWithoutConstraints.resize(3*mNumberOfRelevantControlPoints);
		mLHSWithoutConstraints.clear();
        mRHSWithoutConstraints.clear();                
    }

    // --------------------------------------------------------------------------
    void ComputeLHS()
    {
        std::cout << "\n> Start computing LHS..." << std::endl;           
        boost::timer timer;  
    
        if(isLHSWithoutConstraintsAlreadyComputed)
            for(auto & condition_i : mrReconstructionConditions)
            {
                condition_i->ComputeAndAddLHSContribution( mLHSWithoutConstraints );
                isRHSWithoutConstraintsAlreadyComputed = true;                
            }

        mLHS = mLHSWithoutConstraints;

        for(auto & condition_i : mrReconstructionConstraints)
            condition_i->ComputeAndAddLHSContribution( mLHS );            
        
        std::cout << "> Time needed for computing LHS: " << timer.elapsed() << " s" << std::endl;            
    }

    // --------------------------------------------------------------------------
    void ComputeRHS()
    {
        std::cout << "\n> Start computing RHS..." << std::endl;           
        boost::timer timer;  
    
        if(isRHSWithoutConstraintsAlreadyComputed)
            for(auto & condition_i : mrReconstructionConditions)
            {
                condition_i->ComputeAndAddRHSContribution( mRHSWithoutConstraints );
                isRHSWithoutConstraintsAlreadyComputed = true;
            }

        mRHS = mRHSWithoutConstraints;        

        for(auto & condition_i : mrReconstructionConstraints)
            condition_i->ComputeAndAddRHSContribution( mRHS );              

        std::cout << "> Time needed for computing RHS: " << timer.elapsed() << " s" << std::endl;            
    } 

    // --------------------------------------------------------------------------
    void SolveEquationSystem()
    {
        std::cout << "\n> Start solving reconstruction equation..." << std::endl;           
        boost::timer timer;
        
		Vector control_point_update = ZeroVector(3*mNumberOfRelevantControlPoints);

		// Assign sparse LHS matrix to compressed LHS matrix required by linear solver
        CompressedMatrix compressed_lhs = mLHS;
                
        mpLinearSolver->Solve(compressed_lhs, control_point_update, mRHS);

        mrReconstructionDataBase.UpdateControlPointDisplacements( control_point_update );

        std::cout << "> Time needed for solving reconstruction equation: " << timer.elapsed() << " s" << std::endl;                    
    }

    // --------------------------------------------------------------------------
    void IdentifyControlPointsRelevantForReconstruction()
    {
        std::cout << "\n> Start identifying control points relevant for reconstruction..." << std::endl;           
        boost::timer timer;      

        for(auto & condition_i : mrReconstructionConditions)
            condition_i->FlagControlPointsRelevantForReconstruction();

        for(auto & condition_i : mrReconstructionConstraints)
            condition_i->FlagControlPointsRelevantForReconstruction();              

        std::cout << "> Time needed for identifying control points relevant for reconstruction: " << timer.elapsed() << " s" << std::endl;    
    }      

    // --------------------------------------------------------------------------
    void AssignEquationIdToControlPointsRelevantForReconstruction()
    {
        std::cout << "\n> Start assigning equation ID to control points relevant for reconstruction..." << std::endl;           
        boost::timer timer;      

        int equation_id = 1;
        std::vector<Patch>& patch_vector = mrReconstructionDataBase.GetPatchVector();

        for(auto & patch_i : patch_vector) 
        {
            std::vector<ControlPoint>& control_points = patch_i.GetSurfaceControlPoints();
            for(auto & control_point_i : control_points)
            {
                if(control_point_i.IsRelevantForReconstruction())
                {
                            control_point_i.SetEquationId(equation_id);
                            ++mNumberOfRelevantControlPoints;
                            ++equation_id;
                }
                ++mNumberOfControlPoints;
            }        
        }      

        std::cout << "> Number of control points in total = " << mNumberOfControlPoints << "." << std::endl;
        std::cout << "> Number of control points relevant for mapping = " << mNumberOfRelevantControlPoints << "." << std::endl;  
        std::cout << "> Time needed for assigning equation ID to control points relevant for reconstruction: " << timer.elapsed() << " s" << std::endl;      
    }  

    // ==============================================================================

    /// Turn back information as a string.
    virtual std::string Info() const
    {
		return "CADReconstructionSolver";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream &rOStream) const
    {
		rOStream << "CADReconstructionSolver";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream &rOStream) const
    {
    }

private:

    // ==============================================================================
    // Initialized by class constructor
    // ==============================================================================
    ReconstructionDataBase& mrReconstructionDataBase;    
    std::vector<ReconstructionCondition::Pointer>& mrReconstructionConditions;
    std::vector<ReconstructionCondition::Pointer>& mrReconstructionConstraints;
    CompressedLinearSolverType::Pointer mpLinearSolver;

    // ==============================================================================
    // Member variables
    // ==============================================================================
    int mNumberOfControlPoints = 0;
    int mNumberOfRelevantControlPoints = 0;
    SparseMatrix mLHS;
    SparseMatrix mLHSWithoutConstraints;
    bool isLHSWithoutConstraintsAlreadyComputed = false;    
    Vector mRHS;
    Vector mRHSWithoutConstraints;    
    bool isRHSWithoutConstraintsAlreadyComputed = false;

    /// Assignment operator.
    //      CADReconstructionSolver& operator=(CADReconstructionSolver const& rOther);

    /// Copy constructor.
    //      CADReconstructionSolver(CADReconstructionSolver const& rOther);

}; // Class CADReconstructionSolver
} // namespace Kratos.

#endif // CAD_RECONSTRUCTION_SOLVER_H
