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
   
        if(isLHSWithoutConstraintsComputed==false)
            for(auto & condition_i : mrReconstructionConditions)
            {
                condition_i->ComputeAndAddLHSContribution( mLHSWithoutConstraints );
                isLHSWithoutConstraintsComputed = true;                
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
    
        if(isRHSWithoutConstraintsComputed==false)
            for(auto & condition_i : mrReconstructionConditions)
            {
                condition_i->ComputeAndAddRHSContribution( mRHSWithoutConstraints );
                isRHSWithoutConstraintsComputed = true;
            }

        mRHS = mRHSWithoutConstraints;        

        for(auto & condition_i : mrReconstructionConstraints)
            condition_i->ComputeAndAddRHSContribution( mRHS );              

        std::cout << "> Time needed for computing RHS: " << timer.elapsed() << " s" << std::endl;            
    } 

    // --------------------------------------------------------------------------
    void RegularizeEquationSystem()
    {
        std::cout << "\n> Starting to regularize equation system ..." << std::endl;           
    
        bool is_small_value_on_main_diagonal = false;
        int number_of_small_values_on_main_diagonal = 0;

        for(int i_itr=0; i_itr<mLHS.size1(); i_itr++)
            if(std::abs(mLHS(i_itr,i_itr))<1e-10)
            {
                is_small_value_on_main_diagonal = true;
                number_of_small_values_on_main_diagonal++;

                mLHS(i_itr,i_itr) = 1e-3;
            }

        if(is_small_value_on_main_diagonal)
        {
            std::cout << "> WARNING, number of values on main diagonal < 1e-10: " << number_of_small_values_on_main_diagonal << " !!!!!!!!!!!!!!! " << std::endl;
            std::cout << "> All these values are manually set to 1e-3." << std::endl;
        }

        std::cout << "> Finished regularizing equation system." << std::endl;           
    }   

    // --------------------------------------------------------------------------
    void SolveEquationSystem()
    {
        std::cout << "\n> Start solving reconstruction equation..." << std::endl;           
        boost::timer timer;
        
		Vector control_point_update = ZeroVector(3*mNumberOfRelevantControlPoints);

		// Assign sparse LHS matrix to compressed LHS matrix required by linear solver
        CompressedMatrix compressed_lhs = mLHS;
                
        // KRATOS_WATCH(compressed_lhs.size1())
        // for(size_t i=0; i<compressed_lhs.size1(); i++)
        // {
        //         compressed_lhs(i,i) = 1;
        //         mRHS(i) = 1;
        // }

        mpLinearSolver->Solve(compressed_lhs, control_point_update, mRHS);

        // KRATOS_WATCH(control_point_update[0])

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
    bool isLHSWithoutConstraintsComputed = false;    
    Vector mRHS;
    Vector mRHSWithoutConstraints;    
    bool isRHSWithoutConstraintsComputed = false;

    /// Assignment operator.
    //      CADReconstructionSolver& operator=(CADReconstructionSolver const& rOther);

    /// Copy constructor.
    //      CADReconstructionSolver(CADReconstructionSolver const& rOther);

}; // Class CADReconstructionSolver
} // namespace Kratos.

#endif // CAD_RECONSTRUCTION_SOLVER_H
