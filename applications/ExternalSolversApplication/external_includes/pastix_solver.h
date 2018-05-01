//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

#if !defined(KRATOS_PASTIX_SOLVER )
#define  KRATOS_PASTIX_SOLVER

// External includes
extern "C" {
#include <pastix.h>
}

// Project includes
#include "includes/define.h"
#include "linear_solvers/direct_solver.h"
#include "includes/ublas_interface.h"
#include "includes/kratos_parameters.h"

namespace Kratos
{

template< class TSparseSpaceType, class TDenseSpaceType,
class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType> >
class PastixSolver : public DirectSolver< TSparseSpaceType,
    TDenseSpaceType, TReordererType>
{
public:
    /**
     * Counted pointer of PastixSolver
     */
    KRATOS_CLASS_POINTER_DEFINITION( PastixSolver );

    typedef LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType> BaseType;

    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    typedef typename TSparseSpaceType::VectorType VectorType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;
        
    // Pastix definitions
#ifdef VERSION_PASTIX_6
    #include <spm.h>
    typedef double PastixFloatType;
    typedef pastix_int_t PastixIntegerType;
    typedef pastix_data_t PastisDataType;
#else
    typedef pastix_float_t PastixFloatType;
    typedef pastix_int_t PastixIntegerType;
    typedef pastix_data_t PastisDataType;
#endif

    PastixSolver(Parameters settings)
    {
        Parameters default_settings( R"(
                                    {
                                    "solver_type" : "PastixSolver",
                                    "solution_method": "Direct",
                                    "tolerance":1e-6,
                                    "max_iteration":100,
                                    "gmres_krylov_space_dimension":100,
                                    "ilu_level_of_fill" : 1,
                                    "is_symmetric":false,
                                    "verbosity":0,
                                    "scaling": false,
                                    "block_size": 1,
                                    "use_block_matrices_if_possible" : true
                                }  )" );

        settings.ValidateAndAssignDefaults(default_settings);

        //validate if values are admissible
        std::set<std::string> available_solution_methods = {"Direct","Iterative"};
        if(available_solution_methods.find(settings["solution_method"].GetString()) == available_solution_methods.end())
        {
            KRATOS_ERROR << "trying to choose an inexisting solution method. Options are Direct, Iterative. Current choice is : " << settings["solution_method"].GetString() << std::endl;
        }

        if(settings["solution_method"].GetString() == "Iterative")
            mincomplete = 1;
        else
            mincomplete = 0;

        mTol = settings["tolerance"].GetDouble();
        mmax_it = settings["max_iteration"].GetInt();
        mlevel_of_fill = settings["ilu_level_of_fill"].GetInt();
        mverbosity=settings["verbosity"].GetInt();
        mndof = settings["block_size"].GetInt();


        if(settings["is_symmetric"].GetBool() == false)
                msymmetric = 0;
        else
                msymmetric = 1;


    }
        
    /**
     * Default constructor - uses ILU+GMRES
     * @param NewMaxTolerance tolerance that will be achieved by the iterative solver
     * @param NewMaxIterationsNumber this number represents both the number of iterations AND the size of the krylov space
     * @param level_of_fill of fill that will be used in the ILU
     * @param verbosity a number from 0 (no output) to 2 (maximal output)
     * @param is_symmetric set to True to solve assuming the matrix is symmetric
     */
    PastixSolver(double NewMaxTolerance,
                            int NewMaxIterationsNumber,
                            int level_of_fill,
                            int verbosity,
                            bool is_symmetric)
    {
        std::cout << "setting up pastix for iterative solve " << std::endl;
        mTol = NewMaxTolerance;
        mmax_it = NewMaxIterationsNumber;
        mlevel_of_fill = level_of_fill;
        mincomplete = 1;
        mverbosity=verbosity;

        if(is_symmetric == false)
            msymmetric = 0;
        else
            msymmetric = 1;

        mndof = 1;
    }

    /**
     * Direct Solver
     * @param verbosity a number from 0 (no output) to 2 (maximal output)
     * @param is_symmetric set to True to solve assuming the matrix is symmetric
     */
    PastixSolver(int verbosity, bool is_symmetric)
    {
        KRATOS_INFO("PastixSolver") << "Setting up pastix for direct solve " <<std::endl;
        mTol = -1;
        mmax_it = -1;
        mlevel_of_fill = -1;
        mincomplete = 0;
        mverbosity=verbosity;

        mndof = 1;
        if(is_symmetric == false)
            msymmetric = 0;
        else
            msymmetric = 1;
    }

    /**
     * Destructor
     */
    ~PastixSolver() override {};

    /**
     * Normal solve method.
     * Solves the linear system Ax=b and puts the result on SystemVector& rX.
     * rX is also th initial guess for iterative methods.
     * @param rA System matrix
     * @param rX Solution vector.
     * @param rB Right hand side vector.
     */
    bool Solve(SparseMatrixType& rA, VectorType& rX, VectorType& rB) override
    {
        int state = solvePASTIX(mverbosity, rA.size1(), rA.value_data().size(), rA.value_data().begin(), &(rA.index1_data()[0]), &(rA.index2_data()[0]), &rX[0], &rB[0]
        ,mmax_it,mTol,mincomplete,mlevel_of_fill,mndof,msymmetric);

        return state;
    }

    /**
     * Multi solve method for solving a set of linear systems with same coefficient matrix.
     * Solves the linear system Ax=b and puts the result on SystemVector& rX.
     * rX is also th initial guess for iterative methods.
     * @param rA System matrix
     * @param rX Solution vector.
     * @param rB Right hand side vector.
     */
    bool Solve(SparseMatrixType& rA, DenseMatrixType& rX, DenseMatrixType& rB) override
    {
        bool is_solved = true;
        return is_solved;
    }

    /**
     * Print information about this object.
     */
    void  PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "Pastix solver finished.";
    }

    /**
     * Print object's data.
     */
    void  PrintData(std::ostream& rOStream) const override
    {
    }

    /** Some solvers may require a minimum degree of knowledge of the structure of the matrix. To make an example
     * when solving a mixed u-p problem, it is important to identify the row associated to v and p.
     * another example is the automatic prescription of rotation null-space for smoothed-aggregation solvers
     * which require knowledge on the spatial position of the nodes associated to a given dof.
     * This function tells if the solver requires such data
     */
    bool AdditionalPhysicalDataIsNeeded() override
    {
        return true;
    }

    /** Some solvers may require a minimum degree of knowledge of the structure of the matrix. To make an example
     * when solving a mixed u-p problem, it is important to identify the row associated to v and p.
     * another example is the automatic prescription of rotation null-space for smoothed-aggregation solvers
     * which require knowledge on the spatial position of the nodes associated to a given dof.
     * This function is the place to eventually provide such data
     */
    void ProvideAdditionalData (
        SparseMatrixType& rA,
        VectorType& rX,
        VectorType& rB,
        typename ModelPart::DofsArrayType& rdof_set,
        ModelPart& r_model_part
    ) override
    {
        int old_ndof = -1;
        unsigned int old_node_id = rdof_set.begin()->Id();
        int ndof=0;

        for (ModelPart::DofsArrayType::iterator it = rdof_set.begin(); it!=rdof_set.end(); it++)
        {
            if(it->EquationId() < rA.size1() )
            {
                unsigned int id = it->Id();
                if(id != old_node_id)
                {
                    old_node_id = id;
                    if(old_ndof == -1) old_ndof = ndof;
                    else if(old_ndof != ndof) //if it is different than the block size is 1
                    {
                        old_ndof = -1;
                        break;
                    }

                    ndof=1;
                }
                else
                {
                    ndof++;
                }
            }
        }

        if(old_ndof == -1)
            mndof = 1;
        else
            mndof = ndof;

    //KRATOS_WATCH(mndof);


    }

private:

    double mTol;
    int mmax_it;
    int mincomplete;
    int mlevel_of_fill;
    int mverbosity;
    int mndof;
    int msymmetric;
//     double mDropTol;
//     double mFillTol;
//     double mFillFactor;

    /**
     * Assignment operator.
     */
    PastixSolver& operator=(const PastixSolver& Other);

    /**
     * Copy constructor.
     */
    PastixSolver(const PastixSolver& Other);

#define MPI_COMM_WORLD 0
    void CompRow_to_CompCol(int m, int n, int nnz,
                            double *a, size_t *colind, size_t *rowptr,
                            PastixFloatType **at, PastixIntegerType **rowind, PastixIntegerType **colptr)
    {
        int i, j, col, relpos;
        PastixIntegerType *marker = (PastixIntegerType *)malloc(n*sizeof(PastixIntegerType));
        for(i = 0; i<n; i++)
            marker[i] = 0;

        /* Allocate storage for another copy of the matrix. */
//        *at = (double *) (PastixFloatType *)malloc(nnz);
//        *rowind = (int *) (PastixIntegerType *)malloc(nnz);
//        *colptr = (int *) (PastixIntegerType *)malloc(n+1);

        /* Get counts of each column of A, and set up column pointers */
        for(i = 0; i < m; ++i)
            for(j = rowptr[i]; j < rowptr[i+1]; ++j) ++marker[colind[j]];
        (*colptr)[0] = 0;
        for(j = 0; j < n; ++j)
        {
            (*colptr)[j+1] = (*colptr)[j] + marker[j];
            marker[j] = (*colptr)[j];
        }

        /* Transfer the matrix into the compressed column storage. */
        for(i = 0; i < m; ++i)
        {
            for(j = rowptr[i]; j < rowptr[i+1]; ++j)
            {
                col = colind[j];
                relpos = marker[col];
                (*rowind)[relpos] = i;
                (*at)[relpos] = a[j];
                ++marker[col];
            }
        }

        free(marker);
    }

    int CompressCSCinFortranNumbering(int n, int ndof,
                    PastixFloatType **at, PastixIntegerType **rowind, PastixIntegerType **colptr
                    )
    {
        int nblocks = n/ndof;
        int nnz = (*colptr)[n]-1;
        int nnz_blocks = nnz/(ndof*ndof);
        int i,j,k;
        int v,r,c,ai,aj,nzb,offset,aii,l;

        PastixFloatType *bv      = (PastixFloatType *)malloc(sizeof(PastixFloatType)*(nnz));
        PastixIntegerType *br      = (PastixIntegerType *)malloc(sizeof(PastixIntegerType)*(nnz_blocks));
        PastixIntegerType *bc      = (PastixIntegerType *)malloc(sizeof(PastixIntegerType)*(nblocks+1));

        v = 0;
        r = 0;
        c = 1;

        for(j=0; j<nblocks; j++)
        {
            aj = j*ndof;
            nzb = ((*colptr)[(aj+1)] - (*colptr)[aj])/ndof;

            for(i=0; i<nzb;i++)
            {
                ai = (*colptr)[aj]  - 1 + i*ndof ;
                offset = nzb * ndof;
                for(k = 0; k<ndof; k++)
                {
                    aii = ai+k*offset;
                    for(l = 0; l<ndof; l++)
                    {
                        bv[v] = (*at)[aii+l];
                        v += 1;
                    }
                }
                br[r] = (*rowind)[ai]/ndof + 1;
                r += 1;
            }
            bc[j] = c;
            c += nzb;
        }
        bc[nblocks] = c;

        //free data that is not needed anymore
        free(*at);
        free(*rowind);
        free(*colptr);

        //do pointer swap
        *at = bv;
        *rowind = br;
        *colptr = bc;

        return nblocks;

    }

    int solvePASTIX(int verbosity,int mat_size, int nnz, double* AA, size_t* IA, size_t* JA, double *x, double* b, int m_gmres,
                    double tol, int incomplete, int ilu_level_of_fill, int ndof, int symmetric )
    {
        PastisDataType *pastix_data = NULL;                                                              /* Pointer to a storage structure needed by pastix           */
        PastixIntegerType ncol = mat_size;                                                               /* Size of the matrix                                        */
        PastixIntegerType *rows = (PastixIntegerType *)malloc(sizeof(PastixIntegerType)*(nnz));          /* Indexes of first element of each column in row and values */
        PastixIntegerType *colptr = (PastixIntegerType *)malloc(sizeof(PastixIntegerType)*(mat_size+1)); /* Row of each element of the matrix                         */
        PastixFloatType *values = (PastixFloatType *)malloc(sizeof(PastixFloatType)*(nnz));              /* Value of each element of the matrix                       */
        PastixFloatType *rhs = (PastixFloatType *)malloc(sizeof(PastixFloatType)*mat_size);              /* Right Hand Side                                           */
        PastixIntegerType iparm[IPARM_SIZE];                                                             /* Integer parameters for pastix                             */
        double dparm[DPARM_SIZE];                                                                        /* Floating parameters for pastix                            */
        PastixIntegerType *perm = (PastixIntegerType *)malloc((ncol+1)*sizeof(PastixIntegerType));       /* Permutation tabular                                       */
        PastixIntegerType *invp = (PastixIntegerType *)malloc((ncol+1)*sizeof(PastixIntegerType));       /* Reverse permutation tabular                               */
    #ifdef _OPENMP
        PastixIntegerType nbthread = omp_get_max_threads();                                              /* Number of thread wanted by user                           */
    #else
        PastixIntegerType nbthread = 1;                                                                  /* Number of thread wanted by user                           */
    #endif

        PastixIntegerType verbosemode = verbosity;                                                       /* Level of verbose mode (0, 1, 2)                           */
    #ifdef VERSION_PASTIX_6
        int ordering = PastixOrderScotch;                                                                /* Ordering to use                                           */
    #else
        int ordering = API_ORDER_SCOTCH;                                                                 /* Ordering to use                                           */
    #endif
        PastixIntegerType nbrhs = 1;
    //     int incomplete = 1;                                                                              /* Indicate if we want to use incomplete factorisation       */
        int level_of_fill = ilu_level_of_fill; //6;                                                      /* Level of fill for incomplete factorisation                */
        int amalgamation = 25;                                                                           /* Level of amalgamation for Kass                            */
    #ifndef VERSION_PASTIX_6
    //     int             ooc = 2000;                                                                      /* OOC limit (Mo/percent depending on compilation options)   */
    #endif

    /*    memset(colptr,0,(mat_size+1)*sizeof(PastixIntegerType));
        memset(rows,0,(nnz)*sizeof(PastixIntegerType));*/

        // Compute the transpose
        CompRow_to_CompCol(mat_size, mat_size, nnz,AA, JA, IA,&values, &rows, &colptr);

        int i;
    /*    FILE *fp_columns;
        fp_columns=fopen("columns.txt", "w");
        fprintf(fp_columns,"%d \n",(mat_size+1));
        for(i=0; i<mat_size+1; i++)
            fprintf(fp_columns,"%d \n",colptr[i]);
        fclose(fp_columns);

        FILE *fp_rows;
        fp_rows=fopen("rows.txt", "w");
        fprintf(fp_rows,"%d \n",(nnz));
        for(i=0; i<nnz; i++)
            fprintf(fp_rows,"%d \n",rows[i]);
        fclose(fp_rows);

        FILE *fp_values;
        fp_values=fopen("values.txt", "w");
        fprintf(fp_values,"%d \n",(nnz));
        for(i=0; i<nnz; i++)
            fprintf(fp_values,"%e \n",values[i]);
        fclose(fp_values);*/

        //exit(1);

//         for(i=0; i<nnz; i++)
//           rows[i] = rows[i]+1;
//
//         for(i=0; i<mat_size; i++)
//           colptr[i] = colptr[i]+1;

    #ifdef VERSION_PASTIX_6
//         /* Check the sparse matrix */ // NOTE: This will require transform the current matrix to Pastix sparse matrix
//         pastix_spm_t *AA2;
//         spmPrintInfo( AA, stdout );
//
//         AA2 = spmCheckAndCorrect( AA );
//         if ( AA2 != AA ) {
//             spmExit( AA );
//             free( AA );
//             AA = AA2;
//         }
//         /**
//          * Generate a Fake values array if needed for the numerical part
//          */
//         if ( AA->flttype == PastixPattern ) {
//             spmGenFakeValues( AA );
//         }
    #else
        /**
         * Matrix needs :
         *    - to be in fortran numbering
         *    - to have only the lower triangular part in symmetric case
         *    - to have a graph with a symmetric structure in unsymmetric case
         */
        PastixIntegerType mat_type;
        if(symmetric == 0)
            mat_type = API_SYM_NO;
        else
            mat_type = API_SYM_YES;
        iparm[IPARM_MATRIX_VERIFICATION] = API_YES;
        if(NO_ERR != pastix_checkMatrix(0, verbosemode,
                                        mat_type,
                                        API_YES,
                                        ncol, &colptr, &rows, &values, NULL, 1)) {
            return 1;
        }
    #endif

        // Copy to block format if needed
        if(ndof > 1)
            CompressCSCinFortranNumbering(mat_size,ndof,&values, &rows, &colptr);

        //copy b to the solution. It will be overwritten
        for(i = 0; i < mat_size; i++)
            x[i] = b[i];

        /*******************************************/
        /* Initialize parameters to default values */
        /*******************************************/
    #ifdef VERSION_PASTIX_6
        iparm[IPARM_MODIFY_PARAMETER] = 0;
    #else
        iparm[IPARM_MODIFY_PARAMETER] = API_NO;
    #endif
        pastix(&pastix_data, MPI_COMM_WORLD,
            ncol, colptr, rows, values,
            perm, invp, x, 1, iparm, dparm);

        /*******************************************/
        /*       Customize some parameters         */
        /*******************************************/
        iparm[IPARM_THREAD_NBR] = nbthread;
        iparm[IPARM_VERBOSE] = verbosemode;
        iparm[IPARM_ORDERING] = ordering;
        iparm[IPARM_DOF_NBR] = ndof;
    #ifdef VERSION_PASTIX_6
        switch(symmetric)
        {
        case 1:
            iparm[IPARM_FACTORIZATION] = PastixFactLDLT;
            break;
        default:
            iparm[IPARM_FACTORIZATION] = PastixFactLU;
        }

        iparm[IPARM_INCOMPLETE] = incomplete;

    //    iparm[IPARM_FREE_CSCUSER] = 0;

        if(incomplete == 1) {
            iparm[IPARM_REFINEMENT] = PastixRefineGMRES;
            iparm[IPARM_GMRES_IM] = m_gmres;
            dparm[DPARM_EPSILON_REFINEMENT] = tol;
            iparm[IPARM_LEVEL_OF_FILL] = level_of_fill;
            iparm[IPARM_AMALGAMATION_LVLBLAS] = amalgamation;
            iparm[IPARM_AMALGAMATION_LVLCBLK] = amalgamation;
        }

        iparm[IPARM_START_TASK] = PastixTaskOrdering;
        iparm[IPARM_END_TASK] = PastixTaskClean;
    #else
        iparm[IPARM_SYM] = mat_type;
        switch(mat_type)
        {
        case API_SYM_YES:
            iparm[IPARM_FACTORIZATION] = API_FACT_LDLT;
            break;
        case API_SYM_HER:
            iparm[IPARM_FACTORIZATION] = API_FACT_LDLH;
            break;
        default:
            iparm[IPARM_FACTORIZATION] = API_FACT_LU;
        }
        iparm[IPARM_MATRIX_VERIFICATION] = API_NO;
        if(incomplete == 1)
            iparm[IPARM_INCOMPLETE]          = API_YES;
        else if(incomplete == 0)
            iparm[IPARM_INCOMPLETE]          = API_NO;
        else
            KRATOS_WARNING("PastixSolver") << "Incomplete flag should either be 0 (direct solve) or 1 for ILU solve" << std::endl;

    //     iparm[IPARM_OOC_LIMIT]           = ooc;

    //    iparm[IPARM_FREE_CSCUSER] = API_CSC_PRESERVE;

        if(iparm[IPARM_INCOMPLETE]  == API_YES) {
            iparm[IPARM_REFINEMENT] = API_RAF_GMRES;
            iparm[IPARM_GMRES_IM] = m_gmres;
            dparm[DPARM_EPSILON_REFINEMENT] = tol;
            iparm[IPARM_LEVEL_OF_FILL]       = level_of_fill;
            iparm[IPARM_AMALGAMATION_LEVEL]  = amalgamation;
        }

        iparm[IPARM_RHS_MAKING]          = API_RHS_B;
        iparm[IPARM_START_TASK]          = API_TASK_ORDERING;
        iparm[IPARM_END_TASK]            = API_TASK_CLEAN;
    #endif
        /*******************************************/
        /*           Call pastix                   */
        /*******************************************/
        //note that we pass "x" instead of the rhs
        pastix(&pastix_data, 0,
            ncol/ndof, colptr, rows, values,
            perm, invp, x, nbrhs, iparm, dparm);

//         PRINT_RHS("SOL", rhs, ncol, mpid, iparm[IPARM_VERBOSE]);
//         CHECK_SOL(rhs, rhssaved, ncol, mpid);

        free(colptr);
        free(rows);
        free(values);
        free(perm);
        free(invp);
        free(rhs);
//         free(rhssaved);
//         free(type);
//         free(rhstype);

        return EXIT_SUCCESS;
    }
#undef MPI_COMM_WORLD
}; // Class PastixSolver


/**
 * input stream function
 */
template<class TSparseSpaceType, class TDenseSpaceType,class TReordererType>
inline std::istream& operator >> (std::istream& rIStream, PastixSolver< TSparseSpaceType,
                                  TDenseSpaceType, TReordererType>& rThis)
{
    return rIStream;
}

/**
 * output stream function
 */
template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const PastixSolver<TSparseSpaceType,
                                  TDenseSpaceType, TReordererType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

//#undef MPI_COMM_WORLD

}  // namespace Kratos.


#endif // KRATOS_PASTIX_SOLVER  defined
