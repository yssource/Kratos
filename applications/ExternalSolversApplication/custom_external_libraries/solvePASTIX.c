/*-----------------------------------------------------------------*
 * main test driver for the ARMS2 preconditioner for
 * Matrices in the COO/Harwell Boeing format
 *-----------------------------------------------------------------*
 * Yousef Saad - Aug. 2005.                                        *
 *                                                                 *
 * Report bugs / send comments to: saad@cs.umn.edu                 *
 *-----------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "pastix.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef VERSION_PASTIX_6
typedef pastix_float_t PastixFloatType;
typedef pastix_int_t PastixIntegerType;
typedef pastix_data_t PastisDataType;
#else
typedef pastix_float_t PastixFloatType;
typedef pastix_int_t PastixIntegerType;
typedef pastix_data_t PastisDataType;
#endif

#define MPI_COMM_WORLD 0

void CompRow_to_CompCol(int m, int n, int nnz,
                        double *a, size_t *colind, size_t *rowptr,
                        PastixFloatType **at, PastixIntegerType **rowind, PastixIntegerType **colptr)
{
    register int i, j, col, relpos;
    PastixIntegerType *marker = (PastixIntegerType *)malloc(n*sizeof(PastixIntegerType));
    for(i = 0; i<n; i++)
        marker[i] = 0;

    /* Allocate storage for another copy of the matrix. */
//    *at = (double *) (PastixFloatType *)malloc(nnz);
//    *rowind = (int *) (PastixIntegerType *)malloc(nnz);
//    *colptr = (int *) (PastixIntegerType *)malloc(n+1);

//printf("111\n");

    /* Get counts of each column of A, and set up column pointers */
    for(i = 0; i < m; ++i)
        for(j = rowptr[i]; j < rowptr[i+1]; ++j) ++marker[colind[j]];
    (*colptr)[0] = 0;
    for(j = 0; j < n; ++j)
    {
        (*colptr)[j+1] = (*colptr)[j] + marker[j];
        marker[j] = (*colptr)[j];
    }
//printf("111\n");

    /* Transfer the matrix into the compressed column storage. */
    for(i = 0; i < m; ++i)
    {
        for(j = rowptr[i]; j < rowptr[i+1]; ++j)
        {
            col = colind[j];
            relpos = marker[col];
//printf("col %d \n",col);
//printf("relpos %d \n",relpos);
//printf("j %d \n",j);
//printf("i %d \n",i);
            (*rowind)[relpos] = i;
//printf("(*rowind)[relpos]  %d \n",(*rowind)[relpos] );
            (*at)[relpos] = a[j];
            ++marker[col];
        }
    }
//printf("111\n");

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


//int solvePASTIX(int echo_level,int mat_size, int nnz, double* AA, size_t* IA, size_t* JA, double *x, double* b)
int solvePASTIX(int verbosity,int mat_size, int nnz, double* AA, size_t* IA, size_t* JA, double *x, double* b, int m_gmres, 
                double tol, int incomplete, int ilu_level_of_fill, int ndof, int symmetric )
{
    PastisDataType  *pastix_data = NULL; /* Pointer to a storage structure needed by pastix           */
    PastixIntegerType    ncol = mat_size;               /* Size of the matrix                                        */
    PastixIntegerType   *rows      = (PastixIntegerType *)malloc(sizeof(PastixIntegerType)*(nnz));  /* Indexes of first element of each column in row and values */
    PastixIntegerType   *colptr        = (PastixIntegerType *)malloc(sizeof(PastixIntegerType)*(mat_size+1));   /* Row of each element of the matrix                         */
    PastixFloatType *values      = (PastixFloatType *)malloc(sizeof(PastixFloatType)*(nnz));   /* Value of each element of the matrix                       */
    PastixFloatType *rhs         = (PastixFloatType *)malloc(sizeof(PastixFloatType)*mat_size);  /* right hand side                                           */
    PastixIntegerType    iparm[IPARM_SIZE];  /* integer parameters for pastix                             */
    double          dparm[DPARM_SIZE];  /* floating parameters for pastix                            */
    PastixIntegerType   *perm        = (PastixIntegerType *)malloc((ncol+1)*sizeof(PastixIntegerType)); /* Permutation tabular                                       */
    PastixIntegerType   *invp        = (PastixIntegerType *)malloc((ncol+1)*sizeof(PastixIntegerType)); /* Reverse permutation tabular                               */
#ifdef _OPENMP
    PastixIntegerType             nbthread = omp_get_max_threads();           /* Number of thread wanted by user                           */
#else
        PastixIntegerType             nbthread = 1;           /* Number of thread wanted by user                           */
#endif
        
    PastixIntegerType             verbosemode = verbosity;        /* Level of verbose mode (0, 1, 2)                           */
     int             ordering = API_ORDER_SCOTCH;           /* Ordering to use                                           */
    PastixIntegerType             nbrhs = 1;
    //int             incomplete = 1;         /* Indicate if we want to use incomplete factorisation       */
    int             level_of_fill = ilu_level_of_fill; //6;      /* Level of fill for incomplete factorisation                */
    int             amalgamation = 25;       /* Level of amalgamation for Kass                            */
    //int             ooc = 2000;                /* OOC limit (Mo/percent depending on compilation options)   */
    PastixIntegerType    mat_type;
//    int j;
//        long            i;
//        double norme1, norme2;
    int i;
//    printf("aaa\n");

    if(symmetric == 0)
         mat_type = API_SYM_NO;
    else
         mat_type = API_SYM_YES;

/*    memset(colptr,0,(mat_size+1)*sizeof(PastixIntegerType));
    memset(rows,0,(nnz)*sizeof(PastixIntegerType));*/

    //compute the transpose
       CompRow_to_CompCol(mat_size, mat_size, nnz,AA, JA, IA,&values, &rows, &colptr);


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

//     for(i=0; i<nnz; i++)
//       rows[i] = rows[i]+1;
//
//     for(i=0; i<mat_size; i++)
//       colptr[i] = colptr[i]+1;
//    /*
//     * Matrix needs :
//     *    - to be in fortran numbering
//     *    - to have only the lower triangular part in symmetric case
//     *    - to have a graph with a symmetric structure in unsymmetric case
//     */
    iparm[IPARM_MATRIX_VERIFICATION] = API_YES;
    if(NO_ERR != pastix_checkMatrix(0, verbosemode,
                                    mat_type,
                                    API_YES,
                                    ncol, &colptr, &rows, &values, NULL, 1))
        return 1;

    //copy to block format if needed
    if(ndof > 1)
        CompressCSCinFortranNumbering(mat_size,ndof,&values, &rows, &colptr);

    //copy b to the solution. It will be overwritten
    for(i = 0; i < mat_size; i++)
        x[i] = b[i];
//    printf("bbb\n")    ;
    /*******************************************/
    /* Initialize parameters to default values */
    /*******************************************/
    iparm[IPARM_MODIFY_PARAMETER] = API_NO;
    pastix(&pastix_data, MPI_COMM_WORLD,
           ncol, colptr, rows, values,
           perm, invp, x, 1, iparm, dparm);
//    printf("ccc\n")    ;

    /*******************************************/
    /*       Customize some parameters         */
    /*******************************************/
    iparm[IPARM_THREAD_NBR] = nbthread;
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
    iparm[IPARM_VERBOSE]             = verbosemode;
     iparm[IPARM_ORDERING]            = ordering;
    iparm[IPARM_DOF_NBR] = ndof;
    if(incomplete == 1)
        iparm[IPARM_INCOMPLETE]          = API_YES;
    else if(incomplete == 0)
        iparm[IPARM_INCOMPLETE]          = API_NO;
    else
        printf("incomplete flag should either be 0 (direct solve) or 1 for ILU solve");

    //iparm[IPARM_OOC_LIMIT]           = ooc;

  //  iparm[IPARM_FREE_CSCUSER] = API_CSC_PRESERVE;

    if(iparm[IPARM_INCOMPLETE]  == API_YES)
    {
        iparm[IPARM_REFINEMENT] = API_RAF_GMRES;
        iparm[IPARM_GMRES_IM] = m_gmres;
        dparm[DPARM_EPSILON_REFINEMENT] = tol;
        iparm[IPARM_LEVEL_OF_FILL]       = level_of_fill;
        iparm[IPARM_AMALGAMATION_LEVEL]  = amalgamation;
    }

    iparm[IPARM_RHS_MAKING]          = API_RHS_B;
    iparm[IPARM_START_TASK]          = API_TASK_ORDERING;
    iparm[IPARM_END_TASK]            = API_TASK_CLEAN;




    /*******************************************/
    /*           Call pastix                   */
    /*******************************************/
    //note that we pass "x" instead of the rhs
    pastix(&pastix_data, 0,
           ncol/ndof, colptr, rows, values,
           perm, invp, x, nbrhs, iparm, dparm);

//        PRINT_RHS("SOL", rhs, ncol, mpid, iparm[IPARM_VERBOSE]);
//        CHECK_SOL(rhs, rhssaved, ncol, mpid);




    free(colptr);
    free(rows);
    free(values);
    free(perm);
    free(invp);
    free(rhs);
//        free(rhssaved);
//        free(type);
//        free(rhstype);

    return EXIT_SUCCESS;
}
#undef MPI_COMM_WORLD
