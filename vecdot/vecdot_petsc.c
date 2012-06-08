
/* Program usage:  mpiexec -n <procs> ex2 [-help] [all PETSc options] */ 

static char help[] = "Solves a sparse matrix-vector multiplication and a vector vector multiplication.\n\
Input parameters include:\n\
  -m <mesh_x>       : number of mesh points in x-direction\n\
  -n <mesh_n>       : number of mesh points in y-direction\n\n";

/*T
   Concepts: Sparse matrix-vector multiplication SpMV
   Processors: n
T*/

#include <petscmat.h>

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **args)
{
  Vec            x,y;  /* approx solution, RHS, exact solution */
  PetscInt       m = 8,n = 7;
  PetscErrorCode ierr;
  PetscScalar    v;
#if defined(PETSC_USE_LOG)
  PetscLogStage  stage;
#endif

  PetscInitialize(&argc,&args,(char *)0,help);
  ierr = PetscOptionsGetInt(PETSC_NULL,"-m",&m,PETSC_NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(PETSC_NULL,"-n",&n,PETSC_NULL);CHKERRQ(ierr);

  /* 
     Create parallel vectors.
      - We form 1 vector from scratch and then duplicate as needed.
      - When using VecCreate(), VecSetSizes and VecSetFromOptions()
        in this example, we specify only the
        vector's global dimension; the parallel partitioning is determined
        at runtime. 
      - When solving a linear system, the vectors and matrices MUST
        be partitioned accordingly.  PETSc automatically generates
        appropriately partitioned matrices and vectors when MatCreate()
        and VecCreate() are used with the same communicator.  
      - The user can alternatively specify the local vector and matrix
        dimensions when more sophisticated partitioning is needed
        (replacing the PETSC_DECIDE argument in the VecSetSizes() statement
        below).
  */
  ierr = VecCreate(PETSC_COMM_WORLD,&x);CHKERRQ(ierr);
  ierr = VecSetSizes(x,PETSC_DECIDE,m*n);CHKERRQ(ierr);
  ierr = VecSetFromOptions(x);CHKERRQ(ierr);
  ierr = VecDuplicate(x,&y);CHKERRQ(ierr);
  ierr = VecSetFromOptions(y);CHKERRQ(ierr);

  /* 
     Set exact solution; then compute right-hand-side vector.
     By default we use the vector x with all elements of 1.0;
  */
  ierr = VecSet(x,1.0);CHKERRQ(ierr);
  ierr = VecSet(y,1.9);CHKERRQ(ierr);
  MPI_Pcontrol( 1,"VecDot: x*y");
  ierr = VecDot(x,y,&v);CHKERRQ(ierr);
  MPI_Pcontrol( -1,"VecDot: x*y");

  /*
     Free work space.  All PETSc objects should be destroyed when they
     are no longer needed.
  */
  ierr = VecDestroy(&x);CHKERRQ(ierr);  ierr = VecDestroy(&y);CHKERRQ(ierr);

  /*
     Always call PetscFinalize() before exiting a program.  This routine
       - finalizes the PETSc libraries as well as MPI
       - provides summary and diagnostic information if certain runtime
         options are chosen (e.g., -log_summary). 
  */
  ierr = PetscFinalize();
  return 0;
}
