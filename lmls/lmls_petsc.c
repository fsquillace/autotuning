
static char help[] = "Reads several matrices (M11,M12,L11,L22,L21) obtained from the discretization of the problem from the neutronic diffusion equation and computes the operation y=Ax where A = L11^{-1}*(M11+M12*L22^{-1}*L21).\n";

#include <petsc.h>
#include <petscksp.h>

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **args)
{
  const int iM11=0, iM12=1, iL11=2, iL22=3, iL21=4, ix=5;

  Mat            M11,M12,L11,L22,L21;  /* matrix */
  Vec            x,y;                  /* input and output vectors */
  Vec            omg1,omg2,omg3,omg4;  /* temporary vectors for the operation y=Ax */
  KSP            ksp;                  /* linear solver context */
  PetscViewer    fd[5];                /* viewer */
  char           file[6][PETSC_MAX_PATH_LEN];   /* input file name */

  PetscErrorCode ierr;
  PetscInt       M, N;                 /* number of rows and columns of the GLOBAL matrices (they should be the same) */
  PetscInt       m, n;                 /* number of rows and columns of the LOCAL matrices */
  PetscInt       istart, iend;         /* ownership row range of the process using GLOBAL indexes */
  PetscScalar    one=1.0;
  PetscMPIInt    rank, size;
  PetscBool      flg[6];

  PetscInitialize(&argc,&args,(char *)0,help);

  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);


  PetscPrintf(PETSC_COMM_WORLD, "Number of processes size=%d\n", size);



  // ************ READ MATRICES AND VECTOR x FROM INPUT *****


  // M11 and M12
  ierr = PetscOptionsGetString(PETSC_NULL,"-m11",file[iM11],PETSC_MAX_PATH_LEN,&flg[iM11]);CHKERRQ(ierr);
  ierr = PetscOptionsGetString(PETSC_NULL,"-m12",file[iM12],PETSC_MAX_PATH_LEN,&flg[iM12]);CHKERRQ(ierr);

  // L11 L22 and L21
  ierr = PetscOptionsGetString(PETSC_NULL,"-l11",file[iL11],PETSC_MAX_PATH_LEN,&flg[iL11]);CHKERRQ(ierr);
  ierr = PetscOptionsGetString(PETSC_NULL,"-l22",file[iL22],PETSC_MAX_PATH_LEN,&flg[iL22]);CHKERRQ(ierr);
  ierr = PetscOptionsGetString(PETSC_NULL,"-l21",file[iL21],PETSC_MAX_PATH_LEN,&flg[iL21]);CHKERRQ(ierr);

  // All of the matrix have to be defined by the user.
  // If the user don't specify none of them, it will generate laplacian matrices.
  if (flg[iM11] && flg[iM12] && flg[iL11] && flg[iL22] && flg[iL21]){

  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,file[iM11],FILE_MODE_READ,\
          &fd[iM11]);CHKERRQ(ierr);
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,file[iM12],FILE_MODE_READ,\
          &fd[iM12]);CHKERRQ(ierr);
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,file[iL11],FILE_MODE_READ,\
          &fd[iL11]);CHKERRQ(ierr);
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,file[iL22],FILE_MODE_READ,\
          &fd[iL22]);CHKERRQ(ierr);
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,file[iL21],FILE_MODE_READ,\
         &fd[iL21]);CHKERRQ(ierr);

    // Load the matrix and vector; then destroy the viewer.

    // M11 and M12
   ierr = MatCreate(PETSC_COMM_WORLD,&M11);CHKERRQ(ierr);
    ierr = MatSetType(M11,MATAIJ);CHKERRQ(ierr);
   ierr = MatSetFromOptions(M11);CHKERRQ(ierr);
   ierr = MatLoad(M11,fd[iM11]);CHKERRQ(ierr);
   ierr = MatCreate(PETSC_COMM_WORLD,&M12);CHKERRQ(ierr);
    ierr = MatSetType(M12,MATAIJ);CHKERRQ(ierr);
   ierr = MatSetFromOptions(M12);CHKERRQ(ierr);
   ierr = MatLoad(M12,fd[iM12]);CHKERRQ(ierr);

    // L11 L22 and L21
   ierr = MatCreate(PETSC_COMM_WORLD,&L11);CHKERRQ(ierr);
    ierr = MatSetType(L11,MATAIJ);CHKERRQ(ierr);
   ierr = MatSetFromOptions(L11);CHKERRQ(ierr);
   ierr = MatLoad(L11,fd[iL11]);CHKERRQ(ierr);
   ierr = MatCreate(PETSC_COMM_WORLD,&L22);CHKERRQ(ierr);
    ierr = MatSetType(L22,MATAIJ);CHKERRQ(ierr);
   ierr = MatSetFromOptions(L22);CHKERRQ(ierr);
   ierr = MatLoad(L22,fd[iL22]);CHKERRQ(ierr);
   ierr = MatCreate(PETSC_COMM_WORLD,&L21);CHKERRQ(ierr);
    ierr = MatSetType(L21,MATAIJ);CHKERRQ(ierr);
   ierr = MatSetFromOptions(L21);CHKERRQ(ierr);
   ierr = MatLoad(L21,fd[iL21]);CHKERRQ(ierr);

    ierr = PetscViewerDestroy(&fd[iM11]);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&fd[iM12]);CHKERRQ(ierr);

    ierr = PetscViewerDestroy(&fd[iL11]);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&fd[iL22]);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&fd[iL21]);CHKERRQ(ierr);
  }
  else if(!flg[iM11] && !flg[iM12] && !flg[iL11] && !flg[iL22] && !flg[iL21]){
  // ******************* CREATING FAKE MATRICES *****************
  PetscInt       i,col[3];
  M = N = 100;
  PetscScalar    value[3];

  ierr = MatCreate(PETSC_COMM_WORLD,&M11);CHKERRQ(ierr);
  ierr = MatSetSizes(M11,PETSC_DECIDE,PETSC_DECIDE,M,N);CHKERRQ(ierr);
  ierr = MatSetFromOptions(M11);CHKERRQ(ierr);
  ierr = MatCreate(PETSC_COMM_WORLD,&M12);CHKERRQ(ierr);
  ierr = MatSetSizes(M12,PETSC_DECIDE,PETSC_DECIDE,M,N);CHKERRQ(ierr);
  ierr = MatSetFromOptions(M12);CHKERRQ(ierr);

  ierr = MatCreate(PETSC_COMM_WORLD,&L11);CHKERRQ(ierr);
  ierr = MatSetSizes(L11,PETSC_DECIDE,PETSC_DECIDE,M,N);CHKERRQ(ierr);
  ierr = MatSetFromOptions(L11);CHKERRQ(ierr);
  ierr = MatCreate(PETSC_COMM_WORLD,&L22);CHKERRQ(ierr);
  ierr = MatSetSizes(L22,PETSC_DECIDE,PETSC_DECIDE,M,N);CHKERRQ(ierr);
  ierr = MatSetFromOptions(L22);CHKERRQ(ierr);
  ierr = MatCreate(PETSC_COMM_WORLD,&L21);CHKERRQ(ierr);
  ierr = MatSetSizes(L21,PETSC_DECIDE,PETSC_DECIDE,M,N);CHKERRQ(ierr);
  ierr = MatSetFromOptions(L21);CHKERRQ(ierr);

  // Set values for them
  value[0] = -1.0; value[1] = 2.0; value[2] = -1.0;
  for (i=1; i<M-1; i++) {
    col[0] = i-1; col[1] = i; col[2] = i+1;
    ierr = MatSetValues(M11,1,&i,3,col,value,INSERT_VALUES);CHKERRQ(ierr);
    ierr = MatSetValues(M12,1,&i,3,col,value,INSERT_VALUES);CHKERRQ(ierr);

    ierr = MatSetValues(L11,1,&i,3,col,value,INSERT_VALUES);CHKERRQ(ierr);
    ierr = MatSetValues(L22,1,&i,3,col,value,INSERT_VALUES);CHKERRQ(ierr);
    ierr = MatSetValues(L21,1,&i,3,col,value,INSERT_VALUES);CHKERRQ(ierr);
  }
  i = N - 1; col[0] = N - 2; col[1] = N - 1;
  ierr = MatSetValues(M11,1,&i,2,col,value,INSERT_VALUES);CHKERRQ(ierr);
  ierr = MatSetValues(M12,1,&i,2,col,value,INSERT_VALUES);CHKERRQ(ierr);

  ierr = MatSetValues(L11,1,&i,2,col,value,INSERT_VALUES);CHKERRQ(ierr);
  ierr = MatSetValues(L22,1,&i,2,col,value,INSERT_VALUES);CHKERRQ(ierr);
  ierr = MatSetValues(L21,1,&i,2,col,value,INSERT_VALUES);CHKERRQ(ierr);

  i = 0; col[0] = 0; col[1] = 1; value[0] = 2.0; value[1] = -1.0;
  ierr = MatSetValues(M11,1,&i,2,col,value,INSERT_VALUES);CHKERRQ(ierr);
  ierr = MatSetValues(M12,1,&i,2,col,value,INSERT_VALUES);CHKERRQ(ierr);

  ierr = MatSetValues(L11,1,&i,2,col,value,INSERT_VALUES);CHKERRQ(ierr);
  ierr = MatSetValues(L22,1,&i,2,col,value,INSERT_VALUES);CHKERRQ(ierr);
  ierr = MatSetValues(L21,1,&i,2,col,value,INSERT_VALUES);CHKERRQ(ierr);

  ierr = MatAssemblyBegin(M11,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(M11,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyBegin(M12,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(M12,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  ierr = MatAssemblyBegin(L11,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(L11,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyBegin(L22,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(L22,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyBegin(L21,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(L21,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  // ******************** END CREATING FAKE MATRICES *************

  }
  else{
    SETERRQ(PETSC_COMM_WORLD,1,"You must either indicate the ascii file for each matrix M11 M12 L11 L22 L21 using the option -m11 -m12 .... or without any of these options.");
    PetscFinalize();
    return 1;
  }


  // Get some information about the partitioning of the matrix
  ierr = MatGetSize(M11,&M,&N);CHKERRQ(ierr);
  printf("Global dimension of the matrix M11 M=%d N=%d\n",M,N);
  ierr = MatGetLocalSize(M11,&m,&n);
  printf("Local dimension of the matrix M11 m=%d n=%d\n",m,n);
  ierr = MatGetOwnershipRange(M11,&istart,&iend);
  printf("Ownership range of the rows for process %d istart=%d iend=%d\n",rank,istart,iend);


  // Read vector x from input. If it's not specified by the user, the vector x will be a unitary vector.
  ierr = PetscOptionsGetString(PETSC_NULL,"-x",file[ix],PETSC_MAX_PATH_LEN,&flg[ix]);CHKERRQ(ierr);
  if(!flg[ix]){
    ierr = VecCreate(PETSC_COMM_WORLD,&x);CHKERRQ(ierr);
    ierr = VecSetSizes(x,PETSC_DECIDE,M);CHKERRQ(ierr);
    ierr = VecSetFromOptions(x);CHKERRQ(ierr);
    ierr = VecSet(x,one);CHKERRQ(ierr);
    /*ierr = VecAssemblyBegin(x);CHKERRQ(ierr);*/
    /*ierr = VecAssemblyEnd(x);CHKERRQ(ierr);*/
  }
  else{
          ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,file[ix],FILE_MODE_READ,&fd[ix]);CHKERRQ(ierr);
          ierr = VecLoad(x,fd[ix]);CHKERRQ(ierr);

  }
  ierr = PetscObjectSetName((PetscObject) x, "The input vector");CHKERRQ(ierr);

  // ************ END READ MATRICES AND VECTOR x FROM INPUT *****

  // Create the temporary vectors and y
  ierr = VecDuplicate(x,&y);CHKERRQ(ierr);
  ierr = VecDuplicate(x,&omg1);CHKERRQ(ierr);
  ierr = VecDuplicate(x,&omg2);CHKERRQ(ierr);
  ierr = VecDuplicate(x,&omg3);CHKERRQ(ierr);
  ierr = VecDuplicate(x,&omg4);CHKERRQ(ierr);


  // ****************** COMPUTE y=Ax *******************

  // Set the Krylov object
  ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);

   //  Set operators. Here the matrix that defines the linear system
   //  also serves as the preconditioning matrix.
  ierr = KSPSetOperators(ksp,L22,L22,DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);

  //  Set runtime options, e.g.,
  //      -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
  //  These options will override those specified above as long as
  //  KSPSetFromOptions() is called _after_ any other customization
  // routines.
  ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);


  // Multiplication y = A*x

  // omg1 = M11*x
  // omg2 = L21*x
  // L22*omg3 = omg2
  // omg4 = omg1 + M12*omg3
  // L11*y = omg4

  ierr = MatMult(M11, x, omg1);CHKERRQ(ierr);

  ierr = MatMult(L21, x, omg2);CHKERRQ(ierr);


  ierr = KSPSolve(ksp,omg2,omg3);CHKERRQ(ierr);


  ierr = MatMult(M12, omg3, omg4);CHKERRQ(ierr);
  ierr = VecAXPY(omg4, one, omg1);CHKERRQ(ierr);


  ierr = KSPSetOperators(ksp,L11,L11,DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);
  ierr = KSPSolve(ksp,omg4,y);CHKERRQ(ierr);

  // ****************** END COMPUTE y=Ax **************************

  /*ierr = VecView(y, PETSC_VIEWER_STDOUT_SELF);CHKERRQ(ierr);*/
  PetscBool test = PETSC_FALSE;
  ierr = PetscOptionsGetBool(PETSC_NULL,"-test",&test,PETSC_NULL);CHKERRQ(ierr);
  if(test){
 
  // ******************** TESTING *********************************


    // the testing doesn't work if the number of process are more than one
    // because the type of the matrices must be different from matmpiaij. Let's try matseqaij
  Mat L11_inv, L22_inv, I;
  Mat A;
  Mat M11_d;
  Vec y2;

  PetscInt       i;
  PetscScalar    val;

  // Create identity matrix
  ierr = MatCreate(PETSC_COMM_WORLD,&I);CHKERRQ(ierr);
  ierr = MatSetType(I, MATDENSE);
  ierr = MatSetSizes(I,PETSC_DECIDE,PETSC_DECIDE,M,N);CHKERRQ(ierr);
  ierr = MatSetFromOptions(I);CHKERRQ(ierr);
  val = 1.0;
  for (i=0; i<M; i++) {
    ierr = MatSetValues(I,1,&i,1,&i,&val,INSERT_VALUES);CHKERRQ(ierr);
  }
  ierr = MatAssemblyBegin(I,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(I,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  // Create L11_inv
  ierr = MatCreate(PETSC_COMM_WORLD,&L11_inv);CHKERRQ(ierr);
  ierr = MatSetType(L11_inv, MATDENSE);
  ierr = MatSetSizes(L11_inv,PETSC_DECIDE,PETSC_DECIDE,M,N);CHKERRQ(ierr);
  ierr = MatSetFromOptions(L11_inv);CHKERRQ(ierr);
  ierr = MatAssemblyBegin(L11_inv,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(L11_inv,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  // Create L22_inv
  ierr = MatCreate(PETSC_COMM_WORLD,&L22_inv);CHKERRQ(ierr);
  ierr = MatSetType(L22_inv, MATDENSE);
  ierr = MatSetSizes(L22_inv,PETSC_DECIDE,PETSC_DECIDE,M,N);CHKERRQ(ierr);
  ierr = MatSetFromOptions(L22_inv);CHKERRQ(ierr);
  ierr = MatAssemblyBegin(L22_inv,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(L22_inv,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  // Create A
  ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
  ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,M,N);CHKERRQ(ierr);
  ierr = MatSetFromOptions(A);CHKERRQ(ierr);
  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  // Define M11_d
  ierr = MatCreate(PETSC_COMM_WORLD,&M11_d);CHKERRQ(ierr);
  ierr = MatSetType(M11_d, MATDENSE);
  ierr = MatSetSizes(M11_d,PETSC_DECIDE,PETSC_DECIDE,M,N);CHKERRQ(ierr);
  ierr = MatSetFromOptions(M11_d);CHKERRQ(ierr);
  ierr = MatAssemblyBegin(M11_d,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(M11_d,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);


  // Calculate A = L11^{-1}*(M11 + M12*L22^{-1}*L21)
  IS perm, iperm;
  MatFactorInfo info;
  ierr = MatGetOrdering(L11,MATORDERINGNATURAL,&perm,&iperm);CHKERRQ(ierr);
  ierr = MatFactorInfoInitialize(&info); CHKERRQ(ierr);
  ierr = MatLUFactor(L11, perm, iperm, &info); CHKERRQ(ierr);

  ierr = MatMatSolve(L11, I, L11_inv);CHKERRQ(ierr);
  // TODO try to convert L11_inv to be sparse such as matseqaij

//  ierr = MatView(L11_inv, PETSC_VIEWER_STDOUT_SELF);CHKERRQ(ierr);
  ierr = MatGetOrdering(L22,MATORDERINGNATURAL,&perm,&iperm);CHKERRQ(ierr);
  ierr = MatFactorInfoInitialize(&info); CHKERRQ(ierr);
  ierr = MatLUFactor(L22, perm, iperm, &info); CHKERRQ(ierr);

  ierr = MatMatSolve(L22, I, L22_inv);CHKERRQ(ierr);



  ierr = MatMatMult(M12, L22_inv, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &A);CHKERRQ(ierr);
  ierr = MatMatMult(A, L21, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &A);CHKERRQ(ierr);

  ierr = MatConvert(M11, MATDENSE, MAT_INITIAL_MATRIX, &M11_d);
  ierr = MatAXPY(A,1.0,M11_d,DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);

  ierr = MatMatMult(L11_inv, A, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &A);CHKERRQ(ierr);

  ierr = VecDuplicate(x,&y2);CHKERRQ(ierr);
  ierr = MatMult(A, x, y2);CHKERRQ(ierr);

  /*ierr = VecView(y2, PETSC_VIEWER_STDOUT_SELF);CHKERRQ(ierr);*/


  // Check the error
  PetscReal norm;
  ierr = VecAXPY(y2,-1.0,y);CHKERRQ(ierr);
  ierr  = VecNorm(y2,NORM_2,&norm);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Norm of error %A\n",
                     norm);CHKERRQ(ierr);

  // ******************** END TESTING *****************************
  }


  PetscFinalize();
  return 0;
}
