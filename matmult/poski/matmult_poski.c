#include <stdio.h>
#include <poski.h>


int main(int argc, char **argv)
{

#ifdef VERBOSE
	printf("\n######################################################\n");
	printf("# A First Example:\n");
	printf("#   - This example computes SpMV (y=alpha*Ax + beta*y)\n");
	printf("######################################################\n");
#endif
        int m=2, n=2;    // Dimensions of the laplacian matrix
        int nnz;         // Non zeros in the laplacian matrix
        int i,j,Ii,J;

        for(i=0; i<argc;i++){
            if(strcmp(argv[i],"--m")==0){
                m = atoi(argv[i+1]);
                i++;
            }
            else if(strcmp(argv[i],"--n")==0){
                n = atoi(argv[i+1]);
                i++;
            }
        }

        // Calculates the number of non zeros values in
        // laplacian matrix for CSR format.
        nnz=0;
        for (Ii=0; Ii<n*m; Ii++) {
            i = Ii/n; j = Ii - i*n;
            if (i>0){
                nnz++;
            }
            if (i<m-1){
                nnz++;
            }
            if (j>0){
                nnz++;
            }
            if (j<n-1){
                nnz++;
            }
            nnz++;
        }

        printf("n=%d m=%d nnz=%d Problem size=%d\n",n, m, nnz, m*n);

	int Aptr[(m*n)+1];
	int Aind[nnz];
	double Aval[nnz];

        double v;
        int ptr=0;

        for (Ii=0; Ii<n*m; Ii++) {
            Aptr[Ii] = ptr;

            i = Ii/n; j = Ii - i*n;

            if (i>0){
                J = Ii - n;
                Aval[ptr] = -1.0;
                Aind[ptr] = J;
                ptr++;
            }
            if (j>0){
                J = Ii - 1;
                Aval[ptr] = -1.0;
                Aind[ptr] = J;
                ptr++;
            }

            Aval[ptr] = 4.0;
            Aind[ptr] = Ii;
            ptr++;

            if (j<n-1){
                J = Ii + 1;
                Aval[ptr] = -1.0;
                Aind[ptr] = J;
                ptr++;
            }
            if (i<m-1){
                J = Ii + n;
                Aval[ptr] = -1.0;
                Aind[ptr] = J;
                ptr++;
            }
        }

        Aptr[n*m]=ptr;

        /*printf("Aptr\n");*/
        /*for(i=0; i<n*m+1;i++){*/
            /*printf("%d ",Aptr[i]);*/
        /*}*/
        /*printf("\n");*/

        /*printf("Aval\n");*/
        /*for(i=0; i<nnz;i++){*/
            /*printf("%f ",Aval[i]);*/
        /*}*/
        /*printf("\n");*/

        /*printf("Aind\n");*/
        /*for(i=0; i<nnz;i++){*/
            /*printf("%d ",Aind[i]);*/
        /*}*/
        /*printf("\n");*/


        double x[n*m];
        double y[n*m];
        for(i=0; i<m*n; i++){
            x[i]=1.0;
            y[i]=1.0;
        }
        double alpha=1,beta=0;

        poski_mat_t A_tunable;
        poski_vec_t x_view, y_view;

        poski_Init();
        poski_threadarg_t *threadargs = poski_InitThreads();	

        A_tunable = poski_CreateMatCSR(Aptr, Aind, Aval,\
                n*m, m*n, nnz, SHARE_INPUTMAT, \
                threadargs, NULL, 2, INDEX_ZERO_BASED, MAT_GENERAL);

        x_view = poski_CreateVecView(x, m*n, STRIDE_UNIT, NULL);
        y_view = poski_CreateVecView(y, m*n, STRIDE_UNIT, NULL);

        poski_MatMult(A_tunable, OP_NORMAL, alpha, x_view, beta, y_view);

        poski_DestroyThreads(threadargs);
        poski_DestroyVec(x_view);
        poski_DestroyVec(y_view);
        poski_DestroyMat(A_tunable);

        poski_Close();

#ifdef VERBOSE
        printf(" Given matrix A:\n");
        poski_report_sparse_CSR_va(Aptr, Aind, Aval, m*n, m*n, nnz);

        printf(" Given vectors:\n");
        printf("\t+ x = [ ");
        for(i=0; i<m*n; i++){
            printf("%.4f; ", x[i]);
        }
        printf("]\n");

        printf("\t+ y = [ ");
        for(i=0; i<m*n; i++){
            printf("%.4f; ", y[i]);
        }
        printf("]\n");
        printf("######################################################\n\n");
#endif

	return 0;
}

