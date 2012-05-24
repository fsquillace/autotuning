
module load tau
module load java

export TAU_MAKEFILE=${TAULIBDIR}/Makefile.tau-papi-mpi-pdt-pgi
export TAU_COMM_MATRIX=1
export TAU_CALLPATH=1
export TAU_CALLPATH_DEPTH=100
export TAU_METRICS=TIME:PAPI_FP_INS:PAPI_L2_DCM:PAPI_L1_DCM
