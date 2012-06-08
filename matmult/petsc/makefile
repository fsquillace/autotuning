CC = cc
CFLAGS	         =
CPPFLAGS         =

############################# PETSC RULES #####################

ifneq ($(PETSC_DIR),)
  include ${PETSC_DIR}/conf/variables
  include ${PETSC_DIR}/conf/rules
endif


SRC_PETSC_MATMULT = matmult_petsc

${SRC_PETSC_MATMULT}: ${SRC_PETSC_MATMULT}.o
	-${CLINKER} -o ${SRC_PETSC_MATMULT} ${SRC_PETSC_MATMULT}.o  ${PETSC_KSP_LIB}
	${RM} ${SRC_PETSC_MATMULT}.o


######################## IPM RULES ################################

SRC_PETSC_IPM_MATMULT = matmult_petsc_ipm

${SRC_PETSC_IPM_MATMULT}: ${SRC_PETSC_MATMULT}.o
	-${CLINKER} -o ${SRC_PETSC_IPM_MATMULT} ${SRC_PETSC_MATMULT}.o  ${PETSC_KSP_LIB} ${IPM}
	${RM} ${SRC_PETSC_MATMULT}.o


##########################  TAU RULES ####################################

TAU_PCC              = tau_cc.sh -tau_makefile=$(TAU_MAKEFILE) $(TAUOPTIONS)
CFLAGS= -I${TAUROOTDIR}/include
#TAULIBS=-L${TAULIBDIR} -lTauMpi-papi-mpi-pdt-pgi

CLINKER          = ${PCC}

         ################### AUTOMATIC ##################

SRC_TAU_PETSC_MATMULT = matmult_petsc_tau_auto

${SRC_TAU_PETSC_MATMULT}: ${SRC_PETSC_MATMULT}.o
	-${TAU_PCC} -o ${SRC_TAU_PETSC_MATMULT} ${SRC_PETSC_MATMULT}.o  ${PETSC_KSP_LIB}
	${RM} ${SRC_PETSC_MATMULT}.o

        ################### MANUAL ########################
SRC_TAU_MAN_PETSC_MATMULT = matmult_petsc_tau_man


${SRC_TAU_MAN_PETSC_MATMULT}: ${SRC_PETSC_MATMULT}.o
	-${TAU_PCC} -o ${SRC_TAU_MAN_PETSC_MATMULT} ${SRC_PETSC_MATMULT}.o  ${PETSC_KSP_LIB}
	${RM} ${SRC_PETSC_MATMULT}.o


################## CLEAN #######################

#clean: $(SRC_PETSC_MATMULT)
	#rm -rf ${SRC_PETSC_MATMULT}

