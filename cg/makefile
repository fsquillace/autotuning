CC = cc
CFLAGS	         =
CPPFLAGS         =


############################# PETSC RULES #####################

ifneq ($(PETSC_DIR),)
  include ${PETSC_DIR}/conf/variables
  include ${PETSC_DIR}/conf/rules
endif

SRC_PETSC_CG = cg

${SRC_PETSC_CG}: ${SRC_PETSC_CG}.o
	-${CLINKER} -o ${SRC_PETSC_CG} ${SRC_PETSC_CG}.o  ${PETSC_KSP_LIB}
	${RM} ${SRC_PETSC_CG}.o

########################### IPM RULES ###################

SRC_PETSC_IPM_CG = cg_ipm

${SRC_PETSC_IPM_CG}: ${SRC_PETSC_CG}.o
	-${CLINKER} -o ${SRC_PETSC_IPM_CG} ${SRC_PETSC_CG}.o  ${PETSC_KSP_LIB} ${IPM}
	${RM} ${SRC_PETSC_CG}.o


##########################  TAU/PETSC RULES ####################################

TAU_PCC              = tau_cc.sh -tau_makefile=$(TAU_MAKEFILE) $(TAUOPTIONS)
CLINKER          = ${PCC}
CFLAGS= -I${TAUROOTDIR}/include
#TAULIBS=-L${TAULIBDIR} -lTauMpi-papi-mpi-pdt-pgi



          ################### AUTOMATIC ##################

          ################### MANUAL ########################
SRC_TAU_MAN_PETSC_CG = cg_tau_man


${SRC_TAU_MAN_PETSC_CG}: ${SRC_PETSC_CG}.o
	-${TAU_PCC} -o ${SRC_TAU_MAN_PETSC_CG} ${SRC_PETSC_CG}.o  ${PETSC_KSP_LIB}
	${RM} ${SRC_PETSC_CG}.o



