CC = cc
CFLAGS	         =
CPPFLAGS         =


############################# PETSC RULES #####################

ifneq ($(PETSC_DIR),)
  include ${PETSC_DIR}/conf/variables
  include ${PETSC_DIR}/conf/rules
endif

SRC_PETSC_LMLS = lmls_petsc

${SRC_PETSC_LMLS}: ${SRC_PETSC_LMLS}.o
	-${CLINKER} -o ${SRC_PETSC_MATMULT} ${SRC_PETSC_MATMULT}.o  ${PETSC_KSP_LIB}
	${RM} ${SRC_PETSC_LMLS}.o

clean: $(SRC_PETSC_LMLS)
	rm -rf ${SRC_PETSC_LMLS}


