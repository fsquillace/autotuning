
CFLAGS	         =
CPPFLAGS         =

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

lmls: lmls.o
	-${CLINKER} -o lmls lmls.o  ${PETSC_KSP_LIB}
	${RM} lmls.o

