CC = cc
CFLAGS	         =
CPPFLAGS         =


############################ PAPI INFO #######################
LIBS_PAPI = $(PAPI_POST_LINK_OPTS)
INCS_PAPI = $(PAPI_INCLUDE_OPTS)

SRC_PAPI = papi_info

papi_info: $(SRC_PAPI).c
	$(CC) $(LIBS_PAPI) $(INCS_PAPI) -o $(SRC_PAPI) $(SRC_PAPI).c

clean_papi:
	rm –rf $(SRC_PAPI)


############################# PETSC RULES #####################

ifneq ($(PETSC_DIR),)
  include ${PETSC_DIR}/conf/variables
  include ${PETSC_DIR}/conf/rules
endif

SRC_PETSC_LMLS = lmls_petsc
SRC_PETSC_MATMULT = matmult_petsc

${SRC_PETSC_LMLS}: ${SRC_PETSC_LMLS}.o
	-${CLINKER} -o ${SRC_PETSC_MATMULT} ${SRC_PETSC_MATMULT}.o  ${PETSC_KSP_LIB}
	${RM} ${SRC_PETSC_LMLS}.o

${SRC_PETSC_MATMULT}: ${SRC_PETSC_MATMULT}.o
	-${CLINKER} -o ${SRC_PETSC_MATMULT} ${SRC_PETSC_MATMULT}.o  ${PETSC_KSP_LIB}
	${RM} ${SRC_PETSC_MATMULT}.o

clean_matmult_petsc: $(SRC_PETSC_MATMULT)
	rm -rf ${SRC_PETSC_MATMULT}

#########################  TAU/PETSC RULES ####################################

SRC_TAU_PETSC_MATMULT = matmult_tau_petsc

TAU_PCC              = tau_cc.sh -tau_makefile=$(TAU_MAKEFILE) $(TAUOPTIONS)
#CLINKER          = ${PCC}

${SRC_TAU_PETSC_MATMULT}: ${SRC_PETSC_MATMULT}.o
	-${TAU_PCC} -o ${SRC_TAU_PETSC_MATMULT} ${SRC_PETSC_MATMULT}.o  ${PETSC_KSP_LIB}
	${RM} ${SRC_PETSC_MATMULT}.o

clean_matmult_tau_petsc: $(SRC_TAU_PETSC_MATMULT)
	rm -rf ${SRC_TAU_PETSC_MATMULT}

########################## POSKI RULES ####################################
#Location of OSKI
OSKIDIR = /opt/oski
POSKIDIR = /opt/poski

OSKIINCS = $(OSKIDIR)/include
OSKILIBS = $(OSKIDIR)/lib/oski

#Location of pOSKI
POSKILIB = $(POSKIDIR)/lib
POSKIINC = $(POSKIDIR)/include/poski

#OSKI link flags
OSKILIBS_SHARED = -I$(OSKIINCS) -Wl,-rpath -Wl,$(OSKILIBS) -L$(OSKILIBS) `cat $(OSKILIBS)/site-modules-shared.txt` -loski

#pOSKI link flags
POSKILIBS_SHARED = -I$(POSKIINC) -Wl,-rpath -Wl,$(POSKILIB) -L$(POSKILIB) -lposki

#library link
LDFLAGS_SHARED = $(OSKILIBS_SHARED) $(POSKILIBS_SHARED) $(LDFLAGS) -lm

DEF = #-DVERBOSE

#C compiler & flags (gcc: -fopenmp, icc: -openmp)
ifeq ($(POSKICC),)
POSKICC=${CC}
endif
ifeq ($(POSKICC),icc)
    POSKICFLAGS = -g -O3 -pthread -openmp
else
    POSKICFLAGS = -g -O3 -pthread -fopenmp
endif

SRC_POSKI=matmult_poski

${SRC_POSKI}: ${SRC_POSKI}.o
	$(POSKICC) $(POSKICFLAGS) -o $@ ${SRC_POSKI}.o $(LDFLAGS_SHARED)

${SRC_POSKI}.o:
	$(POSKICC) $(POSKICFLAGS) -o ${SRC_POSKI}.o -c ${SRC_POSKI}.c $(LDFLAGS_SHARED) $(DEF)

clean_poski:
	rm -rf ${SRC_POSKI} ${SRC_POSKI}.o core*~ 

#eof

