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
	CFLAGS = -g -O3 -pthread -openmp 
else
	CFLAGS = -g -O3 -pthread -fopenmp
endif

SRC=matmult_poski

all: ${SRC}-shared

${SRC}-shared: ${SRC}.o
	$(POSKICC) $(CFLAGS) -o $@ ${SRC}.o $(LDFLAGS_SHARED)

.c.o:
	$(POSKICC) $(CFLAGS) -o $@ -c $< $(LDFLAGS_SHARED) $(DEF)

clean:
	rm -rf ${SRC}-shared ${SRC}.o core*~ 

#eof


