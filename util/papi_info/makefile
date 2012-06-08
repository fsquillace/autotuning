CC = cc
CFLAG = –g –O3 –pthread –openmp

LIBS = $(PAPI_POST_LINK_OPTS)
INCS = $(PAPI_INCLUDE_OPTS)

FILE = papi_info

papi_info: $(FILE).c
	$(CC) $(LIBS) $(INCS) -o $(FILE) $(FILE).c

clean:
	rm –rf $(FILE)
