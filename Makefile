TARGET	= libmatmult.so
LIBSRCS	=  matmult_lib.c matmult_mnk.c matmult_mkn.c matmult_nkm.c matmult_nmk.c matmult_kmn.c matmult_knm.c
LIBOBJS	= matmult_lib.o matmult_mnk.o matmult_mkn.o matmult_nkm.o matmult_nmk.o matmult_kmn.o matmult_knm.o

OPT	= -g # -O3 -fopt-info
PIC	= -fPIC

CC	= gcc
CFLAGS= $(OPT) $(PIC) $(XOPTS)

SOFLAGS = -shared 
XLIBS	= -L/usr/lib64/atlas -lsatlas

$(TARGET): $(LIBOBJS)
	$(CC) -o $@ $(SOFLAGS) $(LIBOBJS) $(XLIBS)

clean:
	@/bin/rm -f core core.* $(LIBOBJS) TARGET= libmatmult.so
