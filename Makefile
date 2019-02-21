include make.inc

.SUFFIXES: .o .c

DEFS =

OBJS =  SRC/PAL_coo_to_csr.o \
	SRC/PAL_csr_matvec.o \
	SRC/PAL_effective_W.o \
	SRC/PAL_effective_W_QR.o \
        SRC/PAL_free_Matrix.o \
	SRC/PAL_free_RT.o \
	SRC/PAL_geneig.o \
	SRC/PAL_lineig.o \
	SRC/PAL_op_A_solve.o \
	SRC/PAL_op_B_matvec.o \
	SRC/PAL_op_S_build.o \
	SRC/PAL_pade_approx.o \
	SRC/PAL_pade_coeff.o \
	SRC/PAL_read_matrix.o \
	SRC/PAL_residual.o \
	SRC/PAL_SuperLU_dgstrf.o \
	SRC/PAL_SuperLU_dgstrs.o \
	SRC/PAL_SuperLU_free.o \
	SRC/PAL_V_kron_a.o \
	SRC/PAL_V_kron_b.o \
	SRC/PAL_wtime.o \
	SRC/main.o

.PHONY: all clean

all: main_PAL.x

main_PAL.x: SRC/main.o $(OBJS)  
	$(LOADER) $(LDFLAGS) -o $@ $^ $(LIBS)

clean:
	$(RM) $(OBJS) *.x

.c.o:
	$(CC) $(DEFS) $(CFLAGS) $(INCS) -o $@ -c $<
