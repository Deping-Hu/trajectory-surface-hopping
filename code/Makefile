#
# Paths to directories commands
#
DIR = ./
INCDIR = ${DIR}/include
#
# UNIX commands
#
CP    = /bin/cp -f
RM    = /bin/rm -f
MKEXE = /bin/chmod 755
AR    = /usr/bin/ar
#
# Libraries
#
ifdef LIB
LIB0 = ${LIB}
else
LIB = -L/opt/intel/composer_xe_2013.1.117/mkl/lib/intel64 -lmkl_lapack95_lp64 -lmkl_blas95_lp64 -lmkl_core -lmkl_sequential -lmkl_intel_lp64
endif
#
# Compilers
#
F90C = ifort
F90FLAGS = -stand f03 -check all -traceback -warn all -fstack-protector -assume protect_parens -implicitnone

RUNAR = ar rv
#
# Executable, source and objects
#
EXE  = main.exe
ARCHI = ${DIR}/${EXE:.exe=.a}
MAIN = ${EXE:.exe=.f90}


F90SRC = sub_init_random_seed.o sub_read_parameter.o sub_read_hamiton.o sub_gradient.o sub_coupling.o sub_ele_motion_fun.o sub_ele_rk4.o sub_hop.o sub_nuc_motion_fun.o sub_nuc_rk4.o sub_read_ini_trj.o sub_read_parameter.o sub_vel_adj.o sub_write_trj.o

#
# Dependencies
#
${ARCHI}: ${FSRC} ${F90SRC}
	${AR} rv ${ARCHI} *.o

${EXE}: ${ARCHI} ${MAIN}
	${AR} rv ${ARCHI} *.o
	${F90C}  -o  $@ ${F90FLAGS} ${MAIN} ${LIB} ${ARCHI} ${LIB}  

install: ${EXE} 
	${RM} *.o

clean:
	${RM} *.o *.mod *genmod.f90

veryclean:
	${RM} ${EXE} *.o ${ARCHI}

#
# Suffix rules
#
.SUFFIXES: $(SUFFIXES) .f90
.f90.o:
	${F90C} -c ${F90FLAGS} $*.f90
