all: rbc.out
	./rbc.out
#	gprof rbc.out gmon.out > profile.dat

rbc.out: rbc.f90 solve.f90 tools.f90 basic.f90 def_types.f90
	gfortran rbc.f90 solve.f90 basic.f90 def_types.f90 tools.f90 -llapack -lblas -p -g -o rbc.out
# 	gfortran lapack_prb.f90 -llapack -lblas -o lapack_prb
	#gfortran lapack_prb.o  -L$HOME/libf77/$ARCH -llapack -lblas
    
clean:
	rm rbc.out profile.dat
