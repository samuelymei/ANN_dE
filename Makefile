.SUFFIXES: 
.SUFFIXES: .f90 .o

FC = ifort
FFLAGS = -CB
INCLUDE = 
LIBS = 

EXE = ANN_dE.x

MODULES = precision_m.mod random_m.mod ann_m.mod atom_m.mod molecule_m.mod

OBJS = precision_m.o lib.o random.o ann.o atom.o molecule.o activationfunc.o main.o gdmin.o dfpmin.o lnsrch.o

all:	${EXE}

$(EXE):$(OBJS) ${MODULES}
	$(FC) -o $@ $(FFLAGS) $(OBJS) $(LIBS)

%.o %.mod:%.f90
	$(FC) -c $(FFLAGS) $(INCLUDE) $<

include .depend

depend .depend:
	makedepf90 *.f90 > .depend

clean:
	/bin/rm -f $(EXE) $(OBJS) ${MODULES}

