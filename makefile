FC=gfortran

FFLAGS=-02 
SOURCE=modified_biharmonic.f
TARGET=modified_biharmonic

OBJS=\
d1mach.o\
dasum.o\
daxpy.o\
dcfft.o\
dcopy.o\
ddot.o\
dgeco.o\
dgefa.o\
dgmres.o\
dhels.o\
dheqr.o\
dnrm2.o\
dorth.o\
dpigmr.o\
drlcal.o\
dscal.o\
dxlcal.o\
idamax.o\
isdgmr.o\
matplot.o\
prini.o\
quad2.o\
ribesl.f\
rkbesl.f

modified_biharomonic: ${OBJS} modified_biharmonic.f
	${FC} ${FFLAGS} modified_biharmonic.f -o ${TARGET} \
	  ${OBJS}

d1mach.o: d1mach.f
	${FC} ${FFLAGS} -c d1mach.f

dasum.o: dasum.f
	${FC} ${FFLAGS} -c dasum.f

daxpy.o: daxpy.f
	${FC} ${FFLAGS} -c daxpy.f

dcfft.o: dcfft.f
	${FC} ${FFLAGS} -c dcfft.f

dcopy.o: dcopy.f
	${FC} ${FFLAGS} -c dcopy.f

ddot.o: ddot.f
	${FC} ${FFLAGS} -c ddot.f

dgeco.o: dgeco.f
	${FC} ${FFLAGS} -c dgeco.f

dgefa.o: dgefa.f
	${FC} ${FFLAGS} -c dgefa.f

dgmres.o: dgmres.f
	${FC} ${FFLAGS} -c dgmres.f

dhels.o: dhels.f
	${FC} ${FFLAGS} -c dhels.f

dheqr.o: dheqr.f
	${FC} ${FFLAGS} -c dheqr.f

dnrm2.o: dnrm2.f
	${FC} ${FFLAGS} -c dnrm2.f

dorth.o: dorth.f
	${FC} ${FFLAGS} -c dorth.f

dpigmr.o: dpigmr.f
	${FC} ${FFLAGS} -c dpigmr.f

drlcal.o: drlcal.f
	${FC} ${FFLAGS} -c drlcal.f

dscal.o: dscal.f
	${FC} ${FFLAGS} -c dscal.f

dxlcal.o: dxlcal.f
	${FC} ${FFLAGS} -c dxlcal.f

idamax.o: idamax.f
	${FC} ${FFLAGS} -c idamax.f

isdgmr.o: isdgmr.f
	${FC} ${FFLAGS} -c isdgmr.f

matplot.o: matplot.f
	${FC} ${FFLAGS} -c matplot.f

prini.o: prini.f
	${FC} ${FFLAGS} -c prini.f

quad2.o: quad2.f
	${FC} ${FFLAGS} -c quad2.f

ribesl.o: ribesl.f
	${FC} ${FFLAGS} -c ribesl.f

rkbesl.o: ribesl.f
	${FC} ${FFLAGS} -c ribesl.f

clean:
	rm *.o
	rm modified_biharmonic
