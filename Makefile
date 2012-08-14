
FC = gfortran
FFLAGS = 

RUN = cea2
SRC = cea2.f
INC = cea.inc


$(RUN): $(SRC) $(INC)
	$(FC) $(FFLAGS) -o $@ $(SRC)


thermo.lib: $(RUN) thermo.inp
	echo thermo | ./$(RUN)


trans.lib: $(RUN) trans.inp
	echo trans | ./$(RUN)


clean:
	$(RM) $(RUN)
	$(RM) thermo.lib thermo.out
	$(RM) trans.lib trans.out
	$(RM) cea2.out cea2.plt
