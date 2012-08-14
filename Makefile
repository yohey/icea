
FC = gfortran
FFLAGS = 

RUN = cea2
SRC = cea2.f
INC = cea.inc


$(RUN): $(SRC) $(INC)
	$(FC) $(FFLAGS) -o $@ $(SRC)


check: check-thermo check-trans check-cea2
	@echo
	@echo 'Congratulations! All tests passed.'

check-thermo: $(RUN) thermo.inp
	echo thermo | ./$(RUN)
	diff test/thermo.lib thermo.lib
	diff test/thermo.out thermo.out

check-trans: $(RUN) trans.inp
	echo trans | ./$(RUN)
	diff test/trans.lib trans.lib
	diff test/trans.out trans.out

check-cea2: $(RUN) cea2.inp
	echo cea2 | ./$(RUN)
	diff test/cea2.out cea2.out
	diff test/cea2.plt cea2.plt


clean:
	$(RM) $(RUN)
	$(RM) thermo.lib thermo.out
	$(RM) trans.lib trans.out
	$(RM) cea2.out cea2.plt
