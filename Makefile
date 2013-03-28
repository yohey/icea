
FC = gfortran
FFLAGS = 
DFLAGS = -O0 -std=f2003 -fbounds-check -Wall -ffpe-trap=invalid,zero,overflow -g -fbacktrace

RUN = cea2
SRC = cea.f90 cea2.f90


.PHONY: test clean


$(RUN): $(SRC)
	$(FC) $(FFLAGS) -o $@ $^


debug: $(SRC)
	$(FC) $(DFLAGS) -o $@ $^


test: test-thermo test-trans test-cea2
	@echo
	@echo 'Congratulations! All tests passed.'

test-thermo: $(RUN) thermo.inp
	echo thermo | ./$(RUN)
	diff test/thermo.lib thermo.lib
	diff test/thermo.out thermo.out

test-trans: $(RUN) trans.inp
	echo trans | ./$(RUN)
	diff test/trans.lib trans.lib
	diff test/trans.out trans.out

test-cea2: $(RUN) cea2.inp
	echo cea2 | ./$(RUN)
	diff test/cea2.out cea2.out
	diff test/cea2.plt cea2.plt


clean:
	$(RM) $(RUN)
	$(RM) debug
	$(RM) cea.mod
	$(RM) thermo.lib thermo.out
	$(RM) trans.lib trans.out
	$(RM) cea2.out cea2.plt
