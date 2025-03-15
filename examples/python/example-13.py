#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pycea import CEA_Problem


prob = CEA_Problem()

prob.set_problem(mode = 'rocket', name = 'Example 13', equilibrium = True)
prob.set_output_options(SI = False, trace_tol = 1e-10)

prob.add_reactant('fuel', 'N2H4(L)', ratio = 0.8, T = 298.15, rho = 0.874)
prob.add_reactant('fuel', 'Be(a)',   ratio = 0.2, T = 298.15, rho = 1.431)
prob.add_reactant('oxyd', 'H2O2(L)', ratio = 1.0, T = 298.15)

prob.set_chamber_pressures([3000], unit = 'legacy-psi')
prob.set_mixture_ratios([67], type_ = '%fuel')
prob.set_pressure_ratios([3, 10, 30, 300])

prob.set_insert_species(['BeO(L)'])

prob.set_legacy_mode(True)

prob.run(out_filename = 'example-13.out')
