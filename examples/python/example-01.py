#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from icea import CEA_Problem


prob = CEA_Problem()

prob.set_problem(mode = 'tp', name = 'Example 1', mole_ratios = True)

prob.set_output_options(SI = False)

prob.add_reactant('fuel', 'H2',  ratio = 1.0)
prob.add_reactant('oxyd', 'Air', ratio = 1.0)

prob.set_chamber_pressures([1.0, 0.1, 0.01], unit = 'atm')
prob.set_chamber_temperatures([3000, 2000])
prob.set_mixture_ratios([1.0, 1.5], type_ = 'eq.ratio')

prob.set_only_species(["Ar",   "C",    "CO",   "CO2",  "H",    "H2",   "H2O",  "HNO",  "HO2",  "HNO2", "HNO3",
                       "N",    "NH",   "NO",   "N2",   "N2O3", "O",    "O2",   "OH",   "O3"])

prob.set_legacy_mode(True)

prob.run(out_filename = 'example-01.out')
