#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from icea import CEA_Problem


prob = CEA_Problem()

prob.set_problem(mode = 'tv', name = 'Example 2')

prob.set_output_options(SI = False, transport = True)

prob.add_reactant('fuel', 'H2',  ratio = 100)
prob.add_reactant('oxyd', 'Air', ratio = 100)

prob.set_chamber_temperatures([3000])
prob.set_chamber_densities([9.1864e-5, 8.0877e-6, 6.6054e-7], unit = 'g/cc')
prob.set_mixture_ratios([1.0], type_ = 'eq.ratio')

prob.set_only_species(["Ar",   "C",    "CO",   "CO2",  "H",    "H2",   "H2O",  "HNO",  "HO2",  "HNO2", "HNO3",
                       "N",    "NH",   "NO",   "N2",   "N2O3", "O",    "O2",   "OH",   "O3"])

prob.set_legacy_mode(True)

prob.run(out_filename = 'example-02.out')
