#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from icea import CEA_Problem


prob = CEA_Problem()

prob.set_problem(mode = 'deton', name = 'Example 6')

prob.set_output_options(SI = False, transport = True)

prob.add_reactant('oxyd', 'O2', ratio = 100, T = 298.15)
prob.add_reactant('fuel', 'H2', ratio = 100, T = 298.15)

prob.set_chamber_temperatures([298.15, 500])
prob.set_chamber_pressures([1, 20], unit = 'bar')
prob.set_mixture_ratios([1.0], type_ = 'eq.ratio')

prob.set_legacy_mode(True)

prob.run(out_filename = 'example-06.out')
