#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from icea import CEA_Problem


prob = CEA_Problem()

prob.set_problem(mode = 'tp', name = 'Example 14', mole_ratios = True)

prob.set_output_options(SI = True, debug_points = [4])

prob.add_reactant('name', 'H2(L)', 'H 2', ratio = 100, rho =  874)
prob.add_reactant('name', 'O2(L)', 'O 2', ratio =  60, rho = 1431)

prob.set_chamber_pressures([0.05], unit = 'atm')
prob.set_chamber_temperatures([1000, 500, 350, 305, 304.3, 304.2, 304, 300])

prob.set_legacy_mode(True)

prob.run(out_filename = 'example-14.out')
