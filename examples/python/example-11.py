#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pycea import CEA_Problem


prob = CEA_Problem()

prob.set_problem(mode = 'rocket', name = 'Example 11', equilibrium = True, mole_ratios = True, ions = True)
prob.set_output_options(SI = True, transport = True)

prob.add_reactant('fuel', 'Li(cr)', ratio = 1.0000, T = 298.15)
prob.add_reactant('oxyd', 'F2(L)',  ratio = 0.5556, T =  85.02)

prob.set_chamber_pressures([1000], unit = 'legacy-psi')
prob.set_pressure_ratios([68.0457])
prob.set_subsonic_area_ratios([10])
prob.set_supersonic_area_ratios([10, 20, 100])

prob.set_legacy_mode(True)

prob.run(out_filename = 'example-11.out')
