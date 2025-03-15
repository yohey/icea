#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pycea import CEA_Problem


prob = CEA_Problem()

prob.set_problem(mode = 'rocket', name = 'Example 8', equilibrium = True)

prob.set_output_options(SI = True)

prob.add_reactant('fuel', 'H2(L)', ratio = 100, T = 20.27)
prob.add_reactant('oxyd', 'O2(L)', ratio = 100, T = 90.17)

prob.set_chamber_pressures([53.3172])
prob.set_mixture_ratios([5.55157])
prob.set_pressure_ratios([10, 100, 1000])
prob.set_subsonic_area_ratios([1.58])
prob.set_supersonic_area_ratios([25, 50, 75])

prob.set_legacy_mode(True)

prob.run(out_filename = 'example-08.out')
