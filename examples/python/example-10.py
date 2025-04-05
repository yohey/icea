#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pycea import CEA_Problem


prob = CEA_Problem()

prob.set_problem(mode = 'rocket', name = 'Example 10', equilibrium = False)
prob.set_output_options(SI = True, short = True)

prob.add_reactant('fuel', 'H2(L)', ratio = 100, T = 20.27)
prob.add_reactant('oxyd', 'O2(L)', ratio = 100, T = 90.17)

prob.set_finite_area_combustor(mass_flow_ratio = 1333.9)
prob.set_chamber_pressures([5.33172])
prob.set_mixture_ratios([5.55157])
prob.set_pressure_ratios([10, 100, 1000])
prob.set_supersonic_area_ratios([25, 50, 75])

prob.set_legacy_mode(True)

prob.run(out_filename = 'example-10.out')
