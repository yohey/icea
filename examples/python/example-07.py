#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from icea import CEA_Problem


prob = CEA_Problem()

prob.set_problem(mode = 'shock', name = 'Example 7', mole_ratios = True, equilibrium = True, frozen = True, incident = True)

prob.add_reactant('name', 'H2', ratio = 0.05, T = 300)
prob.add_reactant('name', 'O2', ratio = 0.05, T = 300)
prob.add_reactant('name', 'Ar', ratio = 0.90, T = 300)

prob.set_chamber_pressures([10, 20], unit = 'mmHg')

prob.set_initial_velocities([1000, 1100, 1200, 1250, 1300, 1350, 1400])

prob.set_legacy_mode(True)

prob.run(out_filename = 'example-07.out')
