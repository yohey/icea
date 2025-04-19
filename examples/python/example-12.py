#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from icea import CEA_Problem


prob = CEA_Problem()

prob.set_problem(mode = 'rocket', name = 'Example 12', equilibrium = True, frozen = True, frozen_at_throat = True)
prob.set_output_options(SI = True, mass_fractions = True)

prob.add_reactant('fuel', 'CH6N2(L)', rho =  874)
prob.add_reactant('oxyd', 'N2O4(L)',  rho = 1431)

prob.set_chamber_pressures([1000], unit = 'legacy-psi')
prob.set_mixture_ratios([2.5])
prob.set_pressure_ratios([68.0457])
prob.set_supersonic_area_ratios([5, 10, 25, 50, 75, 100, 150, 200])

prob.set_only_species(['CO', 'CO2', 'H',   'HNO', 'HNO2',   'HO2', 'H2', 'H2O', 'H2O2',
                       'N',  'NO',  'NO2', 'N2',  'N2O',    'O',   'OH', 'O2',  'HCO',
                       'NH', 'CH4', 'NH2', 'NH3', 'H2O(L)', 'C(gr)'])

prob.set_legacy_mode(True)

prob.run(out_filename = 'example-12.out')
