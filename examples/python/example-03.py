#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from icea import CEA_Problem


prob = CEA_Problem()

prob.set_problem(mode = 'hp', name = 'Example 3')

prob.set_output_options(SI = True, trace_tol = 1e-15)

prob.add_reactant('oxyd', 'Air',             ratio = 1.0, T = 700.00)
prob.add_reactant('fuel', 'C7H8(L)',         ratio = 0.4, T = 298.15)
prob.add_reactant('fuel', 'C8H18(L),n-octa', ratio = 0.6, T = 298.15)

prob.set_chamber_pressures([100, 10, 1], unit = 'bar')

prob.set_mixture_ratios([17.0])

prob.set_omit_species(['CCN',             'CNC',             'C2N2',            'C2O',
                       'C3H4,allene',     'C3H4,propyne',    'C3H4,cyclo-',     'C3',
                       'C3H5,allyl',      'C3H6,propylene',  'C3H6,cyclo-',     'C3H3,propargyl',
                       'C3H6O',           'C3H7,n-propyl',   'C3H7,i-propyl',   'Jet-A(g)',
                       'C3O2',            'C4',              'C4H2',            'C3H8O,2propanol',
                       'C4H4,1,3-cyclo-', 'C4H6,butadiene',  'C4H6,2-butyne',   'C3H8O,1propanol',
                       'C4H8,tr2-butene', 'C4H8,isobutene',  'C4H8,cyclo-',     'C4H6,cyclo-',
                       '(CH3COOH)2',      'C4H9,n-butyl',    'C4H9,i-butyl',    'C4H8,1-butene',
                       'C4H9,s-butyl',    'C4H9,t-butyl',    'C4H10,isobutane', 'C4H8,cis2-buten',
                       'C4H10,n-butane',  'C4N2',            'C5',              'C3H8',
                       'C5H6,1,3cyclo-',  'C5H8,cyclo-',     'C5H10,1-pentene', 'C10H21,n-decyl',
                       'C5H10,cyclo-',    'C5H11,pentyl',    'C5H11,t-pentyl',  'C12H10,biphenyl',
                       'C5H12,n-pentane', 'C5H12,i-pentane', 'CH3C(CH3)2CH3',   'C12H9,o-bipheny',
                       'C6H6',            'C6H5OH,phenol',   'C6H10,cyclo-',    'C6H2',
                       'C6H12,1-hexene',  'C6H12,cyclo-',    'C6H13,n-hexyl',   'C6H5,phenyl',
                       'C7H7,benzyl',     'C7H8',            'C7H8O,cresol-mx', 'C6H5O,phenoxy',
                       'C7H14,1-heptene', 'C7H15,n-heptyl',  'C7H16,n-heptane', 'C10H8,azulene',
                       'C8H8,styrene',    'C8H10,ethylbenz', 'C8H16,1-octene',  'C10H8,napthlene',
                       'C8H17,n-octyl',   'C8H18,isooctane', 'C8H18,n-octane',  'C9H19,n-nonyl',
                       'Jet-A(L)',        'C6H6(L)',         'H2O(s)',          'H2O(L)'])

prob.set_legacy_mode(True)

prob.run(out_filename = 'example-03.out')
