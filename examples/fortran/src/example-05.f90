subroutine run_example_05()
  use cea
  implicit none

  type(CEA_Problem):: prob

  call prob%set_problem(mode = 'hp', name = 'Example 5')
  call prob%set_output_options(SI = .false.)

  call prob%add_reactant('name', 'NH4CLO4(I)',  ratio = 72.06d0, T = 298.15d0)
  call prob%add_reactant('name', 'CHOS-Binder', formula = 'C 1 H 1.86955 O .031256 S .008415', &
                                                ratio = 18.58d0, T = 298.15d0, h = -2999.082d0, h_unit = 'cal/mol')
  call prob%add_reactant('name', 'AL(cr)',      ratio =  9.00d0, T = 298.15d0)
  call prob%add_reactant('name', 'MgO(cr)',     ratio =  0.20d0, T = 298.15d0)
  call prob%add_reactant('name', 'H2O(L)',      ratio =  0.16d0, T = 298.15d0)

  call prob%set_chamber_pressures([500d0, 250d0, 125d0, 50d0, 5d0], unit = 'legacy-psi')

  call prob%set_omit_species([character(15):: &
       'COOH',  'C2', 'C2H', 'CHCO,ketyl',  'C2H2,vinylidene',  'CH2CO,ketene',  'C2H3,vinyl', &
       'CH3CO,acetyl',  'C2H4O,ethylen-o', 'CH3CHO,ethanal',  'CH3COOH',   '(HCOOH)2', &
       'C2H5',            'C2H6',            'CH3N2CH3',        'CH3OCH3', &
       'C2H5OH',          'CCN',             'CNC',             'C2N2', &
       'C2O',             'C3',              'C3H3,propargyl',  'C3H4,allene', &
       'C3H4,propyne',    'C3H4,cyclo-',     'C3H5,allyl',      'C3H6,propylene', &
       'C3H6,cyclo-',     'C3H6O',           'C3H7,n-propyl',   'C3H7,i-propyl', &
       'C3H8',            'C3H8O,1propanol', 'C3H8O,2propanol', 'C3O2', &
       'C4',              'C4H2',            'C4H4,1,3-cyclo-', 'C4H6,butadiene', &
       'C4H6,2-butyne',   'C4H6,cyclo-',     'C4H8,1-butene',   'C4H8,cis2-buten', &
       'C4H8,tr2-butene', 'C4H8,isobutene',  'C4H8,cyclo-',     '(CH3COOH)2', &
       'C4H9,n-butyl',    'C4H9,i-butyl',    'C4H9,s-butyl',    'C4H9,t-butyl', &
       'C4H10,isobutane', 'C4H10,n-butane',  'C4N2',            'C5', &
       'C5H6,1,3cyclo-',  'C5H8,cyclo-',     'C5H10,1-pentene', 'C5H10,cyclo-', &
       'C5H11,pentyl',    'C5H11,t-pentyl',  'C5H12,n-pentane', 'C5H12,i-pentane', &
       'CH3C(CH3)2CH3',   'C6H2',            'C6H5,phenyl',     'C6H5O,phenoxy', &
       'C6H6',            'C6H5OH,phenol',   'C6H10,cyclo-',    'C6H12,1-hexene', &
       'C6H12,cyclo-',    'C6H13,n-hexyl',    'C7H7,benzyl',    'C7H8', &
       'C7H8O,cresol-mx', 'C7H14,1-heptene', 'C7H15,n-heptyl',  'C7H16,n-heptane', &
       'C8H8,styrene',    'C8H10,ethylbenz', 'C8H16,1-octene',  'C8H17,n-octyl', &
       'C8H18,isooctane', 'C8H18,n-octane',  'C9H19,n-nonyl',   'C10H8,naphthale', &
       'C10H21,n-decyl',  'C12H9,o-bipheny', 'C12H10,biphenyl', 'Jet-A(g)', &
       'HNCO',   'HNO',  'HNO2',   'HNO3',   'HCCN',    'HCHO,formaldehy',  'HCOOH', &
       'NH',     'NH2',  'NH2OH',  'NCN',    'N2H2',  'NH2NO2',   'N2H4',  'H2O2', &
       '(HCOOH)2',   'C6H6(L)',  'C7H8(L)',  'C8H18(L),n-octa',  'Jet-A(L)',  'H2O(s)', 'H2O(L)'])

  call prob%set_legacy_mode(.true.)

  call prob%run(out_filename = 'example-05.out')

  return
end subroutine run_example_05
