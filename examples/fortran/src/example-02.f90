subroutine run_example_02()
  use cea
  implicit none

  type(CEA_Problem):: prob

  call prob%set_problem(mode = 'tv', name = 'Example 2')
  call prob%set_output_options(SI = .false., transport = .true.)

  call prob%add_reactant('fuel', 'H2',  ratio = 100d0)
  call prob%add_reactant('oxyd', 'Air', ratio = 100d0)

  call prob%set_chamber_temperatures([3000d0])
  call prob%set_chamber_densities([9.1864d-05, 8.0877d-06, 6.6054d-07], unit = 'g/cc')
  call prob%set_mixture_ratios([1.0d0], type = 'eq.ratio')

  call prob%set_only_species(['Ar  ', 'C   ', 'CO  ', 'CO2 ', 'H   ', 'H2  ', 'H2O ', 'HNO ', 'HO2 ', 'HNO2', 'HNO3', &
                              'N   ', 'NH  ', 'NO  ', 'N2  ', 'N2O3', 'O   ', 'O2  ', 'OH  ', 'O3  '])

  call prob%set_legacy_mode(.true.)

  call prob%run(out_filename = 'example-02.out')

  return
end subroutine run_example_02
