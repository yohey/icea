module mod_cea
  implicit none

  type:: CEA_Problem
     ! Information used in variable output format
     character(4):: Fmt(30) = [character(4):: '(1X', ',A15', ',', 'F9.', '0,', 'F9.', '0,', 'F9.', '0,', &
                               'F9.', '0,', 'F9.', '0,', 'F9.', '0,', 'F9.', '0,', 'F9.', '0,', &
                               'F9.', '0,', 'F9.', '0,', 'F9.', '0,', 'F9.', '0,', 'F9.', '0', ')']
  end type CEA_Problem

  ! Atomic Valences
  integer, parameter:: Valence(100) = [1, 1, 0, 1, 2, 3, 4, 0, -2, -1, 0, 1, 2, 3, 4, 5, &
                                       4, -1, 0, 1, 2, 3, 4, 5, 3, 2, 3, 2, 2, 2, 2, 3, 4, 3, 4, &
                                       -1, 0, 1, 2, 3, 4, 5, 6, 7, 3, 3, 2, 1, 2, 3, 4, 3, 4, &
                                       -1, 0, 1, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, &
                                       4, 5, 6, 7, 4, 4, 4, 3, 2, 1, 2, 3, 2, -1, 0, 1, 2, 3, 4, &
                                       5, 6, 5, 4, 3, 3, 3, 3, 3]

end module mod_cea
