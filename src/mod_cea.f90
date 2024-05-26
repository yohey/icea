module mod_cea
  implicit none

  type:: CEA_Problem
     ! INFORMATION USED IN VARIABLE OUTPUT FORMAT
     character(4):: Fmt(30) = [character(4):: '(1X', ',A15', ',', 'F9.', '0,', 'F9.', '0,', 'F9.', '0,', &
                               'F9.', '0,', 'F9.', '0,', 'F9.', '0,', 'F9.', '0,', 'F9.', '0,', &
                               'F9.', '0,', 'F9.', '0,', 'F9.', '0,', 'F9.', '0,', 'F9.', '0', ')']
  end type CEA_Problem

end module mod_cea
