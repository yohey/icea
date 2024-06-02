!***********************************************************************
!                            cea.include
!***********************************************************************
!
!  The following parameters set the maximum dimensions for many variables
!    They are defined in part 2, page 39 of the manuals, NASA RP-1311.
!    The variable Ncol set the number of columns in the output.  It may
!    be increased for wider paper or smaller fonts.
!
module mod_legacy_cea
  use mod_cea
  implicit none

  real(8):: Enn, Ennl, Enlsav, Ensave, Sumn
  real(8):: Deln(maxNgc), Enln(maxNgc), Sln(maxNgc)
  real(8):: En(maxNgc, Ncol)

  ! for subroutine CPHS and ALLCON
  real(8):: cx(7) = [0d0, 0d0, 1d0, 0.5d0, 0.6666666666666667d0, 0.75d0, 0.8d0]
  real(8):: hcx(7) = [0d0, 0d0, 1d0, 0d0, 0d0, 0d0, 0d0]
  real(8):: scx(7)

end module mod_legacy_cea
