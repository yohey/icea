module mod_types
  use mod_constants
  implicit none

  integer, parameter:: maxNgc = 600
  integer, parameter:: maxNc  = 300
  integer, parameter:: Ncol   = 8
  integer, parameter:: maxMat = 50
  integer, parameter:: maxTr  = 40
  integer, parameter:: maxR   = 24
  integer, parameter:: maxEl  = 20
  integer, parameter:: maxNg  = 500
  integer, parameter:: maxMix = 52
  integer, parameter:: maxT   = 51
  integer, parameter:: maxPv  = 26

  ! The following parameters set the input/output unit numbers.
  ! These numbers are also defined in the manual, part 2 p. 39,
  ! and may be adjusted as desired.
  integer:: IOINP
  integer:: IOOUT
  !***********************************************************************

  type:: TransportProperty
     integer:: j(2)
     real(8):: data(6, 3, 2)
  end type TransportProperty


  type:: CEA_Point
     logical:: Debug
     real(8):: En(maxNgc)
     real(8), pointer:: B0(:)
     real(8), pointer:: B0p(:, :)
     real(8):: Cpr, Dlvpt, Dlvtp, Gammas, Hsum
     real(8):: Ppp, Ssum, Totn, Ttt, Vlm, Wm
     real(8):: Sonvel, Spim, Vmoc
     ! ROCKET
     real(8):: AeAt, App
     ! SHCK
     real(8):: U1, Mach1
     ! TRANP
     real(8):: Coneql, Confro, Cpeql, Cpfro, Preql, Prfro, Vis
  end type CEA_Point


  type:: CEA_Problem
     integer:: io_log = 0

     logical:: invalid_case = .false.

     character(MAX_FILENAME):: filename_trans_lib = 'trans.lib'
     character(MAX_FILENAME):: filename_thermo_lib = 'thermo.lib'

     type(TransportProperty), allocatable:: transport_properties(:)
     type(CEA_Point), pointer:: points(:, :)

     integer:: iOF !< index of O/F
     integer:: ipt !< index of point (new)
     integer:: max_points

     character(15):: ensert(20)

     ! for subroutine SETEN
     integer:: Isv, lsave, lsav
     real(8):: Enlsav, Ensave, Tsave

     real(8):: Enn, Ennl, Sumn
     real(8):: Deln(maxNgc), Enln(maxNgc), Sln(maxNgc)

     integer:: Iplt, Nc, Ng, Ngp1, Nlm, Nplt, Nof, Nomit, Nonly, Np, Npr, Ngc, Nsert, Nspr, Nspx, Nt
     integer:: Jcond(45), Jx(maxEl), Nfla(maxR), Ifz(maxNc)

     real(8):: Cpmix, Wmix, Bcheck
     real(8):: Am(2), Hpp(2), Vmin(2), Vpls(2), Wp(2), Oxf(maxMix), P(maxPv), Rh(2), T(maxT), V(maxPv)

     integer:: Imat, Iq1, Jliq, Jsol, Msing

     logical:: Convg, Detdbg, Detn, Eql, Gonly, Hp, Ions, Massf, Moles
     logical:: Pderiv, Shock, Short, SIunit, Sp, Tp, Trnspt, Vol

     real(8):: Eqrat, Hsub0, Oxfl, Pp, R, Size, S0, Tln, Tm, Trace, Tt, Viscns, Vv
     real(8):: Atwt(maxEl), X(maxMat)
     real(8):: A(maxEl, maxNgc), G(maxMat, maxMat+1)

     character(2):: Elmt(maxEl), Ratom(maxR, 12)
     character(8):: Fox(maxR)
     character(10):: Thdate
     character(15):: Case, Energy(maxR), Omit(0:maxNgc), Pltvar(20), Prod(0:maxNgc), Rname(maxR)

     real(8):: Pltout(500, 20)

     integer:: Nreac
     integer:: Jray(maxR)
     real(8):: Dens(maxR), Enth(maxR), Pecwt(maxR), Rmw(maxR), Rtemp(maxR), Rnum(maxR, 12)

     real(8):: Cpsum
     real(8):: Cft(9, maxNc), Coef(9, maxNg, 3), Temp(2, maxNc)
     real(8):: Cp(maxNgc), H0(maxNgc), Mu(maxNgc), Mw(maxNgc), S(maxNgc), Tg(4)

     integer:: Iopt, Isup, Nfz, Npp, Nsub, Nsup
     logical:: Area, Debugf, Fac, Froz, Page1, Rkt
     real(8):: Acat, Awt, Cstr, Tcest, Ma
     real(8):: Pcp(2*Ncol), Subar(13), Supar(13)

     integer:: Nsk
     logical:: Incdeq, Incdfz, Refleq, Reflfz, Shkdbg
     real(8):: A1, Gamma1

     integer:: Nm, Nr, Ntape
     integer:: Ind(maxTr), Jcm(maxEl)
     real(8):: Cprr(maxTr), Con(maxTr), Wmol(maxTr), Xs(maxTr)
     real(8):: Eta(maxTr, maxTr), Stc(maxTr, maxTr)

     ! Information used in variable output format
     character(4):: fmt(30) = [character(4):: '(1X', ',A15', ',', 'F9.', '0,', 'F9.', '0,', 'F9.', '0,', &
                               'F9.', '0,', 'F9.', '0,', 'F9.', '0,', 'F9.', '0,', 'F9.', '0,', &
                               'F9.', '0,', 'F9.', '0,', 'F9.', '0,', 'F9.', '0,', 'F9.', '0', ')']
  end type CEA_Problem

contains

  subroutine init_case(cea)
    implicit none

    type(CEA_Problem), intent(inout):: cea

    integer:: i, j

    ! temporary
    cea%max_points = maxT * maxPv

    allocate(cea%points(maxMix, cea%max_points))

    do concurrent (i = 1:maxMix, j = 1:cea%max_points)
       cea%points(i, j)%Debug = .false.
       cea%points(i, j)%En(:) = 0
       cea%points(i, j)%U1 = 0
       cea%points(i, j)%Mach1 = 0
    end do

    do concurrent (i = 1:maxMix)
       allocate(cea%points(i, 1)%B0(maxEl))
       allocate(cea%points(i, 1)%B0p(maxEl, 2))

       do concurrent (j = 2:cea%max_points)
          cea%points(i, j)%B0 => cea%points(i, 1)%B0
          cea%points(i, j)%B0p => cea%points(i, 1)%B0p
       end do
    end do

    cea%ensert(:) = ""

    cea%iOF = 0

    return
  end subroutine init_case

end module mod_types
