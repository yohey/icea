module mod_cea
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
  integer, parameter:: ioinp =  7
  integer, parameter:: ioout =  8
  integer, parameter:: iosch = 13
  integer, parameter:: iothm = 14
  integer, parameter:: ioplt = 15
  integer, parameter:: iotrn = 18
  !***********************************************************************

  type:: CEA_Problem
     real(8):: Eqrat, Hsub0, Oxfl, Pp, R, Size, S0, Tln, Tm, Trace, Tt, Viscns, Vv
     real(8):: Atwt(maxEl), B0(maxEl), X(maxMat)
     real(8):: A(maxEl, maxNgc), G(maxMat, maxMat+1)

     character(2):: Elmt(maxEl), Ratom(maxR, 12)
     character(8):: Fox(maxR)
     character(10):: Thdate
     character(15):: Case, Energy(maxR), Omit(0:maxNgc), Pltvar(20), Prod(0:maxNgc), Rname(maxR)

     real(8):: Cpr(Ncol), Dlvpt(Ncol), Dlvtp(Ncol), Gammas(Ncol), Hsum(Ncol)
     real(8):: Ppp(Ncol), Ssum(Ncol), Totn(Ncol), Ttt(Ncol), Vlm(Ncol), Wm(Ncol)
     real(8):: Pltout(500, 20)

     integer:: Nreac
     integer:: Jray(maxR)
     real(8):: Dens(maxR), Enth(maxR), Pecwt(maxR), Rmw(maxR), Rtemp(maxR), Rnum(maxR, 12)

     real(8):: Cpsum
     real(8):: Cft(maxNc, 9), Coef(maxNg, 9, 3), Temp(2, maxNc)
     real(8):: Cp(maxNgc), H0(maxNgc), Mu(maxNgc), Mw(maxNgc), S(maxNgc), Tg(4)

     integer:: Iopt, Isup, Nfz, Npp, Nsub, Nsup
     logical:: Area, Debugf, Fac, Froz, Page1, Rkt
     real(8):: Acat, Awt, Cstr, Tcest, Ma
     real(8):: Aeat(Ncol), App(Ncol), Pcp(2*Ncol), Sonvel(Ncol), Spim(Ncol), Subar(13), Supar(13), Vmoc(Ncol)

     integer:: Nsk
     logical:: Incdeq, Incdfz, Refleq, Reflfz, Shkdbg
     real(8):: U1(Ncol), Mach1(Ncol), A1, Gamma1

     integer:: Nm, Nr, Ntape
     integer:: Ind(maxTr), Jcm(maxEl)
     real(8):: Cprr(maxTr), Con(maxTr), Wmol(maxTr), Xs(maxTr)
     real(8):: Eta(maxTr, maxTr), Stc(maxTr, maxTr)

     real(8):: Coneql(Ncol), Confro(Ncol), Cpeql(Ncol), Cpfro(Ncol), Preql(Ncol), Prfro(Ncol), Vis(Ncol)

     ! Information used in variable output format
     character(4):: fmt(30) = [character(4):: '(1X', ',A15', ',', 'F9.', '0,', 'F9.', '0,', 'F9.', '0,', &
                               'F9.', '0,', 'F9.', '0,', 'F9.', '0,', 'F9.', '0,', 'F9.', '0,', &
                               'F9.', '0,', 'F9.', '0,', 'F9.', '0,', 'F9.', '0,', 'F9.', '0', ')']
  end type CEA_Problem

  !***********************************************************************
  ! Fundamental constants from:  Cohen, E. Richard & Taylor, Barry N.,
  ! The 1986 codata recommended values of the fundamental physical
  ! constants, J. Phys. Chem. Ref. Data, vol. 17, No. 4, 1988, pp. 1795--1803.
  !***********************************************************************
  real(8), parameter:: R0 = 8314.51d0
  real(8), parameter:: pi = 3.14159265d0
  real(8), parameter:: Avgdr = 6.0221367d0
  real(8), parameter:: Boltz = 1.380658d0

  ! Atomic Symbols
  character(2), parameter:: atomic_symbol(100) = ['H ', 'D ', 'HE', 'LI', 'BE', 'B ', 'C ', 'N ', 'O ', 'F ', &
                                                  'NE', 'NA', 'MG', 'AL', 'SI', 'P ', 'S ', 'CL', 'AR', 'K ', 'CA', 'SC', &
                                                  'TI', 'V ', 'CR', 'MN', 'FE', 'CO', 'NI', 'CU', 'ZN', 'GA', 'GE', 'AS', &
                                                  'SE', 'BR', 'KR', 'RB', 'SR', 'Y ', 'ZR', 'NB', 'MO', 'TC', 'RU', 'RH', &
                                                  'PD', 'AG', 'CD', 'IN', 'SN', 'SB', 'TE', 'I ', 'XE', 'CS', 'BA', 'LA', &
                                                  'CE', 'PR', 'ND', 'PM', 'SM', 'EU', 'GD', 'TB', 'DY', 'HO', 'ER', 'TM', &
                                                  'YB', 'LU', 'HF', 'TA', 'W ', 'RE', 'OS', 'IR', 'PT', 'AU', 'HG', 'TL', &
                                                  'PB', 'BI', 'PO', 'AT', 'RN', 'FR', 'RA', 'AC', 'TH', 'PA', 'U ', 'NP', &
                                                  'PU', 'AM', 'CM', 'BK', 'CF', 'ES']

  ! Atomic Maxx - Coplen, T. B., Atomic Weights of the Elements 1999.
  !               J. Phys. Chem. Ref. Data, vol. 30, No. 3, 2001, pp. 701--712.
  real(8), parameter:: atomic_mass(100) = [  1.007940d0,    2.014102d0,   4.002602d0,   6.941000d0,   9.0121820d0, &
                                            10.811000d0,   12.010700d0,  14.006700d0,  15.999400d0,  18.9984032d0, &
                                            20.179700d0,   22.989770d0,  24.305000d0,  26.981538d0,  28.0855000d0, &
                                            30.973761d0,   32.065000d0,  35.453000d0,  39.948000d0,  39.0983000d0, &
                                            40.078000d0,   44.955910d0,  47.867000d0, &
                                            50.941500d0,   51.996100d0,  54.938049d0,  55.845000d0,  58.9332000d0,  58.693400d0, &
                                            63.546000d0,   65.390000d0,  69.723000d0, &
                                            72.640000d0,   74.921600d0,  78.960000d0,  79.904000d0, &
                                            83.800000d0,   85.467800d0,  87.620000d0,  88.905850d0, &
                                            91.224000d0,   92.906380d0,  95.940000d0,  97.907200d0, &
                                           101.070000d0,  102.905500d0, 106.420000d0, 107.868200d0, &
                                           112.411000d0,  114.818000d0, 118.710000d0, &
                                           121.760000d0,  127.600000d0, 126.904470d0, &
                                           131.293000d0,  132.905450d0, 137.327000d0, 138.905500d0, &
                                           140.116000d0,  140.907650d0, 144.912700d0, 145.000000d0, &
                                           150.360000d0,  151.964000d0, 157.250000d0, 158.925340d0, &
                                           162.500000d0,  164.930320d0, 167.259000d0, 168.934210d0, &
                                           173.040000d0,  174.967000d0, 178.490000d0, &
                                           180.947900d0,  183.840000d0, 186.207000d0, &
                                           190.230000d0,  192.217000d0, 195.078000d0, 196.966550d0, &
                                           200.590000d0,  204.383300d0, 207.200000d0, 208.980380d0, 208.9824000d0, 209.987100d0, &
                                           222.017600d0,  223.019700d0, 226.025400d0, 227.027800d0, &
                                           232.038100d0,  231.035880d0, 238.028910d0, 237.048200d0, &
                                           244.064200d0,  243.061400d0, 247.070300d0, 247.070300d0, &
                                           251.058700d0,  252.083000d0]

  ! Atomic Valences
  integer, parameter:: atomic_valence(100) = [1, 1, 0, 1, 2, 3, 4, 0, -2, -1, 0, 1, 2, 3, 4, 5, &
                                              4, -1, 0, 1, 2, 3, 4, 5, 3, 2, 3, 2, 2, 2, 2, 3, 4, 3, 4, &
                                              -1, 0, 1, 2, 3, 4, 5, 6, 7, 3, 3, 2, 1, 2, 3, 4, 3, 4, &
                                              -1, 0, 1, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, &
                                              4, 5, 6, 7, 4, 4, 4, 3, 2, 1, 2, 3, 2, -1, 0, 1, 2, 3, 4, &
                                              5, 6, 5, 4, 3, 3, 3, 3, 3]

end module mod_cea
