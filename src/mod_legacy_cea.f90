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

!  The following parameters set the input/output unit numbers.  These
!    numbers are also defined in the manual, part2 p39, and may be 
!    adjusted as desired.

  integer, parameter:: ioinp =  7
  integer, parameter:: ioout =  8
  integer, parameter:: iosch = 13
  integer, parameter:: iothm = 14
  integer, parameter:: ioplt = 15
  integer, parameter:: iotrn = 18
!***********************************************************************
  real(8):: Enn, Ennl, Enlsav, Ensave, Sumn
  real(8):: Deln(maxNgc), Enln(maxNgc), Sln(maxNgc)
  real(8):: En(maxNgc, Ncol)

  integer:: Ip, Iplt, It, Nc, Ng, Ngp1, Nlm, Nplt, Nof, Nomit, Nonly, Np, Npr, Npt, &
       Ngc, Nsert, Nspr, Nspx, Nt
  integer:: Jcond(45), Jx(maxEl), Nfla(maxR), Ifz(maxNc)

  real(8):: Cpmix, Wmix, Bcheck
  real(8):: Am(2), Hpp(2), Vmin(2), Vpls(2), Wp(2), Oxf(maxMix), &
       P(maxPv), Rh(2), T(maxT), V(maxPv)
  real(8):: B0p(maxEl, 2)

  integer:: Imat, Iq1, Isv, Jliq, Jsol, Lsave, Msing

  logical:: Convg, Debug(Ncol), Detdbg, Detn, Eql, Gonly, Hp, Ions, Massf, &
       Moles, Newr, Pderiv, Shock, Short, Siunit, Sp, Tp, Trnspt, Vol

  real(8):: Eqrat, Hsub0, Oxfl, Pp, R, Size, S0, Tln, Tm, &
       Trace, Tt, Viscns, Vv
  real(8):: Atwt(maxEl), B0(maxEl), X(maxMat)
  real(8):: A(maxEl, maxNgc), G(maxMat, maxMat+1)

  character(2):: Elmt(maxEl), Ratom(maxR, 12)
  character(8):: Fox(maxR)
  character(10):: Thdate
  character(15):: Case, Energy(maxR), Omit(0:maxNgc), Pltvar(20), &
       Prod(0:maxNgc), Rname(maxR)
  character(200):: Pfile

  real(8):: Cpr(Ncol), Dlvpt(Ncol), Dlvtp(Ncol), Gammas(Ncol), Hsum(Ncol), &
       Ppp(Ncol), Ssum(Ncol), Totn(Ncol), Ttt(Ncol), Vlm(Ncol), &
       Wm(Ncol)
  real(8):: Pltout(500, 20)

  integer:: Nreac
  integer:: Jray(maxR)
  real(8):: Dens(maxR), Enth(maxR), Pecwt(maxR), Rmw(maxR), Rtemp(maxR)
  real(8):: Rnum(maxR, 12)

  real(8):: Cpsum
  real(8):: Cft(maxNc, 9), Coef(maxNg, 9, 3), Temp(2, maxNc)
  real(8):: Cp(maxNgc), H0(maxNgc), Mu(maxNgc), Mw(maxNgc), S(maxNgc), Tg(4)

  integer:: Iopt, Isup, Nfz, Npp, Nsub, Nsup
  logical:: Area, Debugf, Fac, Froz, Page1, Rkt
  real(8):: Acat, Awt, Cstr, Tcest, Ma
  real(8):: Aeat(Ncol), App(Ncol), Pcp(2*Ncol), Sonvel(Ncol), Spim(Ncol), &
       Subar(13), Supar(13), Vmoc(Ncol)

  integer:: Nsk
  logical:: Incdeq, Incdfz, Refleq, Reflfz, Shkdbg
  real(8):: U1(Ncol), Mach1(Ncol), A1, Gamma1

  integer:: Nm, Nr, Ntape
  integer:: Ind(maxTr), Jcm(maxEl)
  real(8):: Cprr(maxTr), Con(maxTr), Wmol(maxTr), Xs(maxTr)
  real(8):: Eta(maxTr, maxTr), Stc(maxTr, maxTr)

  real(8):: Coneql(Ncol), Confro(Ncol), Cpeql(Ncol), Cpfro(Ncol), &
       Preql(Ncol), Prfro(Ncol), Vis(Ncol)

  ! for subroutine CPHS and ALLCON
  real(8):: cx(7) = [0d0, 0d0, 1d0, 0.5d0, 0.6666666666666667d0, 0.75d0, 0.8d0]
  real(8):: hcx(7) = [0d0, 0d0, 1d0, 0d0, 0d0, 0d0, 0d0]
  real(8):: scx(7)


!***********************************************************************
! FUNDAMENTAL CONSTANTS FROM:  COHEN, E.RICHARD & TAYLOR, BARRY N.,
! THE 1986 CODATA RECOMMENDED VALUES OF THE FUNDAMENTAL PHYSICAL
! CONSTANTS, J.PHYS.CHEM.REF.DATA, VOL.17, NO.4, 1988, PP 1795-1803.
!***********************************************************************

  real(8), parameter:: Rr = 8314.51d0, Pi = 3.14159265d0, Avgdr = 6.0221367d0, Boltz = 1.380658d0
! ATOMIC SYMBOLS
  character(2):: Symbol(100)
  data Symbol/'H ', 'D ', 'HE', 'LI', 'BE', 'B ', 'C ', 'N ', 'O ', 'F ', &
       'NE', 'NA', 'MG', 'AL', 'SI', 'P ', 'S ', 'CL', 'AR', 'K ', 'CA', 'SC', &
       'TI', 'V ', 'CR', 'MN', 'FE', 'CO', 'NI', 'CU', 'ZN', 'GA', 'GE', 'AS', &
       'SE', 'BR', 'KR', 'RB', 'SR', 'Y ', 'ZR', 'NB', 'MO', 'TC', 'RU', 'RH', &
       'PD', 'AG', 'CD', 'IN', 'SN', 'SB', 'TE', 'I ', 'XE', 'CS', 'BA', 'LA', &
       'CE', 'PR', 'ND', 'PM', 'SM', 'EU', 'GD', 'TB', 'DY', 'HO', 'ER', 'TM', &
       'YB', 'LU', 'HF', 'TA', 'W ', 'RE', 'OS', 'IR', 'PT', 'AU', 'HG', 'TL', &
       'PB', 'BI', 'PO', 'AT', 'RN', 'FR', 'RA', 'AC', 'TH', 'PA', 'U ', 'NP', &
       'PU', 'AM', 'CM', 'BK', 'CF', 'ES'/

!  ATOMIC WEIGHTS - Coplen, T.B., Atomic Weights of the Elements 1999. 
!     J.Phys.Chem.Ref.Data, vol.30, no.3, 2001, pp.701-712.

  real(8):: Atmwt(100) = (/1.00794d0, 2.014102d0, 4.002602d0, 6.941d0, &
       9.012182d0, 10.811d0, 12.0107d0, 14.0067d0, 15.9994d0, 18.9984032d0, &
       20.1797d0, 22.989770d0, 24.305d0, 26.981538d0, 28.0855d0, 30.973761d0, &
       32.065d0, 35.453d0, 39.948d0, 39.0983d0, 40.078d0, 44.95591d0, &
       47.867d0, 50.9415d0, 51.9961d0, 54.938049d0, &
       55.845d0, 58.933200d0, 58.6934d0, 63.546d0, 65.39d0, 69.723d0, 72.64d0, &
       74.92160d0, 78.96d0, 79.904d0, 83.80d0, 85.4678d0, 87.62d0, 88.90585d0, &
       91.224d0, 92.90638d0, 95.94d0, 97.9072d0, 101.07d0, 102.9055d0, &
       106.42d0, &
       107.8682d0, 112.411d0, 114.818d0, 118.710d0, 121.760d0, 127.6d0, &
       126.90447d0, 131.293d0, 132.90545d0, 137.327d0, 138.9055d0, 140.116d0, &
       140.90765d0, 144.9127d0, 145.d0, 150.36d0, 151.964d0, 157.25d0, &
       158.92534d0, &
       162.50d0, 164.93032d0, 167.259d0, 168.93421d0, 173.04d0, 174.967d0, &
       178.49d0, 180.9479d0, 183.84d0, 186.207d0, 190.23d0, 192.217d0, &
       195.078d0, 196.96655d0, 200.59d0, 204.3833d0, 207.2d0, 208.98038d0, &
       208.9824d0, 209.9871d0, &
       222.0176d0, 223.0197d0, 226.0254d0, 227.0278d0, 232.0381d0, &
       231.03588d0, 238.02891d0, 237.0482d0, 244.0642d0, 243.0614d0, &
       247.0703d0, 247.0703d0, 251.0587d0, 252.083d0/)

end module mod_legacy_cea
