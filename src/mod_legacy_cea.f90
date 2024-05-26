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

  ! for subroutine CPHS and ALLCON
  real(8):: cx(7) = [0d0, 0d0, 1d0, 0.5d0, 0.6666666666666667d0, 0.75d0, 0.8d0]
  real(8):: hcx(7) = [0d0, 0d0, 1d0, 0d0, 0d0, 0d0, 0d0]
  real(8):: scx(7)

end module mod_legacy_cea
