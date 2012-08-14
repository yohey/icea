!***********************************************************************
!                            cea.include  
!***********************************************************************
!
!  The following parameters set the maximum dimensions for many variables
!    They are defined in part 2, page 39 of the manuals, NASA RP-1311.
!    The variable NCOL set the number of columns in the output.  It may
!    be increased for wider paper or smaller fonts.
!
module cea
  implicit none

      INTEGER MAXNGC,MAXNC,NCOL,MAXMAT,MAXTR,MAXR,MAXEL,MAXNG,MAXMIX, &
              MAXT,MAXPV
      PARAMETER (MAXNGC=600)
      PARAMETER (MAXNC=300)
      PARAMETER (NCOL=8)
      PARAMETER (MAXMAT=50)
      PARAMETER (MAXTR=40)
      PARAMETER (MAXR=24)
      PARAMETER (MAXEL=20)
      PARAMETER (MAXNG=500)
      PARAMETER (MAXMIX=52)
      PARAMETER (MAXT=51)
      PARAMETER (MAXPV=26)

!  The following parameters set the input/output unit numbers.  These
!    numbers are also defined in the manual, part2 p39, and may be 
!    adjusted as desired.

      INTEGER IOSCH,IOTHM,IOPLT,IOTRN,IOINP,IOOUT
      PARAMETER (IOINP= 7)
      PARAMETER (IOOUT= 8)
      PARAMETER (IOSCH=13)
      PARAMETER (IOTHM=14)
      PARAMETER (IOPLT=15)
      PARAMETER (IOTRN=18)
!***********************************************************************
      REAL*8 Enn,Ennl,Enlsav,Ensave,Sumn
      REAL*8 Deln(MAXNGC),Enln(MAXNGC),Sln(MAXNGC)
      REAL*8 En(MAXNGC,NCOL)

      INTEGER Ip,Iplt,It,Nc,Ng,Ngp1,Nlm,Nplt,Nof,Nomit,Nonly,Np,Npr,Npt, &
              Ngc,Nsert,Nspr,Nspx,Nt
      INTEGER Jcond(45),Jx(MAXEL),Nfla(MAXR),Ifz(MAXNC)

      REAL*8 Cpmix,Wmix,Bcheck
      REAL*8 Am(2),Hpp(2),Vmin(2),Vpls(2),Wp(2),Atmwt(100),Oxf(MAXMIX), &
             P(MAXPV),Rh(2),T(MAXT),V(MAXPV),Valnce(100)
      REAL*8 B0p(MAXEL,2)

      INTEGER Imat,Iq1,Isv,Jliq,Jsol,Lsave,Msing

      LOGICAL Convg,Debug(NCOL),Detdbg,Detn,Eql,Gonly,Hp,Ions,Massf, &
              Moles,Newr,Pderiv,Shock,Short,Siunit,Sp,Tp,Trnspt,Vol

      REAL*8 Avgdr,Boltz,Eqrat,Hsub0,Oxfl,Pi,Pp,R,Rr,Size,S0,Tln,Tm, &
             Trace,Tt,Viscns,Vv
      REAL*8 Atwt(MAXEL),B0(MAXEL),X(MAXMAT)
      REAL*8 A(MAXEL,MAXNGC),G(MAXMAT,MAXMAT+1)

      CHARACTER*2 Elmt(MAXEL),Ratom(MAXR,12),Symbol(100)
      CHARACTER*4 Fmt(30)
      CHARACTER*8 Fox(MAXR)
      CHARACTER*10 Thdate
      CHARACTER*15 Case,Energy(MAXR),Omit(0:MAXNGC),Pltvar(20), &
                   Prod(0:MAXNGC),Rname(MAXR)
      CHARACTER*200 Pfile

      REAL*8 Cpr(NCOL),Dlvpt(NCOL),Dlvtp(NCOL),Gammas(NCOL),Hsum(NCOL), &
             Ppp(NCOL),Ssum(NCOL),Totn(NCOL),Ttt(NCOL),Vlm(NCOL), &
             Wm(NCOL)
      REAL*8 Pltout(500,20)

      INTEGER Nreac
      INTEGER Jray(MAXR)
      REAL*8 Dens(MAXR),Enth(MAXR),Pecwt(MAXR),Rmw(MAXR),Rtemp(MAXR)
      REAL*8 Rnum(MAXR,12)

      REAL*8 Cpsum
      REAL*8 Cft(MAXNC,9),Coef(MAXNG,9,3),Temp(2,MAXNC)
      REAL*8 Cp(MAXNGC),H0(MAXNGC),Mu(MAXNGC),Mw(MAXNGC),S(MAXNGC),Tg(4)

      INTEGER Iopt,Isup,Nfz,Npp,Nsub,Nsup
      LOGICAL Area,Debugf,Fac,Froz,Page1,Rkt
      REAL*8 Acat,Awt,Cstr,Tcest,Ma
      REAL*8 Aeat(NCOL),App(NCOL),Pcp(2*NCOL),Sonvel(NCOL),Spim(NCOL), &
             Subar(13),Supar(13),Vmoc(NCOL)

      INTEGER Nsk
      LOGICAL Incdeq,Incdfz,Refleq,Reflfz,Shkdbg
      REAL*8 U1(NCOL),Mach1(NCOL),A1,Gamma1

      INTEGER Nm,Nr,Ntape
      INTEGER Ind(MAXTR),Jcm(MAXEL)
      REAL*8 Cprr(MAXTR),Con(MAXTR),Wmol(MAXTR),Xs(MAXTR)
      REAL*8 Eta(MAXTR,MAXTR),Stc(MAXTR,MAXTR)

      REAL*8 Coneql(NCOL),Confro(NCOL),Cpeql(NCOL),Cpfro(NCOL), &
             Preql(NCOL),Prfro(NCOL),Vis(NCOL)


!***********************************************************************
! FUNDAMENTAL CONSTANTS FROM:  COHEN,E.RICHARD & TAYLOR,BARRY N.,
! THE 1986 CODATA RECOMMENDED VALUES OF THE FUNDAMENTAL PHYSICAL
! CONSTANTS, J.PHYS.CHEM.REF.DATA, VOL.17, NO.4, 1988, PP 1795-1803.
!***********************************************************************

      DATA Rr/8314.51D0/,Pi/3.14159265D0/,Avgdr/6.0221367D0/, &
           Boltz/1.380658D0/
! ATOMIC SYMBOLS
      DATA Symbol/'H ','D ','HE','LI','BE','B ','C ','N ','O ','F ', &
           'NE','NA','MG','AL','SI','P ','S ','CL','AR','K ','CA','SC', &
           'TI','V ','CR','MN','FE','CO','NI','CU','ZN','GA','GE','AS', &
           'SE','BR','KR','RB','SR','Y ','ZR','NB','MO','TC','RU','RH', &
           'PD','AG','CD','IN','SN','SB','TE','I ','XE','CS','BA','LA', &
           'CE','PR','ND','PM','SM','EU','GD','TB','DY','HO','ER','TM', &
           'YB','LU','HF','TA','W ','RE','OS','IR','PT','AU','HG','TL', &
           'PB','BI','PO','AT','RN','FR','RA','AC','TH','PA','U ','NP', &
           'PU','AM','CM','BK','CF','ES'/

!  ATOMIC WEIGHTS - Coplen,T.B., Atomic Weights of the Elements 1999. 
!     J.Phys.Chem.Ref.Data, vol.30, no.3, 2001, pp.701-712.

      DATA atmwt/               1.00794D0,2.014102D0,4.002602D0,6.941D0, &
       9.012182D0,10.811D0,12.0107D0,14.0067D0,15.9994D0,18.9984032D0, &
       20.1797D0,22.989770D0,24.305D0,26.981538D0,28.0855D0,30.973761D0, &
       32.065D0,35.453D0,39.948D0,39.0983D0,40.078D0,44.95591D0, &
       47.867D0, 50.9415D0,51.9961D0,54.938049D0, &
       55.845D0,58.933200D0,58.6934D0,63.546D0,65.39D0,69.723D0,72.64D0, &
       74.92160D0,78.96D0,79.904D0,83.80D0,85.4678D0,87.62D0,88.90585D0, &
       91.224D0,92.90638D0,95.94D0,97.9072D0,101.07D0,102.9055D0, &
       106.42D0, &
       107.8682D0,112.411D0,114.818D0,118.710D0, 121.760D0,127.6D0, &
       126.90447D0,131.293D0,132.90545D0,137.327D0,138.9055D0,140.116D0, &
       140.90765D0,144.9127D0,145.D0,150.36D0,151.964D0,157.25D0, &
       158.92534D0, &
       162.50D0,164.93032D0,167.259D0,168.93421D0,173.04D0,174.967D0, &
       178.49D0,180.9479D0,183.84D0,186.207D0,190.23D0,192.217D0, &
       195.078D0,196.96655D0,200.59D0,204.3833D0,207.2D0,208.98038D0, &
       208.9824D0, 209.9871D0, &
       222.0176D0,223.0197D0,226.0254D0,227.0278D0,232.0381D0, &
       231.03588D0,238.02891D0,237.0482D0,244.0642D0,243.0614D0, &
       247.0703D0,247.0703D0,251.0587D0,252.083D0/
! ATOMIC VALENCES
      DATA Valnce/1.,1.,0.,1.,2.,3.,4.,0., - 2., - 1.,0.,1.,2.,3.,4.,5., &
           4., - 1.,0.,1.,2.,3.,4.,5.,3.,2.,3.,2.,2.,2.,2.,3.,4.,3.,4., &
           - 1.,0.,1.,2.,3.,4.,5.,6.,7.,3.,3.,2.,1.,2.,3.,4.,3.,4., &
           - 1.,0.,1.,2.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3.,3., &
           4.,5.,6.,7.,4.,4.,4.,3.,2.,1.,2.,3.,2., - 1.,0.,1.,2.,3.,4., &
           5.,6.,5.,4.,3.,3.,3.,3.,3./
! INFORMATION USED IN VARIABLE OUTPUT FORMAT
      DATA Fmt/'(1X',',A15',',','F9.','0,','F9.','0,','F9.','0,','F9.', &
           '0,','F9.','0,','F9.','0,','F9.','0,','F9.','0,','F9.','0,', &
           'F9.','0,','F9.','0,','F9.','0,','F9.','0',')'/

end module cea
