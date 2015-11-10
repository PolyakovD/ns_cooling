* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * **
*   CALCULATION OF ELECTRON TRANSPORT COEFFICIENTS IN MAGNETIC FIELDS  *
*           for the case  of very strong electron degeneracy           *
*      (thermal averaging is completely neglected in this version).    *
*   Whenever need arbitrary degeneracy, use version `CONDUCT' instead! *
* Last revision: 01.02.2013.                                           *
*   For theoretical background and references see:                     *
*           http://www.ioffe.rssi.ru/astro/conduct/                    *
* Difference from generic CONDEGEN - allowance for nuclear form factor *
*     in the HFB21 (BSk21) nuclear-profile approximation               *
*      Remarks and suggestions are welcome. Please send them to        *
*       Alexander Potekhin <palex@astro.ioffe.ru>                      *
*   This code stems from condegsc013.f (v.31.01.13)                    *
* Differences: - "smooth composition" is replaced by BSk21             *
*    (that is, OYAFORM is replaced by BSk21NUC);                       *
*   - non-used CONDIsc and HLfit8 subroutines are excluded.            *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * **

*   ----------------------   MAIN block   ---------------------------  *
*       This is auxiliary MAIN program for input/output purposes.      *
*                You can change it or write your own.                  *
*                     (It may be commented-out.)                       *
*     Calculations are performed in the subroutine CONDBSk (below).    *
*       Most of internal quantities are in the relativistic units      *
*         (\hbar = m_e = c = 1; default throughout the file)           *
*   -----------------------------------------------------------------  *
C%C      implicit double precision (A-H), double precision (O-Z)
C%C      save
C%C      character KEY
C%C      data BOHR/137.036/,PI/3.14159265/
C%C   20 continue
C%C      write(*,'('' Impurity parameter (effective Z): ''$)')
C%C      read*,Zimp
C%C   25 continue
C%C      write(*,'('' magnetic field [Gauss]: ''$)')
C%C      read*,B
C%C      B12=B/1.d12
C%C   50 continue
C%C      write(*,'('' lg(T[K]): ''$)')
C%C      read*,Tlg
C%C  100 continue
C%C      write(*,'('' lg(rho[g/cm^3]): ''$)')
C%C      read*,RHOlg
C%C*  Call for the central subroutine which calculates the transport
C%C*  coefficients SIGMA,KAPPA,Q (in Relativistic units):           
C%C      call CONDBSk(Tlg,RHOlg,B12,Zimp,
C%C     &  SIGMA,SIGMAT,SIGMAH,CKAPPA,CKAPPAT,CKAPPAH)
C%C* If you want to allow for ion thermal conduction in the approximation
C%C*  of Chugunov & Haensel (2007), then uncomment the next line:
C%CC      call CONDIsc(TEMP,DENSI,Zion,CMI,CMI1,RKAPi)
C%C      RKAPPA=RKAPPA+RKAPi
C%C      RTKAPPA=RTKAPPA+RKAPi
C%C*   ---------------------   OUTPUT:   ------------------------------   *
C%C      write(*,113)
C%C      write(*,111) RHOlg,Tlg,B12,SIGMA,CKAPPA,SIGMAT,CKAPPAT,
C%C     *   SIGMAH,CKAPPAH
C%C      write(*,'('' New density? (Y/N) ''$)')
C%C      read(*,'(A)') KEY
C%C      if (KEY.ne.'n'.and.KEY.ne.'N') goto 100
C%C      write(*,'('' New temperature? (Y/N) ''$)')
C%C      read(*,'(A)') KEY
C%C      if (KEY.ne.'n'.and.KEY.ne.'N') goto 50
C%C      write(*,'('' New B? (Y/N) ''$)')
C%C      read(*,'(A)') KEY
C%C      if (KEY.ne.'n'.and.KEY.ne.'N') goto 25
C%C      write(*,'('' New Zimp? (Y/N) ''$)')
C%C      read(*,'(A)') KEY
C%C      if (KEY.ne.'n'.and.KEY.ne.'N') goto 20
C%C      stop
C%C  111 format(3F8.3,1P,3(2X,2E10.2))
C%C  113 format(29X,'longitudinal',
C%C     *'          transverse            off-diagonal'/
C%C     *'  lg(rho)  lg(T)  B12 ',3('       sigma     kappa'))
C%C      end

* * * * * * * * * * * * *   Block CONDBSk  * * * * * * * * * * * * * * *
*  This subroutine calculates the electron electric and thermal        *
*      conductivity tensors in strongly degenerate matter              *
*   --------------------------------------------------  Version 01.02.13
      subroutine CONDBSk(Tlg,RHOlg,B12,Zimp,
     & SIGMA,SIGMAT,SIGMAH,CKAPPA,CKAPPAT,CKAPPAH)
* Input: Tlg - lg(T[K]), RHOlg=lg(rho[g/cc]), B12=magn.f./10^12 Gauss,
*        Zimp - impurity param.: Zimp^2 = < n_j (Z-Z_j)^2 > / n,
*          where Z_j, n_j are partial charges and densities of impur.
* Output: SIGMA, SIGMAT, SIGMAH - 
*         longitudinal, transverse and off-diagonal electr.conductiv.
*         CKAPPA,CKAPPAT,CKAPPAH - thermal conductivities in CGS units
      implicit double precision (A-H), double precision (O-Z)
      save
      data PI/3.14159265d0/
      data BOHR/137.036/
      data AUM/1822.9/,AUD/15819.4/,DRIP/4.3d11/,KEOS/21/
* NOTATIONS:
*        AUM - atomic mass unit divided by the electron mass
*        AUD - relativistic unit of density in g/cm^3
*        DRIP - neutron drip density in g/cm^3
      data UNISIG/7.763d20/,UNIKAP/2.778d15/
* Rel.units of SIGMA and KAPPA expressed in CGS.
*   -------------------   PRELIMINARY CONVERSION OF UNITS:
      if (RHOlg.lt.5..or.RHOlg.gt.15.) stop'CONDBSk: RHO out of range'
      TEMP=10.**(Tlg-9.)/5.93 ! Temperature in mc^2
      B=B12/44.14 ! B is magnetic field in relativistic units
      RHO=10.**RHOlg ! density in g/cc
      !!!!!!!!
      !write(*, '(e25.12)') RHO
      !!!!!!!!
      call BSk21NUC(RHO,XN,KZ,CMI,CMI1,xnuc) ! get baryon numbers
      Zion=KZ
      xnuct=xnuc ! fiducial nuclear-size parameter for thermal cond.
      DENSI=RHO/(AUD*AUM*CMI1) ! number density of ions
*   -------------------   RESTRICTIONS:   --------------------------   *
      !!!!!!!!
      !write(*, *) TEMP, DENSI, B, Zion, CMI
      !!!!!!!!
      if (TEMP.le.0..or.DENSI.le.0..or.B.lt.0..or.Zion.le.0.
     &  .or.CMI.le.0.) stop'CONDBSk: Non-positive input parameter'
      if (CMI1.lt.CMI) stop'CONDBSk: Incorrect CMI1'
      if (Zion.lt..5) stop'CONDBSk: Too small ion charge'
      if (CMI.lt.1.) stop'CONDBSk: Too small ion mass'
      if (DENSI.gt.1.d6) stop'CONDBSk: Too high density'
*   -----------------   PLASMA PARAMETERS   ------------------------   *
      DENS=DENSI*Zion ! DENS - number density of electrons
      SPHERION=(.75/PI/DENSI)**.3333333 ! Ion sphere radius
      GAMMA=Zion**2/BOHR/TEMP/SPHERION ! Ion coupling parameter
      XSR=(3.*PI**2*DENS)**.3333333 ! special relativity parameter 
*   XSR equals the non-magnetic Fermi momentum
      EF0=dsqrt(1.+XSR**2) ! non-mag.Fermi energy
      Q2e0=4./PI/BOHR*XSR*EF0 ! non-magn.e-screening at strong degen.
      CST=XSR**3/.75/PI*Zion/BOHR**2 ! =4\pi n_i(Ze^2)^2
*   -------------------  Chemical potential   ----------------------   *
      if (XSR**2.gt.4.d2*B) then ! non-quantizing case
         PCL=XSR
         Q2e=Q2e0
      else ! quantizing magn.field
         PM0=XSR**3/1.5/B ! p_F in strongly quantizing field
        if (PM0**2.le.2.*B) then ! strongly quantizing case
           PCL=PM0
           Q2e=B/PI/BOHR/PCL ! e-screening
        else ! weakly quantizing case - find p_F by iteration
           Pmax=PM0
           Pmin=XSR/2.
   22      PCL=(Pmax+Pmin)/2.
           NL=PCL**2/2./B
           SN=PCL
           SM=1./PCL
          do N=1,NL
             PN=dsqrt(PCL**2-2.*B*N) ! =p_n
             SN=SN+2.*PN
             SM=SM+2./PN
          enddo
           D=SN*B/2./PI**2 ! estimate n_e
          if (D.lt.DENS) then ! increase PCL
             Pmin=PCL
          else ! decrease PCL
             Pmax=PCL
          endif
          if (dabs(D-DENS).gt.1.d-4*DENS) goto 22 ! next iteration
           Q2e=B/PI/BOHR*SM
        endif
      endif
*   -------------------   Relaxation times    ----------------------   *
      call COULINsc(PCL,XSR,GAMMA,B,Zion,CMI,Q2e,xnuc,xnuct,
     *   CLeff,CLlong,CLtran,SN,THtoEL)
      E=dsqrt(1.+PCL**2) ! magn.Fermi energy
      TAU=PCL**3/E/4./PI/DENSI/(Zion/BOHR)**2/CLlong/SN
      GYROM=B/E ! gyrofrequency
      TAUt0=PCL**3/E/CST/CLtran*SN
      if (Zimp.gt.0.) call COUL99I(PCL,XSR,GAMMA,B,Q2e, ! incl.impurity
     *   CLeffI,CLlongI,CLtranI,SN)
      TAUlong=TAU*CLlong*Zion**2/(CLlong*Zion**2+CLlongI*Zimp**2)
      TAUt=TAUt0*CLtran*Zion**2/(CLtran*Zion**2+CLtranI*Zimp**2)
      TAUtran=TAUt/(1.+(TAUt*GYROM)**2)
      TAUhall=TAUt*GYROM*TAUtran
* Modification of thermal conductivity: inclusion of THtoEL (21.11.99)
      TAUlongT=TAU*CLlong*Zion**2/
     /  (CLlong*Zion**2*THtoEL+CLlongI*Zimp**2)
      TAUtT=TAUt0*CLtran*Zion**2/(CLtran*Zion**2*THtoEL+CLtranI*Zimp**2)
      TAUtranT=TAUtT/(1.+(TAUt*GYROM)**2)
      TAUhallT=TAUtT*GYROM*TAUtranT
*   ----------------------------------------------------------------   *
*   Longitudinal transport coefficients:
      C=SN*PCL**3/3./PI**2/E/BOHR ! common factor
      RSIGMA=C*TAUlong
      RTSIGMA=C*TAUtran
      RHSIGMA=C*TAUhall
* Find thermal conductivity from the Wiedemann-Franz law:
      CTH=C*PI**2*TEMP/3.*BOHR
      if (B.eq.0.) then ! corrected 12.11.07
         call TAUEESY(XSR,TEMP,TAUEE) ! eff.e-e relax.time
         EECOR=TAUEE/(TAUlongT+TAUEE)
      else
         EECOR=1.
      endif
      RKAPPA=CTH*TAUlongT*EECOR
      RTKAPPA=CTH*TAUtranT*EECOR
      RHKAPPA=CTH*TAUhallT*EECOR
*   -------  CONVERSION TO ORDINARY PHYSICAL (CGSE) UNITS:   -------   *
      SIGMA=RSIGMA*UNISIG ! SIGMA in s^{-1}
      CKAPPA=RKAPPA*UNIKAP ! KAPPA in erg/(K cm s)
      SIGMAT=RTSIGMA*UNISIG
      CKAPPAT=RTKAPPA*UNIKAP
      SIGMAH=RHSIGMA*UNISIG
      CKAPPAH=RHKAPPA*UNIKAP
      return
      end

*  ================   EFFECTIVE COULOMB LOGARITHM  ==================  *
      subroutine COULINsc(PCL,XSR,GAMMA,B,Zion,CMI,Q2e,xnuc,xnuct,
     *   CLeff,CLlong,CLtran,SN,THtoEL)
*                                                       Version 24.02.00
*   Difference from COUL99 - inclusion of xnuc
*   Input: PCL - non-magnetic electron momentum \equiv \sqrt(E^2-1),
*          Q2e - squared electron screening wavenumber (IN REL.UNITS)
*          XSR = p_F/mc - relativity (density) parameter,
*          GAMMA - Coulomb coupling parameter of ions,
*          B - magnetic field,
*          Zion - mean charge of the ion,
*          CMI - mean atomic weight,
*          xnuc - ion radius divided by Wigner-Seitz cell radius
*   Output: CLlong, CLtran - eff.Coulomb log.,
*           SN = N_e(E)/N_0(E) = (3/2)(eB\hbar/c)\sum_{ns} p_n/p_0^3
      implicit double precision (A-H), double precision (O-Z)
      save
      data Uminus1/2.78/,Uminus2/12.973/,AUM/1822.9/,BOHR/137.036/
* Dimensional quantities are in the relativistic units (m_e=\hbar=c=1)
*        Uminus1,Uminus2 - dimensionless frequency moments of phonons
*        AUM - atomic mass unit divided by the electron mass
*        BOHR - radius of the first Bohr orbit in the rel.units
      data PI/3.14159265/
*   ----------------------   Preliminaries   -----------------------   *
      DENS=XSR**3/3./PI**2 ! number density of electrons
      DENSI=DENS/Zion ! number density of ions (rel.)
      SPHERION=(.75/PI/DENSI)**.3333333 ! Ion sphere radius
      Q2icl=3.*GAMMA/SPHERION**2 ! squared Debye screening momentum
      ECL=sqrt(1.+PCL**2) ! Energy
      VCL=PCL/ECL ! Velocity
      PM2=(2.*PCL)**2 ! squared max.momentum transfer
      TRP=Zion/GAMMA*sqrt(CMI*AUM*SPHERION/3./BOHR) ! =T/T_p
      BORNCOR=VCL*Zion*PI/BOHR ! first non-Born correction
*   ---------------------   Non-magnetic fit   ---------------------   *
      C=(1.+.06*GAMMA)*dexp(-dsqrt(GAMMA))
      Q2s=(Q2icl*C+Q2e)*dexp(-BORNCOR) ! eff.scr.wavenumber in rel.un.
      XS=Q2s/PM2 ! eff.screening param.
      R2W=Uminus2/Q2icl*(1.+.3333*BORNCOR)
      XW=R2W*PM2 ! eff. Debye-Waller param.
** Modification WITH FINITE SIZES OF NUCLEI; xnuc=r_{nuc}/a_i
      XW1=14.7327*xnuc**2 ! =4(9\pi/4)^{2/3} x_{nucl}^2 =coeff.at q^2
      XW1=XW1*(1.+.3333*BORNCOR)*(1.+Zion/13.*dsqrt(xnuc))
      CL=COULAN2(XS,XW,VCL,XW1)
      A0=1.683*sqrt(PCL/CMI/Zion) ! zero-vibr.param.(Baiko&Yakovlev95)
      VIBRCOR=exp(-A0/4.*Uminus1*exp(-9.1*TRP)) ! corr.for zero-vibr.
      T0=.19/Zion**.16667 ! critical T/T_p parameter
      G0=TRP/sqrt(TRP**2+T0**2)*(1.+(Zion/125.)**2) ! mod.10.01.99
      GW=G0*VIBRCOR
      CLeff=CL*GW ! 1st FIT (for non-magnetic electrical conductivity)
      G2=TRP/sqrt(.0081+TRP**2)**3
      THtoEL=1.+G2/G0*(1.+BORNCOR*VCL**3)*.0105*(1.-1./Zion)*
     *  (1.d0+xnuct**2*dsqrt(2.d0*Zion))
      if (PCL**2.gt.4.d2*B) then ! Non-magnetic case
         CLlong=CLeff
         CLtran=CLeff
         SN=1.d0
         goto 50
      endif
*   -----------------------   Magnetic fit   -----------------------   *
      ENU=PCL**2/2.d0/B
      NL=ENU
      SN=0.
      do N=0,NL
         PB=dsqrt(ENU-N) ! =p_n/sqrt(2b)
         SN=SN+PB
        if (N.ne.0) SN=SN+PB
      enddo
      SN=SN*1.5d0*B*dsqrt(2.d0*B)/PCL**3
      if (ENU.le.1.d0) then ! Exact calculation     
         Xis=Q2s/2./B ! Screening parameter, scaled magnetically
         ZETA=R2W*2.*B ! magn.scaled exponential coefficient
         Xi=2.*PCL**2/B
         Xsum=Xi+Xis
         Q2M=(EXPINT(Xsum,1)-
     -     dexp(-ZETA*Xi)*EXPINT((1.+ZETA)*Xsum,1))/Xsum
         CLlong=(PCL*VCL/B)**2*Q2M/1.5*GW
         QtranM=(1.+Xsum)*EXPINT(Xsum,0)-1.-dexp(-ZETA*Xi)*
     *      ((1.+(1.+ZETA)*Xsum)*EXPINT((1.+ZETA)*Xsum,0)-1.)
         QtranP=(1.+Xis)*EXPINT(Xis,0)-1.-
     -   ((1.+(1.+ZETA)*Xis)*EXPINT((1.+ZETA)*Xis,0)-1.)
         Q=(ECL**2*QtranP+QtranM)*B/PCL**2 ! Q(E,b)
         CLtran=.375*Q/ECL**2*GW
      else
*   Preliminaries:
         DNU=ENU-NL
         XS1=(dsqrt(XS)+1./(2.+XW/2.))**2
         PN=dsqrt(2.*B*DNU)
         SQB=dsqrt(B)
         X=dmax1(PN/SQB,1.d-10)
*   Longitudinal:
        if (XW.lt..01) then
           EXW=1.
        elseif (XW.gt.50.) then
           EXW=1./XW
        else
           EXW=(1.d0-dexp(-XW))/XW
        endif
         A1=(30.-15.*EXW-(15.-6.*EXW)*VCL**2)/
     /     (30.-10.*EXW-(20.-5.*EXW)*VCL**2)
         Q1=.25*VCL**2/(1.-.667*VCL**2)
         DLT=SQB/PCL*(A1/X-sqrt(X)*(1.5-.5*EXW+Q1)+
     +     (1.-EXW+.75*VCL**2)/(1.+VCL**2)*(X-sqrt(X))/NL)
         Y1=1./(1.+DLT)
         CL0=dlog(1.d0+1.d0/XS1)
         P2=CL0*(.07+.2*EXW)
         Y2=1.5*CL0*(X**3-X/3.)/(NL+.75/(1.+2.*B)**2*X**2)+P2*X
         PY=1.+.06*CL0**2/NL**2
         DT=dsqrt(PY*Y1**2+Y2**2) ! ratio of relax.times
         CLlong=CLeff/DT
*   Transverse:
         DB=1./(1.+.5/B)
         CL1=XS1*CL0
         P1=.8*(1.+CL1)+.2*CL0
         P2=1.42-.1*DB+sqrt(CL1)/3.
         P3=(.68-.13*DB)*CL1**.165
         P4=(.52-.1*DB)*sqrt(sqrt(CL1))
         DLT=SQB/PCL*(P1/X**2*SQB/PCL+P3*alog(NL+0.)/X-
     -     (P2+P4*alog(NL+0.))*sqrt(X))
         CLtran=CLeff*(1.+DLT)
      endif
   50 return
      end

      subroutine COUL99I(PCL,XSR,GAMMA,B,Q2e,
     *   CLeff,CLlong,CLtran,SN) ! IMPURITY
*                                                       Version 29.10.99
*  This is a simplified version of COUL99 for the el.-impurity scattering
*   Input: XSR = p_F/mc - relativity (density) parameter,
*          PCL - non-magnetic electron momentum \equiv \sqrt(E^2-1),
*          GAMMA - Coulomb coupling parameter of ions,
*          B - magnetic field,
*          Q2e - squared electron screening wavenumber
*    (ALL IN THE RELATIVISTIC UNITS) 
*   Output: CLlong, CLtran - eff.Coulomb log.,
*           SN = N_e(E)/N_0(E) = (3/2)(eB\hbar/c)\sum_{ns} p_n/p_0^3
      implicit double precision (A-H), double precision (O-Z)
      save
* Dimensional quantities are in the relativistic units (m_e=\hbar=c=1)
*        BOHR - radius of the first Bohr orbit in the rel.units
      data PI/3.14159265/,XW/1.d99/ ! XW=infinity
*   ----------------------   Preliminaries   -----------------------   *
      DENS=XSR**3/3./PI**2 ! number density of electrons
      ECL=sqrt(1.+PCL**2) ! Energy
      VCL=PCL/ECL ! Velocity
      PM2=(2.*PCL)**2 ! squared max.momentum transfer
*   ---------------------   Non-magnetic fit   ---------------------   *
      C=(1.+.06*GAMMA)*dexp(-dsqrt(GAMMA))
      Q2s=Q2e ! eff.scr.wavenumber in rel.un.
      XS=Q2s/PM2 ! eff.screening param.
      CLeff=COULAN2(XS,XW,VCL,0.d0) ! 1st FIT (non-magn.el.conductivity)
      if (PCL**2.gt.4.d2*B) then ! Non-magnetic case
         CLlong=CLeff
         CLtran=CLeff
         SN=1.d0
         goto 50
      endif
*   -----------------------   Magnetic fit   -----------------------   *
      ENU=PCL**2/2.d0/B
      NL=ENU
      SN=0.
      do N=0,NL
         PB=dsqrt(ENU-N) ! =p_n/sqrt(2b)
         SN=SN+PB
        if (N.ne.0) SN=SN+PB
      enddo
      SN=SN*1.5d0*B*dsqrt(2.d0*B)/PCL**3
      if (ENU.le.1.d0) then ! Exact calculation     
         Xis=Q2s/2./B ! Screening parameter, scaled magnetically
         Xi=2.*PCL**2/B
         Xsum=Xi+Xis
         Q2M=EXPINT(Xsum,1)
         CLlong=(PCL*VCL/B)**2*Q2M/1.5
         QtranM=(1.+Xsum)*EXPINT(Xsum,0)-1.
         QtranP=(1.+Xis)*EXPINT(Xis,0)-1.
         Q=(ECL**2*QtranP+QtranM)*B/PCL**2 ! Q(E,b)
         CLtran=.375*Q/ECL**2
      else
*   Preliminaries:
         DNU=ENU-NL
         XS1=(dsqrt(XS)+1./(2.+XW/2.))**2
         PN=dsqrt(2.*B*DNU)
         SQB=dsqrt(B)
         X=dmax1(PN/SQB,1.d-10)
*   Longitudinal:
        if (XW.lt..01) then
           EXW=1.
        elseif (XW.gt.50.) then
           EXW=1./XW
        else
           EXW=(1.d0-dexp(-XW))/XW
        endif
         A1=(30.-15.*EXW-(15.-6.*EXW)*VCL**2)/
     /     (30.-10.*EXW-(20.-5.*EXW)*VCL**2)
         Q1=.25*VCL**2/(1.-.667*VCL**2)
         DLT=SQB/PCL*(A1/X-sqrt(X)*(1.5-.5*EXW+Q1)+
     +     (1.-EXW+.75*VCL**2)/(1.+VCL**2)*(X-sqrt(X))/NL)
         Y1=1./(1.+DLT)
         CL0=dlog(1.d0+1.d0/XS1)
         P2=CL0*(.07+.2*EXW)
         Y2=1.5*CL0*(X**3-X/3.)/(NL+.75/(1.+2.*B)**2*X**2)+P2*X
         PY=1.+.06*CL0**2/NL**2
         DT=dsqrt(PY*Y1**2+Y2**2) ! ratio of relax.times
         CLlong=CLeff/DT
*   Transverse:
         DB=1./(1.+.5/B)
         CL1=XS1*CL0
         P1=.8*(1.+CL1)+.2*CL0
         P2=1.42-.1*DB+sqrt(CL1)/3.
         P3=(.68-.13*DB)*CL1**.165
         P4=(.52-.1*DB)*sqrt(sqrt(CL1))
         DLT=SQB/PCL*(P1/X**2*SQB/PCL+P3*alog(NL+0.)/X-
     -     (P2+P4*alog(NL+0.))*sqrt(X))
         CLtran=CLeff*(1.+DLT)
      endif
   50 return
      end

      function COULAN2(XS,XW0,V,XW1)
*  ------    Analytic expression for Coulomb logarithm  Version 23.05.00
*   XS=(q_s/2p)^2, where p - momentum, q_s - eff.scr.momentum
*   XW0=u_{-2} (2p/\hbar q_D)^2, where u_{-2}=13, q_D^2=3\Gamma/a_i^2
*   V=p/(mc)
*   XW1=s1*(2p/\hbar)^2, s1 \approx r_{nuc}^2
      implicit double precision (A-H), double precision (O-Z)
      save
      data EPS/1.d-2/,EPS1/1.D-3/,EULER/0.5772156649 d0/
      if(XS.lt.0..or.XW0.lt.0..or.V.lt.0..or.XW1.lt.0.)stop'COULAN2'
      do I=0,1
        if (I.eq.0) then
           XW=XW0+XW1
           B=XS*XW
        else ! to do the 2nd term
          XW=XW1
          B=XS*XW
        endif
        if (I.eq.0.or.KEY.eq.2) then ! 23.05.00: for KEY=2 re-check
* Check applicability of asymptotes:
          if (XW.lt.EPS) then
             KEY=1
             goto 50
          endif
          if (XW.gt.1./EPS.and.B.gt.1./EPS) then
             KEY=2
          elseif (XS.lt.EPS1.and.B.lt.EPS1/(1.+XW)) then
             KEY=3
          else
             KEY=4
          endif
        endif
   50   continue
         EA=dexp(-XW)
         E1=1.-EA
        if (KEY.ne.1) E2=(XW-E1)/XW
        if (KEY.eq.1) then
           CL0=dlog((XS+1.)/XS)
           CL1=.5*XW*(2.-1./(XS+1.)-2.*XS*CL0)
           CL2=.5*XW*(1.5-3.*XS-1./(XS+1.)+3.*XS**2*CL0)
        elseif (KEY.eq.2) then
           CL0=dlog(1.d0+1./XS)
           CL1=(CL0-1.d0/(1.+XS))/2.
           CL2=(2.*XS+1.)/(2.*XS+2.)-XS*CL0
        elseif (KEY.eq.3) then
           CL1=.5*(EA*EXPINT(XW,0)+dlog(XW)+EULER)
           CL2=.5*E2
        elseif (KEY.eq.4) then
           CL0=dlog((XS+1.)/XS)
           EL=EXPINT(B,0)-EXPINT(B+XW,0)*EA
           CL1=.5*(CL0+XS/(XS+1.)*E1-(1.+B)*EL)
           CL2=.5*(E2-XS*XS/(1.+XS)*E1-2.*XS*CL0+XS*(2.+B)*EL)
        else
           stop'COULAN2:invalid KEY'
        endif
        if (I.eq.0) then ! 1st term calculated
           COULAN2=CL1-V**2*CL2
          if (XW1.lt.EPS1) return ! don't calculate the 2nd term
        else ! 2nd term calculated
           COULAN2=COULAN2-(CL1-V**2*CL2)
        endif
      enddo
      return
      end
*   ==================  AUXILIARY SUBROUTINE  ================   *
      function EXPINT(XI,L) ! = e^XI E_{L+1}(XI)
      implicit double precision (A-H), double precision (O-Z)
      save
      data GAMMA/.5772156649D0/,Nrep/21/
      if (XI.ge.1.) then ! continued fraction
         CL=L
         CI=Nrep
         C=0.
         do 11 I=Nrep,1,-1
         C=CI/(XI+C)
         C=(CL+CI)/(1.+C)
   11    CI=CI-1.
         Q0=1./(XI+C)
      else ! power series
         PSI=-GAMMA
         do 21 K=1,L
   21    PSI=PSI+1./K ! Psi(L+1)
         Q0=0.
         CMX=1. ! (-XI)^M/M!
         CL=L
         CM=-1.
         do 22 M=0,Nrep
         CM=CM+1.
        if (M.ne.0) CMX=-CMX*XI/CM ! (-XI)^M/M!
        if (M.ne.L) then ! Actually DQ=-deltaQ0
           DQ=CMX/(CM-CL)
        else
           DQ=CMX*(dlog(XI+1.d-20)-PSI)
        endif
         Q0=Q0-DQ
   22    continue
         Q0=exp(XI)*Q0
      endif
   50 continue
      EXPINT=Q0
      return
      end

*   ----------------------------------------------------

      subroutine TAUEESY(X,TEMP,TAUEE)
      implicit double precision (A-H), double precision (O-Z)
* Relaxation time of electron-electron collisions
*   according to Shternin & Yakovlev (2006),
*   corrected in the nondegenerate regime so as to match Lampe (1968)
* X - rel.parameter, DENS - electron density, TEMP -temperature[rel.un.]
*                                                      Version 17.07.06
      E=dsqrt(1.+X**2)
      V=X/E
      Y=.0963913/TEMP*X*dsqrt(V) ! .096391=2\sqrt{\alpha/\pi}
      C1=.123636+.016234*V**2
      C2=.0762+.05714*V**4
      A=12.2+25.2*V**3
      C=A*dexp(C1/C2)
      YV=Y*V
      CIL=dlog(1.d0+128.56/(37.1*Y+10.83*Y**2+Y**3))*
     * (.1587-.02538/(1.+.0435*Y))/V
      CIT=V**3*dlog(1.d0+C/(A*YV+YV**2))*
     *  (2.404/C+(C2-2.404/C)/(1.+.1*YV))
      CILT=V*dlog(1.d0+C/(A*Y+10.83*YV**2+YV**(8.d0/3.d0)))*
     *  (18.52*V**2/C+(C2-18.52*V**2/C)/(1.+.1558*Y**(1.d0-.75d0*V)))
      FI=CIL+CIT+CILT
      FREQ=.00021381*X*Y*dsqrt(V)*FI ! .00021381=6\alpha^{3/2}/\pi^{5/2}
* Correction for partial degeneracy:
      if (X.gt..001) then
         THETA=TEMP/(E-1.d0) ! degeneracy
      else
         THETA=2.*TEMP/X**2
      endif
      T=25.*THETA
      TAUEE=(1.+T+.4342*dsqrt(THETA)*T**2)/(1.+T**2)/FREQ
      return
      end

      subroutine BSk21NUC(RHO,XN,KZ,CMI,CMI1,xnuc)
*                                                       Version 01.02.13
* Input: RHO
* Output: XN,KZ,CMI,CMI1
* BSk21 composition of the crust
      implicit double precision (A-H), double precision (O-Z)
      save AX,AA,A1,PARxp
      parameter(XND=2.57541d-4,XNC=.0809) ! neutron drip and crust/core
      parameter(SphPI=.238732415d0) ! (4\pi/3)^{-1}
      parameter(C13=1.d0/3.d0)
      dimension AX(9),AA(4),A1(6),PARxp(5)
      data AX/.085,.802,16.35,3.634,2.931,12.0190,1.050E-05,3.187,.611/
      data AA/133.2,9.5,8.69,.003276/
      data A1/132.6,188.5,227.3,.5218,.00656,.00155/
      data PARxp/.1060,3.75,.600,7.72,3.17/
* Determine baryon number density XN (extraction from BSkNofR):
      X=RHO/1.66d15
      F1=(AX(1)*X**AX(2)+AX(3)*X**AX(4))/(1.d0+AX(5)*X)**3
      F2=X/(AX(7)+AX(8)*X**AX(9))
      RLG=dlog10(RHO)
      XF=RLG-AX(6)
      if (X.gt.40.d0) then
         FRM=0.d0
      elseif (X.lt.-40.d0) then
         FRM=1.d0
      else
         FRM=1.d0/(dexp(X)+1.d0)
      endif
      F=FRM*F2+(1.-FRM)*F1
      XN=X/(1.d0+F)
* XN is set.
      if (XN.gt.XNC) stop'BSk21NUC: NS core is reached.'
      if (XN.lt.XND) then ! outer crust
        if (RHO.lt.8.14d6) then
           KZ=26
           CMI=56.
           Rch=3.743 ! rms charge radius [fm]
        elseif (RHO.lt.2.71d8) then
           KZ=28
           CMI=62.
           Rch=3.860
        elseif (RHO.lt.1.33d9) then
           KZ=28
           CMI=64.
           Rch=3.880
        elseif (RHO.lt.1.44d9) then
           KZ=28
           CMI=66.
           Rch=3.898
        elseif (RHO.lt.3.13d9) then
           KZ=36
           CMI=86.
           Rch=4.187
        elseif (RHO.lt.1.14d10) then
           KZ=34
           CMI=84
           Rch=4.144
        elseif (RHO.lt.2.80d10) then
           KZ=32
           CMI=82.
           Rch=4.094
        elseif (RHO.lt.5.35d10) then
           KZ=30
           CMI=80.
           Rch=4.040
        elseif (RHO.lt.6.26d10) then
           KZ=29
           CMI=79.
           Rch=4.008
        elseif (RHO.lt.9.04d10) then
           KZ=28
           CMI=78.
           Rch=3.967
        elseif (RHO.lt.1.28d11) then
           KZ=28
           CMI=80.
           Rch=3.994
        elseif (RHO.lt.2.09d11) then
           KZ=42
           CMI=124.
           Rch=4.579
        elseif (RHO.lt.2.51d11) then
           KZ=40
           CMI=122.
           Rch=4.532
        elseif (RHO.lt.2.96d11) then
           KZ=39
           CMI=121.
           Rch=4.512
        elseif (RHO.lt.3.28d11) then
           KZ=38
           CMI=120.
           Rch=4.497
        elseif (RHO.lt.4.02d11) then
           KZ=38
           CMI=122.
           Rch=4.523
        else
           KZ=38
           CMI=124.
           Rch=4.540
        endif
         CMI1=CMI
         RWS=(SphPI*CMI1/XN)**C13 ! Wigner-Seitz cell radius [fm]
         xnuc=Rch/RWS/dsqrt(.6d0) ! equivalent unform-density radius
      else ! inner crust
         KZ=40
* extraction from FRACRUST:
         X=XN/XND
         Xlg=dlog10(X)
         CMI=(AA(1)+AA(2)*Xlg+AA(3)*Xlg**3)*
     *     dsqrt(dsqrt(dim(1.d0,(AA(4)*X)**2)))
         CMI1=(A1(1)+A1(2)*Xlg+A1(3)*Xlg**2)/(1+(A1(4)*Xlg)**4)*
     *     (1+A1(5)*X)*dim(1.d0,(A1(6)*X)**2)
         !!write(*, *) 'CMI ', CMI
         !!write(*, *) 'CMI1 ', CMI1
* extraction from NUCPAR:
         xnuc=(PARxp(1)+(PARxp(2)*XN)**PARxp(3))/
     /     (1.d0-(PARxp(4)*XN)**PARxp(5))
      endif
      return
      end