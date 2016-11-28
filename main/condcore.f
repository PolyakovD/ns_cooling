* Calculation of conductivities in the NS core
* Uses: sfreduc.f gusakov.f
      subroutine CONDCO(XN,Yn,Yp,Ye,Ymu,T9,Tcp9,Tcn9,EFMn,EFMp,INMED,
     *  CKAPn,CKAPe,CKAPmu)
*                                                       Version 09.06.16
* Stems from CONDCOR v. 26.07.15. Difference - INMED key included
* Neutron conductivity: Baiko, Haensel, Yakovlev, A&A 374, 151 (2001)
* Lepton conductivity: Shternin, Yakovlev, PRD 75, 103004 (2007)
* Input: XN - baryon number density [fm^{-3}]
*        Yn - neutron number fraction
*        Yp - proton number fraction
*        T9 - temperature [GK],
*        Tcp9 - SF critical temperature T_c for protons
*        Tcn9 - T_c for neutrons
*        EFMn, EFMp = m*/m - effective mass factors (n,p)
*        INMED controls the use of in-medium correction factors
* Output: CKAPn - neutron thermal conductivity [erg/(cm s K)]
*         CKAPe - electron thermal conductivity
*         CKAPmu - muon thermal conductivity
      implicit double precision (A-H), double precision (O-Z)
      FUNC1(A,VT2)=A+dsqrt((1.d0-A)**2+VT2)
      FUNC2(A,VT2)=dsqrt(VT2+A**2)-A
      FUNCX(A,VT)=dexp(VT-FUNC2(A,VT**2))
      FUNCXE(A,VT)=dexp(-FUNC2(A,VT**2))
      include'const.inc'
      parameter(UNIKAP=BOLCON/(COMPTON/1.d13)**2*CSPEED) ! erg/(K cm s)
      parameter(COFREQ=3.48d17)
* COFREQ = 64 m_N (k_B*(1 GK))^2 / (5 \hbar^3) [mb^{-1} s^{-1}]
      parameter(COKAP=3.744111283d40*1.2) ! C=1.2: BHY-Sect.2.7
* (\pi^2/3)*k_B^2*(1 GK)*(1 fm^{-3})/m_n [erg/(K cm s^2)]
      parameter(C16=1.d0/6.d0,TINY=1.d-19)
      parameter(CP3=45.d0*zeta3/TwoPI**2) ! 45\zeta(3)/(4\pi^2)
      parameter(CQ2=4.d0/PI/BOHR) ! 4\alpha/\pi
      parameter(CnuTRAN=24.d0*zeta3/PI**3/BOHR**2) ! coef.in SY Eq.(34)
      parameter(CnuLONG=.8d0*PI**2/BOHR**2) ! 4\pi^2\alpha^2/5: SY(35)
      parameter(CnuPRIM=12.d0*18.52d0/PI**3/BOHR**2) !coef.in SY Eq.(36)
      TRI2XN=3.d0*PI**2*XN ! 3\pi^2 n_b [fm^{-3}]
      CKFN=(TRI2XN*Yn)**C13 ! k_{Fn} [fm^{-1}]
      CKFP=(TRI2XN*Yp)**C13 ! k_{Fp} [fm^{-1}]
      CKFE=(TRI2XN*Ye)**C13 ! k_{Fe} [fm^{-1}]
      CKFM=(TRI2XN*Ymu)**C13 ! k_{F\mu} [fm^{-1}]
      CKFPR=CKFP*COMPTON ! k_{Fp} [rel.un.]
      XSRE=CKFE*COMPTON ! k_{Fe} [rel.un.] = p_{Fe}/m_e c
      CKFMR=CKFM*COMPTON ! k_{F\mu} [rel.un.]
      XSRP=CKFPR/DMPE/EFMp ! p_{Fp}/m_p c
      XSRM=CKFMR/RATMUE ! p_{F\mu}/m_\mu c
      GSRP=dsqrt(1.d0+XSRP**2) ! Lorentz factor of protons
      GSRE=dsqrt(1.d0+XSRE**2) ! Lorentz factor of electrons
      GSRM=dsqrt(1.d0+XSRM**2) ! Lorentz factor of muons
* In-vacuum integrals [mb], Eq.(29):
      S0N1=14.57/CKFN/dsqrt(CKFN)*(1.d0-.0788*CKFN+.0883*CKFN**2)/
     /  (1.d0-.1114*CKFN)
      S0N2=7.88/CKFN**2*(1.d0-.2241*CKFN+.2006*CKFN**2)/
     /  (1.d0-.1742*CKFN)
      S0P1=.8007*CKFP/CKFN**2*
     *  (1.d0+31.28*CKFP-.0004285*CKFP**2+26.85*CKFN+.08012*CKFN**2)/
     /  (1.d0-.5898*CKFN+.2368*CKFN**2+.5838*CKFP**2+.884*CKFN*CKFP)
      S0P2=.383*CKFP**4/(CKFN**5*dsqrt(CKFN))*
     *  (1.d0+102.*CKFP+53.91*CKFN)/
     /  (1.d0-.7087*CKFN+.2537*CKFN**2+9.404*CKFP**2-1.589*CKFN*CKFP)
* In-medium correction factors, BHY, Eq.(30):
      if (INMED.eq.0) then ! no in-medium effects
         CKN1=1.d0
         CKN2=1.d0
         CKP1=1.d0
         CKP2=1.d0
      elseif (INMED.eq.1) then ! in-medium effects according to BHY
         U=CKFN-1.665
         CKN1=(.4583+.892*U**2-.5497*U**3-.06205*CKFP+.04022*CKFP**2+
     +     .2122*U*CKFP)/EFMn**2
         U=CKFN-1.556
         CKN2=(.4891+1.111*U**2-.2283*U**3+.01589*CKFP-.02099*CKFP**2+
     +     .2773*U*CKFP)/EFMn**2
         U=CKFN-2.126
         CKP1=(.04377+1.100*U**2+.1180*U**3+.1626*CKFP+.3871*U*CKFP-
     -     .2990*U**4)/EFMp**2
         U=CKFN-2.116
         CKP2=(.0001313+1.248*U**2+.2403*U**3+.3257*CKFP+.5536*U*CKFP-
     -     .3237*U**4+.09786*U**2*CKFP)/EFMp**2
      elseif (INMED.eq.2) then ! approximation to Shternin et al.(2013)
         CKN1=1.7d0
         CKN2=1.7d0
         CKP1=1.7d0
         CKP2=1.7d0
      else
         stop'CONDCO: unknown INMED'
      endif
* Superfluidity correction factors, BHY Eqs.(45)-(48), (50), (51):
      RN1=1.d0
      RN2=1.d0
      RP1=1.d0
      RP2=1.d0
      RC=1.d0
      VTn=0.
      EVTn=1.d0
      if (T9.lt.Tcn9) then ! neutron SF
         call DLGAPS(2,T9/Tcn9,VTn) ! only type-B SF is implemented
         VTn2=VTn**2
         EVTn=dexp(-VTn)
         RN1=FUNC1(.9468d0,.5346d0*VTn2)**3/1.5d0*
     *     FUNCX(.377d0,2.d0*VTn)+
     +     (1.d0+1.351d0*VTn2)**2/3.d0*
     *     FUNCX(.169d0,3.d0*VTn)*EVTn ! BHY Eq.(45a)*exp(2*VTn)
         RN2=.5d0*(FUNC1(.6242d0,.07198d0*VTn2)*
     *     FUNCX(3.6724d0,2.d0*VTn)+
     +     (1.d0+.01211*VTn2)**9*
     *     FUNCX(7.5351d0,3.d0*VTn)*EVTn) ! BHY (45b)*exp(2*VTn)
         RC1=FUNC1(.647d0,.109d0*VTn2)
         RC=RC1*dsqrt(RC1)*FUNCX(1.39d0,VTn) ! BHY Eq.(51)*exp(VTn) 
        if (T9.lt.Tcp9) then ! neutron and proton SF
           call DLGAPS(1,T9/Tcp9,VTp)
           VTp2=VTp**2
           UMAX=FUNC2(1.485d0,dmax1(VTn2,VTp2)) ! u_+
           UMIN=FUNC2(1.485d0,dmin1(VTn2,VTp2)) ! u_-
           UN=FUNC2(1.485d0,VTn2) ! u_n
           UP=FUNC2(1.485d0,VTp2) ! u_p
           E1=dexp(VTn-UMAX-UMIN)
           E2=dexp(VTn-2.d0*UMAX)
           RP1=E1*(.7751d0+.4823d0*UN+.1124d0*UP+
     +         .04991d0*UN**2+.08513d0*UN*UP+.01284d0*UN**2*UP)+
     +       E2*(.249d0+.3539d0*UMAX-.2189d0*UMIN-
     -         .6069d0*UN*UMIN+.7362d0*UP*UMAX) ! BHY (47)
           UMAX=FUNC2(1.761d0,dmax1(VTn2,VTp2)) ! u_+
           UMIN=FUNC2(1.761d0,dmin1(VTn2,VTp2)) ! u_-
           UN=FUNC2(1.761d0,VTn2) ! u_n
           UP=FUNC2(1.761d0,VTp2) ! u_p
           E1=dexp(VTn-UMAX-UMIN)
           E2=dexp(VTn-2.d0*UMAX)
           RP2=E1*(1.1032d0+.8645d0*UN+.2042d0*UP+
     +         .07937d0*UN**2+.1451d0*UN*UP+.01333d0*UN**2*UP)+
     +       E2*(-.1032d0-.234d0*UMAX+.06152*UN*UMAX+
     +         .7533d0*UN*UMIN-1.007d0*UP*UMAX) ! BHY (48)*e^(VTn)
        else ! only neutron SF: BHY (46)*e^(VTn)
           RP1=FUNC1(.4459d0,.03016d0*VTn2)**2*FUNCX(2.1178d0,VTn)
           RP2=FUNC1(.801d0,.04645d0*VTn2)*FUNCX(2.3569d0,VTn)
        endif
      elseif (T9.lt.Tcp9) then ! only proton SF
         call DLGAPS(1,T9/Tcp9,VTp)
         VTp2=VTp**2
         RP1=.5d0*(FUNC1(.3695d0,.01064*VTp2)*FUNCXE(2.4451d0,VTp)+
     +     (1.d0+.1917d0*VTp2)**1.4d0*FUNCXE(4.6627d0,2.d0*VTp))
         RP2=.0436d0*FUNC1(-3.345d0,19.55d0*VTp2)*FUNCXE(2.0247d0,VTp)+
     +     .0654d0*FUNCXE(8.992d0,dsqrt(1.5d0)*VTp)+
     +     .891d0*FUNCXE(9.627d0,3.d0*VTp)
      endif
* Effective in-medium-and-SF-corrected cross sections, Eq.(50):
      SN1=S0N1*CKN1
      SN2=S0N2*CKN2
      SP1=S0P1*CKP1
      SP2=S0P2*CKP2
      SNN=SN2*RN2+3.d0*SN1*(RN1-RN2)
      SNP=SP2*RP2+SP1*(3.d0*RP1-RP2)/2.d0
      FREQnn=COFREQ*EFMn**3*T9**2*SNN ! Eq.(21): \nu_{nn}*e^(2*VTn)
      FREQnp=COFREQ*EFMn*EFMp**2*T9**2*SNP ! \nu_{np}*e^(VTn)
      if (FREQnp.eq.0.) then ! TAURC2=TAU*RC**2
         TAURC2=RC**2/FREQnn
      elseif (FREQnn*EVTn.lt.FREQnp) then
         TAURC2=RC**2*EVTn/(FREQnn*EVTn+FREQnp)
      else
         TAURC2=RC**2/(FREQnn+FREQnp/EVTn)
      endif
      CKAPn=COKAP*XN*T9/EFMn*TAURC2 ! Eq.(50)
* Leptons:
      CKAPe=0.
      CKAPmu=0.
      if (Ye.lt.TINY) return ! no leptons
      QL2=CQ2*(XSRE*GSRE+CKFMR*GSRM*RATMUE+EFMp*DMPE*GSRP*CKFPR) ! (16)
      QT2=CQ2*(XSRE**2+CKFMR**2+CKFPR**2) !q_t^2 [rel.un.]: SY(17)
      TEMR=T9*1.d3/UNITEMP ! T in rel.un.
* SF reduction factors:
      RLTL=1.d0 ! R_l^{\|\perp}
      RPT=1.d0 ! R_p^\perp
      RTOT=1.d0 ! R_{tot}^\perp
      RPAR=1.d0 ! R_p^\|
      if (T9.lt.Tcp9) then ! proton SF reduction factors
         R=(CKFE**2+CKFM**2)/CKFP**2 ! SY(54)
         R1=R+1.d0
         RLTL=R1**C13/(R1**2-.757d0*VTp+.25655*VTp2)**C16 ! SY(61)
         P1=.48d0-.17d0*R
         P3=((1.d0-P1)*CP3/R)**2
         RTOT=P1*dexp(-.14d0*VTp2)+(1.d0-P1)/dsqrt(1.d0+P3*VTp2) !SY(88)
* RPAR is taken from SY Eq.(92) but with 3 typos corrected:
         RPAR=(1.d0+(26.33d0*VTp2+.376d0*VTp2**2)*FUNCXE(3.675d0,VTp)+
     +     .742d0*(FUNCXE(1.673d0,VTp)-1.d0))*FUNCXE(1.361d0,VTp)
      endif
* \nu_{ci}^\perp: SY Eq.(34):
      FTEE=CnuTRAN*TEMR*XSRE**3/QT2/GSRE ! e-e
      FTEP=FTEE*(CKFPR/XSRE)**2 ! e-p
* \nu_{ci}^\|: SY Eq.(35):
      FLEE=CnuLONG*TEMR**2*GSRE**3/XSRE/QL2/dsqrt(QL2) ! e-e(||):SY(35)
      FLEP=FLEE*(GSRP/GSRE)**2 ! e-p(||):SY(35)
      if (T9.lt.Tcp9) then ! apply proton SF reduction factors
* SY (86), (89) + approximate all R^\perp by R_{tot}^\perp:
         FTEE=FTEE*RTOT
         FTEP=FTEP*RTOT
         FLEP=FLEP*RPAR
      endif
* SY Eq.(33):
      FEE=FTEE+FLEE
      FEP=FTEP+FLEP
      if (Ymu.lt.TINY) then ! no muons
         FREQe=FEP+FEE
         TAUe=1.d0/FREQe
         RKAPe=TEMR*XSRE**3/9.d0/GSRE*TAUe ! el.th.cond. in rel.un.
      else
        if (Ymu.gt.Ye) stop'CONDCOR: Ymu > Ye - inapplicable case'
         FTEM=FTEE*(CKFM/CKFE)**2 ! e-\mu collisions: SY(34)
         FLEM=FLEE*(RATMUE*GSRM/GSRE)**2 ! SY(35)
         XGe2mu=XSRM/XSRE/GSRM*GSRE ! change carriers from e to \mu
         FTME=FTEE*XGe2mu ! \mu-e(\perp): SY(34)*RTOT
         FLME=FLEE/XGe2mu ! \mu-e(||): SY(35)
         FTMP=FTEP*XGe2mu ! \mu-p(\perp): SY(34)*RTOT
         FLMP=FLEP/XGe2mu ! \mu-e(||): SY(35)
         F1EM=CnuPRIM*TEMR*(TEMR**2/QT2)**C13*GSRE*CKFMR**2/XSRE/QL2*
     *     RLTL ! \nu'_{e\mu}: SY Eqs.(36), (58)
         F1ME=F1EM/XGe2mu*(CKFE/CKFM)**2 ! \nu'_{\mu e}: SY Eq.(36)*(58)
         FTMM=FTEM*XGe2mu ! \mu-\mu
         FLMM=FLEM/XGe2mu
         FEM=FTEM+FLEM ! SY Eq.(33)
         FMM=FTMM+FLMM
         FMP=FTMP+FLMP
         FME=FTME+FLME
         FREQmu=FMP+FME+FMM ! SY Eq.(6)
         FREQe=FEP+FEE+FEM
         DENOM=FREQe*FREQmu-F1ME*F1EM
         TAUe=(FREQmu-F1EM)/DENOM
         TAUmu=(FREQe-F1ME)/DENOM
         RKAPe=TEMR*XSRE**3/9.d0/GSRE*TAUe ! e th.cond.\kappa_e[rel.un.]
         RKAPmu=TEMR*CKFMR**2*XSRM/9.d0/GSRM*TAUmu ! \kappa_\mu
         CKAPmu=RKAPmu*UNIKAP ! \kappa_\mu in erg/(K cm s)
      endif
      CKAPe=RKAPe*UNIKAP ! \kappa_e in erg/(K cm s)
      return
      end
