* Collection of subroutines for calculation of neutrino emission rates
*   in the neutron star core (without superfluidity).
* Basic references: Yakovlev et al. (2001) [Phys.Rep. 354, 1-155] (YKGH)
*               and Gusakov (2002)[A&A 389, 702-715]
* Uses: sfreduc.f gusakov.f 
* Last change: 08.10.15
      subroutine CORNSF(XN,Yn,Yp,Ye,Ymu,T9,B12,EFMn,EFMp,
     *  QDUe,QDUmu,QMUne,QMUpe,QMUnmu,QMUpmu,QBnn,QBnp,QBpp)
*                                                       Version 26.06.15
* Main neutrino emissivities from the core without superfluidity
* Equation and section numbers refer to Yakovlev et al. (2001) PR 354, 1
* Input: XN - baryon number density in fm^{-3}
*        Yn - neutron number fraction
*        Yp - proton number fraction
*        Ye - electron number fraction
*        Ymu - muon number fraction
*        T9 = temperature / 10^{9} K
*        B12 - magnetic field [TG]
*        EFMp,EFMn - effective proton and neutron mass factors (m*/m)
* Output - neutrino emission rates in erg/(s cm^3):
*        QDUe - electron Durca,
*        QDUmu - muon Durca,
*        QMUne - Murca, neutron branch with electrons,
*        QMUpe - Murca, proton branch with electrons,
*        QMUnmu - Murca, neutron branch with muons,
*        QMUpmu - Murca, proton branch with muons,
*        QBnn, QBnp, QBpp - nucleon (nn, np, pp) bremsstrahlung.
      implicit double precision (A-H), double precision (O-Z)
      parameter (nNu=3) ! number of neutrino flavors
* constants in Eqs.(165)-(167):
      parameter (alphaNN=0.59d0,alphaNP=1.06d0,alphaPP = 0.11d0)
      parameter (betaNN=0.56d0,betaNP=0.66d0,betaPP=0.7d0)
      parameter (CABNN=alphaNN*betaNN*nNU,
     *  CABNP=alphaNP*betaNP*nNU,CABPP=alphaPP*betaPP*nNU)
      include 'const.inc'
      Ye13=Ye**C13
      Ymu13=Ymu**C13
      Yp13=Yp**C13
      Yn13=Yn**C13
      XN13=(XN/.16d0)**C13 ! (n_b/0.16 fm^{-3})^{1/3}
* Neutrino emission due to nucleon Durca process: Eq.(120)
      QDUCOM=4.d27*XN13*EFMn*EFMp*PNMRAT*T9**6 ! common factor
      call RDUMAG(XN,Yn,Yp,Ye,B12,RBDUe)
      QDUe=QDUCOM*Ye13*RBDUe
      call QMURCA(XN,Yn,Yp,Ye,T9,EFMn,EFMp,QMUne,QMUpe)
* Muon branches: Sect.3.4.4
      QMUnmu=0.
      QMUpmu=0.
      QDUmu=0.
      if (Ymu.gt.0.) then
         call QMURCA(XN,Yn,Yp,Ymu,T9,EFMn,EFMp,QMUnmu,QMUpmu)
* Velocity correction (electron velocity is assumed equal to c):
         VFmu=Ymu13/Ye13 ! v_{F,\mu}/c in non-relativistic approximation
         VFmu=VFmu/dsqrt(1.d0+VFmu**2) ! allow for special relativity
         QMUnmu=QMUnmu*VFmu
         call RDUMAG(XN,Yn,Yp,Ymu,B12,RBDUmu)
         QDUmu=QDUCOM*Ymu13*RBDUmu
      endif
* Neutrino emission due to baryon-baryon collisions: Eqs.(165)-(167)
      QBCOM=7.5d19*XN13*T9**8 ! common factor
      QBnn=EFMn**4*Yn13*CABNN*QBCOM ! neutron-neutron brems.
      QBnp=2.d0*(EFMn*EFMp)**2*Yp13*CABNP*QBCOM ! neutron-proton brems.
      QBpp=EFMp**4*Yp13*CABPP*QBCOM ! proton-proton brems.
      return
      end

      subroutine QMURCA(XN,Yn,Yp,Ye,T9,EFMn,EFMp,QMUne,QMUpe)
*                                                       Version 26.06.15
* Neutrino emission due to nucleon Murca process: electron branches.
* Input: XN - baryon number density in fm^{-3}
*        Yn - neutron number fraction
*        Yp - proton number fraction
*        Ye - electron number fraction
*        EFMp,EFMn - effective proton and neutron mass factors (m*/m)
* NB: ASSUMED Yp < Yn AND Ye < Yn.
      implicit double precision (A-H), double precision (O-Z)
      include 'const.inc'
      parameter (CAB=8.1d21*1.13d0*0.68d0) ! constant in Eq.(140)
      Ye13=Ye**C13
      Yp13=Yp**C13
      Yn13=Yn**C13
      XN13=(XN/.16d0)**C13 ! (n_b/0.16 fm^{-3})^{1/3}
      COMMUP=CAB*EFMn*EFMp*XN13*T9**8 ! common factor for two branches
      QMUne=COMMUP*EFMn**2*Yp13 ! n-branch: YKGH:Eq.(140)
* Proton branch:
      DY=3.d0*Yp13-Yn13
      if (DY.gt.Ye13) then ! YKGH:Eq.(142) + Gusakov(2002):Eq.(15)
         QMUpe=COMMUP*EFMp**2*DY/2.d0
      elseif (DY+Ye13.gt.0.) then ! YKGH:Eq.(142)
         QMUpe=COMMUP*EFMp**2*(DY+Ye13)**2/(8.0d0*Ye13)
      else
         QMUpe=0.
      endif
      return
      end

      subroutine COOPER(T9,Tc9,KSF,EFM,XSRx,QPBF)
*                                                       Version 23.06.15
* Neutrino emission due to Cooper pairing of neutrons in NS core:
*   YKGH(2001), Eqs.(236)-(241), with corrections of Leinson (2009,2010)
* Input: T9 - temperature [GK],
*        Tc9 - critical temperature of superfluidity [GK]
*        KSF - type of superfluidity: 1 - ns, 2 - ps, 3 - nt
*               (n/p stands for neutron/proton, s/t for singlet/triplet)
*        EFM=m*/m - effective mass factor
*        XSRx=p_F/mc - relativity factor (Fermi momentum in units of mc)
* Output: QPBF - PBF (pair breaking-formation) emissivity [erg/(cm^3 s)]
      implicit double precision (A-H), double precision (O-Z)
      dimension LJSF(3)
      save
      parameter (nNu=3) ! number of neutrino flavors
      parameter (QPREF=1.17d21*nNu) ! prefactor in Eq.(236)
      data LJSF/1,1,2/ ! 1S0 (ns), 1S0 (ps), 3P2(m_J=0) (nt)
      if (T9.ge.Tc9) then
         QPBF = 0.0d0
        return
      endif
      JSF=LJSF(KSF) ! SF type: 1,2,3 => A(1S0),B(3P2,mJ=0),C(3P2,mJ=2)
      call DLGAPS(JSF,T9/Tc9,VT)
      VT2=VT**2
* Calculate F(v): Eq.(241)
      if (JSF.eq.1) then ! type A superfluidity
         FV=(.602d0+.5942d0*VT2+.288d0*VT2**2)*VT2*
     *      dsqrt(.5547d0+dsqrt(.4453**2+.0113*VT2))*
     *      dexp(2.245d0-dsqrt(4.d0*VT2+2.245**2))
      elseif (JSF.eq.2) then ! type B superfluidity
         FV=(1.204d0+3.733d0*VT2+.3191d0*VT2**2)*VT2/(1.d0+.3511d0*VT2)*
     *      (.7591d0 + dsqrt(.2409d0**2+.3145d0*VT2))*
     *      dexp(.4616d0-dsqrt(4.d0*VT2+.4616d0**2))  
      elseif (JSF.eq.3) then ! type C superfluidity
         FV=(.4013-.043*VT2+.002172*VT2**2)*VT2/
     /      (1.d0-.2018*VT2+.02601*VT2**2-.001477*VT2**3+4.34d-5*VT2**4)
      else
         print*,'JSF =',JSF
         stop'COOPER: unknown type of superfluidity.'
      endif
* Calculate QPBF: Eqs.(236), (238)
      if (KSF.eq.1) then ! singlet neutron
         AN=1.d0
      elseif (KSF.eq.2) then ! singlet proton
         AN=.0064+1.59*(XSRx/EFM)**2*(EFM+.2619) ! Eq.(238), .2619=11/42
      elseif (KSF.eq.3) then ! triplet neutron
         AN=4.17
      else
         print*,'KSF =',KSF
         stop'COOPER: unknown superfluidity type.'
      endif
      QPBF=QPREF*EFM*XSRx*T9**7*AN*FV ! Eq.(236)
* Leinson corrections:
      if (JSF.eq.1) then ! singlet
         QPBF=QPBF*(XSRx/EFM)**2
      else ! triplet
         QPBF=QPBF*.19
      endif
      return
      end

      subroutine RDUMAG(XN,Yn,Yp,Ye,B12,RBDU)
*                                                       Version 26.06.15
*                        small update to conform with const.inc 08.10.15
* Input: XN - baryon number density [fm^{-3}]
*        Yn - neutron number fraction
*        Yp - proton number fraction
*        Ye - electron number fraction
*        B12 - magnetic field [TG]
* Output: RBDU - magnetic switch factor for Durca: YKGH:Eqs.(248)-(250)
      implicit double precision (A-H), double precision (O-Z)
      include 'const.inc'
      parameter(TPI2=3.d0*PI*PI,TINY=1.d-19) ! C23 excluded 08.10.15
      Ye13=Ye**C13
      Yp13=Yp**C13
      Yn13=Yn**C13
      B=B12/UNIB ! m.f.[rel.un.]=B/B_{QED}=\hbar\omega_c/m_e c^2
      PFp=COMPTON*(TPI2*XN*Yp)**C13 ! p_{Fp}[rel.un.]=p_{Fp}/m_e c
      Y=PFp**2/2.d0/dmax1(B,TINY) ! nonrel.estimate of max.Landau num.
      Y23=Y**C23
      X=(Yn13**2-(Yp13+Ye13)**2)/dmax1(Yp13**2,TINY)*Y23
      if (X.lt.-20.d0) then ! Durca is strongly open
         RBDU=1.d0
      elseif (X.lt.0.) then ! Durca is open and oscillates
         XA=dabs(X)
         PHI=(1.211+.4823*XA+.8453*XA**2.533)/(1.d0+1.438*XA**1.209)
         RBDU=1.d0-dcos(PHI)/(.5816+XA)**1.192
      elseif (X.lt.1.d2) then ! smeared switch-off of Durca
         XC=X*(dsqrt(X+.4176)-.04035)
         RBDU=(3.d0*X+6.8d0)/(XC+6.8d0)/(3.d0+3.4641*X)*dexp(-XC/3.d0)
      else ! Durca is strictly forbidden
         RBDU=0.
      endif
      return
      end

      subroutine FLUXOID(XN,Yp,Ye,T9,Tcp9,B12,EFMp,QFLU)
*                                                       Version 26.06.15
* Input: XN - baryon number density in fm^{-3}
*        Yp - proton number fraction
*        Ye - electron number fraction
*        T9 - temperature [GK],
*        Tcp9 - critical temperature of proton superfluidity [GK]
*        B12 - initial uniform magnetic field [TG]
*        EFMp=m*_p/m_p - effective proton mass factor
* Output: QFLU - electron-fluxoid scattering emissivity [erg/(cm^3 s)]
*   according to YKGH, Sect.
      implicit double precision (A-H), double precision (O-Z)
      include 'const.inc'
      parameter(C136=13.d0/6.d0,QuarPI=PI/4.d0)
      if (T9.ge.Tcp9) then
         QFLU=0.
        return
      endif
      XN1=XN/.16d0 ! n_b/n_0
      ALPT=dsqrt(1.d0-(T9/Tcp9)**4) ! alpha_T, YKGH:(253)
      V=1.d0+.00021*B12 ! YKGH:(268)
      G=8.38d-14*T9*B12**2/(XN1*Ye)**C136*(dsqrt(EFMp)/ALPT)**3!gamma
      U=(2.d0*G+1.d0)/(3.d0*G+1.d0) ! YKGH:(266)
      SQRNPM=dsqrt(XN1*Yp/EFMp) ! \sqrt{n_p/n_0*m_p/m*_p}
      T0=.00786*T9/SQRNPM/ALPT ! t_0 after Eq.(265)
      CL0=QuarPI*(.26*T0+.0133)/(T0**2+.25*T0+.0133) ! L_0
      CL=CL0*U*V
      QFLU=2.66d15*B12*T9**6*SQRNPM*ALPT*CL ! YKGH:(263)
      return
      end
