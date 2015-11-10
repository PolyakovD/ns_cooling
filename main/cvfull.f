** Specific heat components for the crust and ocean of a neutron star
* Needs to be linked with: eosmag14.f eos14.f compos21.f bskfit.f
* Last update: 26.01.15
      subroutine CVFULL21(RHO,T6,B12,
     *  XN,Zion,CMI,CMI1,Yn,Yp,Ye,Ymu,EFMp,EFMn,CVtot,CVE,CVI,CVC,CVN)
*                                                       Version 24.01.15
* Composition and heat capacity of neutron-star matter
*   according to BSk21 EoS.
* Input: RHO - mass density in g/cc
*        T6 - temperature in MK
*        B12 - magnetic field in TG
* Output: XN = n_b - baryon number density in fm^{-3}
*         Zion - charge number of a nucleus
*         CMI - effective atomic mass number of a nucleus
*         CMI1 - total number of nucleons per a nucleus
*         Yn,Yp,Ye,Ymu - numbers of free neutrons, protons, electrons,
*                and muons, respectively, per a nucleon
*         EFMp,EFMn - effective proton and neutron mass factors(mP*/mP)
*				,(mN*/mN), respectively 
*         CVtot - total heat capacity per a nucleon
*         CVE,CVI,CVC,CVN - electron, ion (nucleus, proton), Coulomb,
*                and neutron heat capacity parts, resp., per a nucleon
      implicit double precision (A-H), double precision (O-Z)
      parameter(KEOS=21,XNC=.0809)
*         XNC - number density of baryons in fm^{-3} at the crust/core
*              interface (Pearson et al. 2012),
      call BSkNofR(KEOS,RHO,XN)
      if (XN.lt.XNC) then ! crust
         call CRUBSk21(XN,KZ,CMI,CMI1,xnuc,XND)
         Zion=KZ
         Yn=dim(1.d0,CMI/CMI1) ! number of free neutrons per nucleon
         XNN=XN*Yn ! number density of free neutrons, fm^-3
         Yp=0. ! no free protons
         Ye=Zion/CMI1 ! number of electrons per nucleon
         Ymu=0. ! no muons
         call EFFMASS(KEOS,XNN,Yp,EFMp,EFMn)
         call CVCRU21(XN,T6,B12,Zion,CMI,CMI1,EFMn,
     *     CVtot,CVE,CVI,CVC,CVN)
      else ! core
         call FRACORE(KEOS,XN,Ye,Ymu)
         Yp=Ye+Ymu
         Yn=dim(1.d0,Yp)
         call EFFMASS(KEOS,XN,Yp,EFMp,EFMn)
         CMI=1.008*EFMp ! eff.mass of a proton in a.m.u.
         CMI1=1.d0/Yp
         Zion=1.d0 ! proton charge
         CVC=0. ! neglect Coulomb nonideality
         call CVCOR21(XN,T6,B12,Ye,Ymu,EFMp,EFMn,
     *     CVtot,CVN,CVI,CVE,CVmu)
      endif
      return
      end

      subroutine CVCRU21(XN,T6,B12,Zion,CMI,CMI1,EFMn,
     *  CVtot,CVE,CVI,CVC,CVN)
*                                                       Version 24.01.15
* Input:
*        XN = n_b - baryon number density in fm^{-3}
*        T6 - temperature in MK,
*        B12 - magnetic field in TG,
*        Zion and CMI - nuclear charge and mass numbers,
*        CMI1 - total number of baryons per nucleus,
*        EFMn - effective mass of a neutron relative to bare mass
* Output:
*         CVtot - heat capacity per nucleon, div. by Boltzmann const. kB
*         CVE - ideal electron gas component
*         CVI - ideal ion + nonideal ion and ion-electron components
*         CVN - non-superfluid neutron gas component
      implicit double precision (A-H), double precision (O-Z)
      save
      parameter (PI=3.141592653d0,C13=1.d0/3.d0,TINY=1.d-20,
     *  AUM=1822.888d0, ! a.m.u./m_e
     *  RMN=1838.68366d0, ! neutron-electron mass ratio 
     *  BOHRFM=.52917720859d5, ! Bohr radius in fm
     *  UNIB=44.14004916, ! rel.un.of magnetic field [TG]
     *  UN_T6=.3157746d0, ! = Hartree energy/10^6 K = a.u. of T in MK
     *  FINV=137.036, ! inverse fine-structure constant
     *  GAMIMELT=175., ! OCP value of Gamma_i for melting
     *  RSIMELT=140.) ! ion density parameter of quantum melting
      TEMP=T6/UN_T6 ! temperature in a.u.
      DENSB=XN*BOHRFM**3 ! baryon number density in a.u.
      TEMR=TEMP/FINV**2 ! temperature in relativistic units
* (1) ideal electron gas (including relativity and degeneracy)
      DENSI=DENSB/CMI1 ! number density of ions ( = nuclei) [au]
      DENS=DENSI*Zion ! number density of electrons [au]
      DENR=DENS/FINV**3 ! number density of electrons in relativ. units
      B=B12/UNIB ! magnetic field in relativistic units
      call CHEMAG8(B,TEMR,DENR,
     *  DENR1,CHI,FEid,PEid,UEid,SEid,CVE,CHITE,CHIRE,DDENR)
* (2) ideal ion gas (including relativity, degeneracy, and m.f.):
      GYRO=0.d0
      MULTI=1
      if (dabs(Zion-1.d0).lt.TINY) then ! assume protons
         GYRO=5.5857d0 ! = g_p
         MULTI=2
      endif
      RS=(.75d0/PI/DENS)**C13 ! r_s - electron density parameter
      Z13=Zion**C13
      RSI=RS*CMI*Z13**7*AUM ! R_S - ion density parameter
      GAME=1.d0/(RS*TEMP) ! electron Coulomb parameter Gamma_e
      GAMI=Z13**5*GAME   ! effective Gamma_i - ion Coulomb parameter
      if (GAMI.lt.GAMIMELT.or.RSI.lt.RSIMELT) then
         LIQSOL=0 ! liquid regime
         call MAGNION(B,TEMR,DENRI,Zion,CMI,GYRO,MULTI,
     *     FIid,PIid,UIid,SIid,CVI,CHITI,CHIRI) ! ideal-gas + spin parts
      else
         LIQSOL=1 ! solid regime
         ZETI=Zion/(CMI*AUM)*GAMAG/TEMP ! \hbar\omega_{ci}/kT
         call SPINION(GYRO,ZETI,MULTI,Fspin,Uspin,CVspin)
         CVI=1.5d0+CVspin ! 1.5 (formally, non-magn.id.gas) + spin part
      endif
      GAMAG=B*FINV**2 ! magnetic field in a.u.
* (3) Coulomb+xc (nonideal) contributions
      call EOSFIM14(LIQSOL,CMI,Zion,RS,TEMP,GAMAG,0,0,
     *  FC1,UC1,PC1,SC1,CV1,PDT1,PDR1,
     *  FC2,UC2,PC2,SC2,CV2,PDT2,PDR2)
      CVI=CVI/CMI1 ! ideal-gas + spin part of CV / 1 nucleon
      CVC=CV1/CMI1 ! nonideal part of CV / 1 nucleon
      CVtot=CVE*Zion/CMI1+CVI+CVC ! C_V (e+i+Coul.) / 1 nucleon
* (3) neutrons
      if (CMI1.gt.CMI) then ! there are free neutrons
         RMNe=RMN*EFMn ! effective neutron-electron mass ratio
         TEMRn=TEMR/RMNe ! scaling of T for neutrons
         DENRn=DENSB*(1.d0-CMI/CMI1)/(RMNe*FINV)**3 ! scaled n-density
         call CHEMAG8(0.d0,TEMRn,DENRn,
     *     DENR1n,CHIn,FNid,PNid,UNid,SNid,CVN,CHITN,CHIRN,DDENRn)
         CVN=CVN*(1.d0-CMI/CMI1) ! neutron heat capacity per 1 baryon
         CVtot=CVtot+CVN
      else
         CVN=0.
      endif
      return
      end

      subroutine CVCOR21(XN,T6,B12,Ye,Ymu,EFMp,EFMn,
     *     CVtot,CVN,CVP,CVE,CVmu)
*                                                       Version 24.01.15
* Input:
*        XN = n_b - baryon number density in fm^{-3}
*        T6 - temperature in MK,
*        B12 - magnetic field in TG,
*        Zion and CMI - nuclear charge and mass numbers,
*        CMI1 - total number of baryons per nucleus,
*        EFMn - effective mass of a neutron relative to bare mass
* Output:
*         CVtot - heat capacity per nucleon, div. by Boltzmann const. kB
*         CVN - non-superfluid neutron gas component
*         CVP - non-superfluid proton gas component
*         CVE - ideal electron gas component
*         CVmu - ideal muon gas component
      implicit double precision (A-H), double precision (O-Z)
      save
      parameter (RMN=1838.68366d0, ! neutron-electron mass ratio 
     *  RMP=1836.15267d0, ! proton-electron mass ratio 
     *  RMmu=206.768284d0, ! muon-electron mass ratio 
     *  BOHRFM=.52917720859d5, ! Bohr radius in fm
     *  UN_T6=.3157746d0, ! = Hartree energy/10^6 K = a.u. of T in MK
     *  FINV=137.036, ! inverse fine-structure constant
     *  TINY=1.d-20)
      Yp=Ye+Ymu ! proton fraction
      Yn=dim(1.d0,Yp)
      TEMR=T6/UN_T6/FINV**2 ! temperature in relativistic units
      DENRB=XN*(BOHRFM/FINV)**3 ! number density of baryons in rel. un.
* neutrons:
      RMNe=RMN*EFMn ! effective neutron-electron mass ratio
      TEMRn=TEMR/RMNe ! scaling of T for neutrons
      DENRn=DENRB*Yn/RMNe**3 ! scaled neutron density
      call CHEMAG8(0.d0,TEMRn,DENRn,
     *  DENR1n,CHIn,FNid,PNid,UNid,SNid,CVN,CHITN,CHIRN,DDENRn)
      CVN=CVN*Yn ! normalize per baryon
* protons:
      RMPe=RMP*EFMp ! effective proton-electron mass ratio
      TEMRp=TEMR/RMPe ! scaling of T for protons
      DENRp=DENRB*Yp/RMPe**3 ! scaled proton density
      call CHEMAG8(0.d0,TEMRp,DENRp,
     *  DENR1p,CHIp,FPid,PPid,UPid,SPid,CVP,CHITP,CHIRP,DDENRp)
      CVP=CVP*Yp ! normalize per baryon
* electrons:
      DENR=DENRB*Ye ! electron density in rel.un.
      B=B12/44.14 ! magnetic field in rel.un.
      call CHEMAG8(B,TEMR,DENR,
     *  DENR1,CHI,FEid,PEid,UEid,SEid,CVE,CHITE,CHIRE,DDENR)
      CVE=CVE*Ye ! normalize per baryon
* muons:
      CVmu=0.
      if (Ymu.gt.TINY) then
         TEMRmu=TEMR/RMmu ! scaling of T for muons
         DENRmu=DENRB*Ymu/RMmu**3 ! scaled muon density in rel.un.
         Bmu=B/RMmu**2 ! magnetic field in rel.un. scaled for muons
         call CHEMAG8(Bmu,TEMRmu,DENRmu,
     *     DENR1mu,CHImu,Fmuid,Pmuid,Umuid,Smuid,CVmu,CHITmu,CHIRmu,
     *     DDENRmu)
         CVmu=CVmu*Ymu
      endif
* total:
      CVtot=CVN+CVP+CVE+CVmu ! total heat capacity per baryon
      return
      end
