* Collection of subroutines for calculation of superfluid reduction
*   factors to neutrino emission rates.
* Basic references: Yakovlev et al. (2001) [Phys.Rep. 354, 1-155] (YKGH)
*               and Gusakov (2002)[A&A 389, 702-715]
* USES: gusakov.f
* Last change: 02.09.15
* Contains:
*   SFRED - Collection of SF red.factors for a given SF type
*   RBrPPA - Reduction factor for pp- or nn-brems., SF type A
*   RBrNNC - Reduction factor for pp- or nn-brems., SF type C
*   RBrNPA - for np-brems. with non-SF neutrons and SF type A protons
*   REDUP1 - Durca SF reduction factor for a single kind of SF particles
*   REDUP2 - Durca SF reduction factor when both nucleons are superfluid
*   REDCV - Reduction factor for specific heat
*   DLGAPS - SF dim.less gap \Delta(T)/T for a given type of SF and T/Tc
*   CRITEM - "gap function" and crit.temperature Tc for a SF
*   CRITNP - neutron and proton crit. temperatures for a set of 3 SF's
      subroutine SFRED(T9,Tcp9,Tcn9,JSF,RD,RMn,RMp,Rnn,Rnp,Rpp)
*                                                       Version 16.06.15
* Superfluidity reduction factors for neutrino emissivities
* References: Gusakov (2002) [A&A 389, 702-715]; 
*             Yakovlev et al. (2001) [Phys.Rep. 354, 1-155] (YKGH).
* Input: T9 - temperature [GK]
*        Tcp9 - SF critical temperature T_c for protons
*        Tcn9 - T_c for neutrons
*        JSF=1,2,3 for neutron pairing type:
*               A (1S0), B (3P2, m_J=0), C (3P2 ,m_J=2).
* Output: RD - reduction factor for Durca process
*         RMn - reduction factor for neutron branch of Murca process
*         RMp - reduction factor for proton branch of Murca process
*         Rnn - reduction factor for strong collisions of SF neutrons
*         Rnp - SF reduction factor for strong neutron-proton collisions
*         Rpp - reduction factor for strong collisions of SF protons
      implicit double precision (A-H), double precision (O-Z)
      parameter (TINY=1.d-19)
      if (Tcp9.le.T9.and.Tcn9.le.T9) then ! no SF
         RD=1.d0
         RMn=1.d0
         RMp=1.d0
         Rnn=1.d0
         Rpp=1.d0
         Rnp=1.d0
        return
      endif
      VTn=0.
      if (Tcn9.ge.T9) call DLGAPS(JSF,T9/Tcn9,VTn)
      VTp=0.
      if (Tcp9.ge.T9) call DLGAPS(1,T9/Tcp9,VTp) ! proton SF type A only
      if (Tcp9.lt.T9) then ! no proton SF, only neutron SF
         call REDUP1(JSF,VTn,RD)
         call RBrPPA(VTn,Rnn) ! Eq.(228)
         Rpp=1.d0
         call RBrNPA(VTn,Rnp)
      elseif (Tcn9.le.T9) then ! no neutron SF, only proton SF
         call DLGAPS(1,T9/Tcp9,VTp) ! proton SF type A only
         call REDUP1(1,VTp,RD)
         call RBrPPA(VTp,Rpp)
         Rnn=1.d0
         call RBrNPA(VTp,Rnp)
      else ! both neutron and proton SF
         call DLGAPS(JSF,T9/Tcn9,VTn)
         call DLGAPS(1,T9/Tcp9,VTp)
         call REDUP2(JSF,VTn,VTp,RD) ! Durca reduction factor RD
         call RBrPPA(VTp,Rpp)
        if (JSF.eq.3) then ! C-type SF
           call RBrNNC(VTn,Rnn) ! Gusakov (2002) Eq.(60)
        else ! use factor for A-type SF
           call RBrPPA(VTn,Rnn) ! YKGH:Eq.(228)
        endif
         call RBrNPA(VTp,Rnp)
         call REDUP1(1,VTp,R1)
        if (R1.lt.TINY) R1=TINY ! to prevent machine zero
         Rnp=RD/R1*Rnp ! YKGH:Eq.(229) [or Eq.(54) of Gusakov (2002)]
      endif
      call REMUP2(JSF,VTn,VTp,RMn,RMp)
      return
      end

      subroutine RBrPPA(VT,Rpp)
*                                                       Version 11.06.15
* Reduction factor for neutrino emissivity in strong collisions of
*   singlet-superfluid protons: YKGH(220)
* Input: VT - dim.-less gap for proton superfluidity type A
* Output: Rpp - reduction factor for pp-bremsstrahlung
      implicit double precision (A-H), double precision (O-Z)
      C=.1747d0+dsqrt(.8253d0**2+(.007933d0*VT)**2)
      D=.7333d0+dsqrt(.2667d0**2+(.1678d0*VT)**2)
      Rpp=(C**2*dexp(4.228d0-dsqrt(4.228d0**2+4.d0*VT**2))+
     +  D**7*dsqrt(D)*dexp(7.762d0-dsqrt(7.762d0**2+9.d0*VT**2)))/2.0d0
      return
      end

      subroutine RBrNNC(VT,Rnn)
*                                                       Version 16.06.15
* Reduction factor for neutrino emissivity in strong collisions of
*   superfluid neutrons in state 3P2(m_J=2): Gusakov (2002) Eq.(60)
* Input: VT - dim.-less gap for neuton superfluidity type C
* Output: Rnn - reduction factor for nn-bremsstrahlung
      implicit double precision (A-H), double precision (O-Z)
      parameter (P1=.442995d0,P2=.250953d0,P3=.410517d0,P4=7.09171d-4)
      VT2=VT**2
      Rnn=dexp(-P1*VT2/(1.d0+P2*VT2+P4*VT2**2)**P3)
      return
      end

      subroutine RBrNPA(VT,Rnp)
*                                                       Version 16.06.15
* Reduction factor for neutrino emissivity in collisions of 
*   non-superfluid neutrons with singlet-superfluid protons: YKGH(221)
* Input: VT - dim.-less gap for proton superfluidity type A
* Output: Rnp - reduction factor for np-bremsstrahlung
      implicit double precision (A-H), double precision (O-Z)
      A=.9982d0+dsqrt(.0018d0**2+(.3815d0*VT)**2)
      B=.3949d0+dsqrt(.6051d0**2+(.2666d0*VT)**2)
      Rnp=(A*dexp(1.306d0-dsqrt(1.306d0**2+VT**2))+
     +  1.732d0*B**7*dexp(3.303d0-dsqrt(3.303d0**2+4.d0*VT**2)))/
     /  2.732d0
      return
      end

      subroutine REDUP1(JSF,VT,RD)
*                                                       Version 09.06.15
* Durca SF reduction factor for a single kind of SF particles: YKGH(199)
* Input: JSF=1, 2 or 3 for states 1S0, 3P2(m_J=0), 3P2(m_J=2), resp.
*        VT=\Delta/T - dimensionless gap factor
* Output: RD - reduction factor for Durca emissivity
      implicit double precision (A-H), double precision (O-Z)
      parameter(TINY=1.d-19)
      if (VT.lt.TINY) then ! no SF
         RD=1.d0
        return
      endif
      if (JSF.eq.1) then ! type A
         RD=(.2312d0+dsqrt(.7688d0**2+(.1438d0*VT)**2))**5.5d0*
     *     dexp(3.427d0-dsqrt(3.427d0**2+VT**2)) ! R_A^D
      elseif (JSF.eq.2) then ! type B
         RD=(.2546+dsqrt(.7454d0**2+(0.1284d0*VT)**2))**5*
     *     dexp(2.701d0-dsqrt(2.701**2+VT**2)) ! R_B^D
      elseif (JSF.eq.3) then ! type C
         RD=(.5d0+(.0922*VT)**2)/(1.d0+(.1821*VT)**2+(.16736*VT)**4)+
     +     .5d0*dexp(1.d0-dsqrt(1.d0+(.4129*VT)**2)) ! R_C^D
      else
          print*,'JSF =',JSF
         stop'REDUP1: JSF out of range'
      endif
      return
      end

      subroutine REDUP2(JSF,VTn,VTp,RD)
*                                                       Version 11.06.15
* reduction factor for Durca process when both nucleons are superfluid
* Input: JSF=1, 2 or 3 for [nn] states 1S0, 3P2(m_J=0), 3P2(m_J=2)
*      VTn=T/Tcn - ratio of neutron SF gap to temperature
*      VTp=T/Tcp - ratio of proton SF gap to temperature
* Output: RD - Durca SF reduction factor
      implicit double precision (A-H), double precision (O-Z)
      include 'const.inc'
      parameter(PI2=PI*PI)
      parameter(PREK2=PI2*PI2*7.d0/60.d0) ! used in K_2
      parameter(CI0=457.d0*PI2*PI2*PI2/5040.d0) ! I_0 (page 92 of YKGH)
      VTn2=VTn**2
      VTp2=VTp**2
      if (JSF.eq.1) then ! SF type AA
         U=VTn2+VTp2
        if (U.gt.2.5d3) then ! strong SF
           call REDUP1(JSF,VTn,RD1)
           call REDUP1(1,VTp,RD2)
           RD=dmin1(RD1,RD2) ! Eq.(204)
* Modification of AYP 11.06.15
           RDmax=dmax1(RD1,RD2)
           RD=RD*RDmax/(dsqrt(RD*dim(1.d0,RDmax))+dsqrt(RDmax))**2
        else ! moderate SF: Eqs. (202), (203)
           W=VTp2-VTn2
           W2=W**2
           U1=1.8091d0+dsqrt(VTn2+2.2476d0**2)
           U2=1.8091d0+dsqrt(VTp2+2.2476d0**2)
           D=1.52*dsqrt(U1*U2)**3*(U1**2+U2**2)*dexp(-U1-U2)
           PE=(U+.43847d0+dsqrt(W2+8.368d0*U+491.32d0))/2.d0 ! p_e
           PS=(U+dsqrt(W2+5524.8d0*U+6.7737d0))/2.d0 ! p_s
           PR=U+12.421d0
           PQ=dsqrt(W2+16.35d0*U+45.171d0)
           P=(PR+PQ)/2.d0
           Q=(PR-PQ)/2.d0
           PSQ=dsqrt(P) ! \sqrt{p}
           PQSQ=dsqrt(PQ) ! \sqrt{p-q}
* ! \ln[(\sqrt{p}+\sqrt{p-q})/\sqrt{q}]:
           PQLN=dlog((PSQ+PQSQ)/dsqrt(Q))
* K_0, K_1, K_2; S:
           CK0=PQSQ/120.d0*(6.d0*P**2+83.d0*P*Q+16.d0*Q**2)-
     -       PQSQ*Q/8.d0*(4.d0*P+3.d0*Q)*PQLN
           CK1=PI2/6.d0*PQSQ*(P+2.d0*Q)-PI2/2.d0*Q*PSQ*PQLN
           CK2=PREK2*PQSQ
           S=(CK0+CK1+.42232*CK2)/CI0*dsqrt(PI/2.d0)*
     *       dsqrt(dsqrt(PS))*dexp(-dsqrt(PE))
           RD=U/(U+.9163d0)*S+D
        endif
      elseif (JSF.eq.2) then
        if (VTn2+VTp2.gt.25.d0) then ! strong SF
           call REDUP1(JSF,VTn,RD1)
           call REDUP1(1,VTp,RD2)
           RD=dmin1(RD1,RD2) ! Eq.(204)
* Modification of AYP 11.06.15
           RDmax=dmax1(RD1,RD2)
           RD=RD*RDmax/(dsqrt(RD*dim(1.d0,RDmax))+dsqrt(RDmax))**2
        else ! moderate SF: Eq.(205)
           RD=(1.d4-2.839d0*VTp2**2-5.022d0*VTn2**2)/
     /       (1.0d4+757.0d0*VTp2+1494.d0*VTn2+
     +          211.1d0*VTp2*VTn2+.4832d0*(VTp2*VTn2)**2)
        endif
      elseif (JSF.eq.3) then
        if (VTn2+VTp2.gt.1.d2) then ! strong SF
           call REDUP1(JSF,VTn,RD1)
           call REDUP1(1,VTp,RD2)
           RD=dmin1(RD1,RD2) ! Eq.(204)
* Modification of AYP 11.06.15
           RDmax=dmax1(RD1,RD2)
           RD=RD*RDmax/(dsqrt(RD*dim(1.d0,RDmax))+dsqrt(RDmax))**2
        else ! moderate SF: Eq.(206)
           RD=1.d4/(1.d4+793.9d0*VTp2+457.3d0*VTn2+66.07d0*VTp2*VTn2+
     +       2.093d0*VTn2**2+.3112d0*VTp2**3+1.068*VTp2**2*VTn2+
     +       .0153d0*(VTp2*VTn2)**2+.006312*VTp2**3*VTn2)
        endif
      else
          print*,'JSF =',JSF
         stop'REDUP2: JSF out of range'
      endif
      return
      end

      subroutine REDCV(TAU,JSF,RCV)
*                                                       Version 09.06.15
* Superfluid reduction factors for specific heat
* [Yakovlev et al. 1999 Phys.-Usp. 42, 737]
* Input: TAU=T/Tc,
*        JSF=1, 2 or 3 for states 1S0, 3P2(m_J=0), 3P2(m_J=2), resp.
* Output: RCV - reduction factor for heat capacity
*         (Yakovlev et al. 1999 Phys.-Usp., after correction of a typo.)
      implicit double precision (A-H), double precision (O-Z)
      parameter(TINY=1.d-19)
      if (TAU.lt.TINY) stop'REDUCV: non-positive TAU'
      if (TAU.gt.1.d0) then ! non-superfluid
         RCV=1.d0
        return
      endif
      call DLGAPS(JSF,TAU,VT)
      if (JSF.eq.1) then ! type A superfluidity
         RCV=(.4186+dsqrt(1.01405+(.5010*VT)**2))**2.5*
     *     dexp(1.456-dsqrt(2.12+VT**2))
      elseif (JSF.eq.2) then ! type B superfluidity
         RCV=(.6893+dsqrt(.6241+(.2824*VT)**2))**2*
     *     dexp(1.934-dsqrt(3.74+VT**2))
      elseif (JSF.eq.3) then ! type C superfluidity
         RCV=(2.188-(9.537d-5*VT)**2+(.1491*VT)**4)/
     /     (1.d0+(.2846*VT)**2+(.01335*VT)**4+(.1815*VT)**6)
      else
          print*,'JSF =',JSF
         stop'REDUCV: JSF out of range'
      endif
      if (RCV.lt.TINY) RCV=0.
      return
      end

      subroutine DLGAPS(JSF,TAU,VT)
*                                                       Version 09.06.15
* Dimensionless (normalized by kT) superfluid gap factors: YKGH (188).
* Input: JSF=1, 2 or 3 for states 1S0, 3P2(m_J=0), 3P2(m_J=2), resp.,
*        TAU=T/Tc - ratio of temperature to critical temperature
* Output: VT=\Delta/T - dimensionless gap factor
      implicit double precision (A-H), double precision (O-Z)
      if (TAU.le.0.d0) stop'DLGAPS: non-positive T/T_c'
      if (TAU.ge.1.d0) then ! no SF
         VT=0.
        return
      endif
      if (JSF.eq.1) then ! type A superfluidity
         VT=dsqrt(1.d0-TAU)*(1.456d0-.157d0/dsqrt(TAU)+1.764d0/TAU)
      elseif (JSF.eq.2) then ! type B superfluidity
         VT=dsqrt(1.d0-TAU)*(.7893d0+1.188d0/TAU)
      elseif (JSF.eq.3) then ! type C superfluidity
         T4=TAU**4
         VT=dsqrt(1.d0-T4)/TAU*(2.03d0-.4903d0*T4+.1727d0*T4**2)
      else
         print*,'JSF =',JSF
         stop'DLGAPS: unknown type of superfluidity.'
      endif
      return
      end

      subroutine CRITEM(XN,MODSF,Yp,Yn,KSF,DLT,Tc9,XSRx)
*                                                       Version 02.09.15
* Input: XN - number density of baryons [fm^{-3}]
*        MODSF=1,...,26 - gap type from gaps.d; -1 switches SF off.
*        Yp, Yn - proton and neutron number fractions of baryons
* Output: KSF - type of superfluidity: 1 - ns, 2 - ps, 3 - nt, -1 - none
*               (n/p stands for neutron/proton, s/t for singlet/triplet)
*         DLT - "gap function" \Delta(0) [MeV]
*         Tc9 - critical temperature [GK]
*         XSRx - SF baryon Fermi momentum divided by mc (auxiliary output)
      implicit double precision (A-H), double precision (O-Z)
      save
      include 'const.inc'
      parameter(COMPRO=COMPTON/DMPE) ! \hbar/(m_p c) [fm]
      parameter(MAXMOD=26)
      dimension ADLT0(MAXMOD),AK(0:3,MAXMOD),LSF(MAXMOD)
      data KREAD/0/
      if (MODSF.eq.-1) then
         KSF=-1
         DLT=0.
         Tc9=0.
         XSRx=0.
        return
      endif
      if(MODSF.lt.1..or.MODSF.gt.MAXMOD) then
         print*,'MODSF =',MODSF,', MAXMOD =',MAXMOD
        stop'CRITEM: MODSF out of range'
      endif
      if (KREAD.eq.0) then ! read parameters from file
         open(1,file='gaps.d',status='OLD')
        do I=1,7
           read(1,*)
        enddo
        do I=1,MAXMOD
           read(1,*) ADLT0(I),(AK(J,I),J=0,3),LSF(I) ! gap parameters
        enddo
        close(1)
         KREAD=1
      endif
      XNP=XN*Yp
      XNN=XN*Yn
      KSF=LSF(MODSF)
      if (KSF.eq.1.or.KSF.eq.3) then ! neutron SF (type A or B)
         XNX=XNN
      elseif (KSF.eq.2) then ! proton SF (type A)
         XNX=XNP
      else
         stop'unknown superfluidity type KSF'
      endif
      CKF=(3.d0*PI**2*XNX)**.3333333d0 ! Fermi wavenumber [fm^{-1}]
      if (CKF.gt.AK(0,MODSF).and.CKF.lt.AK(2,MODSF)) then
         DK0=(CKF-AK(0,MODSF))**2
         DK2=(CKF-AK(2,MODSF))**2
         DLT=ADLT0(MODSF)*DK0*DK2/(DK0+AK(1,MODSF))/(DK2+AK(3,MODSF))
      else
         DLT=0.
      endif
      if (KSF.lt.3) then ! singlet
         Tc=.5669*DLT
      else ! triplet type B; .1187 = .8416/4\sqrt{\pi};
         Tc=.1187*DLT ! where .8416 comes from Yakovlev et al. (1999)
      endif
      Tc9=11.60452d0*Tc ! conversion from MeV to GK
      XSRx=COMPRO*CKF ! Fermi mom. p_{Fp} / m_p c
      if (KSF.ne.2) XSRx=XSRx*PNMRAT ! neutron i.o. proton
      return
      end

      subroutine CRITNP(XN,MODSF1,MODSF2,MODSF3,Yp,Yn,
     *  JSF,Tcn9,Tcp9,XSRn,XSRp)
*                                                       Version 12.09.15
* Input: XN - number density of baryons [fm^{-3}]
*        MODSF1=1-9, MODSF2=10-18, MODSF3=19-26 - gap types (ns,ps,nt,
*          where n/p stands for neutron/proton, s/t for singlet/triplet)
*          from gaps.d; MODSF=-1 switches the given SF off.
*        Yp, Yn - proton and neutron number fractions of baryons
* Output: JSF=1,2,3 for neutron pairing type:
*               A (1S0), B (3P2, m_J=0), C (3P2 ,m_J=2).
*               (n/p stands for neutron/proton, s/t for singlet/triplet)
*         Tcn9,Tcp9 - neutron and proton critical temperature [GK]
*         XSRn,XSRp - neutron and proton Fermi wavenumbers [rel.un.]
* If Tcn9=0, then set JSF:=-1 and XSRn:=0; if Tcp9=0, then XSRp:=0.
      implicit double precision (A-H), double precision (O-Z)
      if (Yn.gt.0.) then
         call CRITEM(XN,MODSF3,Yp,Yn,KSFn,DLTn,Tcn9,XSRn)
         JSF=2
         call CRITEM(XN,MODSF1,Yp,Yn,KSFn1,DLTn1,Tcn91,XSRn1)
        if (Tcn91.gt.Tcn9) then ! type-A neutron SF
           Tcn9=Tcn91
           JSF=1
        endif
      else
         Tcn9=0.
         XSRn=0.
         JSF=-1
      endif
      if (Yp.gt.0.) then
         call CRITEM(XN,MODSF2,Yp,Yn,KSFp,DLTp,Tcp9,XSRp)
      else
         Tcp9=0.
         XSRp=0.
      endif
      return
      end
