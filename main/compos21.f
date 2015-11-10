* Crust composition and effective nuxcleon masses for BSk21 EoS
      subroutine EFFMASS(KEOS,XN,Yp,EFMp,EFMn)
*                                                       Version 11.04.13
*                                                     Corrected 19.03.14
* Proton and neutron effective masses in the nuclear matter.
* Calculation according to Eq.(A10) of Chamel, Goriely, & Pearson (2009)
*     with parameters from Table IV of Goriely, Chamel, & Pearson (2010)
* The code by A.Y.Potekhin <palex@astro.ioffe.ru>
* Input: KEOS - number of BSK model (19,20,or 21)
*        XN = n_b - baryon number density in fm^{-3}
*        Yp - proton fraction
* Output: EFMp, EFMn - ratios of the proton and neutron effective masses
*        to the masses of an isolated proton and neutron, respectively.
      implicit double precision (A-H), double precision (O-Z)
      save AT1,AT2X2,AT4,AT5,AX1,AX4,AX5
      dimension AT1(3),AT2X2(3),AT4(3),AT5(3),AX1(3),AX4(3),AX5(3)
      parameter (HbarC=197.3269718) ! \hbar*c [Mev*fm]
      parameter (Ep=938.272046,En=939.565379) ! m_p*c^2, m_n*c^2 [MeV]
      data AT1/403.072,438.219,396.131/,
     *  AT2X2/-1055.55,-1147.64,-1390.38/,
     *  AT4/-60.,-100.,-100./
     *  AT5/-90.,-120.,-150./,
     *  AX1/-0.137960,-0.392047,0.0648452/, ! corrected 19.03.14
     *  AX4/-6.,-3.,2./,
     *  AX5/-13.,-11.,-11./
      K=KEOS-18
      Yn=1.d0-Yp
      if (K.eq.1) then
         BETA=1.d0/3.d0
      elseif (K.eq.2) then
         BETA=1.d0/6.d0
      elseif (K.eq.3) then
         BETA=.5d0
      else
        stop'EFFMASS: unknown EOS'
      endif
      GAMMA=1.d0/12.d0
      DFMp=AT1(K)*((1.d0+.5d0*AX1(K))-(.5d0+AX1(K))*Yp)+
     +  AT2X2(K)*(.5d0+Yp)+
     +  AT4(K)*((1.d0+.5d0*AX4(K))-(.5d0+AX4(K))*Yp)*XN**BETA+
     +  AT5(K)*((1.d0+.5d0*AX5(K))+(.5d0+AX5(K))*Yp)*XN**GAMMA
      DFMn=AT1(K)*((1.d0+.5d0*AX1(K))-(.5d0+AX1(K))*Yn)+
     +  AT2X2(K)*(.5d0+Yn)+
     +  AT4(K)*((1.d0+.5d0*AX4(K))-(.5d0+AX4(K))*Yn)*XN**BETA+
     +  AT5(K)*((1.d0+.5d0*AX5(K))+(.5d0+AX5(K))*Yn)*XN**GAMMA
      EFMp=1.d0/(1.d0+.5d0*DFMp*XN*Ep/HbarC**2)
      EFMn=1.d0/(1.d0+.5d0*DFMn*XN*En/HbarC**2)
      return
      end

      subroutine CRUBSk21(XN,KZ,CMI,CMI1,xnuc,XND)
*                                                       Version 24.01.15
* BSk21 composition of the crust
* Stems from COMPOS21 v.15.12.14: XN i.o. RHO on the input.
* Input: XN = n_b - baryon number density in fm^{-3}
* Output: KZ - nuclear charge number,
*         CMI - nuclear mass number,
*         CMI1 - total number of nucleons per nucleus (incl.free n)
*         xnuc=sqrt{5<r^2>/3} - nucl.size parameter (Gnedin et al. 2003)
*         XND - neutron-drip number density of baryons in fm^{-3} 
*              = the inner/outer crust interface (Pearson et al. 2012),
      implicit double precision (A-H), double precision (O-Z)
      parameter (KEOS=21)
* First, calculation for the inner crust:
      KZ=40
      call INCRUST(KEOS,XN,XND,XNC,CMI,CMI1,Zcell,Zclust,
     *  ap,an,Cp,Cn,Dn_p,Dn_n,xnuc,xnucn)
      if (XN.gt.XND) goto 50 ! inner crust indeed
* Outer crust or ocean (Pearson et al. 2011, PRC 83, 065810):
      if (XN.lt.4.90d-9) then
         KZ=26
         CMI=56.
      elseif (XN.lt.1.63d-7) then
         KZ=28
         CMI=62.
      elseif (XN.lt.8.01d-7) then
         KZ=28
         CMI=64.
      elseif (XN.lt.8.68d-7) then
         KZ=28
         CMI=66.
      elseif (XN.lt.1.88d-6) then
         KZ=36
         CMI=86.
      elseif (XN.lt.6.87d-6) then
         KZ=34
         CMI=84
      elseif (XN.lt.1.68d-5) then
         KZ=32
         CMI=82.
      elseif (XN.lt.3.21d-5) then
         KZ=30
         CMI=80.
      elseif (XN.lt.4.39d-5) then
         KZ=29
         CMI=79.
      elseif (XN.lt.5.42d-5) then
         KZ=28
         CMI=78.
      elseif (XN.lt.7.58d-5) then
         KZ=28
         CMI=80.
      elseif (XN.lt.1.25d-4) then
         KZ=42
         CMI=124.
      elseif (XN.lt.1.50d-4) then
         KZ=40
         CMI=122.
      elseif (XN.lt.1.77d-4) then
         KZ=39
         CMI=121.
      elseif (XN.lt.1.96d-4) then
         KZ=38
         CMI=120.
      elseif (XN.lt.2.40d-4) then
         KZ=38
         CMI=122.
      else
         KZ=38
         CMI=124.
      endif
      CMI1=CMI
      XSR=.0100884*(RHO*KZ/CMI1)**.3333333d0
      xnuc=.00155*(CMI/KZ)**.33333*XSR ! i.e. r_nuc=1.15 A^{1/3} fm
   50 continue
      return
      end
