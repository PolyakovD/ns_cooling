** Superfluid gap functions, critical temperatures, reduction factors
* Needs to be linked with: bskfit.f
* Last update: 28.02.15
      subroutine SFCRIT(XN,KEOS,MODSF,KSF,DLT,Tc)
*                                                       Version 28.02.15
* Input: XN - number density of baryons [fm^{-3}]
*        KEOS=19,20,21 - BSk EOS number
*        MODSF=1,...,26 - gap type from gaps.d
* Output: KSF - type of superfluidity: 1 - ns, 2 - ps, 3 - nt
*               (n/p stands for neutron/proton, s/t for singlet/triplet)
*         DLT - gap function [MeV]
*         Tc - critical temperature [K]
      implicit double precision (A-H), double precision (O-Z)
      save
      parameter(MAXMOD=26)
      parameter (PI=3.14159265d0)
      dimension ADLT0(MAXMOD),AK(0:3,MAXMOD),LSF(MAXMOD)
      data KREAD/0/
      if (KEOS.eq.19) then
         XND=2.63464d-4
         XNC=.0885
      elseif (KEOS.eq.20) then
         XND=2.62873d-4
         XNC=.0854
      elseif (KEOS.eq.21) then
         XND=2.57541d-4
         XNC=.0809
      else
         stop'SFCRIT: unknown EOS'
      endif
      if(MODSF.lt.1..or.MODSF.gt.MAXMOD)stop'SFCRIT: MODSF out of range'
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
      if (XN.lt.XNC) call INCRUST(KEOS,XN,XND,XNC,CMI,CMI1,Zcell,Zclust,
     *  ap,an,Cp,Cn,Dn_p,Dn_n,xnuc,xnucn)
      if (XN.lt.XND) then ! outer crust
         XNP=0.
         XNN=0.
      elseif (XN.lt.XNC) then ! inner crust
         XNP=0.
         XNN=(CMI1-CMI)/CMI1*XN
      else ! core
         call FRACORE(KEOS,XN,Ye,Ymu)
         XNP=XN*(Ye+Ymu) ! proton number density [fm^{-3}]
         XNN=XN-XNP ! neutron number density
      endif
      KSF=LSF(MODSF)
      if (KSF.eq.1.or.KSF.eq.3) then
         XNX=XNN
      elseif (KSF.eq.2) then
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
      else ! triplet; .1183 --> .1187 = .8416/4\sqrt{\pi};
         Tc=.1187*DLT ! where .8416 comes from Yakovlev et al. (1999)
      endif
      Tc=1.16045d10*Tc ! conversion from MeV to K
      return
      end

      subroutine REDUCV(TAU,KSF,RCV)
*                                                       Version 28.02.15
* Superfluid reduction factors for specific heat
* [Yakovlev et al. 1999 Phys.-Usp.]
* Input: TAU=T/Tc,
*        KSF - type of superfluidity: 1 - ns, 2 - ps, 3 - nt
*               (n/p stands for neutron/proton, s/t for singlet/triplet)
* NEGATIVE values of KSF may be used to directly input JSF=-KSF:
*        JSF=1, 2 or 3 for states 1S0, 3P2(m_J=0), 3P2(m_J=2), resp.
* Output: RCV - reduction factor for heat capacity
*         (Yakovkev et al. 1999 Phys.-Usp., after correction of a typo.)
      implicit double precision (A-H), double precision (O-Z)
      dimension LJSF(3)
      save
      parameter(TINY=1.d-19)
      data LJSF/1,1,2/ ! 1S0 (ns), 1S0 (ps), 3P2(m_J=0) (nt)
      if (TAU.lt.TINY) stop'REDUCV: non-positive TAU'
      if (TAU.gt.1.d0) then ! non-superfluid
         RCV=1.d0
        return
      endif
      if (KSF.gt.0) then
         JSF=LJSF(KSF)
      else
         JSF=-KSF
      endif
      if (JSF.eq.1) then ! type A superfluidity
         VT=dsqrt(1.d0-TAU)*(1.456-.157/dsqrt(TAU)+1.764/TAU)
         RCV=(.4186+dsqrt(1.01405+(.5010*VT)**2))**2.5*
     *     dexp(1.456-dsqrt(2.12+VT**2))
      elseif (JSF.eq.2) then ! type B superfluidity
         VT=dsqrt(1.d0-TAU)*(.7893+1.188/TAU)
         RCV=(.6893+dsqrt(.6241+(.2824*VT)**2))**2*
     *     dexp(1.934-dsqrt(3.74+VT**2))
      elseif (JSF.eq.3) then ! type C superfluidity
         T4=TAU**4
         VT=dsqrt(1.d0-T4)/TAU*(2.03-.4903*T4+.1727*T4**2)
         RCV=(2.188-(9.537d-5*VT)**2+(.1491*VT)**4)/
     /     (1.d0+(.2846*VT)**2+(.01335*VT)**4+(.1815*VT)**6)
      else
         stop'REDUCV: JSF out of range'
      endif
      if (RCV.lt.TINY) RCV=0.
      return
      end

      subroutine SFCURV(KEOS,MODSF,N,RHOD,TcD)
*                                                       Version 23.04.15
* Input: KEOS=19,20,21 - BSk EOS number
*        MODSF=1,...,26 - gap type from gaps.d
*        N - spatial knots number
*        RHOD - rho distribution
* Output: TcD - critical temperature distribution [K]
      implicit none
	  integer KEOS, MODSF, N
	  real*8 RHOD(*), TcD(*)
	  integer i, KSF
	  real*8 XN, DLT
	  do i = 1, N
		call BSkNofR(KEOS,RHOD(i),XN)
		call SFCRIT(XN,KEOS,MODSF,KSF,DLT,TcD(i))
	  end do
	  return
      end