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
*         (Yakovlev et al. 1999 Phys.-Usp., after correction of a typo.)
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
