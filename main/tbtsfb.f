* Last update: 28.04.15
      subroutine TbTsFb(GS14,Tb9,B12,COSLAT,TS6,FB)
*                                                       Version 28.04.15
* Stems from TsTbBSk v.12.05.13 with a correction of 07.12.14
* Tsurface as function of Tint for a ground-state (BSk21) envelope
* with a relativistic dipole field.
* Fit to 2D simulations with neutrino emission but without convection.
* Input: GS14 (g in 10^{14} cm/s^2)
*        Tb9 (Tb in 10^9 K)
*        B12 (magnetic field in 1e12 G at the magnetic pole)
*        COSLAT (cosine of colatitude to the magnetic pole(polar angle))
* Output: TS6 - surface temperature Ts in 10^6 K
*         FB - thermal flux density Fb in erg/(cm^2 s)
      implicit double precision (A-H), double precision (O-Z)
      save
      B12sq=dsqrt(B12)
      Tb9sq=dsqrt(Tb9)
* polar Ts:
      T0=(15.7*Tb9sq*Tb9+1.36*Tb9)**.3796 ! no neutrinos, + e-e
      P2=.337/(1+.02*B12sq)
      T1=1.13*B12**.119*TB9**P2
      P3=1.d0+.15*B12sq
      Ts=dsqrt(dsqrt(GS14*(T1**4+P3*T0**4))) ! Ts without neutrino limit
      Tsmax=5.2*GS14**.65+.093*dsqrt(GS14)*B12sq ! neutrino limitation
      REDUCT=dsqrt(dsqrt(1.d0+(Ts/Tsmax)**4))
      Tpol6=Ts/REDUCT
* Ratio Tpol/Teq:
      b=B12sq**5/Tb9sq
      rat=1.d0+(1230.*Tb9)**3.35*B12*
     *  dsqrt(1.+2.d0*B12**2)/(B12+450.*Tb9+119.*B12*Tb9)**4+
     +  .0066*b/(1+.00258*b)
      Teq6=Tpol6/rat
* Fraction of the difference Tpol-Teq at given angle:
      b=10.d0*B12/(Tb9sq+.1d0*B12/dsqrt(Tb9sq))
      a=b*Tb9sq/3.d0
      c=a+b+1.d0
      frac=COSLAT**2*c/(1.d0+a*COSLAT+b*COSLAT**2)
* Final result:
      TS6=Teq6+(Tpol6-Teq6)*frac
      TS6=DMIN1(TS6,Tb9*1.d3)
      FB=5.6704d19*(TS6*REDUCT)**4 ! 5.6704e-5 - Stefan-Boltzmann const.
      FB=FB*1.027 ! correction for the envelope thickness at 1.4 M_Sun
      return
      end
