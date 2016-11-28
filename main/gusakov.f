* Last change: 16.06.15
      subroutine REMUP2(JSF,VTn,VTp,RMn,RMp)
*                                                       Version 14.06.15
* Murca SF reduction factors acc.to Gusakov M.E., A&A 389, 702 (2002)
      implicit double precision (A-H), double precision (O-Z)
      save
      parameter (HALFPI=3.14159265359d0/2.d0)
      VTn2=VTn**2
      VTp2=VTp**2
      PHI=datan2(VTn,VTp)
      VT2=VTn2+VTp2
      if (JSF.eq.1) then ! neutron SF type A
         call GUSNA(VTn2,VTp2,PHI,An,Bn,Cn,INDRn)
         call GUSNA(VTp2,VTn2,HALFPI-PHI,Ap,Bp,Cp,INDRp)
      elseif (JSF.eq.2) then ! neutron SF type B
         call GUSNB(VTn2,VTp2,PHI,An,Bn,Cn,INDRn)
         call GUSPB(VTn2,VTp2,PHI,Ap,Bp,Cp,INDRp)
      elseif (JSF.eq.3) then ! neutron SF type C
         call GUSNC(VTn2,VTp2,PHI,An,Bn,Cn,INDRn)
         call GUSPC(VTn2,VTp2,PHI,Ap,Bp,Cp,INDRp)
      else
         print*,'JSF =',JSF
         stop'REMUP2: unknown type of superfluidity.'
      endif
      if (INDRn.eq.4) then
         RMn=Cn*dexp(-An/Bn)
      else
         RMn=dexp(-An*VT2/(1.d0+Bn*VT2)**Cn)
      endif
      if (INDRp.eq.4) then
         RMp=Cp*dexp(-Ap/Bp)
      else
         RMp=dexp(-Ap*VT2/(1.d0+Bp*VT2)**Cp)
      endif
      return
      end

      subroutine GUSNA(VTn2,VTp2,PHI,A,B,C,INDREG)
*                                                       Version 15.06.15
* Coefficients ABC of Gusakov (2002) for R^n in the case of AA-type SF
      implicit double precision (A-H), double precision (O-Z)
      save
      parameter (PI=3.14159265359d0)
      parameter (COEPH=2.d0*PI/.321750554d0) ! coef. in Eq.(A5)
      parameter (PI34=.75d0*PI)
      dimension ANA(18,4)
      data ANA/.257798,.003532,19.57034,.03635,.173561,.039996,
     *  .101014,16.61755,.063353,.101188,.343374,-.135307,
     *  2.404372,1.055914,1.08636,3*0.,
     *    -9.495146,-1.909172,.82025,10.17103,5.874262,.023332,
     *  .003191,201.8576,5.520899,1.257021,-2.367854,1.096571,
     *  .481874,487.429,-.452688,-257.9342,17.83708,0.,
     *    -2.678004,64.33063,-2.736549,.093232, .380818,-.015405,
     *  -16.7934,112.4511,517.5343,.134529,-.174503,-.029008,
     *  1.277903,-25.70616,558.1592,.328108,.642631,.260288,
     *    .26873,.089294,.002913,1.752838d-5,3.047384d-7,.022415,
     *  .001835,5.84941d-7,.00161,9*0./
      CP2=dcos(PHI)**2
      if (VTn2+VTp2.lt.25.d0) then ! Region IV in Fig.1: Eq.(A6)
         INDREG=4
         A=ANA(1,4)*VTn2+ANA(2,4)*VTp2+ANA(3,4)*VTn2*VTp2+
     +     ANA(4,4)*VTn2**3+ANA(5,4)*VTp2**3
         B=1.d0+ANA(6,4)*VTn2+ANA(7,4)*VTp2+ANA(8,4)*VTn2**2
         C=1.d0+ANA(9,4)*VTp2**2
      elseif (VTn2.gt.VTp2) then ! Region I: Eq.(A3)
         INDREG=1
         Y=dsin(PHI+ANA(15,1))**2
         Z=dcos(PHI+ANA(14,1))**2
         A=ANA(1,1)+ANA(2,1)*CP2**2+
     +     ANA(4,1)/(1.d0+ANA(3,1)*CP2)-ANA(5,1)*CP2
         B=ANA(6,1)+ANA(7,1)*CP2**2+
     +     ANA(9,1)/(1.d0+ANA(8,1)*CP2)-ANA(10,1)*CP2
         C=ANA(11,1)-ANA(12,1)/(Y*(1.d0+ANA(13,1)*Z**2)**3)
      elseif (VTn2.gt.VTp2/9.d0) then ! Region II: Eq.(A4)
         INDREG=2
         T=dcos(PHI+ANA(8,2))**2
         Q=dsin(ANA(9,2)*PHI+ANA(10,2))**2
         Y=dsin(PHI+PI34)**2
         A=ANA(1,2)+ANA(2,2)*CP2**2+
     +     ANA(4,2)/(1.d0+ANA(3,2)*CP2)+ANA(5,2)*CP2
         B=ANA(6,2)-ANA(7,2)*Q/(1.d0+ANA(11,2)*T**2)+
     +     ANA(12,2)*Q*T**2
         C=ANA(13,2)+
     +     ((ANA(14,2)+ANA(16,2)/(1.d0+ANA(15,2)*CP2**2))*Y-
     -     ANA(17,2))*Y**2
      else ! Region III: Eq.(A5)
         INDREG=3
         T=dsin(COEPH*PHI)**2
         Y=dsin(PHI)**2
         A=ANA(4,3)+(ANA(1,3)/(1.d0+ANA(2,3)*Y+ANA(3,3)*Y**2)+
     +     ANA(5,3)+ANA(6,3)*T)*Y
         B=ANA(10,3)+(ANA(7,3)/(1.d0+ANA(8,3)*Y+ANA(9,3)*Y**2)+
     +     ANA(11,3)+ANA(12,3)*T)*Y
         C=ANA(16,3)+(ANA(13,3)/(1.d0+ANA(14,3)*Y+ANA(15,3)*Y**2)+
     +      ANA(17,3)+ANA(18,3)*T)*Y
      endif
      return
      end

      subroutine GUSNB(VTn2,VTp2,PHI,A,B,C,INDREG)
*                                                       Version 14.06.15
* Coefficients ABC of Gusakov (2002) for R^n in the case of BA-type SF
      implicit double precision (A-H), double precision (O-Z)
      save
      parameter (PI=3.14159265359d0)
      parameter (PI34=.75d0*PI)
      dimension ANB(17,4)
      data ANB/-.719681,-.024591,.297357,1.260056,.100466,.148464,
     *  .253881,140.3699,.132615,.280765,.375796,-.096843,
     *  3.100942,.275434,.330574,0.,0.,
     *    -6.475443,-1.186294,.591347,6.953996,3.366945,-9.172994,
     *  -2.675793,1.053679,10.38526,7.138369,7*0.,
     *    .316041,-289.2964,2480.961,-268.8219,1984.115,3503.094,
     *  .331551,-.265977,1098.324,65528.01,.0245,.120536,
     *  89.79866,5719.134,285.8473,.402111,16657.19,
     *    .565001,.087929,.006756,1.667194d-4,3.782805d-6,.173165,
     *  1.769413d-5,7.710124d-8,.001695,8*0./
      if (VTn2+VTp2.lt.25.d0) then ! Region IV in Fig.1: Eq.(A10)
         INDREG=4
         A=ANB(1,4)*VTn2+ANB(2,4)*VTp2+ANB(3,4)*VTn2*VTp2+
     +     ANB(4,4)*VTn2**3
         B=dsqrt(1.d0+ANB(6,4)*VTn2+ANB(7,4)*VTp2+
     +     (ANB(5,4)*VTn2+ANB(8,4))*VTn2**3)
         C=1.d0+ANB(9,4)*VTp2**2
      elseif (VTn2.gt.VTp2) then ! Region I: Eq.(A7)
         INDREG=1
         T=dcos(PHI)**2
         Y=dsin(PHI+ANB(15,1))**2
         Z=dcos(PHI+ANB(14,1))**2
         A=ANB(1,1)+ANB(2,1)*T**2+ANB(4,1)/(1.d0+ANB(3,1)*T)-
     -     ANB(5,1)*T
         B=ANB(6,1)+ANB(7,1)*T**2+ANB(9,1)*T/(1.d0+ANB(8,1)*T**2)-
     -     ANB(10,1)*T
         C=ANB(11,1)-ANB(12,1)/(Y*(1.d0+ANB(13,1)*Z**2)**3)
      elseif (VTn2.gt.VTp2/9.d0) then ! Region II: Eq.(A8)
         INDREG=2
         T=dcos(PHI+ANB(8,2))**2
         Z=dcos(PHI)**2
         Q=dsin(ANB(9,2)*PHI+ANB(10,2))**2
         Y=dsin(PHI+PI34)
         A=ANB(1,2)+ANB(2,2)*Z**2+ANB(4,2)/(1.d0+ANB(3,2)*Z)+
     +     ANB(5,2)*Z
         B=.035
         C=ANB(6,2)+ANB(7,2)*Z**2+ANB(9,2)/(1.d0+ANB(8,2)*Z)+
     +     ANB(10,2)*Z
      else ! Region III: Eq.(A9)
         INDREG=3
         PHI2=PHI**2
         PHI3=PHI2*PHI
         PHI4=PHI3*PHI
         A=ANB(12,3)*
     *     (1.d0+ANB(15,3)*PHI2+ANB(16,3)*PHI3+ANB(17,3)*PHI4)/
     /     (1.d0+ANB(13,3)*PHI2+ANB(14,3)*PHI3)
         B=ANB(11,3)+ANB(7,3)/
     /     (1.d0+ANB(8,3)*PHI2+ANB(9,3)*PHI3+ANB(10,3)*PHI4)
         C=ANB(1,3)*
     *     (1.d0+ANB(4,3)*PHI2+ANB(5,3)*PHI3+ANB(6,3)*PHI4)/
     /     (1.d0+ANB(2,3)*PHI2+ANB(3,3)*PHI3)
      endif
      return
      end

      subroutine GUSPB(VTn2,VTp2,PHI,A,B,C,INDREG)
*                                                       Version 14.06.15
* Coefficients ABC of Gusakov (2002) for R^p in the case of BA-type SF
      implicit double precision (A-H), double precision (O-Z)
      save
      dimension APB(17,4)
      data APB/.288203,-.124974,17.39273,.083392,.059046,.028084,
     *  -.01999,28.3721,.244471,-.61047,.023288,.475196,
     *  -.18042,25.51325,.281721,-.080480,-.191637,
     *     0.398261,-0.054952,-0.084964,-0.036240,-0.168712,-0.704750,
     *  -0.066981,1.223731,0.363094,-0.357641,0.869196,-0.364248,
     *  2.668230,-0.765093,-4.198753,2*0.,
     *    0.387542,-195.5462,3032.985,-189.0452,3052.617,442.6031,
     *  0.041901,-0.022201,5608.168,-10761.76,0.064643,0.296253,
     *  106.3387,-75.36126,84.65801,0.530223,-86.76801,
     *    0.272730,0.165858,0.005903,2.555386d-5,2.593057d-7,0.023930,
     *  0.006180,1.289532d-5,0.005368,8*0./
      if (VTn2+VTp2.lt.25.d0) then ! Region IV in Fig.1: Eq.(A13)
         INDREG=4
         A=APB(1,4)*VTp2+APB(2,4)*VTn2+APB(3,4)*VTn2*VTp2+
     +     APB(4,4)*VTp2**3+APB(5,4)*VTn2**3
         B=1.d0+APB(6,4)*VTp2+APB(7,4)*VTn2+APB(8,4)*VTp2**2
         C=1.d0+APB(9,4)*VTn2**2
      elseif (VTn2.gt.VTp2) then ! Region I: Eq.(A11)
         INDREG=1
         T=dcos(PHI)**2
         Y=dsin(PHI+APB(15,1))**2
         A=APB(1,1)+
     *     (APB(2,1)+APB(4,1)/(1.d0+APB(3,1)*T*PHI)**2+
     +       APB(5,1)*T*PHI)*PHI
         B=APB(6,1)+(APB(7,1)+
     +     APB(9,1)/(1.d0+APB(8,1)*T*PHI+APB(10,1)*Y*T)**2+
     +     APB(11,1)*PHI)*PHI
         C=APB(12,1)-APB(13,1)*T+
     +     APB(16,1)/(1.d0+APB(14,1)*T**2)**2+APB(17,1)*T*PHI
      elseif (VTn2.gt.VTp2/9.d0) then ! Region II: Eq.(A12)
         INDREG=2
         Z=dcos(PHI)**2
         A=APB(1,2)+APB(2,2)*Z**2+APB(4,2)/(1.d0+APB(3,2)*Z)+
     +     APB(5,2)*PHI
         B=APB(6,2)+APB(7,2)*Z**2+APB(8,2)/(1.d0+APB(9,2)*Z)+
     +     APB(10,2)*PHI
         C=APB(11,2)+
     +     (APB(13,2)/(1.d0+APB(14,2)*Z)-APB(12,2))*PHI**2+
     +     APB(15,2)*PHI
      else ! Region III: Eq.(A9')
         INDREG=3
         PHI2=PHI**2
         PHI3=PHI2*PHI
         PHI4=PHI3*PHI
         A=APB(12,3)*
     *     (1.d0+APB(15,3)*PHI2+APB(16,3)*PHI3+APB(17,3)*PHI4)/
     /     (1.d0+APB(13,3)*PHI2+APB(14,3)*PHI3)
         B=APB(11,3)+APB(7,3)/
     /     (1.d0+APB(8,3)*PHI2+APB(9,3)*PHI3+APB(10,3)*PHI4)
         C=APB(1,3)*
     *     (1.d0+APB(4,3)*PHI2+APB(5,3)*PHI3+APB(6,3)*PHI4)/
     /     (1.d0+APB(2,3)*PHI2+APB(3,3)*PHI3)
      endif
      return
      end

      subroutine GUSNC(VTn2,VTp2,PHI,A,B,C,INDREG)
*                                                       Version 14.06.15
* Coefficients ABC of Gusakov (2002) for R^n in the case of CA-type SF
      implicit double precision (A-H), double precision (O-Z)
      save
      dimension ANC(19,4)
      data ANC/.897393,-0.045357,.309724,-0.739962,.222597,.032104,
     *  -0.054011,61.73448,.195679,-0.001851,.482581,-0.001637,
     *  -0.685659,1.528415,-0.053834,-0.452426,-0.053502,2*0.,
     *    -3.471368,-0.133540,.143230,3.634659,.496579,.030609,
     *  .005056,.438608,-2.970431,.284703,.898355,-0.036420,
     *  0.407393,-0.058942,.605413,2.851209,-0.800218,1.497718,
     *  1.476375,
     *    .322115,-15.05047,112.9733,-13.79012,128.3156,39.82789,
     *  .164614,49.07699,-3.145006,5132.076,.018737,.100223,
     *  4.055407,390.6242,6.594365,175.7396,441.3965,2*0.,
     *    0.175090,0.088159,3.055763d-3,3.984607d-7,5.591497d-8,.046496,
     *  1.452790d-5,4.505614d-8,1.779724d-3,2.136809d-4,5.365717d-4,
     *  8*0./
      if (VTn2+VTp2.lt.25.d0) then ! Region IV in Fig.1: Eq.(A16)
         INDREG=4
         A=ANC(1,4)*VTn2+ANC(2,4)*VTp2+ANC(3,4)*VTn2*VTp2+
     +     ANC(4,4)*VTn2**3
         B=dsqrt(1.d0+ANC(6,4)*VTn2+ANC(7,4)*VTp2+
     +     (ANC(5,4)*VTn2+ANC(8,4))*VTn2**3)
         C=1.d0+ANC(9,4)*VTp2**2+ANC(10,4)*VTp2+ANC(11,4)*VTn2*VTp2
      elseif (VTn2.gt.VTp2) then ! Region I: Eq.(A14)
         INDREG=1
         T=dcos(PHI)**2
         Z=dcos(PHI+ANC(14,1))**2
         Q=dsin(PHI+ANC(16,1))**2
         S=dsin(PHI+ANC(15,1))**2
         A=ANC(1,1)+ANC(2,1)*T**2+ANC(4,1)/(1.d0+ANC(3,1)*T)-ANC(5,1)*T
         B=ANC(6,1)+
     *     (ANC(7,1)*T+ANC(9,1)/(1.d0+ANC(8,1)*T*Q)**2-ANC(10,1))*T
         C=ANC(11,1)-ANC(12,1)/(S*(1.d0+ANC(13,1)*Z**2)**3)+
     +     ANC(17,1)*S**2
      elseif (VTn2.gt.VTp2/9.d0) then ! Region II: Eq.(A15)
         INDREG=2
         Z=dcos(PHI)**2
         TB=dcos(PHI+ANC(8,2))**2
         TC=dcos(PHI+ANC(15,2))**2
         QB=dsin(ANC(9,2)*PHI+ANC(10,2))**2
         QC=dsin(ANC(16,2)*PHI+ANC(17,2))**2
         A=ANC(1,2)+ANC(2,2)*Z**2+ANC(4,2)/(1.d0+ANC(3,2)*Z)+ANC(5,2)*Z
         B=ANC(6,2)+
     +     (ANC(12,2)*TB**2-ANC(7,2)*QB/(1.d0+ANC(11,2)*TB**2))*QB
         C=ANC(13,2)+
     +     (ANC(19,2)*TC**2-ANC(14,2)/(1.d0+ANC(18,2)*TC**2))*QC
      else ! Region III: Eq.(A9'')
         INDREG=3
         PHI2=PHI**2
         PHI3=PHI2*PHI
         PHI4=PHI3*PHI
         A=ANC(12,3)*
     *     (1.d0+ANC(15,3)*PHI2+ANC(16,3)*PHI3+ANC(17,3)*PHI4)/
     /     (1.d0+ANC(13,3)*PHI2+ANC(14,3)*PHI3)
         B=ANC(11,3)+ANC(7,3)/
     /     (1.d0+ANC(8,3)*PHI2+ANC(9,3)*PHI3+ANC(10,3)*PHI4)
         C=ANC(1,3)*
     *     (1.d0+ANC(4,3)*PHI2+ANC(5,3)*PHI3+ANC(6,3)*PHI4)/
     /     (1.d0+ANC(2,3)*PHI2+ANC(3,3)*PHI3)
      endif
      return
      end

      subroutine GUSPC(VTn2,VTp2,PHI,A,B,C,INDREG)
*                                                       Version 15.06.15
* Coefficients ABC of Gusakov (2002) for R^n in the case of CA-type SF
      implicit double precision (A-H), double precision (O-Z)
      save
      parameter (PI=3.14159265359d0)
      parameter (COEPH=2.d0*PI/.321750554d0) ! coef. in Eq.(A19)
      dimension APC(17,4)
      data APC/.049947,-0.029006,3872.363,.250385,-0.245758,.018241,
     *  .090256,108.8302,1.007326,.061586,.797695,175.5965,
     *  9.306619,-0.551550,1.203014,.096598,-0.441039,
     *    -4.985248,-0.025984,-0.007404,5.294455,-0.201654,.184431,
     *  -0.139729,.415562,2.692073,-0.385832,1.055347,.013667,
     *  -0.509106,-0.267675,.034585,2*0.,
     *    .100241,.005432,-0.748377,.050631,.007900,-0.032915,
     *  -0.000768,.044312,-0.697892,.032534,.080109,.031994,
     *  8.724039,2.982355,-0.062076,2*0.,
     *    0.272905,.058684,2.053694d-3,1.800867d-7,1.911708d-8,.052786,
     *  2.043824d-5,4.458912d-8,1.101541d-3,3.312811d-4,2.682799d-4,
     *  6*0./
      if (VTn2+VTp2.lt.25.d0) then ! Region IV in Fig.1: Eq.(A20)
         INDREG=4
         A=APC(1,4)*VTp2+APC(2,4)*VTn2+APC(3,4)*VTp2*VTn2+
     +     APC(4,4)*VTp2**3
         B=dsqrt(1.d0+APC(6,4)*VTp2+APC(7,4)*VTn2+
     +     (APC(5,4)*VTp2+APC(8,4))*VTp2**3)
         C=1.d0+APC(9,4)*VTn2**2+APC(10,4)*VTn2+APC(11,4)*VTn2*VTp2
      elseif (VTn2.gt.VTp2) then ! Region I: Eq.(A17)
         INDREG=1
         T=dcos(PHI)**2
         Z=dcos(PHI+APC(14,1))**2
         S=dsin(PHI+APC(16,1))**2
         Q=dsin(PHI+APC(15,1))**2
         A=APC(1,1)+APC(2,1)*T**2+APC(4,1)*T/(1.d0+APC(3,1)*T)-
     -     APC(5,1)*T
         B=APC(6,1)+
     *     APC(7,1)*T*Q+APC(9,1)*T/(1.d0+APC(8,1)*T*Q)**2-APC(10,1)*T
         C=APC(11,1)+
     +     (APC(17,1)*S-APC(12,1)/(S*(1.d0+APC(13,1)*Z**2)**3))*T
      elseif (VTn2.gt.VTp2/9.d0) then ! Region II: Eq.(A18)
         INDREG=2
         Z=dcos(PHI)**2
         A=APC(1,2)+APC(2,2)*Z**2+APC(4,2)/(1.d0+APC(3,2)*Z)+
     +     APC(5,2)*PHI
         B=APC(6,2)+APC(7,2)*Z**2+APC(8,2)/(1.d0+APC(9,2)*Z)+
     +     APC(10,2)*PHI
         C=APC(11,2)+APC(12,2)*Z**2+APC(13,2)/(1.d0+APC(14,2)*Z)+
     +     APC(15,2)*PHI
      else ! Region III: Eq.(A19)
         INDREG=3
         T=dsin(COEPH*PHI)**2
         Z=dcos(PHI)**2
         A=APC(1,3)+APC(2,3)*PHI*Z+APC(4,3)/(1.d0+APC(3,3)*Z)+
     +     APC(5,3)*PHI
         B=APC(6,3)+APC(7,3)*PHI*T+APC(8,3)*Z/(1.d0+APC(9,3)*Z**2)+
     +     APC(10,3)*PHI
         C=APC(11,3)+APC(12,3)*PHI*Z+APC(14,3)/(1.d0+APC(13,3)*Z)+
     +     APC(15,3)*PHI
      endif
      return
      end
