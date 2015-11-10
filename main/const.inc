! Fundamental constants
! Last change: 22.05.11
      parameter (PI=3.14159265d0,DECLN=2.302585093d0, &
       eV=1.602176487d-12, & ! erg/eV
       AUM=1822.9, & ! a.m.u./m_e
       AUD=15819.3832d0, & ! relativistic unit of density in g/cc
       SBCREL=0.1644933, & ! Stefan-Boltzmann constant in rel.un.
       UAMU=1.660538782d-24, & ! unified atomic mass unit in grams
       UNIPR=1.4217752544d13, & ! relativistic unit of pressure in Mbar
       UNIB=44.14004916, & ! rel.un.of magnetic field [TG]
       GAMIMELT=175., & ! Gamma_i value at which the melting occurs
       RSIMELT=140.)   ! ion density parameter of quantum melting
! [approximate value, acc.to Jones & Ceperley (1996) PRL 76, 4572]
      parameter(Ryd=13.6057d0, & ! Rydberg energy / eV
       UN_T6=.3157746d0, & ! = 2*(Rydberg energy)/10^6 K = a.u. of T [MK]
       BOHR=137.03599968d0, & ! inverse fine-structure constant
       AURHO=11.2058713d0, & ! (atomic unit of mass density)*AUM in g/cc
       BOHRAD=.52917720859d-8, & ! Bohr radius in cm
       COMPTON=3.8615926459d-11, & ! Compton length over 2\pi in cm
       AUTIME=2.4188843265d-17) ! atomic unit of time in s
      parameter(C13 = 1.d0 / 3.d0, C53 = 5.d0 / 3.d0, C73 = 7.d0 / 3.d0)
      parameter(TwoPI = 2.d0 * PI, FourPI = 4.d0 * PI)
      parameter(UNITEMP = UN_T6 * BOHR**2, & ! rel.un. of T [MK] = 5930
       UNIOP=1.d0/(AUD*COMPTON), & ! rel.un.of opacity [cm^2/g]=1.637d6
       AUNUMDEN=1.d0/BOHRAD**3, & ! a.u. of number density
       UN_B12=BOHR**2/UNIB) ! inverse a.u. of magnetic field [TG]
