c==============================================================================
      real function psr_ne(x, y, z, dmod, ip, lip)
c==============================================================================
c
c     General function to compute the local electron
c     density given the Cartesian coords (x,y,z)
c     in the Galaxy wrt the Galactic centre (0,0,0)
c     in pc/cc for various distance models which are
c     specified by the integer mod. Options for mod
c     are as follows...
c
c                        mod = 0  : LMT model
c                              1  : BWHV model
c                              2  : CWFSR model
c                              3  : TC93 model
c                              4  : NE2001 model
c                              5  : NE=0.03 model
c                              6  : GBC01 model
c
      implicit none
      real rsun
      parameter (rsun = 8.5)
c
c     passed down variables...
c
      real x, y, z
      integer dmod
c
c     local variables
c
      real ne1,ne2,nea,negc,nelism,necN,nevN
      real F1, F2, Fa, Fgc, Flism, FcN, FvN
      integer whicharm,wlism,wLDR,wLHB,wLSB,wLOOPI,hitclump,hitvoid
      integer wvoid,lout
      real negum,a
      real lgn, bgn, dgn, rgn
      logical first
      character ip*(*)
      integer lip
      data first /.true./
      data lgn/260.0/,bgn/0.0/,dgn/0.5/,rgn/0.115/
      real xgn, ygn, zgn, alpha, psr_dist, r, sech2
      common /outputtxt/ lout
      if (first) then
c
c       Calc Position of the Gum nebula
c
        call calc_xyz(lgn,bgn,dgn,xgn,ygn,zgn)
c        if (dmod.eq.0) write(lout,*) 'LMT distance model'
c        if (dmod.eq.1) write(lout,*) 'BWHV distance model'
c        if (dmod.eq.2) write(lout,*) 'CWFSR distance model'
c        if (dmod.eq.3) write(lout,*) 'TC93 distance model'
c        if (dmod.eq.4) write(lout,*) 'NE2001 distance model'
c        if (dmod.eq.5) write(lout,*) 'NE0.03 distance model'
c        if (dmod.eq.6) write(lout,*) 'GBC01 distance model'
	write(lout,*)
        first = .false.
      end if
c
c     Does the line of sight to the pulsar pass
c     through the Gum nebula if so set alpha to
c     1, to take into account the extra cont
c     ribution to n_{e} from this nearby
c     H II region. LMT parameters Dgn=0.5kpc
c     lgn=260, bgn=0, radius=0.115kpc
c
c      if (sqrt((x-xgn)**2.0+(y-ygn)**2.0).lt.rgn) then
      if (psr_dist(x,y,z,xgn,ygn,zgn).lt.rgn) then
        alpha = 1.0
      else
        alpha = 0.0
      end if
      alpha=0.0
c
c     Galactocentric radius
c
      R = sqrt((x ** 2) + (y ** 2))
      if (dmod.eq.0) then
c
c       Calculate n_{e} according to LMT formula
c
        psr_ne = (0.28 * alpha) + ((((0.025 + (0.015 *
     &  exp(- (1.0 * abs(z / 0.070))))) *
     &  (2.0 / (1.0 + (R / rsun))) )))
      else if (dmod.eq.1) then
c
c       Calculate n_{e} according to Bhattacharya et al formula
c
        psr_ne = (0.28 * alpha) + (((( (0.0294*exp(- (1.0 *
     &  abs(z / 0.68)))) + (0.0176 * exp(- (1.0 *
     &  abs(z / 0.070))))) * (2.0 / (1.0 + (R / rsun))))) )
      else if (dmod.eq.2) then
c
c       Calculate n_{e} according to Cordes et al formula
c       DRL added extra Gum nebula term onto published model
c       25/03/93 @ JB
c
        psr_ne = 0.025*exp(-1.0*abs(z))*exp(-1.0*(R**2/400.0)) +
     &  0.2*exp(-1.0*abs(z/0.15))*exp(-1.0*((R-4.0)**2/4.0)) +
     &  (0.28*alpha)
      else if (dmod.eq.3) then
         call density(x,y,z,ne1,ne2,nea,negum,ip,lip)
         psr_ne=ne1+ne2+nea+negum
      else if (dmod.eq.4) then
          write(*,*) "made it into psr_ne.f"
            call density_2001(x,y,z,
     .                ne1,ne2,nea,negc,nelism,necN,nevN,
     .                F1, F2, Fa, Fgc, Flism, FcN, FvN,
     .                whicharm, wlism, wLDR, wLHB, wLSB, wLOOPI,
     .                hitclump, hitvoid, wvoid)
            psr_ne=ne1+ne2+nea+negc+nelism+necN+nevN
      else if (dmod.eq.5) then
c
c     Constant electron density model
c
         psr_ne = 0.03
      else if (dmod.eq.6) then
c
c       Calculate n_{e} according to Gomez et al
c
         psr_ne = 1.77e-2*sech2(R/15.4)*sech2(z/1.10)/sech2(Rsun/15.4)+
     &           1.07e-2*sech2(R/3.6)*sech2(z/0.04)/sech2(Rsun/3.6)
c         a = abs(z)
c         psr_ne =.0203*exp(-1.0*R/30.4)*exp(Rsun/30.4)*exp(-1.0*a/1.07)+
c     &           .0071*exp(-1.0*R/1.5)*exp(Rsun/1.5)*exp(-1.0*a/0.05)
      else
        stop 'Distance model undefined!!!'
      end if
      end
