cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 24 June 2002: added calculations of path lengths through LISM components
c 26  May 2002: modified for NE2001 routines which are cleaned-up
c             versions of development routines.
c Nov 1999 - May 2002: development versions
c 1992-1993: TC93 version

      subroutine dmdsm(
     .   l,b,ndir,dmpsr,dist,limit,sm,smtau,smtheta,smiso,ip,lip)

      implicit none
      real l,b,dmpsr,dist,sm,smtau,smtheta,smiso
      integer ndir
      character*1 limit
      character ip*(*)
      integer lip

c  Computes pulsar distance and scattering measure
c  from model of Galactic electron distribution.

c  Input: real l	galactic longitude in radians
c         real b	galactic latitude in radians
c         integer ndir  >= 0 calculates dist from dmpsr
c                       < 0 for dmpsr from dist
c Input or output:
c	  real dmpsr	(dispersion measure in pc/cm^3)
c         real dist	(distance in kpc)

c  Output:
c	  char*1 limit	(set to '>' if only a lower distance limit can be
c			 given; otherwise set to ' ')
c         sm            (scattering measure, uniform weighting) (kpc/m^{20/3})
c         smtau         (scattering measure, weighting for pulse broadening)
c         smtheta       (scattering measure, weighting for angular broadening
c                        of galactic sources)
c	  smiso 	(scattering measure appropriate for calculating the
c			isoplanatic angle at the source's location


c       parameter(alpha = 11./3.)
c       parameter(pi = 3.14159)
c       parameter(c_sm = (alpha - 3.) / 2. * (2.*pi)**(4.-alpha) )


	real c_sm, c_u, sm_factor
        parameter(c_sm = 0.181)         ! constant in sm definition
        parameter(c_u = 10.16)          ! units conversion for sm
        parameter(sm_factor = c_sm * c_u)


	integer wg1, wg2, wga, wggc, wglism, wgcN, wgvN
        common /modelflags/ wg1, wg2, wga, wggc, wglism, wgcN, wgvN

c parameters of large-scale components (inner+outer+arm components):
        real n1h1,h1,A1,F1,n2,h2,A2,F2,na,ha,wa,Aa,Fa
        common/galparams/n1h1,h1,A1,F1,n2,h2,A2,F2,
     .                na,ha,wa,Aa,Fa

c factors for controlling individual spiral arms:
c       narm:   multiplies electron density (in addition to the`fac'
c                     quantities)
c       warm:   arm width factors that multiply nominal arm width
c       harm:   arm scale height factors
c       farm:   factors that multiply n_e^2 when calculating SM

        integer narmsmax, narmsmax1
        parameter (narmsmax=5, narmsmax1=narmsmax+1)
        real narm, warm, harm, farm
        common/armfactors/
     .     harm(narmsmax),narm(narmsmax),warm(narmsmax),farm(narmsmax)

	real armpaths, armdistances
	common/armpathlengths/ armpaths(narmsmax1),
     .     armdistances(narmsmax1)


	real dx0, dy0, dz0
	common/dxyz/dx0,dy0,dz0

	integer whicharm

c Large scale components:

	real ne1, ne2, nea
	real F1val, F2val, Faval

c Galactic center:

	real negc, Fgc

c LISM:
	real nelism, Flism
        integer wlism, wLDR, wLHB, wLSB, wLOOPI
	real ldr_path, lhb_path, lsb_path, loopI_path
	real ldr_dist, lhb_dist, lsb_dist, loopI_dist
	integer wtemp

c clumps:
        real necN, FcN
	integer hitclump

c voids:
	real nevN, FvN
	integer hitvoid, wvoid

c subroutines needed:
c	density_2001 (and those that it calls) in density.NE2001.f
c       scattering routines in scattering98.f

	real R0, rrmax, zmax, dmax
	data R0/8.5/
c	data rrmax/30.0/		! Max radius for reliable ne
	data rrmax/50.0/		! Max radius for reliable ne
c	data zmax/1.76/			! Max |z|
c	data zmax/5.00/			! Max |z|
	data zmax/25.00/		! Max |z|
        data dmax/50.0/                 ! maximum distance calculated

	logical first
	data first/.true./

	save

c other variables
	real x, y, z, r, rr
	real sl, cl, sb, cb

	real d, dstep, dtest, dstep_pc, dd

	real dm, dmstep
	real sm_sum1, sm_sum2, sm_sum3, sm_sum4, sm_term
        real sm_sum1_last, sm_sum2_last, sm_sum3_last, sm_sum4_last
	integer nstep
	integer ncount
	integer i

	real dm1, dm2, dma, dmgc, dmlism, dmcN, dmvN
	real sm1, sm2, sma, smgc, smlism, smcN, smvN
	real dsm1, dsm2, dsma, dsmgc, dsmlism, dsmcN, dsmvN

	integer wtotal
	real ne


c	open(24,file='fort.24', status='unknown')
c	open(25,file='fort.25', status='unknown')
c	write(25,*) l*180./acos(-1.), b*180./acos(-1.), ' = l, b'
c        write(25,1000)
c 1000 format(
c     .   '    d       x       y       z       ne',
c     .   '      dsm    arm  cl void')


	if(first) then
c initial call to density routine to set variable values
c through read-in of parameter file:
        x = 0.0
        y = R0
        z = 0.0
            call density_2001(x,y,z,
     .        ne1,ne2,nea,negc,nelism,necN,nevN,
     .        F1val, F2val, Faval, Fgc, Flism, FcN, FvN,
     .        whicharm, wlism, wldr, wlhb, wlsb, wloopI,
     .        hitclump, hitvoid, wvoid, ip, lip)
c       write(6,*) 'ne1,ne2,negc,nelism,necN,nevN = ',
c    .              ne1,ne2,negc,nelism,necN,nevN
	first=.false.
        endif

	sl=sin(l)
	cl=cos(l)
	sb=sin(b)
	cb=cos(b)
	limit=' '
c	dstep=0.02			! Step size in kpc
c       dstep = min(h1, h2) / 10.       ! step size in terms of scale heights
c	dstep=0.01
        dstep=0.05
        if(ndir.lt.0) dtest=dist
        if(ndir.ge.0) dtest=dmpsr/(n1h1/h1)   ! approximate test distance
        nstep = dtest / dstep	        ! approximate number of steps
        if(nstep.lt.10) dstep=dtest/10  ! make no steps >= 10

c  Sum until dm is reached (ndir >= 0) or dist is reached (ndir < 0).
c  Guard against too few terms by counting number of terms (ncount) so that
c  routine will work for n_e models with large n_e near the Sun.

    5   continue
 	dstep_pc = 1000.*dstep
	dm=0.0
        sm_sum1 = 0.                    ! sum of C_n^2
        sm_sum2 = 0.                    ! sum of C_n^2 * s
        sm_sum3 = 0.                    ! sum of C_n^2 * s^2
	sm_sum4 = 0.			! sum of C_n^2 * s^{5./3.}

	do i=1,narmsmax1
	  armpaths(i) = 0.
	  armdistances(i) = 0.
	enddo

	dm1 = 0.
	dm2 = 0.
	dma = 0.
	dmgc = 0.
	dmlism = 0.
	dmcN = 0.
	dmvN = 0.

	sm1 = 0.
	sm2 = 0.
	sma = 0.
	smgc = 0.
	smlism = 0.
	smcN = 0.
	smvN = 0.

	ldr_path = 0.
	lhb_path = 0.
	lsb_path = 0.
	loopI_path = 0.

	ldr_dist = 0.
	lhb_dist = 0.
	lsb_dist = 0.
	loopI_dist = 0.

        ncount = 0

	d=-0.5*dstep
	do 10 i=1,99999
          ncount = ncount + 1
	  d=d+dstep			! Distance from Sun in kpc
	  r=d*cb
	  x=r*sl
	  y=R0-r*cl
	  z=d*sb
	  rr=sqrt(x**2 + y**2)		! Galactocentric radius
	  if(ndir.ge.0.and.
     +      (d.gt.dmax.or.abs(z).gt.zmax.or.rr.gt.rrmax)) go to 20

          if(ndir.lt.3) then
            call density_2001(x,y,z,
     .        ne1,ne2,nea,negc,nelism,necN,nevN,
     .        F1val, F2val, Faval, Fgc, Flism, FcN, FvN,
     .        whicharm, wlism, wldr, wlhb, wlsb, wloopI,
     .        hitclump, hitvoid, wvoid,ip,lip)
	  endif

	  if(ndir.ge.3) then
            call density_2001(x+dx0,y+dy0,z+dz0,
     .        ne1,ne2,nea,negc,nelism,necN,nevN,
     .        F1val, F2val, Faval, Fgc, Flism, FcN, FvN,
     .        whicharm, wlism, wldr, wlhb, wlsb, wloopI,
     .        hitclump, hitvoid, wvoid,ip,lip)
	  endif


c wlism = 1 causes the lism component to override smooth Galactic components
c wvoid = 1 overrides everything except clumps
	  ne=
     .       (1.-wglism*wlism)*
     .       (wg1*ne1 +
     .        wg2*ne2 +
     .        wga*nea +
     .        wggc*negc) +
     .        wglism*wlism*nelism
          ne = (1-wgvN*wvoid)*ne + wgvN*wvoid*nevN + wgcN*necN
	  dmstep=dstep_pc*ne
	  dm=dm+dmstep			! Add DM for this step
	  wtotal = (1-wgvN*wvoid)*(1-wglism*wlism)
	  dm1 = dm1 + wtotal*wg1*ne1
	  dm2 = dm2 + wtotal*wg2*ne2
	  dma = dma + wtotal*wga*nea
	  dmgc = dmgc + wtotal*wggc*negc
	  dmlism = dmlism + (1.-wgvN*wvoid)*wglism*wlism*nelism
	  dmcN = dmcN + wgcN*necN
	  dmvN = dmvN + wgvN*wvoid*nevN

c         write(24,"('n:',7f10.6,1x))")
c    .        ne1,ne2,nea,negc,nelism,necN,nevN
c        write(24,"(i2,1x,7(f10.5,1x))")
c    .      wtotal,dm1,dm2,dma,dmgc,dmlism,dmcN,dmvN

c         sm_term =
c    .       (1.-wglism*wlism)*
c    .       (wg1   * F1  * ne1**2 +
c    .        wg2   * F2  * ne2**2 +
c    .        wga   * Fa  * nea**2 +
c    .        wggc  * Fgc * negc**2) +
c    .        wglism*wlism * Flism * nelism**2
c	  sm_clumps = FcN * necN**2
c	  sm_voids  = FvN * nevN**2
c         sm_term = (1-wgvN*wvoid) * sm_term
c    .            + wgvN * wvoid * sm_voids
c    .            + wgcN * sm_clumps


	dsm1 = wtotal*wg1*ne1**2*F1
	dsm2 = wtotal*wg2*ne2**2*F2
	dsma =wtotal*wga*nea**2*Fa
	dsmgc = wtotal*wggc*negc**2*Fgc
	dsmlism = (1.-wgvN*wvoid)*wglism*wlism*nelism**2*Flism
	dsmcN = wgcN*necN**2*FcN
	dsmvN = wgvN*wvoid*nevN**2*FvN

	sm_term = dsm1+dsm2+dsma+dsmgc+dsmlism+dsmcN+dsmvN

	sm1 = sm1 + dsm1
	sm2 = sm2 + dsm2
	sma = sma + dsma
	smgc = smgc + dsmgc
	smlism = smlism + dsmlism
	smcN = smcN + dsmcN
	smvN = smvN + dsmvN


        sm_sum1 = sm_sum1 + sm_term
        sm_sum2 = sm_sum2 + sm_term * d
        sm_sum3 = sm_sum3 + sm_term * d**2
        sm_sum4 = sm_sum4 + sm_term * d**1.66667

c pathlengths through LISM components:
c take into account the weighting hierarchy, LHB:LOOPI:LSB:LDR


	if(wlism .eq. 1) then
	  if(wlhb .eq. 1) then
	    lhb_path = lhb_path + dstep
	    lhb_dist = lhb_dist + d
	  endif
	  if(wloopI .eq. 1) then
    	    wtemp = (1-wLHB)
	    loopI_path = loopI_path + wtemp*dstep
	    loopI_dist = loopI_dist + wtemp*d
	  endif
	  if(wlsb .eq. 1) then
    	    wtemp = (1-wLHB)*(1-wloopI)
	    lsb_path = lsb_path + wtemp*dstep
	    lsb_dist = lsb_dist + wtemp*d
	  endif
	  if(wldr .eq. 1) then
    	    wtemp = (1-wLHB)*(1-wloopI)*(1-wlsb)
	    ldr_path = ldr_path + wtemp*dstep
	    ldr_dist = ldr_dist + wtemp*d
	  endif
	endif

c pathlengths: whicharm = 0,5 (currently).
c 	                  1,4 for the equivalent of the TC93 arms
c                         5   for the local arm
c                         0   means interarm paths

	armpaths(whicharm+1) = armpaths(whicharm+1) + dstep
	armdistances(whicharm+1) = armdistances(whicharm+1) + d

c       write(99,"(2(f8.3,1x), 7f10.6)")
c    .     d, dm, sm_term,  sm_sum1, sm_sum2, sm_sum3,
c    .     sm_sum1_last, sm_sum2_last, sm_sum3_last
	if(ndir.ge.0.and.dm.ge.dmpsr) go to 30	! Reached pulsar's DM?
	if(ndir.lt.0.and.d.ge.dist) go to 40	! Reached pulsar's dist?
	sm_sum1_last = sm_sum1
	sm_sum2_last = sm_sum2
	sm_sum3_last = sm_sum3
	sm_sum4_last = sm_sum4
c        write(25,
c     .     "(4(f7.3,1x),f8.4,1x,e10.3,1x,i1,1x,i4,1x,i2)")
c     .     d,x,y,z,ne,sm_term,whicharm,hitclump,hitvoid
10	continue
	stop 'loop limit'

20	limit='>'			! Only lower limit is possible
	dist=d-0.5*dstep
	go to 999

30	dist=d+0.5*dstep - dstep*(dm-dmpsr)/dmstep  ! Interpolate last step

        if(ncount .lt. 10) then
	  dstep  = dstep / 10.
	  go to 5
	endif
	go to 999

40	dmpsr=dm-dmstep*(d+0.5*dstep-dist)/dstep
        if(ncount .lt. 10) then
	   dstep = dstep / 10.
	   go to 5
	endif

999	continue

c normalize the mean distances:


	if(ldr_path .gt. 0.) then
          ldr_dist = ldr_dist / (ldr_path / dstep)
	endif
	if(lhb_path .gt. 0.) then
          lhb_dist = lhb_dist / (lhb_path / dstep)
	endif
	if(lsb_path .gt. 0.) then
          lsb_dist = lsb_dist / (lsb_path / dstep)
	endif
	if(loopI_path .gt. 0.) then
          loopI_dist = loopI_dist / (loopI_path / dstep)
	endif

        dd = d+0.5*dstep-dist
c subtract dd from armpath for latest arm (or iterarm) at end of LOS
	armpaths(whicharm) = armpaths(whicharm)-dd

	do i=1,narmsmax1
	  armdistances(i) =
     .        armdistances(i) / (max(1.,armpaths(i)/dstep))	! mean distance of arm
	enddo
	dm1 = dm1 * dstep_pc
	dm2 = dm2 * dstep_pc
	dma = dma * dstep_pc
	dmgc = dmgc * dstep_pc
	dmlism = dmlism * dstep_pc
	dmcN = dmcN * dstep_pc
	dmvN = dmvN * dstep_pc

c       dsm = sm_term * (d+0.5*dstep - dist)
c       dsm = sm_term * dd

c       sm_sum2 = sm_sum2 - dsm * d
c       sm_sum3 = sm_sum3 - dsm * d**2
c       sm_sum4 = sm_sum4 - dsm * d**1.67

c       sm_sum1 = sm_sum1 - dsm
c       write(99,*) 'dmdsm: sm_term, sm_sum1, sm_sum1_last = ',
c    .    sm_term, sm_sum1, sm_sum1_last

c	write(6,*) 'dmdsm: dsum1, sm_term = ',
c    .     sm_sum1-sm_sum1_last, sm_term

 	sm_sum1 = sm_sum1 - dd*(sm_sum1-sm_sum1_last)/dstep
	sm_sum2 = sm_sum2 - dd*(sm_sum2-sm_sum2_last)/dstep
	sm_sum3 = sm_sum3 - dd*(sm_sum3-sm_sum3_last)/dstep
	sm_sum4 = sm_sum4 - dd*(sm_sum4-sm_sum4_last)/dstep

c       sm_sum2 = sm_sum2 - dsm * dist
c       sm_sum3 = sm_sum3 - dsm * dist**2
c       sm_sum4 = sm_sum4 - dsm * dist**1.67


        sm = sm_factor * dstep * sm_sum1
        smtau =
     +     6. * sm_factor * dstep * (sm_sum2 / dist - sm_sum3 / dist**2)
        smtheta =
     +     3. * sm_factor * dstep * (sm_sum1 + sm_sum3 / dist**2 -
     +     2. * sm_sum2 / dist)
        smiso = sm_factor * dstep * sm_sum4

	sm1 = sm1 * sm_factor * dstep
	sm2 = sm2 * sm_factor * dstep
	sma = sma * sm_factor * dstep
	smgc = smgc * sm_factor * dstep
	smlism = smlism * sm_factor * dstep
	smcN = smcN * sm_factor * dstep
	smvN = smvN * sm_factor * dstep

c       write(24,*) dm1, dm2, dma, dmgc, dmlism, dmcN, dmvN, dm
c        write(24,"(a,a)") 'LISM path lengths (kpc)',
c     .    ' with weighting hierarchy LHB:LOOPI:LSB:LDR'
c        write(24,"(t15, a)") '  LHB     LoopI     LSB      LDR'
c	write(24, "(t3, a, t15, 4(f6.3, 3x))") 'Length',
c     .         lhb_path, loopI_path, lsb_path, ldr_path
c	write(24, "(t3, a, t15, 4(f6.3, 3x))") 'Mean Dist.',
c     .         lhb_dist, loopI_dist, lsb_dist, ldr_dist


c        write(24,"(a)") 'Fractional contributions to DM:'
c        write(24,"(a,a)")
c     .  '  outer   inner    arms     gc    lism',
c     .  '    clumps  voids       DM'
c        write(24,"(7(f7.3,1x), f10.3)")
c     .              dm1/dm, dm2/dm, dma/dm, dmgc/dm,
c     .              dmlism/dm, dmcN/dm, dmvN/dm, dm
c        write(24,"(a)") 'Fractional contributions to SM:'
c        write(24,"(a,a)")
c     .  '  outer   inner    arms     gc    lism',
c     .  '    clumps  voids       SM'
c        write(24,"(7(f7.3,1x), e10.3)")
c     .              sm1/sm, sm2/sm, sma/sm, smgc/sm,
c     .              smlism/sm, smcN/sm, smvN/sm, sm
c        write(24,"(a)") 'Path lengths through spiral arms:'
c        write(24,"(t1,a,t10, a, t30, a)")
c     .      'Arm','Mean Distance','Path Length    (arm=0 => interarm)'
c	do i=1,narmsmax1
c          write(24,"(i2,t10,f8.3,t30,f8.3)")
c     .       i-1, armdistances(i), armpaths(i)
c	enddo
c	close(24)
c	close(25)
        return
	end

c density.NE2001.f
c final version of NE2001
c returns densities, F parameters and weights of the various components

      SUBROUTINE DENSITY_2001(x,y,z,
     .                ne1,ne2,nea,negc,nelism,necN,nevN,
     .                F1, F2, Fa, Fgc, Flism, FcN, FvN,
     .                whicharm, wlism, wLDR, wLHB, wLSB, wLOOPI,
     .                hitclump, hitvoid, wvoid, ip,lip)

c----------------------------------------------------------------------------
c  Returns seven components of the free electron density of the
c  interstellar medium at Galactic location (x,y,z).
c  Calling arguments:
c  input:
c	x, y, z = galactocentric location (kpc)
c       Right-handed coordinate system
c       x is in l=90 direction
c       y is in l=180 direction
c       The sun is at (x,y,z) = (0,R0,0)
c  output:
c    electron densities in cm^{-3}:
c	ne1:	outer, thick disk
c	ne2:	inner, thin disk (annular in form)
c	nea:	spiral arms
c	negc:   galactic center component
c       nelism: local ISM component
c       necN:   contribution from discrete 'clumps'
c       nevN:   contribution from voids
c    fluctuation parameters (one for each ne component):
c       F1, F2, Fa, Fgc, Flism, FcN, FvN
c    flags:
c       whicharm: which of the 5 arms x,y,z is in (0 for interarm region)
c          wlism: 1 if x,y,z is in any of the four LISM components
c           wLDR: 1 if in LDR, 0 if not
c           wLHB: 1 if in LHB, 0 if not
c           wLSB: 1 if in LSB, 0 if not
c         wLOOPI: 1 if in LoopI, 0 if not
c       (nb: nelism is calculated according to LHB:LOOPI:LSB:LDR)
c       hitclump: clump number that x,y,z is in (0 if none)
c        hitvoid: void number that x,y,z is in (0 if none)
c 25 May 2002
c based on routines from TC93 and test routines from 1999-2002 by JMC.
c----------------------------------------------------------------------------
	implicit none
        real x,y,z
        real ne1,ne2,nea,negc,nelism,necN,nevN
        real F1, F2, Fa, Fgc, Flism, FcN, FvN
        integer whicharm,wlism,wLDR,wLHB,wLSB,wLOOPI,hitclump,hitvoid
        integer wvoid
        character ip*(*)
        integer lip

	integer wg1, wg2, wga, wggc, wglism, wgcN, wgvN
	common /modelflags/ wg1, wg2, wga, wggc, wglism, wgcN, wgvN

	logical first
	data first/.true./
	save

c functions:
	real ne_inner
	real ne_outer
        real ne_arms_log_mod
	real ne_gc
        real ne_lism
c subroutines needed:
c	neclumpN
c	nevoidN

	if(first) then			! get parameters first time through
	   call get_parameters(ip,lip)
	   first=.false.
	endif

	ne1 = ne_outer(x,y,z, F1)
	ne2 = ne_inner(x,y,z, F2)
	nea = ne_arms_log_mod(x,y,z,whicharm,Fa,ip,lip)
        negc = ne_gc(x,y,z, Fgc,ip,lip)
        nelism = ne_lism(x,y,z,Flism,wlism,wldr,wlhb,wlsb,wloopI,ip,lip)
        call neclumpN(x,y,z,necN,FcN,hitclump, ip,lip)
        call nevoidN(x,y,z,nevN,FvN,hitvoid,wvoid, ip, lip)

c       write(21, "(3(f8.2,1x),5(f8.4,1x),i2)")
c    .       x,y,z,ne1,ne2,nea,nelism,negc,
c    .       whicharm

	return
	end

	SUBROUTINE GET_PARAMETERS(ip, lip)
c-----------------------------------------------------------------------
	implicit none
	real rsun
	common /mw/ rsun
	data rsun/8.5/

c control flags for turning components on and off:

	integer wg1, wg2, wga, wggc, wglism, wgcN, wgvN
      character ip*(*)
      integer lip
	common /modelflags/ wg1, wg2, wga, wggc, wglism, wgcN, wgvN

c parameters of large-scale components (inner+outer+arm components):
	real n1h1,h1,A1,F1,n2,h2,A2,F2,na,ha,wa,Aa,Fa
	common/galparams/n1h1,h1,A1,F1,n2,h2,A2,F2,
     .                na,ha,wa,Aa,Fa

c factors for controlling individual spiral arms:
c	narm: 	multiplies electron density (in addition to the`fac'
c		      quantities)
c	warm:	arm width factors that multiply nominal arm width
c	harm:	arm scale height factors
c	farm:	factors that multiply n_e^2 when calculating SM

        character fullname*200

	integer narmsmax
	parameter (narmsmax=5)
	real narm, warm, harm, farm
	common/armfactors/
     .     harm(narmsmax),narm(narmsmax),warm(narmsmax),farm(narmsmax)


	real negc0, Fgc0
        common /gcparms/ negc0, Fgc0

        call mkfile('gal01.inp',fullname,ip,lip)
        open(11,file=fullname,status='old')
c	   open(11,file='gal01.inp',status='old')
	   read(11,*)
	   read(11,*) wg1, wg2, wga, wggc, wglism, wgcN, wgvN
	   read(11,1020) n1h1,h1,A1,F1,n2,h2,A2,F2,
     .          na,ha,wa,Aa,Fa,
     .          narm(1), narm(2), narm(3), narm(4), narm(5),
     .          warm(1), warm(2), warm(3), warm(4), warm(5),
     .          harm(1), harm(2), harm(3), harm(4), harm(5),
     .          farm(1), farm(2), farm(3), farm(4), farm(5)
 1020	   format(7x,f8.0)
	   close(11)
c          write(6,*) 'get_parms: weights: ',
c    .           wg1, wg2, wga, wggc, wglism, wgcN, wgvN
c          write(6,*) 'get_parms: ',
c    .           n1h1,h1,A1,F1,n2,h2,A2,F2,
c    .           na,ha,wa,Aa,Fa,
c    .           narm(1), narm(2), narm(3), narm(4), narm(5),
c    .           warm(1), warm(2), warm(3), warm(4), warm(5),
c    .           harm(1), harm(2), harm(3), harm(4), harm(5),
c    .           farm(1), farm(2), farm(3), farm(4), farm(5)
	return
	end


	REAL FUNCTION NE_ARMS_LOG_MOD(x,y,z, whicharm, Farms, ip, lip)
c-----------------------------------------------------------------------
c  Spiral arms are defined as logarithmic spirals using the
c    parameterization in Wainscoat et al. 1992, ApJS, 83, 111-146.
c  But arms are modified selectively at various places to distort them
c    as needed (08 Aug 2000).
c  Note that arm numbering follows that of TC93 for the four large arms
c (after remapping).
c  The local spiral arm is number 5.
c  06 Apr 02:   removed TC type modifications of arms 2,3 (fac calculations)
c  		and replaced with new versions.  Data for these are hard wired.

	implicit none
	real x, y, z
	integer whicharm
	real Farms

	integer whicharm_spiralmodel

        real n1h1,h1,A1,F1,n2,h2,A2,F2,na,ha,wa,Aa,Fa
        common/galparams/n1h1,h1,A1,F1,n2,h2,A2,F2,
     .                na,ha,wa,Aa,Fa

c see get_parameters for definitions of narm, warm, harm.

        integer narmsmax
        parameter (narmsmax=5)
	real narm, warm, harm, farm
	common/armfactors/
     .     harm(narmsmax),narm(narmsmax),warm(narmsmax),farm(narmsmax)


        real rsun
        common /mw/ rsun

        real rad
        parameter(rad=57.29577 95130 823)

	integer ks
        data ks/3/

	integer NN
	data NN/7/

	integer NNmax
	parameter(NNmax=20)

	integer narms
	parameter(narms=5)

	real aarm(narms), rmin(narms), thmin(narms), extent(narms)

	integer armmap(5)		! for remapping from Wainscoat
	data armmap/1, 3, 4, 2, 5/	! order to TC93 order, which is
					! from GC outwards toward Sun.
 	integer NNj(narms)
	data NNj/20, 20, 20, 20, 20/

	real th1(NNmax,narms),r1(NNmax,narms)

	real arm
	integer kmax, narmpoints, ncoord
        parameter(narmpoints=500, ncoord=2)
	dimension arm(narms,narmpoints,ncoord),kmax(narms)

	real nea
        integer j, k, n, jj

	logical first
	data first /.true./

	real rr
 	real dth, th, r
	real smin, sminmin
	real sq, sq1, sq2, sqmin
	integer kk, kmi, kma, kl
	real emm, ebb, exx, eyy,  test
	real ga, fac
	real thxy
	real arg
	real th3a, th3b, fac3min, test3
	real th2a, th2b, fac2min, test2

        character fullname*200
        character ip*(*)
        integer lip

c function:
	real sech2
	save

	rr=sqrt(x**2 + y**2)
	if(first) then			! Reconstruct spiral arm axes

c read arm parameters:
        call mkfile('ne_arms_log_mod.inp',fullname,ip,lip)
        open(11,file=fullname,status='old')
c        open(11, file='ne_arms_log_mod.inp', status='old')
c       write(6,*) 'ne_arms_log_mod.inp:'
        read(11,*)
        read(11,*)
        do j=1,narms
          read(11,*) aarm(j), rmin(j), thmin(j), extent(j)
c         write(6,*) aarm(j), rmin(j), thmin(j), extent(j)
        enddo
        close(11)

	do j=1,narms			! fill sampling array
	  do n=1,NNj(j)
	   th1(n,j) = thmin(j)+(n-1)*extent(j)/(NNj(j)-1.) 	! rad
	   r1(n,j) = rmin(j)*exp((th1(n,j)-thmin(j))/aarm(j))
	   th1(n,j) = th1(n,j)*rad				! deg
c *** begin sculpting spiral arm 2 == TC arm 3***
           if(armmap(j) .eq. 3) then
 	   if(th1(n,j) .gt. 370. .and. th1(n,j) .le. 410.) then
     	      r1(n,j) = r1(n,j) *
     .           (1. + 0.04*cos((th1(n,j)-390.)*180./(40.*rad)))
c    .           (1. + 0.01*cos((th1(n,j)-390.)*180./(40.*rad)))
	   endif
 	   if(th1(n,j) .gt. 315. .and. th1(n,j) .le. 370.) then
     	      r1(n,j) = r1(n,j) *
     .           (1. - 0.07*cos((th1(n,j)-345.)*180./(55.*rad)))
c    .           (1.0 - 0.08*cos((th1(n,j)-345.)*180./(55.*rad)))
	   endif
 	   if(th1(n,j) .gt. 180. .and. th1(n,j) .le. 315.) then
	      r1(n,j) = r1(n,j) *
c    ,           (1 + 0.13*cos((th1(n,j)-260.)*180./(135.*rad)))
     ,           (1 + 0.16*cos((th1(n,j)-260.)*180./(135.*rad)))
	   endif
	   endif
c *** begin sculpting spiral arm 4 == TC arm 2***
           if(armmap(j) .eq. 2) then
 	   if(th1(n,j) .gt. 290. .and. th1(n,j) .le. 395.) then
     	      r1(n,j) = r1(n,j) *
c    .            1.
     .           (1. - 0.11*cos((th1(n,j)-350.)*180./(105.*rad)))
	   endif
	   endif
c *** end arm sculpting ***
c	   write(6,*) j,n, th1(n,j), r1(n,j)
	  enddo
	enddo

c        open(11,file='log_arms.out', status='unknown')
c        write(11,*) 'arm  n   xa     ya'

	   do 21 j=1,narms
	      dth=5.0/r1(1,j)
	      th=th1(1,j)-0.999*dth
	      call cspline(th1(1,j),r1(1,j),-NNj(j),th,r)
c	      write(6,*) 'doing arm ', j, ' with ', NNj(j), ' points',
c    .               dth
c	      write(6,*) (th1(k,j), r1(k,j), k=1,NNj(j))
	      do 10 k=1,narmpoints-1
		 th=th+dth
		 if(th.gt.th1(NNj(j),j)) go to 20
		 call cspline(th1(1,j),r1(1,j),NNj(j),th,r)
		 arm(j,k,1)=-r*sin(th/rad)
		 arm(j,k,2)= r*cos(th/rad)
c                 write(11,"(1x,i2,1x,i3,1x,2(f7.3,1x))")
c     .              j,k,arm(j,k,1),arm(j,k,2)
 10	      continue
 20	   continue
	   kmax(j)=k
 21        continue
c	   close(11)

	first = .false.
	endif
c
c Get spiral arm component:  30 do loop finds a coarse minimum distance
c from line of sight to arm; 40 do loop finds a fine minimum distance
c from line of sight to arm; line 35 ensures that arm limits are not
c exceeded; linear interpolation beginning at line 41 finds the
c minimum distance from line of sight to arm on a finer scale than gridding
c of arms allows (TJL)

	nea=0.0
        ga = 0.
	whicharm = 0
        whicharm_spiralmodel = 0
        sminmin = 1.e10
	thxy = atan2(-x, y) * rad		! measured ccw from +y axis
						! (different from tc93 theta)
	if(thxy.lt.0.) thxy=thxy+360.
	if(abs(z/ha).lt.10.) then
           do 50 j=1,narms
              jj = armmap(j)
              sqmin=1.e10
              do 30 k=1+ks,kmax(j)-ks,2*ks+1
                 sq=(x-arm(j,k,1))**2 + (y-arm(j,k,2))**2
                 if(sq.lt.sqmin) then
                    sqmin=sq
                    kk=k
                 endif
 30           continue
 35           kmi = max(kk-2*ks, 1)
              kma = min(kk+2*ks, kmax(j))
              do 40 k=kmi,kma
                 sq=(x-arm(j,k,1))**2 + (y-arm(j,k,2))**2
                 if(sq.lt.sqmin) then
                    sqmin=sq
                    kk=k
                 endif
 40           continue
 41           if (kk.gt.1.and.kk.lt.kmax(j)) then
                 sq1 = (x - arm(j,kk-1,1))**2 + (y - arm(j,kk-1,2))**2
                 sq2 = (x - arm(j,kk+1,1))**2 + (y - arm(j,kk+1,2))**2
                 if (sq1.lt.sq2) then
                    kl = kk - 1
                 else
                    kl = kk + 1
                 endif
                 emm = (arm(j,kk,2) - arm(j,kl,2))
     $                /(arm(j,kk,1) - arm(j,kl,1))
                 ebb = arm(j,kk,2) - emm*arm(j,kk,1)
                 exx = (x + emm*y - emm*ebb)/(1.0 + emm**(2))
                 test = (exx - arm(j,kk,1))/(arm(j,kl,1) - arm(j,kk,1))
                 if (test.lt.0.0.or.test.gt.1.0) exx = arm(j,kk,1)
                 eyy = emm*exx + ebb
              else
                 exx = arm(j,kk,1)
                 eyy = arm(j,kk,2)
              endif
              sqmin = (x - exx)**(2) + (y - eyy)**(2)
	      smin=sqrt(sqmin)		! Distance of (x,y,z) from arm axis
c           write(23,"(4(f5.2,1x),i2,1x,3(f8.3,1x))")
c    .        x,y,z,rr,j,exx,eyy,smin
	    if(smin.lt.3*wa) then		! If (x,y,z) is close to this
	      ga=exp(-(smin/(warm(jj)*wa))**2)	! arm, get the arm weighting factor
	      if(smin .lt. sminmin) then
		whicharm_spiralmodel = j
		sminmin = smin
	      endif
	      if(rr.gt.Aa) then 	! Galactocentric radial dependence of arms
		ga=ga*sech2((rr-Aa)/2.0)
c		write(6,*) 'd99a: rr,Aa,sech2() = ',
c                 rr, Aa, sech2((rr-Aa)/2.0)
	      endif


c arm3 reweighting:
	      th3a=320.
	      th3b=390.
	      th3b=370.
	      th3a=290.
	      th3b=363.
	      th3b=363.
	      fac3min=0.0
	      test3 = thxy-th3a
	      if(test3 .lt.0.) test3 = test3+360.
	      if(jj.eq.3
     .            .and. 0. .le. test3
     .            .and. test3 .lt. (th3b-th3a))
     .        then
	        arg=6.2831853*(thxy-th3a)/(th3b-th3a)
c		fac = (3.0 + cos(arg))/4.0
		fac = (1.+fac3min + (1.-fac3min)*cos(arg))/2.
		fac = fac**4.0
c		write(90,*) x, y, thxy, th3a, th3b, test3, fac
		ga=ga*fac
	      endif


c arm2 reweighting:
c    first: as in tc93 (note different definition of theta)
	      th2a=35.
	      th2b=55.
	      test2 = thxy-th2a
	      fac = 1.

	      if(jj.eq.2
     .            .and. 0. .le. test2
     .            .and. test2 .lt. (th2b-th2a))
     .        then
	         fac=1.+ test2/(th2b-th2a)
		 fac = 1.		!**** note turned off
                 ga=ga*fac
	      endif

	      if (jj.eq.2 .and. test2 .gt. (th2b-th2a)) then
	        fac = 2.
		fac = 1.		!**** note turned off
                ga=ga*fac
	      endif
c    second:  weaken the arm in a short range:
	      th2a=340.
	      th2b=370.
c note fix does nothing if fac2min = 1.0
	      fac2min=0.1
	      test2 = thxy-th2a
	      if(test2 .lt.0.) test2 = test2+360.

	      if(jj.eq.2
     .            .and. 0. .le. test2
     .            .and. test2 .lt. (th2b-th2a))
     .        then
	        arg=6.2831853*(thxy-th2a)/(th2b-th2a)
		fac = (1.+fac2min + (1.-fac2min)*cos(arg))/2.
c		fac = fac**3.5
c		write(90,*) x, y, thxy, th2a, th2b, test2, fac
		ga=ga*fac
	      endif

	      nea=nea +
     .            narm(jj)*na*ga*sech2(z/(harm(jj)*ha))   ! Add this arm contribution
	    endif
50	  continue
	endif

        ne_arms_log_mod = nea
        Farms = 0
	if(whicharm_spiralmodel .eq. 0) then
	  whicharm = 0
	else
	  whicharm = armmap(whicharm_spiralmodel)	! remap arm number
	  Farms = Fa * farm(whicharm)
	endif
	return
	end



	REAL FUNCTION NE_OUTER(x,y,z, F_outer)
c-----------------------------------------------------------------------
c Thick disk component:
	implicit none
	real x,y,z, F_outer

	real n1h1,h1,A1,F1,n2,h2,A2,F2,na,ha,wa,Aa,Fa
	common/galparams/n1h1,h1,A1,F1,n2,h2,A2,F2,
     .                na,ha,wa,Aa,Fa


        real pihalf, rad
	parameter(pihalf=3.14159 26535 89793 23846 264/2.0)
	parameter(rad=57.29577 95130 823)

	real rsun
	common /mw/ rsun

	real g1
	real sech2
	real rr, suncos, ne1

	logical first
	data first /.true./
	save

c 	g1=sech2(rr/A1)/sech2(8.5/A1)		! TC93 function
	rr=sqrt(x**2 + y**2)
	suncos = cos(pihalf*rsun/A1)
	if (rr.gt.A1) then
	   g1 = 0.0
	else
	   g1 = cos(pihalf*rr/A1)/suncos
	endif
 	ne1=(n1h1/h1)*g1*sech2(z/h1)
	ne_outer = ne1
        F_outer = F1

	return
	end

	REAL FUNCTION NE_INNER(x,y,z, F_inner)
c-----------------------------------------------------------------------
c Thin disk (inner Galaxy) component:
c (referred to as 'Galactic center component' in circa TC93 density.f)
	implicit none
	real x,y,z, F_inner

	real n1h1,h1,A1,F1,n2,h2,A2,F2,na,ha,wa,Aa,Fa
	common/galparams/n1h1,h1,A1,F1,n2,h2,A2,F2,
     .                na,ha,wa,Aa,Fa
	real g2, rr, rrarg
	real sech2
	real ne2
	save

	g2=0.0
	rr=sqrt(x**2 + y**2)
	rrarg=((rr-A2)/1.8)**2
	if(rrarg.lt.10.0) g2=exp(-rrarg)
	ne2=n2*g2*sech2(z/h2)

	ne_inner = ne2
        F_inner = F2
	return
	end


      REAL FUNCTION NE_GC(x, y, z, F_gc, ip, lip)
c-----------------------------------------------------------------------
c     Determines the contribution of the Galactic center to the free
c     electron density of the interstellar medium at Galactic location
c     (x,y,z).  Combine with `fluctuation' parameter to obtain the
c     scattering measure.
c
c     NOTE: This is for the hyperstrong scattering region in the
c     Galactic center.  It is distinct from the inner Galaxy
c     (component 2) of the TC93 model.
c
c     Origin of coordinate system is at Galactic center; the Sun is at
c     (x,y,z) = (0,R0,0), x is in l=90 direction
c
c     Based on Section 4.3 of Lazio & Cordes (1998, ApJ, 505, 715)
c
c Input:
c REAL X - location in Galaxy [kpc]
c REAL Y - location in Galaxy [kpc]
c REAL Z - location in Galaxy [kpc]
c
c COMMON:
c REAL NEGC0 - nominal central density
c
c PARAMETERS:
c REAL RGC - radial scale length of Galactic center density enhancement
c REAL HGC - z scale height of Galactic center density enhancement
c
c Output:
c REAL NE_GC - Galactic center free electron density contribution [cm^-3]
c-----------------------------------------------------------------------
c
      implicit none
      real x, y, z, F_gc

      real rgc, hgc
      real xgc, ygc, zgc
c     parameter (xgc=-0.010, ygc=0., zgc=-0.020)
c     parameter (rgc=0.145)
c     parameter (hgc=0.026)

      character fullname*200
      character ip*(*)
      integer lip

      real rr, zz
      real arg

      real negc0, Fgc0
      common /gcparms/ negc0, Fgc0

      real n1h1,h1,A1,F1,n2,h2,A2,F2,na,ha,wa,Aa,Fa
      common /galparams/ n1h1,h1,A1,F1,n2,h2,A2,F2,
     .                na,ha,wa,Aa,Fa

      logical first
      data first /.true./
      save

      ne_gc = 0.
      F_gc = 0.

      if(first) then
        call mkfile('ne_gc.inp',fullname,ip,lip)
        open(11,file=fullname,status='old')
c        open(11, file='ne_gc.inp', status='old')
          read(11,*)
          read(11,*) xgc, ygc, zgc
	  read(11,*) rgc
	  read(11,*) hgc
	  read(11,*) negc0
          read(11,*) Fgc0
        close(11)
	first = .false.
      endif


c GALACTOCENTRIC RADIUS

      rr = sqrt( (x-xgc)**2 + (y-ygc)**2)
      if(rr .gt. rgc) return				! truncate at 1/e point

c Z-HEIGHT.

      zz =abs(z-zgc)
      if(zz .gt. hgc) return
      arg = (rr/rgc)**2 + (zz/hgc)**2
      if(arg .le. 1.) then
         ne_gc = negc0
	 F_gc = Fgc0
c        write(21,*) 'ne_gc: rr,zz,arg,ne_gc,F_gc ',
c    .                rr, zz, arg, ne_gc, F_gc
      endif

      return
      end


c
c%%%%%%%%%%%%%%%%%%%%%%%%%  cspline.f  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
	subroutine cspline(x,y,nn,xout,yout)
	integer nnmax
	parameter(nnmax=20)
	real x(nnmax),y(nnmax),y2(nnmax),u(nnmax)
	save

	if(nn .gt. nnmax) then
	  write(6,*)
     .    ' too many points to spline. Change parameter statement'
	  write(6,*)
     .    ' in cspline'
	endif

	n=abs(nn)
	if(nn.lt.0) then
	  y2(1)=0.
	  u(1)=0.
	  do 10 i=2,n-1
	    sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
	    p=sig*y2(i-1)+2.
	    y2(i)=(sig-1.)/p
	    u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))
     +        /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
10        continue
	  qn=0.
	  un=0.
	  y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
	  do 20 k=n-1,1,-1
20	  y2(k)=y2(k)*y2(k+1)+u(k)
	endif

	klo=1
	khi=n
30	if (khi-klo.gt.1) then
	  k=(khi+klo)/2
	  if(x(k).gt.xout)then
	    khi=k
	  else
	    klo=k
	  endif
	goto 30
	endif
	h=x(khi)-x(klo)
	if (h.eq.0.) pause 'bad x input.'
	a=(x(khi)-xout)/h
	b=(xout-x(klo))/h
	yout=a*y(klo)+b*y(khi)+
     +    ((a**3-a)*y2(klo)+(b**3-b)*y2(khi))*(h**2)/6.
	return
	end

	REAL FUNCTION SECH2(z)
c-----------------------------------------------------------------------
	sech2=0.0
	if(abs(z).lt.20.0) sech2=(2.0/(exp(z)+exp(-z)))**2
	return
	end

c23456789012345678901234567890123456789012345678901234567890123456789012
c routines to calculate the electron density for the
c Local Interstellar Medium
c
c JMC 26 August-11 Sep. 2000
c     25 October 2001: modified to change weighting scheme
c                      so that the ranking is LHB: LSB: LDR
c                      (LHB overrides LSB and LDR; LSB overrides LDR)
c     16 November 2001: added Loop I component with weighting scheme
c		        LHB:LOOPI:LSB:LDR
c		        LHB   overides everything,
c			LOOPI overrides LSB and LDR
c			LSB   overrides LDR
c			LISM  overrides general Galaxy
c     20 November 2001: The LOOPI component is truncated below z=0
c
c after discussions with Shami Chatterjee
c the sizes, locations and densities of the LISM components
c are based largely on work by Toscano et al. 1999
c and Bhat et al. 1999

	real function ne_LISM(x,y,z,FLISM,wLISM,
     .                    wldr, wlhb,wlsb,wloopI, ip, lip)
	implicit none
	real x,y,z,FLISM
        integer wLISM

	character fullname*200
        character ip*(*)
        integer lip

	real
     .    aldr,bldr,cldr,xldr,yldr,zldr,thetaldr,neldr0,Fldr,
     .    alsb,blsb,clsb,xlsb,ylsb,zlsb,thetalsb,nelsb0,Flsb,
     .    alhb,blhb,clhb,xlhb,ylhb,zlhb,thetalhb,nelhb0,Flhb,
     .    xlpI,ylpI,zlpI,rlpI,drlpI,nelpI,FlpI,dnelpI,dFlpI
        common/nelismparms/
     .    aldr,bldr,cldr,xldr,yldr,zldr,thetaldr,neldr0,Fldr,
     .    alsb,blsb,clsb,xlsb,ylsb,zlsb,thetalsb,nelsb0,Flsb,
     .    alhb,blhb,clhb,xlhb,ylhb,zlhb,thetalhb,nelhb0,Flhb,
     .    xlpI,ylpI,zlpI,rlpI,drlpI,nelpI,FlpI,dnelpI,dFlpI


        logical first
	data first /.true./

c functions:

	real neLDRQ1
	real neLSB
	real neLHB2
	real neLOOPI

c other variables:

	real nelsbxyz, nelhbxyz, neldrq1xyz, neloopIxyz
	real FLDRQ1r, FLSBr, FLHBr, FLOOPIr		! 'r' for returned value
	integer wLDR, wLSB, wLHB, wLOOPI

	if(first) then					! read parameters for LISM
        call mkfile('nelism.inp',fullname,ip,lip)
        open(11,file=fullname,status='old')
c	  open(11,file='nelism.inp',status='unknown')
	  read(11,*)
	  read(11,*) aldr,bldr,cldr
	  read(11,*) xldr,yldr,zldr
	  read(11,*) thetaldr,neldr0,Fldr

	  read(11,*) alsb,blsb,clsb
	  read(11,*) xlsb,ylsb,zlsb
	  read(11,*) thetalsb,nelsb0,Flsb

	  read(11,*) alhb,blhb,clhb
	  read(11,*) xlhb,ylhb,zlhb
	  read(11,*) thetalhb,nelhb0,Flhb

	  read(11,*) xlpI,ylpI,zlpI
	  read(11,*) rlpI,drlpI
	  read(11,*) nelpI,dnelpI,FlpI,dFlpI

c   	  write(99,*) xlpI,ylpI,zlpI
c 	  write(99,*) rlpI,drlpI
c 	  write(99,*) nelpI,dnelpI,FlpI,dFlpI

	  first=.false.
	endif

	neldrq1xyz = neLDRQ1(x,y,z,FLDRQ1r,wLDR)	! low density region in Q1
	nelsbxyz   = neLSB(x,y,z,FLSBr,wLSB)		! Local Super Bubble
	nelhbxyz = neLHB2(x,y,z,FLHBr,wLHB)		! Local Hot Bubble
	neloopIxyz = neLOOPI(x,y,z,FLOOPIr,wLOOPI)	! Loop I


c weight the terms so that the LHB term overrides the other
c terms (we want the density to be low in the LHB, lower than
c in the other terms.

 	ne_LISM =   (1-wLHB)   *
     .           (
     .             (1-wLOOPI) * (wLSB*nelsbxyz + (1-wLSB)*neldrq1xyz)
     .         +     wLOOPI * neloopIxyz
     .           )
     .         +     wLHB  * nelhbxyz

 	FLISM = (1-wLHB) *
     .          (
     .             (1-wLOOPI) * (wLSB*FLSBr + (1-wLSB)*FLDRQ1r)
     .         +     wLOOPI * FLOOPIr
     .          )
     .         +     wLHB  * FLHBr

c return the maximum weight of any of the terms for
c combining with additional terms external to this routine.

	wLISM = max(wLOOPI, max(wLDR, max(wLSB, wLHB)))

c temporary next 3 lines:
c	ne_LISM = nelhbxyz
c	flism = flhb
c	wlism = wlhb

c	write(97,"(9(f8.3,1x))") wLDR, wLSB, wLHB,
c    .      nelsbxyz, neldrq1xyz,nelhbxyz,
c    .      flsb, fldrq1, flhb


	return
	end

	real function  neLDRQ1(x,y,z,FLDRQ1r,wLDRQ1) 	! Low Density Region in Q1
	implicit none
	real x,y,z,FLDRQ1r
	integer wLDRQ1

c input:
c 	x,y,z = coordinates w.r.t. Galaxy as in TC93, CL00
c output:
c	neLDRQ1 = electron density in local hot bubble that
c	        is modeled as an ellipsoidal trough.
c	FLDRQ1 = fluctuation parameter
c	wLDRQ1  = weight of LDRQ1 component used to combine
c		with other components of electron density.
c		wLDRQ1 =  1  at and inside the annular ridge
c		     <  1  outside the annular ridge
c	             -> 0  far outside the annular ridge
c	e.g. total electron density would be evaluated as
c            ne = (1-wLDRQ1)*ne_other + neLDRQ1
	real
     .    aldr,bldr,cldr,xldr,yldr,zldr,thetaldr,neldr0,Fldr,
     .    alsb,blsb,clsb,xlsb,ylsb,zlsb,thetalsb,nelsb0,Flsb,
     .    alhb,blhb,clhb,xlhb,ylhb,zlhb,thetalhb,nelhb0,Flhb,
     .    xlpI,ylpI,zlpI,rlpI,drlpI,nelpI,FlpI,dnelpI,dFlpI
        common/nelismparms/
     .    aldr,bldr,cldr,xldr,yldr,zldr,thetaldr,neldr0,Fldr,
     .    alsb,blsb,clsb,xlsb,ylsb,zlsb,thetalsb,nelsb0,Flsb,
     .    alhb,blhb,clhb,xlhb,ylhb,zlhb,thetalhb,nelhb0,Flhb,
     .    xlpI,ylpI,zlpI,rlpI,drlpI,nelpI,FlpI,dnelpI,dFlpI


	real aa,bb,cc			! scales of ellipsoidal ridge
	real netrough			! ne of annulus, trough
	real Ftrough			! fluctuation parameters
c	real xldr, yldr, zldr		! center of ellipsoid
	real theta 			! position angle of major axis,
					!    measured from x axis
					!    (x axis points toward l=90)

	real q

	real ap, bp, cp, dp
	real s, c

        real radian
        parameter(radian = 57.29577951)

	logical first
	data first /.true./
	save

c	data aa, bb, cc       	/ 0.6, 0.40, 0.3 / 	! GUESS
c	data xldr, yldr, zldr 	/ 0.6, 7.86, 0. / 	! GUESS
c	data theta		/ -45. / 		! GUESS

c	data netrough  /0.010/		! GUESS
c	data Ftrough   / 2./		! GUESS

	aa=aldr
	bb=bldr
	cc=cldr
	theta=thetaldr
	netrough =neldr0
	Ftrough=Fldr

	if(first) then
	  s = sin(theta/radian)
	  c = cos(theta/radian)
	  ap = (c/aa)**2 + (s/bb)**2
	  bp = (s/aa)**2 + (c/bb)**2
	  cp = 1./cc**2
	  dp =  2.*c*s*(1./aa**2 - 1./bb**2)
	  first = .false.
c	  write(6,*) aa,bb,cc,theta,ap,bp,cp,dp
	endif

	neLDRQ1 = 0.
	wLDRQ1 = 0
	FLDRQ1r = 0.
	q = (x-xldr)**2*ap
     .    + (y-yldr)**2*bp
     .    + (z-zldr)**2*cp
     .    + (x-xldr)*(y-yldr)*dp
	  if(q .le. 1.0) then	! inside
	    neLDRQ1 = netrough
	    FLDRQ1r = Ftrough
	    wLDRQ1 = 1
	  endif

	return
	end


	real function neLSB(x,y,z,FLSBr,wLSB)	! Local Super Bubble
	implicit none
	real x,y,z,FLSBr
	integer wLSB
c input:
c 	x,y,z = coordinates w.r.t. Galaxy as in TC93, CL00
c output:
c	neLSB = electron density in local hot bubble that
c	        is modeled as an ellisoidal trough.
c	FLSB = fluctuation parameter
c	wLSB  = weight of LSB component used to combine
c		with other components of electron density.
c		wLSB =  1  at and inside the annular ridge
c		     <  1  outside the annular ridge
c	             -> 0  far outside the annular ridge
c	e.g. total electron density would be evaluated as
c            ne = (1-wLSB)*ne_other + neLSB

	real
     .    aldr,bldr,cldr,xldr,yldr,zldr,thetaldr,neldr0,Fldr,
     .    alsb,blsb,clsb,xlsb,ylsb,zlsb,thetalsb,nelsb0,Flsb,
     .    alhb,blhb,clhb,xlhb,ylhb,zlhb,thetalhb,nelhb0,Flhb,
     .    xlpI,ylpI,zlpI,rlpI,drlpI,nelpI,FlpI,dnelpI,dFlpI
        common/nelismparms/
     .    aldr,bldr,cldr,xldr,yldr,zldr,thetaldr,neldr0,Fldr,
     .    alsb,blsb,clsb,xlsb,ylsb,zlsb,thetalsb,nelsb0,Flsb,
     .    alhb,blhb,clhb,xlhb,ylhb,zlhb,thetalhb,nelhb0,Flhb,
     .    xlpI,ylpI,zlpI,rlpI,drlpI,nelpI,FlpI,dnelpI,dFlpI


	real aa,bb,cc			! scales of ellipsoidal ridge
	real netrough			! ne of annulus, trough
	real Ftrough			! fluctuation parameters
c	real xlsb, ylsb, zlsb		! center of ellipsoid
	real theta 			! position angle of major axis,
					!    measured from x axis
					!    (x axis points toward l=90)

	real q

	real ap, bp, cp, dp
	real s, c

        real radian
        parameter(radian = 57.29577951)

	logical first
	data first /.true./
	save

c	data aa, bb, cc       	/ 0.6, 0.25, 0.3 / 	! GUESS
c	data xlsb, ylsb, zlsb 	/ -0.7, 9.0, 0. / 	! GUESS
c	data theta		/ 150. / 		! GUESS

c	data netrough  /0.01/		! GUESS
c	data Ftrough   / 1./		! GUESS

	aa=alsb
	bb=blsb
	cc=clsb
	theta=thetalsb
	netrough=nelsb0
	Ftrough=Flsb

	if(first) then
	  s = sin(theta/radian)
	  c = cos(theta/radian)
	  ap = (c/aa)**2 + (s/bb)**2
	  bp = (s/aa)**2 + (c/bb)**2
	  cp = 1./cc**2
	  dp =  2.*c*s*(1./aa**2 - 1./bb**2)
	  first = .false.
c	  write(6,*) aa,bb,cc,theta,ap,bp,cp,dp
	endif

	neLSB = 0.
	wLSB = 0
	FLSBr = 0.
	q = (x-xlsb)**2*ap
     .    + (y-ylsb)**2*bp
     .    + (z-zlsb)**2*cp
     .    + (x-xlsb)*(y-ylsb)*dp
	  if(q .le. 1.0) then	! inside
	    neLSB = netrough
	    FLSBr = Ftrough
	    wLSB = 1
	  endif

	return
	end

	real function neLHB(x,y,z,FLHBr,wLHB)	! Local Hot Bubble
	implicit none
	real x,y,z,FLHBr
	integer wLHB
c input:
c 	x,y,z = coordinates w.r.t. Galaxy as in TC93, CL00
c output:
c	neLHB = electron density in local hot bubble that
c	        is modeled as an ellisoidal trough.
c	FLHB = fluctuation parameter
c	wLHB  = weight of LBH component used to combine
c		with other components of electron density.
c		wLBH =  1  at and inside the annular ridge
c		     <  1  outside the annular ridge
c	             -> 0  far outside the annular ridge
c	e.g. total electron density would be evaluated as
c            ne = (1-wLHB)*ne_other + neLHB

	real
     .    aldr,bldr,cldr,xldr,yldr,zldr,thetaldr,neldr0,Fldr,
     .    alsb,blsb,clsb,xlsb,ylsb,zlsb,thetalsb,nelsb0,Flsb,
     .    alhb,blhb,clhb,xlhb,ylhb,zlhb,thetalhb,nelhb0,Flhb,
     .    xlpI,ylpI,zlpI,rlpI,drlpI,nelpI,FlpI,dnelpI,dFlpI
        common/nelismparms/
     .    aldr,bldr,cldr,xldr,yldr,zldr,thetaldr,neldr0,Fldr,
     .    alsb,blsb,clsb,xlsb,ylsb,zlsb,thetalsb,nelsb0,Flsb,
     .    alhb,blhb,clhb,xlhb,ylhb,zlhb,thetalhb,nelhb0,Flhb,
     .    xlpI,ylpI,zlpI,rlpI,drlpI,nelpI,FlpI,dnelpI,dFlpI


	real aa,bb,cc			! scales of ellipsoidal ridge
	real netrough			! ne of annulus, trough
	real Ftrough			! fluctuation parameters
c	real xlhb, ylhb, zlhb		! center of ellipsoid
	real theta 			! position angle of major axis,
					!    measured from x axis
					!    (x axis points toward l=90)
	real q

	real ap, bp, cp, dp
	real s, c

        real radian
        parameter(radian = 57.29577951)

	logical first
	data first /.true./

c	data aa, bb, cc       	/ 0.15, 0.08, 0.2 /
c	data xlhb, ylhb, zlhb 	/ 0., 8.5, 0. /
c	data theta		/ 135. /

c	data netrough  /0.005/
c	data Ftrough   / 1./

	save

	aa=alhb
	bb=blhb
	cc=clhb
	theta=thetalhb
	netrough=nelhb0
	Ftrough=Flhb


	if(first) then
	  s = sin(theta/radian)
	  c = cos(theta/radian)
	  ap = (c/aa)**2 + (s/bb)**2
	  bp = (s/aa)**2 + (c/bb)**2
	  cp = 1./cc**2
	  dp = 2.*c*s*(1./aa**2 - 1./bb**2)
	  first = .false.
c	  write(6,*) aa,bb,cc,theta,ap,bp,cp,dp
	endif

	neLHB = 0.
	wLHB = 0
	FLHBr = 0.
	q = (x-xlhb)**2*ap
     .    + (y-ylhb)**2*bp
     .    + (z-zlhb)**2*cp
     .    + (x-xlhb)*(y-ylhb)*dp
	  if(q .le. 1.0) then	! inside
	    neLHB = netrough
	    FLHBr = Ftrough
	    wLHB = 1
	  endif

	return
	end



	real function neLHB2(x,y,z,FLHBr,wLHB)	! Local Hot Bubble
c LHB modeled as a cylinder
c the cylinder slants in the y direction vs. z as described by parameter yzslope
c the cylinder cross-sectional size in the 'a' direction (major axis)
c       varies with z, tending to zero at its smallest z point.
	implicit none
	real x,y,z,FLHBr
	integer wLHB
c input:
c 	x,y,z = coordinates w.r.t. Galaxy as in TC93, CL00
c output:
c	neLHB2 = electron density in local hot bubble that
c	        is modeled as an ellisoidal trough.
c	FLHB = fluctuation parameter
c	wLHB  = weight of LBH component used to combine
c		with other components of electron density.
c		wLHB =  1  at and inside the annular ridge
c		     <  1  outside the annular ridge
c	             -> 0  far outside the annular ridge
c	e.g. total electron density would be evaluated as
c            ne = (1-wLHB)*ne_other + neLHB2

	real
     .    aldr,bldr,cldr,xldr,yldr,zldr,thetaldr,neldr0,Fldr,
     .    alsb,blsb,clsb,xlsb,ylsb,zlsb,thetalsb,nelsb0,Flsb,
     .    alhb,blhb,clhb,xlhb,ylhb,zlhb,thetalhb,nelhb0,Flhb,
     .    xlpI,ylpI,zlpI,rlpI,drlpI,nelpI,FlpI,dnelpI,dFlpI
        common/nelismparms/
     .    aldr,bldr,cldr,xldr,yldr,zldr,thetaldr,neldr0,Fldr,
     .    alsb,blsb,clsb,xlsb,ylsb,zlsb,thetalsb,nelsb0,Flsb,
     .    alhb,blhb,clhb,xlhb,ylhb,zlhb,thetalhb,nelhb0,Flhb,
     .    xlpI,ylpI,zlpI,rlpI,drlpI,nelpI,FlpI,dnelpI,dFlpI


	real aa,bb,cc			! scales of ellipsoidal ridge
	real netrough			! ne of annulus, trough
	real Ftrough			! fluctuation parameters
c	real xlhb, ylhb, zlhb		! center of ellipsoid
	real theta 			! slant angle in yz plane of cylinder
					!    measured from z axis
	real qxy, qz

        real radian
        parameter(radian = 57.29577951)

	logical first
	data first /.true./

	real yzslope

	real yaxis
	save

	aa=alhb
	bb=blhb
	cc=clhb
	theta=thetalhb
	netrough=nelhb0
	Ftrough=Flhb


	if(first) then
	  yzslope = tan(theta/radian)
	  first = .false.
	endif

	neLHB2 = 0.
	wLHB = 0
	FLHBr = 0.
	yaxis = ylhb + yzslope*z
c cylinder has cross sectional area = constant for z>0
c area -> 0 for z<0 by letting aa->0 linearly for z<0:
c (0.001 = 1 pc is to avoid divide by zero)
	if(z .le. 0. .and. z .ge. zlhb-clhb) then
	  aa = 0.001 + (alhb-0.001)*(1. - (1./(zlhb-clhb))*z)
	else
	  aa = alhb
	endif
c	write(99, *) x, y, z, aa, bb, cc
	qxy =  ( (x-xlhb)/aa )**2 + ( (y-yaxis)/bb )**2
	qz =  abs(z-zlhb)/cc
	if(qxy .le. 1.0 .and. qz .le. 1.0) then ! inside
	    neLHB2 = netrough
	    FLHBr = Ftrough
	    wLHB = 1
	endif

	return
	end

	real function neLOOPI(x,y,z,FLOOPI,wLOOPI)	! Loop I
c component is a spheroid truncated for z<0.
	implicit none
	real x,y,z,FLOOPI
	integer wLOOPI
c input:
c 	x,y,z = coordinates w.r.t. Galaxy as in TC93, CL00
c output:
c	neLOOPI = electron density in LOOP I that
c	        is modeled as an ellisoidal trough
c		with an enhanced shell
c	FLOOPI = fluctuation parameter
c	wLOOPI  = weight of LOOP I component used to combine
c		with other components of electron density.
c		wLOOPI =  1  at and inside the annular ridge
c		       <  1  outside the annular ridge

	real
     .    aldr,bldr,cldr,xldr,yldr,zldr,thetaldr,neldr0,Fldr,
     .    alsb,blsb,clsb,xlsb,ylsb,zlsb,thetalsb,nelsb0,Flsb,
     .    alhb,blhb,clhb,xlhb,ylhb,zlhb,thetalhb,nelhb0,Flhb,
     .    xlpI,ylpI,zlpI,rlpI,drlpI,nelpI,FlpI,dnelpI,dFlpI
        common/nelismparms/
     .    aldr,bldr,cldr,xldr,yldr,zldr,thetaldr,neldr0,Fldr,
     .    alsb,blsb,clsb,xlsb,ylsb,zlsb,thetalsb,nelsb0,Flsb,
     .    alhb,blhb,clhb,xlhb,ylhb,zlhb,thetalhb,nelhb0,Flhb,
     .    xlpI,ylpI,zlpI,rlpI,drlpI,nelpI,FlpI,dnelpI,dFlpI


	real r
	real a1, a2
	logical first
	data first /.true./
	save

	if(first) then
	  a1 = rlpI
	  a2 = rlpI+drlpI
	  first = .false.
	endif

        if(z .lt. 0.) then
	  neLOOPI = 0.
	  FLOOPI = 0.
	  wLOOPI = 0
          return
        endif
	r = sqrt( (x-xlpI)**2 + (y-ylpI)**2 + (z-zlpI)**2)
        if(r .gt. a2) then 	! outside Loop I
	  neLOOPI = 0.
	  FLOOPI = 0.
	  wLOOPI = 0
	else if(r .le. a1) then	! inside volume
	    neLOOPI= nelpI
	    FLOOPI = FlpI
	    wLOOPI = 1
c           write(99,*) x,y,z, r, neLOOPI, ' inside volume'
	else			! inside boundary shell
	    neLOOPI= dnelpI
	    FLOOPI = dFlpI
	    wLOOPI = 1
c           write(99,*) x,y,z,r, neLOOPI, ' inside shell'
	endif

	return
	end

	subroutine neclumpN(x,y,z,necN,FcN,hitclump, ip, lip)
c returns electron density necN and fluctuation parameter FcN
c at position designated by l,b,d,x,y,z c for a set of
c clumps with parameters read in from file  neclumpN.dat

c input:
c	x,y,z	coordinates	(kpc)  (as in TC93)
c
c output:
c	necN	electron density in clump at (x,y,z)
c	FcN	fluctuation parameter
c 	hitclump = 0:   no clump hit
c		   j>0: j-th clump hit

	implicit none
	real x,y,z,necN,FcN
	integer hitclump

	integer nclumpsmax
	parameter (nclumpsmax=2000)

c	character*15 losname(nclumpsmax)
c	character*1 type(nclumpsmax)
	real lc(nclumpsmax), bc(nclumpsmax), dc(nclumpsmax)
        real xc(nclumpsmax), yc(nclumpsmax), zc(nclumpsmax)
        real nec(nclumpsmax), rc(nclumpsmax), Fc(nclumpsmax)
        integer edge(nclumpsmax)
        character ip*(*)
        integer lip
	integer nclumps
	integer hitclumpflag
	common /clumps/ nclumps, hitclumpflag(nclumpsmax)

c parameters:
c	lc	= galactic longitude of clump center
c	bc	= galactic latitude of clump center
c	(xc,yc,zc) = clump center location (calculated)
c       nec	= internal peak electron density
c	rc	= clump radius at 1/e
c       Fc      = clump fluctuation parameter
c	edge    = 0 => use exponential rolloff out to 5rc
c                 1 => uniform and truncated at 1/e


	real radian
	parameter(radian = 57.29577951)

	real rsun
	parameter (rsun=8.5)

	logical first
	data first/.true./

	real slc, clc, sbc, cbc
	real rgalc

	character fullname*200

	real arg

	integer luclump
	data luclump/11/
	integer j
	integer clumpflag

	save

c first time through, read input clump parameters and calculate
c LOS quantities.
c lc,bc = Galactic coordinates (deg)
c   nec = clump electron density (cm^{-3})
c    Fc = fluctuation parameter
c    dc = clump distance from Earth (kpc)
c    rc = clump radius (kpc)
c  edge = 0,1  0=> Gaussian, 1=> Gaussian w/ hard edge at e^{-1}
c  type = LOS type (P pulsar, G other Galactic, X extragalactic
c losname = useful name
	if(first) then 		!read clump parameters
	  j=1
c	  write(6,*) 'reading neclumpN.NE2001.dat'

        call mkfile('neclumpN.NE2001.dat',fullname,ip,lip)
        open(luclump,file=fullname,status='old')
c	  open(luclump, file='neclumpN.NE2001.dat', status='old')
	  read(luclump,*)				! label line
    5     read(luclump,*,end=99) clumpflag,lc(j),bc(j),nec(j),Fc(j),
     .           dc(j),rc(j),edge(j)
          if(clumpflag .eq. 0) then
	    slc = sin(lc(j)/radian)
	    clc = cos(lc(j)/radian)
	    sbc = sin(bc(j)/radian)
	    cbc = cos(bc(j)/radian)
	    rgalc = dc(j)*cbc
	    xc(j) = rgalc*slc
	    yc(j) = rsun-rgalc*clc
	    zc(j) = dc(j)*sbc
c	  write(6,"(a15,1x,8(f8.3,1x))")
c    .           losname(j),lc(j),bc(j),dc(j),
c    .           nec(j),Fc(j),xc(j),yc(j),zc(j)
	    j=j+1
	  endif
	  go to 5
   99     continue
	  first = .false.
	  nclumps = j-1
	  close(luclump)
	endif

	necN = 0.
	hitclump = 0
	FcN = 0.
	do j=1,nclumps
	  arg =
     .      ((x-xc(j))**2 + (y-yc(j))**2 + (z-zc(j))**2) / rc(j)**2
	  if(edge(j) .eq. 0 .and. arg .lt. 5.) then
	    necN = necN + nec(j) * exp(-arg)
            FcN = Fc(j)
	    hitclump = j
            hitclumpflag(j) = 1
	  endif
	  if(edge(j) .eq. 1 .and. arg .le. 1.) then
c    	    necN = necN + nec(j) * exp(-arg)
	    necN = necN + nec(j)
            FcN = Fc(j)
	    hitclump = j
            hitclumpflag(j) = 1
	  endif
	enddo


	return
	end




	subroutine nevoidN(x,y,z,nevN,FvN,hitvoid,wvoid,ip,lip)
c returns electron density nevN and fluctuation parameter FvN
c at position designated by l,b,d,x,y,z c for a set of
c voids with parameters read in from file  nevoidN.dat

c input:
c	x,y,z	coordinates	(kpc)  (as in TC93)
c
c output:
c	nevN	electron density in void at (x,y,z)
c	FvN	fluctuation parameter
c 	hitvoid =   0:   no void hit
c		  j>0:   j-th void hit
c	wvoid = 0,1:	 void weight

	implicit none
	real x,y,z,nevN,FvN
	integer hitvoid, wvoid

	integer nvoidsmax
	parameter (nvoidsmax=2000)

      character ip*(*)
      integer lip

c	character*12 losname(nvoidsmax)
	real lv(nvoidsmax), bv(nvoidsmax), dv(nvoidsmax)
        real nev(nvoidsmax), Fv(nvoidsmax)
        real aav(nvoidsmax), bbv(nvoidsmax), ccv(nvoidsmax)
        real thvy(nvoidsmax), thvz(nvoidsmax)

        real xv(nvoidsmax), yv(nvoidsmax), zv(nvoidsmax)
	real c1(nvoidsmax), s1(nvoidsmax), c2(nvoidsmax), s2(nvoidsmax),
     .       cc12(nvoidsmax), ss12(nvoidsmax),
     .       cs21(nvoidsmax), cs12(nvoidsmax)
        integer edge(nvoidsmax)
	integer nvoids
        integer hitvoidflag
	common /voids/ nvoids, hitvoidflag(nvoidsmax)

c parameters:
c	lv	= galactic longitude of void center
c	bv	= galactic latitude of void center
c	dv	= distance from Sun of void center
c	(xv,yv,zv) = void center location (calculated)
c       nev	= internal peak electron density
c       Fv      = void fluctuation parameter
c	aav	= void major axis at 1/e
c	bbv	= void minor axis at 1/e
c	ccv	= void minor axis at 1/e
c	thvy	= rotation axis of void about y axis
c	thvz	= rotation axis of void about z axis
c	edge    = 0 => use exponential rolloff out to 5rc
c                 1 => uniform and truncated at 1/e


	real radian
	parameter(radian = 57.29577951)

	real rsun
	parameter (rsun=8.5)

	character fullname*200

	logical first
	data first/.true./

	real slc, clc, sbc, cbc
	real rgalc
	real q
	real dx, dy, dz
	real th1, th2

	integer luvoid
	data luvoid/11/
	integer j
	integer voidflag

	save

c first time through, calculate xc, yc, zc

	if(first) then 		!read void parameters
	  j=1
c	  write(6,*) 'reading nevoidN.dat.clean'
        call mkfile('nevoidN.NE2001.dat',fullname,ip,lip)
        open(luvoid,file=fullname,status='old')
c	  open(luvoid, file='nevoidN.NE2001.dat', status='old')
	  read(luvoid,*)				! label line
    5     read(luvoid,*,end=99) voidflag,
     .      lv(j),bv(j),dv(j),				! deg, deg, kpc
     .      nev(j),Fv(j),				! cm^{-3}, dimensionless
     .      aav(j),bbv(j),ccv(j),thvy(j), thvz(j),  	! kpc,kpc,kpc,deg,deg
     .      edge(j)					! 0 or 1
          if(voidflag .eq. 0) then
	    slc = sin(lv(j)/radian)
	    clc = cos(lv(j)/radian)
	    sbc = sin(bv(j)/radian)
	    cbc = cos(bv(j)/radian)
	    rgalc = dv(j)*cbc
	    xv(j) = rgalc*slc
	    yv(j) = rsun-rgalc*clc
	    zv(j) = dv(j)*sbc
            th1=thvy(j)
            th2=thvz(j)
            s1(j) = sin(th1/radian)
            c1(j) = cos(th1/radian)
            s2(j) = sin(th2/radian)
            c2(j) = cos(th2/radian)
            cc12(j) = c1(j)*c2(j)
            ss12(j) = s1(j)*s2(j)
            cs21(j) = c2(j)*s1(j)
            cs12(j) = c1(j)*s2(j)

c	  write(6,"(a12,1x,13(f7.3,1x))")
c    .           losname(j),lv(j),bv(j),dv(j),
c    .           nev(j),Fv(j),xv(j),yv(j),zv(j),
c    .           aav(j), bbv(j), ccv(j),
c    .           th1, th2
	    j=j+1
	  endif
	  go to 5
   99     continue
	  first = .false.
	  nvoids = j-1
	  close(luvoid)
	endif


	nevN = 0.
	FvN = 0.
	hitvoid = 0
	wvoid = 0
c note rotation matrix in the 'q = ' statement below
c corresponds to \Lambda_z\Lambda_y
c where \Lambda_y = rotation around y axis
c       \Lambda_z = rotation around z axis
c defined as
c \Lambda_y =  c1  0  s1
c               0  1   0
c             -s1  0  c1

c \Lambda_z =  c2 s2   0
c             -s2 c2   0
c               0  0   1
c =>
c \Lambda_z\Lambda_y =  c1*c2   s2   s1*c2
c                      -s2*c1   c2  -s1*s2
c                         -s1    0      c1

c so the rotation is around the y axis first, then the z axis
	do j=1,nvoids
	  dx = x-xv(j)
	  dy = y-yv(j)
	  dz = z-zv(j)
	  q = ( cc12(j)*dx + s2(j)*dy + cs21(j)*dz)**2 / aav(j)**2
     .      + (-cs12(j)*dx + c2(j)*dy - ss12(j)*dz)**2 / bbv(j)**2
     .      + (  -s1(j)*dx         +      c1(j)*dz)**2 / ccv(j)**2
	  if(edge(j) .eq. 0 .and. q .lt. 3.) then
	    nevN = nev(j) * exp(-q)
            FvN = Fv(j)
	    hitvoid = j
            hitvoidflag(j)=1
	  endif
	  if(edge(j) .eq. 1 .and. q .le. 1.) then
	    nevN = nev(j)
            FvN = Fv(j)
	    hitvoid = j
            hitvoidflag(j)=1
	  endif
	enddo

	if(hitvoid .ne. 0) wvoid = 1


	return
	end




c this version from 18 March 1998 uses revised coefficients
c that are consistent with Cordes \& Rickett (1998, ApJ, submitted)
c modifications:
c	28 March 2001: added FUNCTION TRANSITION_FREQUENCY

      REAL FUNCTION TAUISS(d, sm, nu)
c
c calculates the pulse broadening time in ms
c from distance, scattering measure, and radio frequency
c
c input:      d = pulsar distance       (kpc)
c            sm = scattering measure    (kpc m^{-20/3})
c            nu = radio frequency       (GHz)
c output: tauss = pulse broadening time (ms)
c
      implicit none
      real d, sm,  nu
      tauiss = 1000. * (sm / 292.)**1.2 * d * nu**(-4.4)
      end
c
c
      REAL FUNCTION SCINTBW(d, sm, nu)
c
c calculates the scintillation bandwidth in kHz
c from distance, scattering measure, and radio frequency
c
c input:        d = pulsar distance       (kpc)
c              sm = scattering measure    (kpc m^{-20/3})
c              nu = radio frequency       (GHz)
c output: scintbw = scintillation bandwidth (kHz)
c
      implicit none
      real d, sm, nu
      real c1
      parameter(c1=1.16)		! for uniform, Kolmogorov medium
      real tauiss
      tauiss = 1000. * (sm / 292.)**1.2 * d * nu**(-4.4)	! ms
      scintbw = c1 / (2. * 3.14159 * tauiss)			! kHz
      end

      REAL FUNCTION SCINTIME(sm, nu, vperp)
c
c calculates the scintillation speed for given distance, galactic
c longitude and latitude, frequency, and transverse velocity
c
c input:   sm = scattering measure	(kpc m^{-20/3})
c          nu = radio frequency 	(GHz)
c       vperp = psr transverse speed  	(km/s)
c
c output: scintime = scintillation time (sec)
c
c usage: should be called with sm = smtau for appropriate
c        line of sight weighting
c reference: eqn (46) of Cordes & Lazio 1991, ApJ, 376, 123.
c
      implicit none
      real sm, nu, vperp
c nb: formerly, the coeff. in the following line was 2.3 from
c     Cordes & Lazio (1991)
      scintime = 3.3 * nu**1.2 * sm**(-0.6) * (100./vperp)
      end


      REAL FUNCTION SPECBROAD(sm, nu, vperp)
c
c calculates the bandwdith of spectral broadening
c for given scattering measure, , frequency, and transverse velocity

c input:   sm = scattering measure	(kpc m^{-20/3})
c          nu = radio frequency 	(GHz)
c       vperp = psr transverse speed  	(km/s)

c output: specbroad = spectral broadening bandwidth (Hz)
c
c usage: should be called with sm = smtau for appropriate
c        line of sight weighting
c reference: eqn (47) of Cordes & Lazio 1991, ApJ, 376, 123.
c
      implicit none
      real sm, nu, vperp
c nb: the coeff. in the following line is 0.14 Hz  from Cordes & Lazio (1991)
c it is changed to 0.097 to conform with FUNCTION SCINTIME and
c a new calculation consistent with Cordes & Rickett (1998)

      specbroad = 0.097 * nu**(-1.2) * sm**0.6 * (vperp/100.)	! Hz
      end


      REAL FUNCTION THETA_XGAL(sm, nu)
c
c calculates angular broadening for an extragalactic
c source of plane waves
c
c sm = scattering measure
c nu = radio frequency
c theta_xgal = angular broadening FWHM (mas)
c
      implicit none
      real sm, nu
      theta_xgal = 128. * sm**0.6 * nu**(-2.2)
      end
c
      REAL FUNCTION THETA_GAL(sm, nu)
c
c calculates angular broadening for a galactic
c source of spherical waves
c
c sm = scattering measure
c nu = radio frequency
c theta_gal = angular broadening FWHM (mas)
c
      implicit none
      real sm, nu
      theta_gal = 71. * sm**0.6 * nu**(-2.2)
      end
c
      FUNCTION EM (sm)
c
c units of sm are kpc m^{-20/3}
c units of em are pc cm^{-6}
c
c calculates the emission measure from the scattering measure
c using an assumed outer scale and spectral index of the
c wavenumber spectrum.
c
c for a wavenumber spectrum P_n(q) = q^{-alpha} from q_0 to q_1
c the mean square electron density is
c
c <n_e^2> =~  4pi*[C_n^2 / (alpha - 3) ] * q_0^{3 - alpha)
c
c ( an approximate form that assumes (q_0 / q_1)^{3-alpha} >> 1.
c
c Jim Cordes 18 Dec 1989
c
      data router /1./	! outer scale = 1 pc
      data pc/3.086e+18/
      data alpha/3.6666667/
      data pi/3.14159/
c
      em = sm *
     1         ( (4. * pi * 1000.) / (alpha - 3.) ) *
     2         (router*pc / (2. * 3.14159) )**(alpha-3.) *
     3         (0.01) ** (20./3.)
c
      return
      end

      REAL FUNCTION THETA_ISO(smiso, nu)
      real smiso, nu

c smiso in (kpc m^{-20/3}) x kpc^{5/3}
c    nu in GHz
c returns the isoplanatic angle in microarcsec
c 12 October 1998
c JMC


c \theta_{iso} = \delta r_s / d
c              = \left [
c	         (\lambda r_e)^2 f_{\alpha} SM_{iso}
c		 \right ]^{1/\alpha}
c where \alpha = 5/3 for Kolmogorov case.
c NB SM_{iso} = \int_0^d ds s^{\alpha} \cnsq
c    so SM_{iso} does not have the units of scattering
c    measure, but rather units of SM x Length^{\alpha}
c
c f_{\alpha} = 8\pi^2 \Gamma(1-\alpha/2) / [\alpha 2^{\alpha} \Gamma(1+\alpha/2)]
c for \alpha = 5/3, f_{\alpha}= 88.3
c
c     real r_e
c     parameter(r_e = 2.82e-13)			!cm
c     real kpc
c     parameter(kpc = 3.086e21)			!cm
c     real falpha
c     parameter(falpha=88.3)

      theta_log_radian =
     .    13.287  				! 0.6*log10(30cm*r_e)
     .  + 1.2 * alog10(nu)
     .  - 1.1676				! 0.6*log10(f_alpha)
     .  - 0.6 * alog10(smiso)
     .  - 34.383				! 1.6 * alog10(kpc)
     .  + 8.					! -(20/3)*log(100)
      theta_log_microarcsec =
     .      theta_log_radian + 11.314425	! 11.314425=alog10(microarsec/rad)
      theta_iso = 10.**theta_log_microarcsec
      return
      end



      REAL FUNCTION THETA_ISO_TEST(smiso, nu)
      real smiso, nu

c smiso in (kpc m^{-20/3}) x kpc^{5/3}
c    nu in GHz
c returns the isoplanatic angle in microarcsec
c 12 October 1998
c JMC


c \theta_{iso} = \delta r_s / d
c              = \left [
c	         (\lambda r_e)^2 f_{\alpha} SM_{iso}
c		 \right ]^{1/\alpha}
c where \alpha = 5/3 for Kolmogorov case.
c NB SM_{iso} = \int_0^d ds s^{\alpha} \cnsq
c    so SM_{iso} does not have the units of scattering
c    measure, but rather units of SM x Length^{\alpha}
c
c f_{\alpha} = 8\pi^2 \Gamma(1-\alpha/2) / [\alpha 2^{\alpha} \Gamma(1+\alpha/2)]
c for \alpha = 5/3, f_{\alpha}= 88.3
c
c     real r_e
c     parameter(r_e = 2.82e-13)			!cm
c     real kpc
c     parameter(kpc = 3.086e21)			!cm
c     real falpha
c     parameter(falpha=88.3)

      theta_log_radian =
     .    13.287  				! 0.6*log10(30cm*r_e)
     .  + 1.2 * alog10(nu)
     .  - 1.1676				! 0.6*log10(f_alpha)
     .  - 0.6 * alog10(smiso)
     .  - 34.383				! 1.6 * alog10(kpc)
     .  + 8.					! -(20/3)*log(100)
      theta_log_microarcsec =
     .      theta_log_radian + 11.314425	! 11.314425=alog10(microarsec/rad)
      theta_iso_test = 10.**theta_log_microarcsec
c     write(6,*) 'smiso, nu = ', smiso, nu
c     write(6,*) 'theta_log_radian = ', theta_log_radian
c     write(6,*) 'theta_log_microarcsec = ', theta_log_microarcsec
c     write(6,*) 'theta_iso = ', theta_iso_test
      return
      end

      REAL FUNCTION TRANSITION_FREQUENCY(sm, smtau, smtheta, dintegrate)
      implicit none
      real sm, smtau, smtheta, dintegrate
c returns the transition frequency between weak and strong scattering
c 28 March 2001
c JMC

c input:
c (all sm values in (kpc m^{-20/3}))
c	    sm = int[\cnsq]
c        smtau = int[(s/D)(1-s/D) \cnsq]
c      smtheta = int[ (1-s/D) \cnsq]
c   dintegrate = distance used to integrate \cnsq (kpc)
c output:
c   transition_frequency = GHz given by
c 	\nu_t  = 318 GHz \xi^{10/17}  SM^{6/17} D_{eff}^{5/17}
c   where
c        D_{eff} = effective path length through medium
c        D_{eff} = \int_0^dintegrate ds s \cnsq / \int_0^dintegrate ds  \cnsq
c
c   Note we can calculate D_{eff} using
c        D_{eff} = dintegrate * (sm - smtau/6 - smtau/3) / sm
c
      real deff
      real xi
      parameter(xi= 0.3989)   ! (2.*pi)^{-1/2} = fresnel scale definition factor
      real coefficient
      parameter(coefficient=318.)	     ! GHz; see NE2001 paper
      deff = (dintegrate*(sm - smtau/6. - smtheta/3.)) / sm
      transition_frequency =
     .     coefficient * xi**(10./17.) * sm**(6./17.) * deff**(5./17.)
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine mkfile(filename,fullname,ip,lip)

      implicit none

      character ip*(*)
      integer lip
      character path*200, filename*(*), fullname*(*)
      integer length,lpth
c      call getpath(path,lpth)
c      path=path(1:lpth)//'/lookuptables'
      path = ip(1:lip)//'/lookuptables'
      call makefilename(path,filename,fullname,length)

      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine makeFilename(path,filename,fullname,length)

      implicit none

      character path*(*), filename*(*),fullname*(*)

      integer length

      call stripBlanks(path,fullname,length)

ccccc THE FOLLOWING LINE HAS TO BE CHANGED IN A NON UNIX-VERSION !!!

      fullname=fullname(1:length)//'/'//filename

      call stripBlanks(fullname,fullname,length)

      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine stripBlanks(instring, outstring, length)

      implicit none

      character instring*(*), outstring*(*), blank*1

      integer length, icount

      data blank/' '/

      icount=1
 100  continue
         if (instring(icount:icount).ne.blank) goto 200
         if (icount.eq.len(instring)) then
            length=1
            outstring=instring
            return
	 else
            icount=icount+1
            goto 100
	 endif
 200  continue

      outstring=instring(icount:)

      length=len(outstring)
 300  continue
         if (outstring(length:length).ne.blank) goto 400
         if (length.gt.1) then
            length=length-1
            goto 300
	 endif
 400  continue

      end
