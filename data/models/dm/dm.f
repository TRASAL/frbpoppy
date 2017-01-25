c==============================================================================
      real function dm(dkpc,l,b,dmod,sm,ip,lip)
c==============================================================================
c
c     Routine to integrate numerically from the sun to a pulsar given its
c     distance (kpc), galactic longitude l and latitude b (degrees). Uses
c     a model for the free electron density distribution in the Galaxy
c     chosen by the integer mod. The DM is returned in usual units (cm-3 pc)
c
      implicit none
c
c     passed down variables
c
      real dkpc, l, b
      integer dmod,lout
c
c     local variables...
c
      real dstep, dmint, maxd, locne, x, y, z, dist, sm
c
c     functions
c
      real psr_ne
c
c     Added ip and lip (inpath and length of inpath) to remove getpath
c     dependence!
c
      character ip*(*)
      integer lip

c
c     define parameters & initial values
c
      parameter(dstep=0.05)        ! integrating step in kpc
      parameter(maxd=30.0)         ! max distance to go for kpc
      character*1 limit
      real s1,s2,s3,s4,dr
      parameter(dr=0.0174532925199)
      common /outputtxt/ lout
      logical first
      data first/.true./
      save
      if (dmod.eq.4) then
c          if (first) write(lout,*)
	      first = .false.
         call dmdsm(l*dr,b*dr,-1,dm,dkpc,limit,s1,s2,s3,s4,ip,lip)
         sm = s2 ! return scattering measure
         return
      endif
      sm = 0.0
      dmint = 0.0                  ! integrated dm
c
c     Main loop, integrate out to given dm
c
      do dist = dstep, maxd, dstep
c
c        calculate local position
c
         call calc_xyz(l,b,dist,x,y,z)
c
c        local electron density
c
         locne = psr_ne(x, y, z, dmod, ip, lip)
c
c        integrated dm
c
         dmint = dmint + locne * dstep * 1000.0
c
c        integrated far enough?
c
         if (dist.gt.dkpc) goto 999
      end do
 999  dm=dmint
      end
