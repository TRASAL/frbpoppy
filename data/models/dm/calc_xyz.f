c==============================================================================
      subroutine calc_xyz(l, b, d, x, y, z)
c==============================================================================
c
c     calculates the cartesian xyz positions in the Galaxy
c     given galactic latitude and longitude and distance.
c                       DRL 30/1/92
c
c     Last Change 93/05/20 DRL @ JB Changed to sun at (-rsun,0,0)
c
      implicit none
c
c     passed down variables: in l & b in degrees d in kpc
c                           out x,y,z in kpc
c
      real l, b, d, x, y, z
c
c     local variables
c
      real cl, cb, sl, sb, lold, bold, rsun, rad
      parameter(rsun=8.5,rad=57.29577951)
      logical first
      data first/.true./
c
c     work out cos & sin if changed from last call
c
      if (l.ne.lold.or.b.ne.bold.or.first) then
        cl=cos(l/rad)
        sl=sin(l/rad)
        cb=cos(b/rad)
        sb=sin(b/rad)
        first=.false.
        lold=l
        bold=b
      endif
c
c     calculate (x,y,z)
c
c      x = d*cb*cl - rsun
c      y = d*cb*sl
c      z = d*sb

	x=d*cb*sl
	y=rsun-d*cb*cl
	z=d*sb

      end
c==============================================================================
