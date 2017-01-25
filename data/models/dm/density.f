      subroutine density(x,y,z,ne1,ne2,nea,negum,ip,lip)

c  Returns four components of the free electron density of the 
c  interstellar medium at Galactic location (x,y,z).  Add them together 
c  to get the total density.  Combine terms with `fluctuation'
c  parameters to get the scattering measure.  
c
c  Origin of coordinate system is at Galactic center; 
c  x,y,z are in kpc; 
c  the sun os at (x,y,z) = (0,R0,0), x is in l=90 direction;
c  electron densities in cm^-3.

C	ne1:	outer zone
C	ne2:	galactic center zone
C	nea:	spiral arms
C	negum:	Gum Nebula

C  Spiral arm locations are contained in array arm(j,k,i), where
C    j=1,4 distinguishes the four independent arms;
C    k=1,kmax(j) selects numerical spiral arm coordinates;
C    i=1,2 for coordinate=x,y
C  The arm locations were digitized and tabulated by hand; loop "do 20"
C  re-creates close approximations to curves fitting the arm axes.

      real n1h1,n2,na,ngum,ne1,ne2,nea,negum
      real th1(7,4),r1(7,4)
      logical first
      save
      common/params/n1h1,h1,A1,F1,n2,h2,A2,F2,na,ha,wa,Aa,Fa,Fg
      common/armcom/arm(4,500,2),kmax(4)
      character*120 path
      integer ldir,lun
      character ip*(*)
      integer lip
      data first/.true./,ks/3/,NN/7/
      data rad/57.2957795/,R0/8.5/
      data th1/164.,200.,240.,280.,290.,315.,330.,
     +            63.,120.,160.,200.,220.,250.,288.,
     +            52.,120.,170.,180.,200.,220.,252.,
     +            20.,70.,100.,160.,180.,200.,223./
	data r1/3.1,3.3,3.9,4.6,4.7,5.1,5.1,
     +          3.3,4.0,4.2,5.0,5.7,6.4,7.2,
     +          4.3,5.5,5.7,6.1,7.2,7.8,8.4,
     +          5.2,6.2,6.9,8.5,9.1,10.0,10.6/

	data xgum/-0.492/,ygum/8.587/		! Location of Gum Nebula
	data rgum/0.182/			! Radius of Gum nebula
	data ngum/0.20/				! n_e inside Gum Nebula

        Fg=0.                           ! Fluctuation parameter for Gum
	if(first) then			! Reconstruct spiral arm axes
c	call getpath(path,ldir)
        call glun(lun)
	do 20 j=1,4
	dth=5.0/r1(1,j)
	th=th1(1,j)-0.999*dth
	call cspline(th1(1,j),r1(1,j),-NN,th,r)
	do 10 k=1,499
	th=th+dth
	if(th.gt.th1(7,j)) go to 20
	call cspline(th1(1,j),r1(1,j),NN,th,r)
	arm(j,k,1)=(R0/7.46)*r*sin(th/rad)
	arm(j,k,2)= -(R0/7.46)*r*cos(th/rad)
10	continue
20	kmax(j)=k
        open(lun,file=ip(1:lip)//'/lookuptables/gal.dat')
	read(lun,1020) n1h1,h1,A1,F1,n2,h2,A2,F2,na,ha,wa,Aa,Fa
1020	format(6x,f8.0)
	close(lun)
	first=.false.
	endif

	rr=sqrt(x**2 + y**2)		! Galactocentric radius
 	g1=sech2(rr/A1)/sech2(8.5/A1)
 	ne1=(n1h1/h1)*g1*sech2(z/h1)	! Outer component

	g2=0.0
	rrarg=((rr-A2)/1.8)**2
	if(rrarg.lt.10.0) g2=exp(-rrarg)
	ne2=n2*g2*sech2(z/h2)		! Galactic center component

	nea=0.
	if(abs(z/ha).lt.3) then		! Get spiral arm component
	  do 50 j=1,4
	    sqmin=1.e10
	    do 30 k=1+ks,kmax(j)-ks,2*ks+1 ! Find min dist to this arm
	      sq=(x-arm(j,k,1))**2 + (y-arm(j,k,2))**2
	      if(sq.lt.sqmin) then
	        sqmin=sq
		kk=k
	      endif
30	    continue
	    do 40 k=kk-ks,kk+ks
	      sq=(x-arm(j,k,1))**2 + (y-arm(j,k,2))**2
	      if(sq.lt.sqmin) then
	        sqmin=sq
		kk=k
	      endif
40	    continue
	    smin=sqrt(sqmin)		! Distance of (x,y,z) from arm axis
	    if(smin.lt.3*wa) then	! If (x,y,z) is close to this
	      ga=exp(-(smin/wa)**2)	! arm, get the arm weighting factor 
	      if(rr.gt.Aa) ga=ga*sech2((rr-Aa)/2.0) ! Radial dependence of arms
	      if(j.eq.2.and.k.ge.101) then
		fac=2.0
		if(k.le.115) fac=1.0 + (k-101)/14.0
		ga=ga*fac
	      endif
	      if(j.eq.3.and.k.ge.60.and.k.le.94) then
	        th=6.2831853*(k-60.0)/34.0
		fac=(3.0+cos(th))/4.0
		ga=ga*fac
	      endif
	      nea=nea + na*ga*sech2(z/ha) ! Add this arm contribution
	    endif
50	  continue
	endif

	negum=0.0			! Gum Nebula component
	r=sqrt((x-xgum)**2 + (y-ygum)**2 + z**2)
	if(r.lt.2*rgum) then
	  negum=ngum
	  ag=0.7
	  r2=ag*rgum
	  if(r.gt.r2) negum=ngum*exp(-((r-r2)/(rgum-r2))**2)
	endif

	return
	end
