c	Test matrel subroutines
	parameter (ncmax=151)
	parameter (ncmaz=74)
c	Next line is the number of interior GLL points in 1D
	parameter (lmax=4)
c	N is the total number of GLL points in 1D including the points +-1
	parameter (N=lmax+2)
c	Next line is the maximum number of GLL points on the x-side of
c	the global grid times the maximum number on the z-side.
	parameter (nptmax=(N*ncmax-(ncmax-1))*(N*ncmaz-(ncmaz-1)))
	parameter (N2=N*N)
	parameter (kmax=ncmax*ncmaz*N2*N2*9)
	character*80 recfil,disfil,progfil
c***
	real*8 phi,dphi,ddphi
	real*8 phi2,dphi2x,dphi2z,dphixx,dphixz,dphizz,wt2
	parameter (lmax2=(lmax+2)*(lmax+2))
	common/basisfns/phi(lmax+2,lmax+2),dphi(lmax+2,lmax+2),ddphi(lmax+2,lmax+2),
     &	phi2(lmax2,lmax2),dphi2x(lmax2,lmax2),dphi2z(lmax2,lmax2),
     &	dphixx(lmax2,lmax2),dphixz(lmax2,lmax2),dphizz(lmax2,lmax2),
     &	wt2(lmax2)
c	phi(i,j) has phi_i (x_j)
c	dphi(i,j) has [d(phi_i)/dx]|x=x_j
	real*8 wt,x
	common/glpts/wt(lmax+2),x(lmax+2)
c***
	common/grdpar/NX,NZ
	common/struc3/dx(ncmax),dz(ncmaz)
	real*8 vp,vs,rho,qb,qk
	common/struc4/vp(ncmax,ncmaz,N*N),vs(ncmax,ncmaz,N*N),
     &	rho(ncmax,ncmaz,N*N),qb(ncmax,ncmaz,N*N),qk(ncmax,ncmaz,N*N)
	common/topo/dztopo(ncmaz,ncmax+1)
c	xg and zg have the (x,z) coordinates at each of the global gridpoints
c	assuming dimensions of 2x2 for each cell (x and z each run from -1 to +1).
	common/glbxy/xg(nptmax),zg(nptmax)
	common/glbgrd/igrd(ncmax,ncmaz,N*N)
	real*8 xgdble,phdble,cdelt,sdelt,sfhi,cfhi,delta,sphi,cphi,phirec,fhi,cpsi,spsi,psi
	real*8 dlat,plat,plon,phiref
c	real*8 pi,rad
	common/sphpar/dlat,plat,plon,phiref
	common/minfo/minc,mmax
cseis
	common/sinfo/nleno,dt,corper
c--
c***
cseis
	parameter (maxom=681)
	  parameter (maxphi=261)
	complex*16 d,dsrec,dsrecp,ds,dsverp,dsvrpp,dsvrp2,dsvert,dsvrtp,dsvrts
	common/datavec/d(3*nptmax,4)
	parameter (maxrec=300)
	real*8 slat,slon
	real*8 fac1
	real*8 reclat(maxrec),reclon(maxrec),recz(maxrec)
c	dimension ncell(maxrec)
	dimension ncxs(maxrec),nczs(maxrec)
	dimension xrecd(maxrec),xrec(maxrec),zrec(maxrec)
	dimension phir(maxrec),deltar(maxrec)
	dimension dsverp(3*nptmax),dsvrpp(3*nptmax),dsvrp2(3*nptmax)
	common/vertpr/dsvert(3*nptmax,maxom),dsvrtp(3*nptmax,maxom),dsvrts(3*nptmax,maxom)
	common/msumr/dsrec(3*N*N*maxrec,maxom),dsrecp(3*N*N*maxrec,maxom)
	common/msum/ds(maxphi,3*ncmax*N,maxom)
c--
c***
cUSED        parameter (nzmax = 32800000, nmax = 215000)
        parameter (nzmax = kmax, nmax = 3*(N*ncmax-(ncmax-1))*(N*ncmaz-(ncmaz-1))+1)
cUSED	integer*8 Ap,Ai,jAi
cUSED	dimension Ap(nmax),Ai(nzmax),jAi(kmax)
	integer*8 Ap
	integer(8), ALLOCATABLE :: Ai(:)
	integer(8), ALLOCATABLE :: jAi(:)
	complex*16, ALLOCATABLE :: AA0(:) , AA1(:) , AA2(:)
	dimension Ap(nmax)
C***
cseis
	real*8 phiv,dphivx,dphivz
	complex*16 w0,wval,ui
	complex*16 facp,facm
	dimension wval(maxom)
	parameter (iommax=3*maxom/2)
	complex*8 xo
	common/spec/xo(iommax),xt(1024)
	common/invtyp/ityp
	dimension xtx(1024),xty(1024),xtz(1024)
	dimension xtcx(1024),xtcy(1024),xtcz(1024)
	dimension xtSx(1024),xtSy(1024),xtSz(1024)
	dimension xtPx(1024),xtPy(1024),xtPz(1024)
	dimension xtxx(1024),xtxy(1024),xtxz(1024)
	dimension xtyy(1024),xtyz(1024),xtzz(1024)
c--
	real*8 ky
	parameter (bigr=6371.)
c
	pi=3.1415926535
        twopi=2.*pi
	rad=360./twopi
	ui=dcmplx(0.,1.)
	call init
	nmat=(NX+1)*(NZ+1)*N2*N2*9
	allocate ( Ai(nmat) )
	allocate ( jAi(nmat) )
	allocate ( AA0(nmat) )
	allocate ( AA1(nmat) )
	allocate ( AA2(nmat) )
		write(6,*)'nmat=',nmat
		write(6,*)'AA size of Ai=',size(Ai)
		write(6,*)'AA size of jAi=',size(jAi)
		write(6,*)'AA Ai(20)=',Ai(20)
		write(6,*)'AA jAi(20)=',jAi(20)
	write(6,*)'out of init: nleno,dt,corper=',nleno,dt,corper
	write(6,*)'plat,plon=',plat,plon
cUMF
	write(6,*)'Determine matrix ordering? (yes=1)'
	write(6,*)'(no, and do inverse for all outputs FT>=2)'
	write(6,*)'(no, and do inverse for only CURL+DIV and S+P outputs FT=3)'
	write(6,*)'(no, and do not do FT=any other value)'
	write(6,*)'(no, and write out only seis2pt5d.outxyz_vertp, seis2pt5d.outxyz_rec, seis2pt5d.outstrains_rec =4)'
	read(5,*) imatr
	if(imatr.eq.1) then
	write(6,*)'calling sporder'
	call sporder(nmat,Ap,Ai,jAi)
	stop
	else
	write(6,*)'Reading in previously calculated matrix ordering'
	open(2,file='sporder.out',form='unformatted')
	read(2) Ap
	read(2) Ai
	read(2) jAi
	close(2)
	endif
c		do j=1,5
c		write(6,*)'Ap(',j,')=',Ap(j)
c		enddo
c		write(6,*)'Ap(500)=',Ap(500),'Ap(501)=',Ap(501)
c		ncol=3*(N*NX-(NX-1))*(N*NZ-(NZ-1))
c		write(6,*)'Ap(',ncol-1,')=',Ap(ncol-1),'Ap(',ncol,')=',Ap(ncol)
c		write(6,*)'Ai(column 5)=', (Ai(i), i=Ap(5)+1,Ap(5+1))
c		write(6,*)'Ai(column 500)=', (Ai(i), i=Ap(500)+1,Ap(500+1))
c		write(6,*)'Ai(column',ncol-1,')=', (Ai(i), i=Ap(ncol-1)+1,Ap(ncol))
c		pause
c
	write(6,*)'line B: plat,plon=',plat,plon
c---
c	Read in lat,lon,dep of receivers
	open(2,file='receivers-latlondep.txt',status='old',err=25)
	read(2,*) nrec
	do nr=1,nrec
	read(2,*) alatin,alonin,depin
		write(6,*) alatin,alonin,depin
	reclat(nr)=(pi/2.d0-alatin/rad)
	reclon(nr)=alonin/rad
	recz(nr)=depin
		write(6,*)'reclat,reclon(',nr,')=',reclat(nr),reclon(nr),'recz=',recz(nr)
		write(6,*)'------------------'
	if(recz(nr).gt.0.d0) recz(nr)=0.d0
	enddo
	close(2)
	write(6,*)'line C'
c
	icell=0
	  do nr=1,nrec
	slat=reclat(nr)
	slon=reclon(nr)
		write(6,*)'receiver colat,lon=',slat,slon
		write(6,*)'pole colat,lon=plat,plon',plat,plon
	zs=recz(nr)
c*	Determine angular distance (rad) and azimuth (rad CCW from due N)
c	of the receiver from the pole of the spherical coordinate system.
c	the azimuth is referenced to the azimuth of the origin (theta=dlat,z=0)
c	of the 2D grid.
        cdelt=dcos(slat)*dcos(plat)+dsin(slat)*dsin(plat)*dcos(plon-slon)
        delta=dacos(cdelt)
	deltar(nr)=delta
	sdelt=dsin(delta)
        sphi=dsin(slat)*dsin(plon-slon)/sdelt
        cphi=(dcos(slat)-dcos(plat)*cdelt)/(dsin(plat)*sdelt)
	phirec=datan2(sphi,cphi)-phiref
	phir(nr)=phirec
c--
		write(6,*)'RECEIVER: phi_receiver=',rad*phirec
c*	The theta-coord of the receiver is the angular distance from (pclat,plon)
	xs=real(delta)
		write(6,*)'RECEIVER: theta-coordinate x[receiver]=',xs
		write(6,*)'RECEIVER: z-coordinate z[receiver]=',zs
	xrecd(nr)=xs
c	Determine which cell number the receiver is located in
	ncxs(nr)=0
	nczs(nr)=0
	do ncz=1,NZ
	do ncx=1,NX
c	lower left corner
	ilx=1
	ilz=1
	il=N*(ilz-1)+ilx
	ig=igrd(ncx,ncz,il)
	x1=xg(ig)
	z1=zg(ig)
c		write(6,*)'RECEIVER: cell #',ncx,ncz,'dimensionless lower left corner=',x1,z1
c	upper right corner
	ilx=N
	ilz=N
	il=N*(ilz-1)+ilx
	ig=igrd(ncx,ncz,il)
	x2=xg(ig)
	z2=zg(ig)
	tx=(xs-x1)*(xs-x2)
	tz=(zs-z1)*(zs-z2)
	  if(tx.le.0.0.and.tz.le.0.0) then
	  ncxs(nr)=ncx
	  nczs(nr)=ncz
	  xrec(nr)=(2./dx(ncx))*(xs-x1)-1.
	  zrec(nr)=(2./dz(ncz))*(zs-z1)-1.
c		write(6,*)'xs,x1=',xs,x1
c		write(6,*)'zs,z1=',zs,z1
		write(6,*)'xrec,zrec=',xrec(nr),zrec(nr)
c	  ilc=N*((N/2)-1)+(N/2)
c	  igc=igrd(ncx,ncz,ilc)
c		write(6,*)'RECEIVER #',nr,': gridpoint of cell center=',igc
	  endif
	enddo
	enddo
	write(6,*)'RECEIVER #',nr,': Cell #s of receiver =(',ncxs(nr),nczs(nr),')'
	if(ncxs(nr).eq.0.or.nczs(nr).eq.0) then
	write(6,*)'RECEIVER #',nr,': Receiver location lies outside model domain'
	stop
	endif
	write(6,*)'Local dimensionless receiver coords=(',xrec(nr),zrec(nr),')'
cc	See if this is a new cell.  If so, keep track of the first nr-value of the new cell.

c	if(icell.eq.0) then
c	icell=1
c	ncell(icell)=nr
c	else

c	inew=1
c	do ic0=1,icell-1
c	nr0=ncell(ic0)
c	if(ncxs(nr).eq.ncxs(nr0).and.nczs(nr).eq.nczs(nr0)) inew=0
c	enddo
c	if(inew.eq.1) then
c	icell=icell+1
c	ncell(icell)=nr
c	endif

c	endif

	  enddo
c	icellt=icell
c	if(i.ne.9999) stop
c---
cseis
	ommin=twopi*(1./dt)/real(nleno)
	obeta=2.0*ommin
c	Anticipating corner periods in WAVE1 that are >= 6 x (1/dt),
c	do only nleno/3 frequencies (thus anticipated value of jmax1 in
c	INVFT, called below, will not exceed nleno/3).
	  if((nleno/3).gt.maxom) then
	  write(6,*)'nleno/3=',nleno/3,' gt maxom=',maxom
	  stop
	  endif
	do iom=1,nleno/3
	wval(iom)=dble(iom-1)*dble(ommin) - ui*dble(obeta)
	enddo
c--
	write(6,*)'minc,mmax=',minc,mmax
cTE	Reset mmax to do just ky=0 and ky=minc
c	minc=2000
c	mmax=minc
c--
	write(6,*)'Number of wavenumber components=',1+2*(mmax/minc)
	write(6,*)'File for storing displacement spectra at receivers?'
	read(5,10) recfil
	write(6,*)'File for storing displacement spectra at surface?'
	read(5,10) disfil
	write(6,*)'Progress file?'
	read(5,10) progfil
	write(6,*)'Beginning and ending fraction of frequency indices?'
	read(5,*) f1,f2
	iom1=int(f1*real(nleno/3))+1
	iom2=int(f2*real(nleno/3))
	write(6,*)'Doing frequency indices ',iom1,'to ',iom2
cTE
c	Continue with the inverse Fourier transforms only if imatr>=2
c	If so, it is assumed that
c	the displacement spectra are in 'displ-j-m', regardless of
c	where earlier partial iom-loops may have stored these spectra.
cOLD	if(iom1.eq.1.and.iom2.eq.(nleno/3)) go to 20
	if(imatr.ge.2) then
	write(6,*)'Will calculate inverse FTs after one call to SOURCE'
	w0=wval(1)
	ky=0.d0
	call source(w0,ky)
	go to 20
	endif
c--
cseis
	do iom=iom1,iom2
cOLD	do iom=1,nleno/3
	w0=wval(iom)
c	Zero out the dsverp-array
	do ig=1,(N*NX-(NX-1))*(N*NZ-(NZ-1))
	i1=3*ig-2
	i2=i1+1
	i3=i2+1
	dsverp(i1)=0.d0
	dsverp(i2)=0.d0
	dsverp(i3)=0.d0
	dsvrpp(i1)=0.d0
	dsvrpp(i2)=0.d0
	dsvrpp(i3)=0.d0
	dsvrp2(i1)=0.d0
	dsvrp2(i2)=0.d0
	dsvrp2(i3)=0.d0
	enddo
c--
c	Start loop over azimuthal order number
c		write(6,*)'minc,mmax+minc,minc=',minc,mmax+minc,minc
c		if(i.ne.9999) stop
	  do m=minc,mmax+minc,minc
cTE	should be +
c--
	ky=dble(m-minc)
	write(6,*)'MLOOP: m=',m,'out of',mmax+minc,'with increment',minc
	write(6,*)'MLOOP: ky=',ky
c - - -
	open(2,file=progfil,status='old',access='append')
	write(2,*)'Doing iom=',iom,'out of',nleno/3
	write(2,*) 'MLOOP: m=',m,'out of',mmax+minc,'with increment',minc
	close(2)
c - - -
	call source(w0,ky)
	write(6,*)'entering matrel, iom=',iom,'out of',nleno/3
	write(6,*)'w0=',w0,'ky=',ky
cOLD	write(6,*)'KU,KL=',KU,KL
	call matrel(w0,ky,nmat,Ap,Ai,jAi,AA0,AA1,AA2)
	write(6,*)'out of matrel'
c**	Store values of Laplace-transformed and phi-transformed displacement
c	at the cells containing receivers.
	if(iom.eq.iom1.and.m.eq.minc) then
	open(4,file=recfil,form='unformatted')
	else
	open(4,file=recfil,form='unformatted',status='old',access='append')
	endif
c	do icell=1,icellt
c	nr=ncell(icell)
	do nr=1,nrec
	ncz=nczs(nr)
	do ilz=1,N
	ncx=ncxs(nr)
	do ilx=1,N
	il=N*(ilz-1)+ilx
	ig=igrd(ncx,ncz,il)
	i1=3*ig-2
	i2=i1+1
	i3=i2+1
c	x- and z-displacements are symmetric wrt ky for mxx,myy,mzz,mxz (or fx,fz) components
c	and cos(ky*phisrc) dependence (i.e. those in d(_,1) output)
c	x- and z-displacements are antisymmetric wrt ky for mxy,myz (or fy) components
c	and cos(ky*phisrc) dependence (i.e. those in d(_,2) output)
c	y-displacement is antisymmetric wrt ky for mxx,myy,mzz,mxz (or fx,fz) components
c	and cos(ky*phisrc) dependence (i.e. those in d(_,1) output)
c	y-displacement is symmetric wrt ky for mxy,myz (or fy) components
c	and cos(ky*phisrc) dependence (i.e. those in d(_,2) output)
c	The situation is reversed for the components with
c	-ui*sin(ky*phisrc) dependence (i.e., those with d(_,3) and d(_,4) output)
	write(4) d(i1,1)+d(i1,2)+d(i1,3)+d(i1,4),d(i2,1)+d(i2,2)+d(i2,3)+d(i2,4),
     &	d(i3,1)+d(i3,2)+d(i3,3)+d(i3,4)
	if(m.gt.minc) write(4) d(i1,1)-d(i1,2)-d(i1,3)+d(i1,4),-d(i2,1)+d(i2,2)+d(i2,3)-d(i2,4),
     &	d(i3,1)-d(i3,2)-d(i3,3)+d(i3,4)
	enddo
	enddo
	enddo
	close(4)
c**	Store values of Laplace-transformed and phi-transformed displacement
c	at the free surface only.  (Could do all points, but this might generate
c	a very large displ-j-m file.)
	if(iom.eq.iom1.and.m.eq.minc) then
	open(4,file=disfil,form='unformatted')
	else
	open(4,file=disfil,form='unformatted',status='old',access='append')
	endif
cOLD	do ig=1,(N*NX-(NX-1))*(N*NZ-(NZ-1))
	ncz=NZ
	ilz=N
	do ncx=1,NX
	do ilx=1,N
	il=N*(ilz-1)+ilx
	ig=igrd(ncx,ncz,il)
	i1=3*ig-2
	i2=i1+1
	i3=i2+1
c	x- and z-displacements are symmetric wrt ky for mxx,myy,mzz,mxz components
c	and cos(ky*phisrc) dependence (i.e. those in d(_,1) output)
c	x- and z-displacements are antisymmetric wrt ky for mxy,myz components
c	and cos(ky*phisrc) dependence (i.e. those in d(_,2) output)
c	y-displacement is antisymmetric wrt ky for mxx,myy,mzz,mxz components
c	and cos(ky*phisrc) dependence (i.e. those in d(_,1) output)
c	y-displacement is symmetric wrt ky for mxy,myz components
c	and cos(ky*phisrc) dependence (i.e. those in d(_,2) output)
c	The situation is reversed for the components with
c	-ui*sin(ky*phisrc) dependence (i.e., those with d(_,3) and d(_,4) output)
	write(4) d(i1,1)+d(i1,2)+d(i1,3)+d(i1,4),d(i2,1)+d(i2,2)+d(i2,3)+d(i2,4),
     &	d(i3,1)+d(i3,2)+d(i3,3)+d(i3,4)
	if(m.gt.minc) write(4) d(i1,1)-d(i1,2)-d(i1,3)+d(i1,4),-d(i2,1)+d(i2,2)+d(i2,3)-d(i2,4),
     &	d(i3,1)-d(i3,2)-d(i3,3)+d(i3,4)
cOLD	if(m.gt.minc) write(4) d(i1,1)-d(i1,2),-d(i2,1)+d(i2,2),d(i3,1)-d(i3,2)
c	if(ig.eq.1.or.ig.eq.2308) then
c	write(6,*)'After matrel ig=',ig,'s=',s,'m=',m,'ky=',ky
c	write(6,*)'After matrel d(',i1,1,')=',d(i1,1)
c	write(6,*)'After matrel d(',i2,1,')=',d(i2,1)
c	write(6,*)'After matrel d(',i3,1,')=',d(i3,1)
c	write(6,*)'After matrel d(',i1,2,')=',d(i1,2)
c	write(6,*)'After matrel d(',i2,2,')=',d(i2,2)
c	write(6,*)'After matrel d(',i3,2,')=',d(i3,2)
c	write(6,*)'After matrel d(',i1,3,')=',d(i1,3)
c	write(6,*)'After matrel d(',i2,3,')=',d(i2,3)
c	write(6,*)'After matrel d(',i3,3,')=',d(i3,3)
c	write(6,*)'After matrel d(',i1,4,')=',d(i1,4)
c	write(6,*)'After matrel d(',i2,4,')=',d(i2,4)
c	write(6,*)'After matrel d(',i3,4,')=',d(i3,4)
c	endif
	enddo
	enddo	
	close(4)
c	Update dsverp and dsvrpp arrays with sum over wavenumber at
c	presumed phi=0.
	do ig=1,(N*NX-(NX-1))*(N*NZ-(NZ-1))
	i1=3*ig-2
	i2=i1+1
	i3=i2+1
	facp=dble(minc)/twopi
	facm=facp
	dsverp(i1)=dsverp(i1)+(d(i1,1)+d(i1,2)+d(i1,3)+d(i1,4))*facp
	dsverp(i2)=dsverp(i2)+(d(i2,1)+d(i2,2)+d(i2,3)+d(i2,4))*facp
	dsverp(i3)=dsverp(i3)+(d(i3,1)+d(i3,2)+d(i3,3)+d(i3,4))*facp
	dsvrpp(i1)=dsvrpp(i1)+(d(i1,1)+d(i1,2)+d(i1,3)+d(i1,4))*facp*ui*dble(m-minc)
	dsvrpp(i2)=dsvrpp(i2)+(d(i2,1)+d(i2,2)+d(i2,3)+d(i2,4))*facp*ui*dble(m-minc)
	dsvrpp(i3)=dsvrpp(i3)+(d(i3,1)+d(i3,2)+d(i3,3)+d(i3,4))*facp*ui*dble(m-minc)
	dsvrp2(i1)=dsvrp2(i1)-(d(i1,1)+d(i1,2)+d(i1,3)+d(i1,4))*facp*dble(m-minc)**2
	dsvrp2(i2)=dsvrp2(i2)-(d(i2,1)+d(i2,2)+d(i2,3)+d(i2,4))*facp*dble(m-minc)**2
	dsvrp2(i3)=dsvrp2(i3)-(d(i3,1)+d(i3,2)+d(i3,3)+d(i3,4))*facp*dble(m-minc)**2
	if(m.gt.minc) then
	dsverp(i1)=dsverp(i1)+(d(i1,1)-d(i1,2)-d(i1,3)+d(i1,4))*facm
	dsverp(i2)=dsverp(i2)+(-d(i2,1)+d(i2,2)+d(i2,3)-d(i2,4))*facm
	dsverp(i3)=dsverp(i3)+(d(i3,1)-d(i3,2)-d(i3,3)+d(i3,4))*facm
	dsvrpp(i1)=dsvrpp(i1)-(d(i1,1)-d(i1,2)-d(i1,3)+d(i1,4))*facm*ui*dble(m-minc)
	dsvrpp(i2)=dsvrpp(i2)-(-d(i2,1)+d(i2,2)+d(i2,3)-d(i2,4))*facm*ui*dble(m-minc)
	dsvrpp(i3)=dsvrpp(i3)-(d(i3,1)-d(i3,2)-d(i3,3)+d(i3,4))*facm*ui*dble(m-minc)
	dsvrp2(i1)=dsvrp2(i1)-(d(i1,1)-d(i1,2)-d(i1,3)+d(i1,4))*facm*dble(m-minc)**2
	dsvrp2(i2)=dsvrp2(i2)-(-d(i2,1)+d(i2,2)+d(i2,3)-d(i2,4))*facm*dble(m-minc)**2
	dsvrp2(i3)=dsvrp2(i3)-(d(i3,1)-d(i3,2)-d(i3,3)+d(i3,4))*facm*dble(m-minc)**2
	endif
c		if(ig.ge.2467.and.ig.le.2472) then
c		write(6,*)'DSVERP: d(',i1,'1)=',d(i1,1),'d(',i2,'1)=',d(i2,1),'d(',i3,'1)=',d(i3,1),'d(',i4,'1)=',d(i4,1)
c		write(6,*)'DSVERP: facp=',facp,'dble(m-minc)=',dble(m-minc)
c		write(6,*)'DSVERP - - - - - - - - - - - -'
c		endif
	enddo
	  enddo
c	Write out latest dsverp-array and dsvrpp-array
	if(iom.eq.iom1) then
	open(2,file='vertp-j',form='unformatted')
	else
	open(2,file='vertp-j',form='unformatted',status='old',access='append')
	endif
	do ig=1,(N*NX-(NX-1))*(N*NZ-(NZ-1))
	i1=3*ig-2
	i2=i1+1
	i3=i2+1
	write(2) dsverp(i1),dsverp(i2),dsverp(i3)
	write(2) dsvrpp(i1),dsvrpp(i2),dsvrpp(i3)
	write(2) dsvrp2(i1),dsvrp2(i2),dsvrp2(i3)
c		write(6,*)'VERTP-J writeout -- ig=',ig
c		write(6,*)'VERTP-J writeout ',dsverp(i1),dsverp(i2),dsverp(i3)
c		write(6,*)'VERTP-J writeout ',dsvrpp(i1),dsvrpp(i2),dsvrpp(i3)
c		write(6,*)'VERTP-J writeout ',dsvrp2(i1),dsvrp2(i2),dsvrp2(i3)
c		write(6,*)'------------------'
	enddo
	close(2)

	enddo
	close(4)
	close(8)
	write(6,*)'line C: plat,plon=',plat,plon
	if(imatr.lt.2) stop
c****************
20	continue
	write(6,*)'after 20 continue'
	deallocate (Ai)
	deallocate (jAi)
	deallocate (AA0)
	deallocate (AA1)
	deallocate (AA2)
c*-*-*-*
c        jmax1=nleno/2
c	jmax1=int(2.*real(nleno)*dt/corper)
        jmax1=nleno/3
	write(6,*)'SEIS2PT5D: jmax1=',jmax1
c	Write out the time derivative of the displacement field and strain field.
	ityp=1
		write(6,*)'calling readp'
	  call readp(nleno)
		write(6,*)'out of readp'
cSAVESPACE
	if(imatr.eq.3) go to 21
c--
	open(2,file='seis2pt5d.outxyz_vertp',form='unformatted')

	do ncx=1,NX
	do ncz=1,NZ
	do ilz=1,N
	do ilx=1,N
	l=N*(ilz-1)+ilx
	ig=igrd(ncx,ncz,l)
cOLD	do ig=1,(N*NX-(NX-1))*(N*NZ-(NZ-1))
	i1=3*ig-2
	i2=i1+1
	i3=i2+1
	xgdble=dble(xg(ig))
	fac1=1.d0/dble(bigr+zg(ig))
c	First x-displacement
	  do iom=1,jmax1
	  xo(iom)=dsvert(i1,iom)
c		write(6,*)'ncx,ncy,ilz,ilx=',ncx,ncy,ilz,ilx,'xg(',ig,')=',xg(ig),'fac1=',fac1
c		write(6,*)'x-displ: dsvert(',i1,iom,')=',dsvert(i1,iom)
c		write(6,*)'----------------'
	  enddo
		write(6,*)'VERTP: At surface: entering invft with args',nleno,dt,obeta,corper
	call invft(nleno,dt,obeta,corper)
	do it=1,1024
c		write(6,*)'Out of invft: xt(',it,')=',xt(it)
	xtx(it)=xt(it)
	enddo
c	Next y-displacement
	  do iom=1,jmax1
	  xo(iom)=dsvert(i2,iom)
	  enddo
		write(6,*)'At surface: entering invft with args',nleno,dt,obeta,corper
	call invft(nleno,dt,obeta,corper)
	do it=1,1024
	xty(it)=xt(it)
	enddo
c	Next z-displacement
	  do iom=1,jmax1
	  xo(iom)=dsvert(i3,iom)
	  enddo
		write(6,*)'At surface: entering invft with args',nleno,dt,obeta,corper
	call invft(nleno,dt,obeta,corper)
	do it=1,1024
	xtz(it)=xt(it)
	enddo

c	Write out in GMT format using distance from symmetry pole and depth.
	phival=0.
	xgdble=dble(xg(ig))
	phdble=phiref+dble(phival)
		write(6,*)'xgdble,phdble=',xgdble,phdble
	cdelt=dcos(xgdble)*dcos(plat)+dsin(xgdble)*dsin(plat)*dcos(phdble)
        delta=dacos(cdelt)
	sdelt=dsin(delta)
        sfhi=dsin(plat)*dsin(phdble)/sdelt
        cfhi=(dcos(plat)-dcos(xgdble)*cdelt)/(dsin(xgdble)*sdelt)
	fhi=datan2(sfhi,cfhi)
        spsi=dsin(xgdble)*dsin(phdble)/sdelt
	cpsi=(dcos(xgdble)-dcos(plat)*cdelt)/(dsin(plat)*sdelt)
	psi=datan2(spsi,cpsi)
		write(6,*)'delta,fhi,psi=',delta,fhi,psi,'plat,plon=',plat,plon

	do it=1,1024
	t=dt*(real(nleno)/2048.)*real(it)
c	Full wavefield displacements
	dspx=xtx(it)
	dspy=xty(it)
	dspz=xtz(it)
c		write(2,*) t,real(plon-psi)*rad,90.-rad*real(delta)
c		write(2,*) (1.e+4)*(dspy*real(cfhi)-dspx*real(sfhi)),(1.e+4)*(-dspy*real(sfhi)-dspx*real(cfhi)),
c     &  (1.e+4)*real(dspz)
c	write(2,*) t,xg(ig),zg(ig),
c* * * 	Units: KU=F; F has units of (10^20 N m)/(10^3 m)^2 (subroutine source);
c	K has units of 10^10 Pa.  Hence U has units of 10^4 m.
c* * *
	write(2) t,real(plon-psi)*rad,90.-rad*real(delta),zg(ig),
     &	(1.e+4)*(dspy*real(cfhi)-dspx*real(sfhi)),(1.e+4)*(-dspy*real(sfhi)-dspx*real(cfhi)),
     &	(1.e+4)*real(dspz)
	enddo

	enddo
	enddo
	enddo
	enddo
	close(2)
cSAVESPACE
	if(imatr.eq.4) go to 23
	open(4,file='seis2pt5d.outstrains_vertp',form='unformatted')

	do ncx=1,NX
	do ncz=1,NZ
	do ilz=1,N
	do ilx=1,N
	l=N*(ilz-1)+ilx
	ig=igrd(ncx,ncz,l)
cOLD	do ig=1,(N*NX-(NX-1))*(N*NZ-(NZ-1))
	i1=3*ig-2
	i2=i1+1
	i3=i2+1
	xgdble=dble(xg(ig))
	fac1=1.d0/dble(bigr+zg(ig))

c	Next tensor strain component exx=d(x-displacement)/dx
	do iom=1,jmax1
	xo(iom)=0.d0
	enddo
c	For strains at gridpoint #l in cell (ncx,ncz), use Lagrangian interpolation.
	do ikz=1,N
	do ikx=1,N
	k=N*(ikz-1)+ikx
	igk=igrd(ncx,ncz,k)
	k1=3*igk-2
	k2=k1+1
	k3=k2+1
	do iom=1,jmax1
	xo(iom)=xo(iom) + dsvert(k1,iom)*dphi2x(k,l)*dble(2./dx(ncx)) * fac1
c		if(iom.eq.2) then
c		write(6,*)'REC LOOP: dsvert(',k1,iom,')=',dsvert(k1,iom),'dphi2x(',k,l,')=',dphi2x(k,l)
c		write(6,*)'REC LOOP: dx(',ncx,')=',dx(ncx),'fac1=',fac1
c		write(6,*)'REC LOOP: latest xo(',iom,')=',xo(iom)
c		write(6,*)'- - - - - - - - - - - - - -'
c		endif
	enddo
	enddo
	enddo
c		write(6,*)'rec loop: entering invft with nleno,dt,obeta,corper=',nleno,dt,obeta,corper
c		write(6,*)'xo=',(xo(iom), iom=1,jmax1)
	call invft(nleno,dt,obeta,corper)
	do it=1,1024
c		write(6,*)'after invft: xt(',it,')=',xt(it)
	xtxx(it)=xt(it)
	enddo

c	Next tensor strain component exy=0.5*[d(x-displacement)/dy
c	+ d(y-displacement)/dx]
	do iom=1,jmax1
	xo(iom)=0.5d0*dsvrtp(i1,iom)/dsin(xgdble) * fac1
	enddo
c	For strains at gridpoint #l in cell (ncx,ncz), use Lagrangian interpolation.
	do ikz=1,N
	do ikx=1,N
	k=N*(ikz-1)+ikx
	igk=igrd(ncx,ncz,k)
	k1=3*igk-2
	k2=k1+1
	k3=k2+1
	do iom=1,jmax1
     	xo(iom)=xo(iom) + 0.5d0*dsvert(k2,iom)*dphi2x(k,l)*dble(2./dx(ncx)) * fac1
	enddo
	enddo
	enddo
c		write(6,*)'rec loop: entering invft with nleno,dt,obeta,corper=',nleno,dt,obeta,corper
c		write(6,*)'xo=',(xo(iom), iom=1,jmax1)
	call invft(nleno,dt,obeta,corper)
	do it=1,1024
c		write(6,*)'after invft: xt(',it,')=',xt(it)
	xtxy(it)=xt(it)
	enddo

c	Next tensor strain component exz=0.5*[d(x-displacement)/dz
c	+ d(z-displacement)/dx]
	do iom=1,jmax1
	xo(iom)=0.d0
	enddo
c	For strains at gridpoint #l in cell (ncx,ncz), use Lagrangian interpolation.
cTOPO
c	Determine topographic stretching factor at (xs,zs) in cell (ncx,ncz).
c	thfac1=x' of notes; x goes from -1 to 1, x' goes from 0 to 1 across the cell.
	thfac1=(xs-x(1))/(x(N)-x(1))
c	thfac2=1-x' of notes.
	thfac2=(x(N)-xs)/(x(N)-x(1))
	if(ncz.gt.1) then
	rfac1=1.+(dztopo(ncz,ncx+1)-dztopo(ncz-1,ncx+1))/dz(ncz)
	rfac2=1.+(dztopo(ncz,ncx)-dztopo(ncz-1,ncx))/dz(ncz)
	endif
	if(ncz.eq.1) then
	rfac1=1.+dztopo(ncz,ncx+1)/dz(ncz)
	rfac2=1.+dztopo(ncz,ncx)/dz(ncz)
	endif
	sjacm1=thfac1*rfac1+thfac2*rfac2
c--
	do ikz=1,N
	do ikx=1,N
	k=N*(ikz-1)+ikx
	igk=igrd(ncx,ncz,k)
	k1=3*igk-2
	k2=k1+1
	k3=k2+1
	do iom=1,jmax1
	xo(iom)=xo(iom) + 0.5d0*dsvert(k1,iom)*dphi2z(k,l)*dble(2./(sjacm1*dz(ncz)))
     &	+ 0.5d0*dsvert(k3,iom)*dphi2x(k,l)*dble(2./dx(ncx)) * fac1
	enddo
	enddo
	enddo
c		write(6,*)'rec loop: entering invft with nleno,dt,obeta,corper=',nleno,dt,obeta,corper
c		write(6,*)'xo=',(xo(iom), iom=1,jmax1)
	call invft(nleno,dt,obeta,corper)
	do it=1,1024
c		write(6,*)'after invft: xt(',it,')=',xt(it)
	xtxz(it)=xt(it)
	enddo

c	Next tensor strain component eyy=d(y-displacement)/dy
	do iom=1,jmax1
	xo(iom)=dsvrtp(i2,iom)/dsin(xgdble) * fac1
	enddo
c		write(6,*)'rec loop: entering invft with nleno,dt,obeta,corper=',nleno,dt,obeta,corper
c		write(6,*)'xo=',(xo(iom), iom=1,jmax1)
	call invft(nleno,dt,obeta,corper)
	do it=1,1024
c		write(6,*)'after invft: xt(',it,')=',xt(it)
	xtyy(it)=xt(it)
	enddo

c	Next tensor strain component eyz=0.5*[d(y-displacement)/dz
c	+ d(z-displacement)/dy]
	do iom=1,jmax1
	xo(iom)=0.5d0*dsvrtp(i3,iom)/dsin(xgdble) * fac1
	enddo
c	For strains at gridpoint #l in cell (ncx,ncz), use Lagrangian interpolation.
cTOPO
c	Determine topographic stretching factor at (xs,zs) in cell (ncx,ncz).
c	thfac1=x' of notes; x goes from -1 to 1, x' goes from 0 to 1 across the cell.
	thfac1=(xs-x(1))/(x(N)-x(1))
c	thfac2=1-x' of notes.
	thfac2=(x(N)-xs)/(x(N)-x(1))
	if(ncz.gt.1) then
	rfac1=1.+(dztopo(ncz,ncx+1)-dztopo(ncz-1,ncx+1))/dz(ncz)
	rfac2=1.+(dztopo(ncz,ncx)-dztopo(ncz-1,ncx))/dz(ncz)
	endif
	if(ncz.eq.1) then
	rfac1=1.+dztopo(ncz,ncx+1)/dz(ncz)
	rfac2=1.+dztopo(ncz,ncx)/dz(ncz)
	endif
	sjacm1=thfac1*rfac1+thfac2*rfac2
c--
	do ikz=1,N
	do ikx=1,N
	k=N*(ikz-1)+ikx
	igk=igrd(ncx,ncz,k)
	k1=3*igk-2
	k2=k1+1
	k3=k2+1
	do iom=1,jmax1
	xo(iom)=xo(iom) + 0.5d0*dsvert(k2,iom)*dphi2z(k,l)*dble(2./(sjacm1*dz(ncz)))
	enddo
	enddo
	enddo
c		write(6,*)'rec loop: entering invft with nleno,dt,obeta,corper=',nleno,dt,obeta,corper
c		write(6,*)'xo=',(xo(iom), iom=1,jmax1)
	call invft(nleno,dt,obeta,corper)
	do it=1,1024
c		write(6,*)'after invft: xt(',it,')=',xt(it)
	xtyz(it)=xt(it)
	enddo

c	Next tensor strain component ezz=d(z-displacement)/dz
	do iom=1,jmax1
	xo(iom)=0.d0
	enddo
c	For strains at gridpoint #l in cell (ncx,ncz), use Lagrangian interpolation.
cTOPO
c	Determine topographic stretching factor at (xs,zs) in cell (ncx,ncz).
c	thfac1=x' of notes; x goes from -1 to 1, x' goes from 0 to 1 across the cell.
	thfac1=(xs-x(1))/(x(N)-x(1))
c	thfac2=1-x' of notes.
	thfac2=(x(N)-xs)/(x(N)-x(1))
	if(ncz.gt.1) then
	rfac1=1.+(dztopo(ncz,ncx+1)-dztopo(ncz-1,ncx+1))/dz(ncz)
	rfac2=1.+(dztopo(ncz,ncx)-dztopo(ncz-1,ncx))/dz(ncz)
	endif
	if(ncz.eq.1) then
	rfac1=1.+dztopo(ncz,ncx+1)/dz(ncz)
	rfac2=1.+dztopo(ncz,ncx)/dz(ncz)
	endif
	sjacm1=thfac1*rfac1+thfac2*rfac2
c--
	do ikz=1,N
	do ikx=1,N
	k=N*(ikz-1)+ikx
	igk=igrd(ncx,ncz,k)
	k1=3*igk-2
	k2=k1+1
	k3=k2+1
	do iom=1,jmax1
	xo(iom)=xo(iom) + dsvert(k3,iom)*dphi2z(k,l)*dble(2./(sjacm1*dz(ncz)))
	enddo
	enddo
	enddo
c		write(6,*)'rec loop: entering invft with nleno,dt,obeta,corper=',nleno,dt,obeta,corper
c		write(6,*)'xo=',(xo(iom), iom=1,jmax1)
	call invft(nleno,dt,obeta,corper)
	do it=1,1024
c		write(6,*)'after invft: xt(',it,')=',xt(it)
	xtzz(it)=xt(it)
	enddo


c	Write out in GMT format using distance from symmetry pole and depth.
	phival=0.
	xgdble=dble(xg(ig))
	phdble=phiref+dble(phival)
		write(6,*)'xgdble,phdble=',xgdble,phdble
	cdelt=dcos(xgdble)*dcos(plat)+dsin(xgdble)*dsin(plat)*dcos(phdble)
        delta=dacos(cdelt)
	sdelt=dsin(delta)
        sfhi=dsin(plat)*dsin(phdble)/sdelt
        cfhi=(dcos(plat)-dcos(xgdble)*cdelt)/(dsin(xgdble)*sdelt)
	fhi=datan2(sfhi,cfhi)
        spsi=dsin(xgdble)*dsin(phdble)/sdelt
	cpsi=(dcos(xgdble)-dcos(plat)*cdelt)/(dsin(plat)*sdelt)
	psi=datan2(spsi,cpsi)
		write(6,*)'delta,fhi,psi=',delta,fhi,psi,'plat,plon=',plat,plon

	do it=1,1024
	t=dt*(real(nleno)/2048.)*real(it)
c	Rotate strain components into primed coord system: 
c	xhat' = -sfhi*xhat + cfhi*yhat
c	yhat' = -cfhi*xhat - sfhi*yhat
c	zhat = zhat
	expxp = xtxx(it)*sfhi**2 - 2.*sfhi*cfhi*xtxy(it) + xtyy(it)*cfhi**2
	eypyp = xtxx(it)*cfhi**2 + 2.*sfhi*cfhi*xtxy(it) + xtyy(it)*sfhi**2
	expyp = (xtxx(it)-xtyy(it))*sfhi*cfhi + xtxy(it)*(sfhi**2-cfhi**2)
	expzp = -xtxz(it)*sfhi + xtyz(it)*cfhi
	eypzp = -xtxz(it)*cfhi - xtyz(it)*sfhi
	ezpzp = xtzz(it)
c* * *	Units: Original U=F/K displacement has units of 10^4 m; division by 10^3 meters
c	yields a dimensionless factor of 10.  Hence multiplication by 1.e+7
c* * *	to yield strain in units of microstrain.
	if(imatr.ne.4) write(4) t,real(plon-psi)*rad,90.-rad*real(delta),zg(ig),(1.e+7)*expxp,
     &	(1.e+7)*eypyp,(1.e+7)*expyp,(1.e+7)*expzp,(1.e+7)*eypzp,(1.e+7)*ezpzp
	enddo

	enddo
	enddo
	enddo
	enddo
	close(4)
c*-*-*-*
21	continue
	write(6,*)'after do 21'
c	Write out the time derivative of the (vector) curl and divergence of the displacement field.
	ityp=1
	open(4,file='seis2pt5d.outCURL+DIVxyz_vertp')
	do ncx=1,NX
	do ncz=1,NZ
	do ilz=1,N
	do ilx=1,N
	l=N*(ilz-1)+ilx
	ig=igrd(ncx,ncz,l)
cOLD	do ig=1,(N*NX-(NX-1))*(N*NZ-(NZ-1))
	i1=3*ig-2
	i2=i1+1
	i3=i2+1
	xgdble=dble(xg(ig))
	fac1=1.d0/dble(bigr+zg(ig))

c	Next x-component of curl of displacement =[- d(y-displacement)/dz
c       + d(z-displacement)/dy]
	do iom=1,jmax1
	xo(iom)=dsvrtp(i3,iom)/dsin(xgdble) * fac1
	enddo
c	For strains at gridpoint #l in cell (ncx,ncz), use Lagrangian interpolation.
cTOPO
c	Determine topographic stretching factor at (xs,zs) in cell (ncx,ncz).
c	thfac1=x' of notes; x goes from -1 to 1, x' goes from 0 to 1 across the cell.
	thfac1=(xs-x(1))/(x(N)-x(1))
c	thfac2=1-x' of notes.
	thfac2=(x(N)-xs)/(x(N)-x(1))
	if(ncz.gt.1) then
	rfac1=1.+(dztopo(ncz,ncx+1)-dztopo(ncz-1,ncx+1))/dz(ncz)
	rfac2=1.+(dztopo(ncz,ncx)-dztopo(ncz-1,ncx))/dz(ncz)
	endif
	if(ncz.eq.1) then
	rfac1=1.+dztopo(ncz,ncx+1)/dz(ncz)
	rfac2=1.+dztopo(ncz,ncx)/dz(ncz)
	endif
	sjacm1=thfac1*rfac1+thfac2*rfac2
c--
	do ikz=1,N
	do ikx=1,N
	k=N*(ikz-1)+ikx
	igk=igrd(ncx,ncz,k)
	k1=3*igk-2
	k2=k1+1
	k3=k2+1
	do iom=1,jmax1
	xo(iom)=xo(iom) - dsvert(k2,iom)*dphi2z(k,l)*dble(2./(sjacm1*dz(ncz)))

	enddo
	enddo
	enddo
c		write(6,*)'rec loop: entering invft with nleno,dt,obeta,corper=',nleno,dt,obeta,corper
c		write(6,*)'xo=',(xo(iom), iom=1,jmax1)
	call invft(nleno,dt,obeta,corper)
	do it=1,1024
c		write(6,*)'after invft: xt(',it,')=',xt(it)
	xtcx(it)=xt(it)
	enddo

c	Next y-component of S-wavefield (curl of curl of displacement) =[d(x-displacement)/dz
c       - d(z-displacement)/dx]
	do iom=1,jmax1
        xo(iom)=0.
	enddo
c	For strains at gridpoint #l in cell (ncx,ncz), use Lagrangian interpolation.
cTOPO
c	Determine topographic stretching factor at (xs,zs) in cell (ncx,ncz).
c	thfac1=x' of notes; x goes from -1 to 1, x' goes from 0 to 1 across the cell.
	thfac1=(xs-x(1))/(x(N)-x(1))
c	thfac2=1-x' of notes.
	thfac2=(x(N)-xs)/(x(N)-x(1))
	if(ncz.gt.1) then
	rfac1=1.+(dztopo(ncz,ncx+1)-dztopo(ncz-1,ncx+1))/dz(ncz)
	rfac2=1.+(dztopo(ncz,ncx)-dztopo(ncz-1,ncx))/dz(ncz)
	endif
	if(ncz.eq.1) then
	rfac1=1.+dztopo(ncz,ncx+1)/dz(ncz)
	rfac2=1.+dztopo(ncz,ncx)/dz(ncz)
	endif
	sjacm1=thfac1*rfac1+thfac2*rfac2
c--
	do ikz=1,N
	do ikx=1,N
	k=N*(ikz-1)+ikx
	igk=igrd(ncx,ncz,k)
	k1=3*igk-2
	k2=k1+1
	k3=k2+1
	do iom=1,jmax1
	xo(iom)=xo(iom) + dsvert(k1,iom)*dphi2z(k,l)*dble(2./(sjacm1*dz(ncz)))
     &	- dsvert(k3,iom)*dphi2x(k,l)*dble(2./dx(ncx)) * fac1
	enddo
	enddo
	enddo
c		write(6,*)'rec loop: entering invft with nleno,dt,obeta,corper=',nleno,dt,obeta,corper
c		write(6,*)'xo=',(xo(iom), iom=1,jmax1)
	call invft(nleno,dt,obeta,corper)
	do it=1,1024
c		write(6,*)'after invft: xt(',it,')=',xt(it)
	xtcy(it)=xt(it)
	enddo

c	Next z-component of curl of displacement =[d(y-displacement)/dx
c       - d(x-displacement)/dy]
	do iom=1,jmax1
	xo(iom)=-dsvrtp(i1,iom)/dsin(xgdble) * fac1
	enddo
c	For strains at gridpoint #l in cell (ncx,ncz), use Lagrangian interpolation.
	do ikz=1,N
	do ikx=1,N
	k=N*(ikz-1)+ikx
	igk=igrd(ncx,ncz,k)
	k1=3*igk-2
	k2=k1+1
	k3=k2+1
	do iom=1,jmax1
     	xo(iom)=xo(iom) + dsvert(k2,iom)*dphi2x(k,l)*dble(2./dx(ncx)) * fac1
	enddo
	enddo
	enddo
c		write(6,*)'rec loop: entering invft with nleno,dt,obeta,corper=',nleno,dt,obeta,corper
c		write(6,*)'xo=',(xo(iom), iom=1,jmax1)
	call invft(nleno,dt,obeta,corper)
	do it=1,1024
c		write(6,*)'after invft: xt(',it,')=',xt(it)
	xtcz(it)=xt(it)
	enddo

c	Next tensor strain component exx=d(x-displacement)/dx
	do iom=1,jmax1
	xo(iom)=0.d0
	enddo
c	For strains at gridpoint #l in cell (ncx,ncz), use Lagrangian interpolation.
	do ikz=1,N
	do ikx=1,N
	k=N*(ikz-1)+ikx
	igk=igrd(ncx,ncz,k)
	k1=3*igk-2
	k2=k1+1
	k3=k2+1
	do iom=1,jmax1
	xo(iom)=xo(iom) + dsvert(k1,iom)*dphi2x(k,l)*dble(2./dx(ncx)) * fac1
c		if(iom.eq.2) then
c		write(6,*)'REC LOOP: dsvert(',k1,iom,')=',dsvert(k1,iom),'dphi2x(',k,l,')=',dphi2x(k,l)
c		write(6,*)'REC LOOP: dx(',ncx,')=',dx(ncx),'fac1=',fac1
c		write(6,*)'REC LOOP: latest xo(',iom,')=',xo(iom)
c		write(6,*)'- - - - - - - - - - - - - -'
c		endif
	enddo
	enddo
	enddo
c		write(6,*)'rec loop: entering invft with nleno,dt,obeta,corper=',nleno,dt,obeta,corper
c		write(6,*)'xo=',(xo(iom), iom=1,jmax1)
	call invft(nleno,dt,obeta,corper)
	do it=1,1024
c		write(6,*)'after invft: xt(',it,')=',xt(it)
	xtxx(it)=xt(it)
	enddo

c	Next tensor strain component eyy=d(y-displacement)/dy
	do iom=1,jmax1
	xo(iom)=dsvrtp(i2,iom)/dsin(xgdble) * fac1
	enddo
c		write(6,*)'rec loop: entering invft with nleno,dt,obeta,corper=',nleno,dt,obeta,corper
c		write(6,*)'xo=',(xo(iom), iom=1,jmax1)
	call invft(nleno,dt,obeta,corper)
	do it=1,1024
c		write(6,*)'after invft: xt(',it,')=',xt(it)
	xtyy(it)=xt(it)
	enddo

c	Next tensor strain component ezz=d(z-displacement)/dz
	do iom=1,jmax1
	xo(iom)=0.d0
	enddo
c	For strains at gridpoint #l in cell (ncx,ncz), use Lagrangian interpolation.
cTOPO
c	Determine topographic stretching factor at (xs,zs) in cell (ncx,ncz).
c	thfac1=x' of notes; x goes from -1 to 1, x' goes from 0 to 1 across the cell.
	thfac1=(xs-x(1))/(x(N)-x(1))
c	thfac2=1-x' of notes.
	thfac2=(x(N)-xs)/(x(N)-x(1))
	if(ncz.gt.1) then
	rfac1=1.+(dztopo(ncz,ncx+1)-dztopo(ncz-1,ncx+1))/dz(ncz)
	rfac2=1.+(dztopo(ncz,ncx)-dztopo(ncz-1,ncx))/dz(ncz)
	endif
	if(ncz.eq.1) then
	rfac1=1.+dztopo(ncz,ncx+1)/dz(ncz)
	rfac2=1.+dztopo(ncz,ncx)/dz(ncz)
	endif
	sjacm1=thfac1*rfac1+thfac2*rfac2
c--
	do ikz=1,N
	do ikx=1,N
	k=N*(ikz-1)+ikx
	igk=igrd(ncx,ncz,k)
	k1=3*igk-2
	k2=k1+1
	k3=k2+1
	do iom=1,jmax1
	xo(iom)=xo(iom) + dsvert(k3,iom)*dphi2z(k,l)*dble(2./(sjacm1*dz(ncz)))
	enddo
	enddo
	enddo
c		write(6,*)'rec loop: entering invft with nleno,dt,obeta,corper=',nleno,dt,obeta,corper
c		write(6,*)'xo=',(xo(iom), iom=1,jmax1)
	call invft(nleno,dt,obeta,corper)
	do it=1,1024
c		write(6,*)'after invft: xt(',it,')=',xt(it)
	xtzz(it)=xt(it)
	enddo

c	Write out in GMT format using distance from symmetry pole and depth.
	phival=0.
	xgdble=dble(xg(ig))
	phdble=phiref+dble(phival)
	cdelt=dcos(xgdble)*dcos(plat)+dsin(xgdble)*dsin(plat)*dcos(phdble)
        delta=dacos(cdelt)
	sdelt=dsin(delta)
        sfhi=dsin(plat)*dsin(phdble)/sdelt
        cfhi=(dcos(plat)-dcos(xgdble)*cdelt)/(dsin(xgdble)*sdelt)
	fhi=datan2(sfhi,cfhi)
        spsi=dsin(xgdble)*dsin(phdble)/sdelt
	cpsi=(dcos(xgdble)-dcos(plat)*cdelt)/(dsin(plat)*sdelt)
	psi=datan2(spsi,cpsi)
c		write(2,*)'delta,fhi,psi=',delta,fhi,psi,'plat,plon=',plat,plon

	do it=1,1024
	t=dt*(real(nleno)/2048.)*real(it)
c	time derivative of the vector curl of displacement wavefield, followed by
c	time derivative of the divergence of the displacement wavefield
	dspx=xtcx(it)
	dspy=xtcy(it)
	dspz=xtcz(it)
	write(4,*) t,real(plon-psi)*rad,90.-rad*real(delta),zg(ig),
     &	(1.e+7)*(dspy*real(cfhi)-dspx*real(sfhi)),(1.e+7)*(-dspy*real(sfhi)-dspx*real(cfhi)),
     &	(1.e+7)*real(dspz),(1.e+7)*real(xtxx(it)+xtyy(it)+xtzz(it))
	enddo

	enddo
	enddo
	enddo
	enddo
	close(4)
c*-*-*-*
c	Write out the acceleration wavefield associated with P or S-wave propagation
	ityp=0
	open(4,file='seis2pt5d.outS+Pxyz_vertp')
	do ncx=1,NX
	do ncz=1,NZ
	do ilz=1,N
	do ilx=1,N
	l=N*(ilz-1)+ilx
	ig=igrd(ncx,ncz,l)
cOLD	do ig=1,(N*NX-(NX-1))*(N*NZ-(NZ-1))
	i1=3*ig-2
	i2=i1+1
	i3=i2+1
	xgdble=dble(xg(ig))
	fac1=1.d0/dble(bigr+zg(ig))

c	Next x-component of S-wavefield (curl of curl of displacement)
	do iom=1,jmax1
        xo(iom) = -dsvrts(k1,iom) * (fac1/dsin(xgdble))**2
	enddo
c	For strains at gridpoint #l in cell (ncx,ncz), use Lagrangian interpolation.
cTOPO
c	Determine topographic stretching factor at (xs,zs) in cell (ncx,ncz).
c	thfac1=x' of notes; x goes from -1 to 1, x' goes from 0 to 1 across the cell.
	thfac1=(xs-x(1))/(x(N)-x(1))
c	thfac2=1-x' of notes.
	thfac2=(x(N)-xs)/(x(N)-x(1))
	if(ncz.gt.1) then
	rfac1=1.+(dztopo(ncz,ncx+1)-dztopo(ncz-1,ncx+1))/dz(ncz)
	rfac2=1.+(dztopo(ncz,ncx)-dztopo(ncz-1,ncx))/dz(ncz)
	endif
	if(ncz.eq.1) then
	rfac1=1.+dztopo(ncz,ncx+1)/dz(ncz)
	rfac2=1.+dztopo(ncz,ncx)/dz(ncz)
	endif
	sjacm1=thfac1*rfac1+thfac2*rfac2
c--
	do ikz=1,N
	do ikx=1,N
	k=N*(ikz-1)+ikx
	igk=igrd(ncx,ncz,k)
	k1=3*igk-2
	k2=k1+1
	k3=k2+1
	do iom=1,jmax1
	xo(iom)=xo(iom) - dsvert(k1,iom)*dphizz(k,l)*dble(2./(sjacm1*dz(ncz)))**2
     &			+ dsvrtp(k2,iom) * (fac1/dsin(xgdble)) * dphi2x(k,l)*dble(2./dx(ncx)) * fac1
     &			+ dsvert(k3,iom) * dphixz(k,l) * dble(2./dx(ncx)) * fac1
     &			  * dble(2./(sjacm1*dz(ncz)))

	enddo
	enddo
	enddo
c		write(6,*)'rec loop: entering invft with nleno,dt,obeta,corper=',nleno,dt,obeta,corper
c		write(6,*)'xo=',(xo(iom), iom=1,jmax1)
	call invft(nleno,dt,obeta,corper)
	do it=1,1024
c		write(6,*)'after invft: xt(',it,')=',xt(it)
	xtSx(it)=xt(it)
	enddo

c	Next y-component of S-wavefield (curl of curl of displacement)
	do iom=1,jmax1
        xo(iom) = 0.
	enddo
c	For strains at gridpoint #l in cell (ncx,ncz), use Lagrangian interpolation.
cTOPO
c	Determine topographic stretching factor at (xs,zs) in cell (ncx,ncz).
c	thfac1=x' of notes; x goes from -1 to 1, x' goes from 0 to 1 across the cell.
	thfac1=(xs-x(1))/(x(N)-x(1))
c	thfac2=1-x' of notes.
	thfac2=(x(N)-xs)/(x(N)-x(1))
	if(ncz.gt.1) then
	rfac1=1.+(dztopo(ncz,ncx+1)-dztopo(ncz-1,ncx+1))/dz(ncz)
	rfac2=1.+(dztopo(ncz,ncx)-dztopo(ncz-1,ncx))/dz(ncz)
	endif
	if(ncz.eq.1) then
	rfac1=1.+dztopo(ncz,ncx+1)/dz(ncz)
	rfac2=1.+dztopo(ncz,ncx)/dz(ncz)
	endif
	sjacm1=thfac1*rfac1+thfac2*rfac2
c--
	do ikz=1,N
	do ikx=1,N
	k=N*(ikz-1)+ikx
	igk=igrd(ncx,ncz,k)
	k1=3*igk-2
	k2=k1+1
	k3=k2+1
	do iom=1,jmax1
	xo(iom)=xo(iom) + dsvrtp(k1,iom) * (fac1/dsin(xgdble)) * dphi2x(k,l)*dble(2./dx(ncx)) * fac1 
     &			- dsvert(k2,iom)*dphizz(k,l) * dble(2./(sjacm1*dz(ncz)))**2
     &			- dsvert(k2,iom)*dphixx(k,l) * (dble(2./dx(ncx)) * fac1)**2
     &			+ dsvrtp(k3,iom) * (fac1/dsin(xgdble)) * dphi2z(k,l)
     &			  * dble(2./(sjacm1*dz(ncz)))
	enddo
	enddo
	enddo
c		write(6,*)'rec loop: entering invft with nleno,dt,obeta,corper=',nleno,dt,obeta,corper
c		write(6,*)'xo=',(xo(iom), iom=1,jmax1)
	call invft(nleno,dt,obeta,corper)
	do it=1,1024
c		write(6,*)'after invft: xt(',it,')=',xt(it)
	xtSy(it)=xt(it)
	enddo

c	Next z-component of S-wavefield (curl of curl of displacement)
	do iom=1,jmax1
        xo(iom) = -dsvrts(k3,iom) * (fac1/dsin(xgdble))**2
	enddo
c	For strains at gridpoint #l in cell (ncx,ncz), use Lagrangian interpolation.
cTOPO
c	Determine topographic stretching factor at (xs,zs) in cell (ncx,ncz).
c	thfac1=x' of notes; x goes from -1 to 1, x' goes from 0 to 1 across the cell.
	thfac1=(xs-x(1))/(x(N)-x(1))
c	thfac2=1-x' of notes.
	thfac2=(x(N)-xs)/(x(N)-x(1))
	if(ncz.gt.1) then
	rfac1=1.+(dztopo(ncz,ncx+1)-dztopo(ncz-1,ncx+1))/dz(ncz)
	rfac2=1.+(dztopo(ncz,ncx)-dztopo(ncz-1,ncx))/dz(ncz)
	endif
	if(ncz.eq.1) then
	rfac1=1.+dztopo(ncz,ncx+1)/dz(ncz)
	rfac2=1.+dztopo(ncz,ncx)/dz(ncz)
	endif
	sjacm1=thfac1*rfac1+thfac2*rfac2
c--
	do ikz=1,N
	do ikx=1,N
	k=N*(ikz-1)+ikx
	igk=igrd(ncx,ncz,k)
	k1=3*igk-2
	k2=k1+1
	k3=k2+1
	do iom=1,jmax1
	xo(iom)=xo(iom) + dsvert(k1,iom) * dphixz(k,l) * dble(2./dx(ncx)) * fac1
     &                    * dble(2./(sjacm1*dz(ncz)))
     &			+ dsvrtp(k2,iom) * (fac1/dsin(xgdble)) * dphi2z(k,l)
     &			  * dble(2./(sjacm1*dz(ncz)))
     &			- dsvert(k3,iom)*dphixx(k,l) * (dble(2./dx(ncx)) * fac1)**2
	enddo
	enddo
	enddo
c		write(6,*)'rec loop: entering invft with nleno,dt,obeta,corper=',nleno,dt,obeta,corper
c		write(6,*)'xo=',(xo(iom), iom=1,jmax1)
	call invft(nleno,dt,obeta,corper)
	do it=1,1024
c		write(6,*)'after invft: xt(',it,')=',xt(it)
	xtSz(it)=xt(it)
	enddo

c	Next x-component of P-wavefield (grad of div of displacement)
	do iom=1,jmax1
        xo(iom) = 0.
	enddo
c	For strains at gridpoint #l in cell (ncx,ncz), use Lagrangian interpolation.
cTOPO
c	Determine topographic stretching factor at (xs,zs) in cell (ncx,ncz).
c	thfac1=x' of notes; x goes from -1 to 1, x' goes from 0 to 1 across the cell.
	thfac1=(xs-x(1))/(x(N)-x(1))
c	thfac2=1-x' of notes.
	thfac2=(x(N)-xs)/(x(N)-x(1))
	if(ncz.gt.1) then
	rfac1=1.+(dztopo(ncz,ncx+1)-dztopo(ncz-1,ncx+1))/dz(ncz)
	rfac2=1.+(dztopo(ncz,ncx)-dztopo(ncz-1,ncx))/dz(ncz)
	endif
	if(ncz.eq.1) then
	rfac1=1.+dztopo(ncz,ncx+1)/dz(ncz)
	rfac2=1.+dztopo(ncz,ncx)/dz(ncz)
	endif
	sjacm1=thfac1*rfac1+thfac2*rfac2
c--
	do ikz=1,N
	do ikx=1,N
	k=N*(ikz-1)+ikx
	igk=igrd(ncx,ncz,k)
	k1=3*igk-2
	k2=k1+1
	k3=k2+1
	do iom=1,jmax1
	xo(iom)=xo(iom) + dsvert(k1,iom)*dphixx(k,l) * (dble(2./dx(ncx)) * fac1)**2
     &			+ dsvrtp(k2,iom) * (fac1/dsin(xgdble)) * dphi2x(k,l)*dble(2./dx(ncx)) * fac1
     &			+ dsvert(k3,iom) * dphixz(k,l) * dble(2./dx(ncx)) * fac1
     &			  * dble(2./(sjacm1*dz(ncz)))

	enddo
	enddo
	enddo
c		write(6,*)'rec loop: entering invft with nleno,dt,obeta,corper=',nleno,dt,obeta,corper
c		write(6,*)'xo=',(xo(iom), iom=1,jmax1)
	call invft(nleno,dt,obeta,corper)
	do it=1,1024
c		write(6,*)'after invft: xt(',it,')=',xt(it)
	xtPx(it)=xt(it)
	enddo

c	Next y-component of P-wavefield (grad of div of displacement)
	do iom=1,jmax1
        xo(iom) = dsvrts(k2,iom) * (fac1/dsin(xgdble))**2
	enddo
c	For strains at gridpoint #l in cell (ncx,ncz), use Lagrangian interpolation.
cTOPO
c	Determine topographic stretching factor at (xs,zs) in cell (ncx,ncz).
c	thfac1=x' of notes; x goes from -1 to 1, x' goes from 0 to 1 across the cell.
	thfac1=(xs-x(1))/(x(N)-x(1))
c	thfac2=1-x' of notes.
	thfac2=(x(N)-xs)/(x(N)-x(1))
	if(ncz.gt.1) then
	rfac1=1.+(dztopo(ncz,ncx+1)-dztopo(ncz-1,ncx+1))/dz(ncz)
	rfac2=1.+(dztopo(ncz,ncx)-dztopo(ncz-1,ncx))/dz(ncz)
	endif
	if(ncz.eq.1) then
	rfac1=1.+dztopo(ncz,ncx+1)/dz(ncz)
	rfac2=1.+dztopo(ncz,ncx)/dz(ncz)
	endif
	sjacm1=thfac1*rfac1+thfac2*rfac2
c--
	do ikz=1,N
	do ikx=1,N
	k=N*(ikz-1)+ikx
	igk=igrd(ncx,ncz,k)
	k1=3*igk-2
	k2=k1+1
	k3=k2+1
	do iom=1,jmax1
	xo(iom)=xo(iom) + dsvrtp(k1,iom) * (fac1/dsin(xgdble)) * dphi2x(k,l)*dble(2./dx(ncx)) * fac1 
     &			+ dsvrtp(k3,iom) * (fac1/dsin(xgdble)) * dphi2z(k,l)
     &			  * dble(2./(sjacm1*dz(ncz)))
	enddo
	enddo
	enddo
c		write(6,*)'rec loop: entering invft with nleno,dt,obeta,corper=',nleno,dt,obeta,corper
c		write(6,*)'xo=',(xo(iom), iom=1,jmax1)
	call invft(nleno,dt,obeta,corper)
	do it=1,1024
c		write(6,*)'after invft: xt(',it,')=',xt(it)
	xtPy(it)=xt(it)
	enddo

c	Next z-component of P-wavefield (grad of div of displacement)
	do iom=1,jmax1
        xo(iom) = 0.
	enddo
c	For strains at gridpoint #l in cell (ncx,ncz), use Lagrangian interpolation.
cTOPO
c	Determine topographic stretching factor at (xs,zs) in cell (ncx,ncz).
c	thfac1=x' of notes; x goes from -1 to 1, x' goes from 0 to 1 across the cell.
	thfac1=(xs-x(1))/(x(N)-x(1))
c	thfac2=1-x' of notes.
	thfac2=(x(N)-xs)/(x(N)-x(1))
	if(ncz.gt.1) then
	rfac1=1.+(dztopo(ncz,ncx+1)-dztopo(ncz-1,ncx+1))/dz(ncz)
	rfac2=1.+(dztopo(ncz,ncx)-dztopo(ncz-1,ncx))/dz(ncz)
	endif
	if(ncz.eq.1) then
	rfac1=1.+dztopo(ncz,ncx+1)/dz(ncz)
	rfac2=1.+dztopo(ncz,ncx)/dz(ncz)
	endif
	sjacm1=thfac1*rfac1+thfac2*rfac2
c--
	do ikz=1,N
	do ikx=1,N
	k=N*(ikz-1)+ikx
	igk=igrd(ncx,ncz,k)
	k1=3*igk-2
	k2=k1+1
	k3=k2+1
	do iom=1,jmax1
	xo(iom)=xo(iom) + dsvert(k1,iom) * dphixz(k,l) * dble(2./dx(ncx)) * fac1
     &                    * dble(2./(sjacm1*dz(ncz)))
     &			+ dsvrtp(k2,iom) * (fac1/dsin(xgdble)) * dphi2z(k,l)
     &			  * dble(2./(sjacm1*dz(ncz)))
     &			+ dsvert(k3,iom)*dphizz(k,l) * dble(2./(sjacm1*dz(ncz)))**2
	enddo
	enddo
	enddo
c		write(6,*)'rec loop: entering invft with nleno,dt,obeta,corper=',nleno,dt,obeta,corper
c		write(6,*)'xo=',(xo(iom), iom=1,jmax1)
	call invft(nleno,dt,obeta,corper)
	do it=1,1024
c		write(6,*)'after invft: xt(',it,')=',xt(it)
	xtPz(it)=xt(it)
	enddo

c	Write out in GMT format using distance from symmetry pole and depth.
	phival=0.
	xgdble=dble(xg(ig))
	phdble=phiref+dble(phival)
	cdelt=dcos(xgdble)*dcos(plat)+dsin(xgdble)*dsin(plat)*dcos(phdble)
        delta=dacos(cdelt)
	sdelt=dsin(delta)
        sfhi=dsin(plat)*dsin(phdble)/sdelt
        cfhi=(dcos(plat)-dcos(xgdble)*cdelt)/(dsin(xgdble)*sdelt)
	fhi=datan2(sfhi,cfhi)
        spsi=dsin(xgdble)*dsin(phdble)/sdelt
	cpsi=(dcos(xgdble)-dcos(plat)*cdelt)/(dsin(plat)*sdelt)
	psi=datan2(spsi,cpsi)
c		write(2,*)'delta,fhi,psi=',delta,fhi,psi,'plat,plon=',plat,plon

	do it=1,1024
	t=dt*(real(nleno)/2048.)*real(it)
c	S-wavefield displacements; note multiplication by the negative squared S wavespeed
c	to convert to acceleration
	dspx=xtSx(it)*(-vs(ncx,ncz,l)**2)
	dspy=xtSy(it)*(-vs(ncx,ncz,l)**2)
	dspz=xtSz(it)*(-vs(ncx,ncz,l)**2)
c	P-wavefield displacements; note multiplication by the squared P wavespeed
c	to convert to acceleration
	dspx2=xtPx(it)*(vp(ncx,ncz,l)**2)
	dspy2=xtPy(it)*(vp(ncx,ncz,l)**2)
	dspz2=xtPz(it)*(vp(ncx,ncz,l)**2)
	write(4,*) t,real(plon-psi)*rad,90.-rad*real(delta),zg(ig),
     &	(1.e+4)*(dspy*real(cfhi)-dspx*real(sfhi)),(1.e+4)*(-dspy*real(sfhi)-dspx*real(cfhi)),
     &	(1.e+4)*real(dspz),
     &	(1.e+4)*(dspy2*real(cfhi)-dspx2*real(sfhi)),(1.e+4)*(-dspy2*real(sfhi)-dspx2*real(cfhi)),
     &	(1.e+4)*real(dspz2)
	enddo

	enddo
	enddo
	enddo
	enddo
	close(4)
c*-*-*-*
c21	continue
23	continue
	ityp=1
	call addmr(nleno,nrec,phir)
	open(2,file='seis2pt5d.outxyz_rec')
	open(4,file='seis2pt5d.outstrains_rec')
c*
		do nr=1,nrec
	xs=xrec(nr)
	zs=zrec(nr)
	fac1=1.d0/(dble(bigr)+recz(nr))
c	First x-displacement
	do iom=1,jmax1
	xo(iom)=0.d0
	enddo
c	For displacements at (xs,zs) in cell (ncx,ncz), use Lagrangian interpolation.
c	ncz=nczs(nr)
	ic=N*N*(nr-1)
	do ilz=1,N
c	ncx=ncxs(nr)
	do ilx=1,N
	ic=ic+1
	i1=3*ic-2
	i2=i1+1
	i3=i2+1
	do iom=1,jmax1
	xo(iom)=xo(iom)+dsrec(i1,iom)*phiv(ilx,ilz,xs,zs)
c		if(iom.eq.1) then
c		write(6,*)'dsrec(',i1,iom,')=',dsrec(i1,iom)
c		write(6,*)'phiv(',ilx,ilz,xs,zs,')=',phiv(ilx,ilz,xs,zs)
c		write(6,*)'latest xo(',iom,')=',xo(iom)
c		write(6,*)'------'
c		endif
	enddo
	enddo
	enddo
c		write(6,*)'rec loop: entering invft with nleno,dt,obeta,corper=',nleno,dt,obeta,corper
c		write(6,*)'xo=',(xo(iom), iom=1,jmax1)
	call invft(nleno,dt,obeta,corper)
	do it=1,1024
c		write(6,*)'after invft: xt(',it,')=',xt(it)
	xtx(it)=xt(it)
	enddo
c	Next y-displacement
	do iom=1,jmax1
	xo(iom)=0.d0
	enddo
c	For displacements at (xs,zs) in cell (ncx,ncz), use Lagrangian interpolation.
c	ncz=nczs(nr)
	ic=N*N*(nr-1)
	do ilz=1,N
c	ncx=ncxs(nr)
	do ilx=1,N
	ic=ic+1
	i1=3*ic-2
	i2=i1+1
	i3=i2+1
	do iom=1,jmax1
	xo(iom)=xo(iom)+dsrec(i2,iom)*phiv(ilx,ilz,xs,zs)
	enddo
	enddo
	enddo
	call invft(nleno,dt,obeta,corper)
	do it=1,1024
	xty(it)=xt(it)
	enddo
c	Next z-displacement
	do iom=1,jmax1
	xo(iom)=0.d0
	enddo
c	For displacements at (xs,zs) in cell (ncx,ncz), use Lagrangian interpolation.
c	ncz=nczs(nr)
	ic=N*N*(nr-1)
	do ilz=1,N
c	ncx=ncxs(nr)
	do ilx=1,N
	ic=ic+1
	i1=3*ic-2
	i2=i1+1
	i3=i2+1
	do iom=1,jmax1
	xo(iom)=xo(iom)+dsrec(i3,iom)*phiv(ilx,ilz,xs,zs)
	enddo
	enddo
	enddo
	call invft(nleno,dt,obeta,corper)
	do it=1,1024
	xtz(it)=xt(it)
	enddo
c	Next tensor strain component exx=d(x-displacement)/dx
	do iom=1,jmax1
	xo(iom)=0.d0
	enddo
c	For strains at (xs,zs) in cell (ncx,ncz), use Lagrangian interpolation.
	ncx=ncxs(nr)
	ncz=nczs(nr)
	ic=N*N*(nr-1)
	do ilz=1,N
	do ilx=1,N
	ic=ic+1
	i1=3*ic-2
	i2=i1+1
	i3=i2+1
	do iom=1,jmax1
	xo(iom)=xo(iom) + dsrec(i1,iom)*dphivx(ilx,ilz,xs,zs)*dble(2./dx(ncx)) * fac1
	enddo
	enddo
	enddo
c		write(6,*)'rec loop: entering invft with nleno,dt,obeta,corper=',nleno,dt,obeta,corper
c		write(6,*)'xo=',(xo(iom), iom=1,jmax1)
	call invft(nleno,dt,obeta,corper)
	do it=1,1024
c		write(6,*)'after invft: xt(',it,')=',xt(it)
	xtxx(it)=xt(it)
	enddo

c	Next tensor strain component exy=0.5*[d(x-displacement)/dy
c	+ d(y-displacement)/dx]
	do iom=1,jmax1
	xo(iom)=0.d0
	enddo
c	For strains at (xs,zs) in cell (ncx,ncz), use Lagrangian interpolation.
	ncx=ncxs(nr)
	ncz=nczs(nr)
	ic=N*N*(nr-1)
	do ilz=1,N
	do ilx=1,N
	ic=ic+1
	i1=3*ic-2
	i2=i1+1
	i3=i2+1
	do iom=1,jmax1
	xo(iom)=xo(iom) + 0.5d0*dsrecp(i1,iom)*phiv(ilx,ilz,xs,zs)/dble(sin(xrecd(nr))) * fac1
     &	+ 0.5d0*dsrec(i2,iom)*dphivx(ilx,ilz,xs,zs)*dble(2./dx(ncx)) * fac1
	enddo
	enddo
	enddo
c		write(6,*)'rec loop: entering invft with nleno,dt,obeta,corper=',nleno,dt,obeta,corper
c		write(6,*)'xo=',(xo(iom), iom=1,jmax1)
	call invft(nleno,dt,obeta,corper)
	do it=1,1024
c		write(6,*)'after invft: xt(',it,')=',xt(it)
	xtxy(it)=xt(it)
	enddo

c	Next tensor strain component exz=0.5*[d(x-displacement)/dz
c	+ d(z-displacement)/dx]
	do iom=1,jmax1
	xo(iom)=0.d0
	enddo
c	For strains at (xs,zs) in cell (ncx,ncz), use Lagrangian interpolation.
	ncx=ncxs(nr)
	ncz=nczs(nr)
cTOPO
c	Determine topographic stretching factor at (xs,zs) in cell (ncx,ncz).
c	thfac1=x' of notes; x goes from -1 to 1, x' goes from 0 to 1 across the cell.
	thfac1=(xs-x(1))/(x(N)-x(1))
c	thfac2=1-x' of notes.
	thfac2=(x(N)-xs)/(x(N)-x(1))
	if(ncz.gt.1) then
	rfac1=1.+(dztopo(ncz,ncx+1)-dztopo(ncz-1,ncx+1))/dz(ncz)
	rfac2=1.+(dztopo(ncz,ncx)-dztopo(ncz-1,ncx))/dz(ncz)
	endif
	if(ncz.eq.1) then
	rfac1=1.+dztopo(ncz,ncx+1)/dz(ncz)
	rfac2=1.+dztopo(ncz,ncx)/dz(ncz)
	endif
	sjacm1=thfac1*rfac1+thfac2*rfac2
c--
	ic=N*N*(nr-1)
	do ilz=1,N
	do ilx=1,N
	ic=ic+1
	i1=3*ic-2
	i2=i1+1
	i3=i2+1
	do iom=1,jmax1
	xo(iom)=xo(iom) + 0.5d0*dsrec(i1,iom)*dphivz(ilx,ilz,xs,zs)*dble(2./(sjacm1*dz(ncz)))
     &	+ 0.5d0*dsrec(i3,iom)*dphivx(ilx,ilz,xs,zs)*dble(2./dx(ncx)) * fac1
	enddo
	enddo
	enddo
c		write(6,*)'rec loop: entering invft with nleno,dt,obeta,corper=',nleno,dt,obeta,corper
c		write(6,*)'xo=',(xo(iom), iom=1,jmax1)
	call invft(nleno,dt,obeta,corper)
	do it=1,1024
c		write(6,*)'after invft: xt(',it,')=',xt(it)
	xtxz(it)=xt(it)
	enddo

c	Next tensor strain component eyy=d(y-displacement)/dy
	do iom=1,jmax1
	xo(iom)=0.d0
	enddo
c	For strains at (xs,zs) in cell (ncx,ncz), use Lagrangian interpolation.
	ncx=ncxs(nr)
	ncz=nczs(nr)
	ic=N*N*(nr-1)
	do ilz=1,N
	do ilx=1,N
	ic=ic+1
	i1=3*ic-2
	i2=i1+1
	i3=i2+1
	do iom=1,jmax1
	xo(iom)=xo(iom) + dsrecp(i2,iom)*phiv(ilx,ilz,xs,zs)/dble(sin(xrecd(nr))) * fac1 
	enddo
	enddo
	enddo
c		write(6,*)'rec loop: entering invft with nleno,dt,obeta,corper=',nleno,dt,obeta,corper
c		write(6,*)'xo=',(xo(iom), iom=1,jmax1)
	call invft(nleno,dt,obeta,corper)
	do it=1,1024
c		write(6,*)'after invft: xt(',it,')=',xt(it)
	xtyy(it)=xt(it)
	enddo

c	Next tensor strain component eyz=0.5*[d(y-displacement)/dz
c	+ d(z-displacement)/dy]
	do iom=1,jmax1
	xo(iom)=0.d0
	enddo
c	For strains at (xs,zs) in cell (ncx,ncz), use Lagrangian interpolation.
	ncx=ncxs(nr)
	ncz=nczs(nr)
cTOPO
c	Determine topographic stretching factor at (xs,zs) in cell (ncx,ncz).
c	thfac1=x' of notes; x goes from -1 to 1, x' goes from 0 to 1 across the cell.
	thfac1=(xs-x(1))/(x(N)-x(1))
c	thfac2=1-x' of notes.
	thfac2=(x(N)-xs)/(x(N)-x(1))
	if(ncz.gt.1) then
	rfac1=1.+(dztopo(ncz,ncx+1)-dztopo(ncz-1,ncx+1))/dz(ncz)
	rfac2=1.+(dztopo(ncz,ncx)-dztopo(ncz-1,ncx))/dz(ncz)
	endif
	if(ncz.eq.1) then
	rfac1=1.+dztopo(ncz,ncx+1)/dz(ncz)
	rfac2=1.+dztopo(ncz,ncx)/dz(ncz)
	endif
	sjacm1=thfac1*rfac1+thfac2*rfac2
c--
	ic=N*N*(nr-1)
	do ilz=1,N
	do ilx=1,N
	ic=ic+1
	i1=3*ic-2
	i2=i1+1
	i3=i2+1
	do iom=1,jmax1
	xo(iom)=xo(iom) + 0.5d0*dsrec(i2,iom)*dphivz(ilx,ilz,xs,zs)*dble(2./(sjacm1*dz(ncz)))
     &	+ 0.5d0*dsrecp(i3,iom)*phiv(ilx,ilz,xs,zs)/dble(sin(xrecd(nr))) * fac1
	enddo
	enddo
	enddo
c		write(6,*)'rec loop: entering invft with nleno,dt,obeta,corper=',nleno,dt,obeta,corper
c		write(6,*)'xo=',(xo(iom), iom=1,jmax1)
	call invft(nleno,dt,obeta,corper)
	do it=1,1024
c		write(6,*)'after invft: xt(',it,')=',xt(it)
	xtyz(it)=xt(it)
	enddo

c	Next tensor strain component ezz=d(z-displacement)/dz
	do iom=1,jmax1
	xo(iom)=0.d0
	enddo
c	For strains at (xs,zs) in cell (ncx,ncz), use Lagrangian interpolation.
	ncx=ncxs(nr)
	ncz=nczs(nr)
cTOPO
c	Determine topographic stretching factor at (xs,zs) in cell (ncx,ncz).
c	thfac1=x' of notes; x goes from -1 to 1, x' goes from 0 to 1 across the cell.
	thfac1=(xs-x(1))/(x(N)-x(1))
c	thfac2=1-x' of notes.
	thfac2=(x(N)-xs)/(x(N)-x(1))
	if(ncz.gt.1) then
	rfac1=1.+(dztopo(ncz,ncx+1)-dztopo(ncz-1,ncx+1))/dz(ncz)
	rfac2=1.+(dztopo(ncz,ncx)-dztopo(ncz-1,ncx))/dz(ncz)
	endif
	if(ncz.eq.1) then
	rfac1=1.+dztopo(ncz,ncx+1)/dz(ncz)
	rfac2=1.+dztopo(ncz,ncx)/dz(ncz)
	endif
	sjacm1=thfac1*rfac1+thfac2*rfac2
c--
	ic=N*N*(nr-1)
	do ilz=1,N
	do ilx=1,N
	ic=ic+1
	i1=3*ic-2
	i2=i1+1
	i3=i2+1
	do iom=1,jmax1
	xo(iom)=xo(iom) + dsrec(i3,iom)*dphivz(ilx,ilz,xs,zs)*dble(2./(sjacm1*dz(ncz)))
	enddo
	enddo
	enddo
c		write(6,*)'rec loop: entering invft with nleno,dt,obeta,corper=',nleno,dt,obeta,corper
c		write(6,*)'xo=',(xo(iom), iom=1,jmax1)
	call invft(nleno,dt,obeta,corper)
	do it=1,1024
c		write(6,*)'after invft: xt(',it,')=',xt(it)
	xtzz(it)=xt(it)
	enddo

c	Write out in GMT format using geographic longitude and latitude.
	xgdble=deltar(nr)
	phirec=phir(nr)
	phdble=phiref+dble(phirec)
	cdelt=dcos(xgdble)*dcos(plat)+dsin(xgdble)*dsin(plat)*dcos(phdble)
        delta=dacos(cdelt)
	sdelt=dsin(delta)
        sfhi=dsin(plat)*dsin(phdble)/sdelt
        cfhi=(dcos(plat)-dcos(xgdble)*cdelt)/(dsin(xgdble)*sdelt)
	fhi=datan2(sfhi,cfhi)
        spsi=dsin(xgdble)*dsin(phdble)/sdelt
	cpsi=(dcos(xgdble)-dcos(plat)*cdelt)/(dsin(plat)*sdelt)
	psi=datan2(spsi,cpsi)
c		write(2,*)'delta,fhi,psi=',delta,fhi,psi,'plat,plon=',plat,plon

	do it=1,1024
	t=dt*(real(nleno)/2048.)*real(it)
	dspx=xtx(it)
	dspy=xty(it)
	dspz=xtz(it)
	write(2,*) t,real(plon-psi)*rad,90.-rad*real(delta),
     &	(1.e+4)*(dspy*real(cfhi)-dspx*real(sfhi)),(1.e+4)*(-dspy*real(sfhi)-dspx*real(cfhi)),
     &	(1.e+4)*real(dspz)
c	Rotate strain components into primed coord system: 
c	xhat' = -sfhi*xhat + cfhi*yhat
c	yhat' = -cfhi*xhat - sfhi*yhat
c	zhat = zhat
	expxp = xtxx(it)*sfhi**2 - 2.*sfhi*cfhi*xtxy(it) + xtyy(it)*cfhi**2
	eypyp = xtxx(it)*cfhi**2 + 2.*sfhi*cfhi*xtxy(it) + xtyy(it)*sfhi**2
	expyp = (xtxx(it)-xtyy(it))*sfhi*cfhi + xtxy(it)*(sfhi**2-cfhi**2)
	expzp = -xtxz(it)*sfhi + xtyz(it)*cfhi
	eypzp = -xtxz(it)*cfhi - xtyz(it)*sfhi
	ezpzp = xtzz(it)
	write(4,*) t,real(plon-psi)*rad,90.-rad*real(delta),(1.e+7)*expxp,
     &	(1.e+7)*eypyp,(1.e+7)*expyp,(1.e+7)*expzp,(1.e+7)*eypzp,(1.e+7)*ezpzp
	enddo
		enddo
	close(2)
	close(4)
c*-*-*-*
25	continue
c       Use the needed frequency samples for the inverse FT.
c       jmax1=int((2.*pi*2./corper)/real(ommin))
        jmax1=nleno/3
c
cSAVESPACE
	if(imatr.eq.3.or.imatr.eq.4) go to 22
c	Evaluate seismic displacements at Earth's surface at longitudes ranging from
c	-0.30 to +0.30 x the spatial periodicity, i.e.
c	-0.30*real(twopi)/(real(minc)) to
c	+0.30*real(twopi)/(real(minc))
	ityp=1
cTE
c	minc=320
c--
	  pscal=0.50
	  nphi=259
	  nphio2=nphi/2
	write(6,*)'calling addm with pscal=',pscal,'nphi=',nphi
	  call addm(nleno,pscal,nphi)
	write(6,*)'out of addm'
	open(2,file='seis2pt5d.outxyz')
	    do iphi=1,nphi
	  phival=pscal*real(twopi)/(real(minc))*real(iphi-(nphio2+1))/real(nphio2)
		write(6,*)'iphi=',iphi,'out of',nphi,'phival=',phival
c	write(6,*)'line D: iphi=',iphi,'plat,plon=',plat,plon
cOLD	write(6,*)'calling addm with iphi=',iphi,'out of 121'
cOLD	call addm(nleno,phival)
cOLD	write(6,*)'out of addm'
c	write(6,*)'line E: plat,plon=',plat,plon
	ic=0
	ncz=NZ
	do ncx=1,NX
	ilz=N
	do ilx=1,N
	il=N*(ilz-1)+ilx
	ig=igrd(ncx,ncz,il)
	ic=ic+1
	i1=3*ic-2
	i2=i1+1
	i3=i2+1
c	First x-displacement
	  do iom=1,jmax1
	  xo(iom)=ds(iphi,i1,iom)
	  enddo
		write(6,*)'At surface: entering invft with args',nleno,dt,obeta,corper
	call invft(nleno,dt,obeta,corper)
	do it=1,1024
	xtx(it)=xt(it)
	enddo
c	Next y-displacement
	  do iom=1,jmax1
	  xo(iom)=ds(iphi,i2,iom)
	  enddo
		write(6,*)'At surface: entering invft with args',nleno,dt,obeta,corper
	call invft(nleno,dt,obeta,corper)
	do it=1,1024
	xty(it)=xt(it)
	enddo
c	Next z-displacement
	  do iom=1,jmax1
	  xo(iom)=ds(iphi,i3,iom)
	  enddo
		write(6,*)'At surface: entering invft with args',nleno,dt,obeta,corper
	call invft(nleno,dt,obeta,corper)
	do it=1,1024
	xtz(it)=xt(it)
	enddo

c	Write out in GMT format using geographic longitude and latitude.
	xgdble=dble(xg(ig))
	phdble=phiref+dble(phival)
	cdelt=dcos(xgdble)*dcos(plat)+dsin(xgdble)*dsin(plat)*dcos(phdble)
        delta=dacos(cdelt)
	sdelt=dsin(delta)
        sfhi=dsin(plat)*dsin(phdble)/sdelt
        cfhi=(dcos(plat)-dcos(xgdble)*cdelt)/(dsin(xgdble)*sdelt)
	fhi=datan2(sfhi,cfhi)
        spsi=dsin(xgdble)*dsin(phdble)/sdelt
	cpsi=(dcos(xgdble)-dcos(plat)*cdelt)/(dsin(plat)*sdelt)
	psi=datan2(spsi,cpsi)
c		write(2,*)'delta,fhi,psi=',delta,fhi,psi,'plat,plon=',plat,plon

	do it=1,1024
	t=dt*(real(nleno)/2048.)*real(it)
	dspx=xtx(it)
	dspy=xty(it)
	dspz=xtz(it)
c		write(2,*) t,real(plon-psi)*rad,90.-rad*real(delta)
c		write(2,*) (1.e+4)*(dspy*real(cfhi)-dspx*real(sfhi)),(1.e+4)*(-dspy*real(sfhi)-dspx*real(cfhi)),
c     &  (1.e+4)*real(dspz)
	write(2,*) t,real(plon-psi)*rad,90.-rad*real(delta),
     &	(1.e+4)*(dspy*real(cfhi)-dspx*real(sfhi)),(1.e+4)*(-dspy*real(sfhi)-dspx*real(cfhi)),
     &	(1.e+4)*real(dspz)
	enddo

	enddo
	enddo
	    enddo
	close(2)
22	continue

10	format(a80)
	end
