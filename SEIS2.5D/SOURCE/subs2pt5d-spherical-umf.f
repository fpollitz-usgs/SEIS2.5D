	subroutine sporder(nmat,Ap,Ai,jAi)
c	For the matrix G to be constructed by subroutine MATREL,
c	characterize it in terms of row-counter matrices Ap and Ai
c	used in UMFPACK.  
c	Use the same loops as in subroutine MATREL to determine this
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
	common/grdpar/NX,NZ
	common/glbgrd/igrd(ncmax,ncmaz,N*N)
c	Next line for kmax should really be kmax=nptmax*N2*9
cc	so we assume that one of the dimensions, e.g. NZ, is <= 4/9 of ncmax.
	parameter (kmax=ncmax*ncmaz*N2*N2*9)
cUSED	dimension irow(kmax),icol(kmax)
	integer(4), ALLOCATABLE :: irow(:)
	integer(4), ALLOCATABLE :: icol(:)
	integer*8, ALLOCATABLE :: Ai0(: , :)
	integer, ALLOCATABLE :: X(:)
	integer, ALLOCATABLE :: Y(:)
	integer*8 nzmax,nmax
cUSED        parameter (nzmax = 32800000, nmax = 215000)
        parameter (nzmax = kmax, nmax = 3*(N*ncmax-(ncmax-1))*(N*ncmaz-(ncmaz-1))+1)
cUSED	integer*8 Ap,Ai,jAi
cUSED	dimension Ap(nmax),Ai(nzmax),jAi(kmax)
	integer*8 Ap
	dimension Ap(nmax)
	integer*8 Ai(nmat),jAi(nmat)
	dimension nrow(nmax-1)
c
	allocate ( irow(nmat) )
	allocate ( icol(nmat) )
	k=0
c--	Following lines follow structure of subroutine MATREL
	do ncz=1,NZ
	do ncx=1,NX
	do ilz=1,N
	do ilx=1,N
	il=N*(ilz-1)+ilx
	i1=3*igrd(ncx,ncz,il)-2
	i2=i1+1
	i3=i2+1
c	i1,i2,i3 are the row numbers corresponding to the normal equation
c	resulting from multiplication of the momentum eqn with phi_[il]
c	in cell (ncx,ncz) and integration over that cell.
c	i1,i2,i3 refer to the x,y,z components.
c	[il] is equivalent to index [j] of notes.
	jl=0
10	jl=jl+1
	if(jl.gt.N2) go to 20
c	jl is the index for the expansion coefficients of displacement
c	[jl] is equivalent to index [k] of notes.
	j1=3*igrd(ncx,ncz,jl)-2
	j2=j1+1
	j3=j2+1
c	j1,j2,j3 are the column numbers, corresponding to the
c	a_jl, b-jl, c_jl expansion coeff of displacement in cell (ncx,ncz)
	if(ilz.eq.1.and.ncz.eq.1) go to 25
cNEW
	if(ilx.eq.1.and.ncx.eq.1) go to 30
	if(ilx.eq.N.and.ncx.eq.NX) go to 30
c--
c	Do the elements (i1,j1),(i1,j2),(i1,j3),(i2,j1),(i2,j2),(i2,j3),
c	(i3,j1),(i3,j2),(i3,j3) in that order.
	k=k+1
c		if(k.gt.kmax) write(6,*)'k=',k,'kmax=',kmax,'nmax=',nmax
	irow(k)=i1
	icol(k)=j1
	k=k+1
	irow(k)=i1
	icol(k)=j2
	k=k+1
	irow(k)=i1
	icol(k)=j3
	k=k+1
	irow(k)=i2
	icol(k)=j1
	k=k+1
	irow(k)=i2
	icol(k)=j2
	k=k+1
	irow(k)=i2
	icol(k)=j3
	k=k+1
	irow(k)=i3
	icol(k)=j1
	k=k+1
	irow(k)=i3
	icol(k)=j2
	k=k+1
	irow(k)=i3
	icol(k)=j3
	go to 10
25	continue
	if(i1.eq.j1) then
	k=k+1
	irow(k)=i1
	icol(k)=j1
	endif
	if(i2.eq.j2) then 
	k=k+1
	irow(k)=i2
	icol(k)=j2
	endif
	if(i3.eq.j3) then
	k=k+1
	irow(k)=i3
	icol(k)=j3
	endif
	go to 10
30	continue
c	Do the elements (i1,j1),(i2,j2),(i3,j3) in that order
	k=k+1
	irow(k)=i1
	icol(k)=j1
	k=k+1
	irow(k)=i2
	icol(k)=j2
	k=k+1
	irow(k)=i3
	icol(k)=j3
	go to 10
20	continue
	enddo
	enddo
	enddo
	enddo	
	ktot=k
	write(6,*)'ktot=',ktot,'should be <= kmax=',kmax,'nmax=',nmax

c	Assume 0-based Ap and Ai as in UMFPACK conventions.
	ncol=3*(N*NX-(NX-1))*(N*NZ-(NZ-1))
		write(6,*)'ncol=',ncol
c	do a first pass to determine how many rows per column there are 
c	(including duplicate irow's for a given icol) and the maximum #irow's for any column
	do j=1,ncol
	nrow(j)=0
	enddo
	maxr=0
	do k=1,ktot
	ic=icol(k)
	nrow(ic)=nrow(ic)+1
	if(nrow(ic).gt.maxr) maxr=nrow(ic)
	enddo
	write(6,*)'maxinum number of row entries for any column=',maxr
	write(6,*)'Allocating Ai0'
c-----------------------------
	allocate ( Ai0(ncol,maxr) )
	allocate ( X(maxr) )
	allocate ( Y(maxr) )
c-----------------------------
	write(6,*)'Construct a preliminary Ai0 array that has the row values arranged in column indices'
	write(6,*)'including redundant entries'
c	Construct a preliminary Ai0 array that has the row values arranged in column indices
c	including redundant entries
	do j=1,ncol
	nrow(j)=0
	enddo
	do k=1,ktot
	ic=icol(k)
	nrow(ic)=nrow(ic)+1
	l=nrow(ic)
	Ai0(ic,l)=irow(k)-1
	enddo
	write(6,*)'Construct the Ai array from Ai0 while also removing the redundant entries'
	write(6,*)'Also construct Ap on the fly'
c	Construct the Ai array from Ai0 while also removing the redundanct entries
c	Also construct Ap on the fly
	l1=0
	Ap(1)=0
	do ic=1,ncol
		k=ic-10000*(ic/10000)
		if(k.eq.0) write(6,*)'Doing ic=',ic,'out of',ncol
	do l=1,nrow(ic)
c	See if row ir (ir-1 in zero-based approach) is already represented in the column ic tally
	  ival=0
	  if(l.gt.1) then
	  do l0=1,l-1
	  if(Ai0(ic,l0).eq.Ai0(ic,l)) ival=1
	  enddo
	  endif
	if(ival.eq.0) then
	l1=l1+1
c	Ai(l1)=Ai0(ic,l)
	lind=l1-Ap(ic)
	X(lind)=Ai0(ic,l)
c		write(6,*)'ic=',ic,'l1=',l1,'Ap(',ic,')=',Ap(ic),'X(',lind,')=',X(lind)
	endif
	enddo
	Ap(ic+1)=l1
c	Sort the Ai-values currently in X so they are in increasing order
	  N0=Ap(ic+1)-Ap(ic)
c	  write(6,*)'ic=',ic
c	  write(6,*)'entering SORT: X=',(X(lind), lind=1,N0)
	  call SORT(X,maxr,N0,Y)
c	  write(6,*)'out of SORT: Y=',(Y(lind), lind=1,N0)
	  do l=Ap(ic)+1,Ap(ic+1)
	  lind=l-Ap(ic)
	  Ai(l)=Y(lind)
	  enddo
	enddo

c	Assign an element of Ai to every G-matrix update #k
c	The jAi-index will be used to assign values to the AA array in MATREL
	do k=1,ktot
	ir=irow(k)
	ic=icol(k)
	do j=Ap(ic)+1,Ap(ic+1)
	if(Ai(j).eq.(ir-1)) jAi(k)=j
	enddo
c		write(6,*)'Ai element # associated with k=',k,'=',jAi(k)
	enddo
c	Write out Ap,Ai,jAi
	open(2,file='sporder.out',form='unformatted')
	write(2) Ap
	write(2) Ai
	write(2) jAi
	close(2)
	deallocate (irow)
	deallocate (icol)
	deallocate (Ai0)
	deallocate (X)
	deallocate (Y)

	return
	end

	subroutine globgrd(dlat)
c	Determine global grid numbers as a function of
c	(ncx,ncz), cell number indices
c	and local grid number il = N*(ilz-1) + ilx,
c	where (ilx,ilz) are local grid number indices
c
c	OUTPUT
c	igrd(ncx,ncz,il) array of global grid numbers
c
	parameter (ncmax=151)
	parameter (ncmaz=74)
c	Next line is the number of interior GLL points in 1D
	parameter (lmax=4)
c	N is the total number of GLL points in 1D including the points +-1
	parameter (N=lmax+2)
c	Next line is the maximum number of GLL points on the x-side of
c	the global grid times the maximum number on the z-side.
	parameter (nptmax=(N*ncmax-(ncmax-1))*(N*ncmaz-(ncmaz-1)))
	common/struc3/dx(ncmax),dz(ncmaz)
c	xg and zg have the (theta[radians],z[km]) coordinates at each of the global gridpoints
c	assuming dimensions of (dx)x(dz) for each cell 
c	The origin of the 2D grid is (theta=dlat,z=0)
c	(dimensionless x and z each run from -1 to +1).
	common/glbxy/xg(nptmax),zg(nptmax)
c
	common/grdpar/NX,NZ
	common/glbgrd/igrd(ncmax,ncmaz,N*N)
	real*8 dlat
	real*8 wt,x
	common/glpts/wt(lmax+2),x(lmax+2)
c*
	sumzn=0.
	do iz=1,NZ
	sumzn=sumzn+dz(iz)/2.
	enddo
c*
	ig=0
	do ncz=1,NZ
	do ncx=1,NX
	do ilz=1,N
	do ilx=1,N
	il=N*(ilz-1)+ilx
	inew=1
c
c	When ilx=1 and ncx>1, then reset ilx-->N and ncx-->ncx-1
	if(ilx.eq.1.and.ncx.gt.1) then
	il1=N
	nc1=ncx-1
	il0=N*(ilz-1)+il1
	igrd(ncx,ncz,il)=igrd(nc1,ncz,il0)
c		write(6,*)'ilx=1: igrd(',ncx,ncz,il,')=',igrd(ncx,ncz,il)
	inew=0
	endif
c
c	When ilz=1 and ncz>1, then reset ilz-->N and ncz-->ncz-1
	if(ilz.eq.1.and.ncz.gt.1) then
	il2=N
	nc2=ncz-1
	il0=N*(il2-1)+ilx
	igrd(ncx,ncz,il)=igrd(ncx,nc2,il0)
c		write(6,*)'ilz=1: igrd(',ncx,ncz,il,')=',igrd(ncx,ncz,il)
	inew=0
	endif
c
	if(inew.eq.1) then
	ig=ig+1
	igrd(ncx,ncz,il)=ig
		write(6,*)'igrd(',ncx,ncz,il,')=',igrd(ncx,ncz,il)
	sumx=real(dlat)/2.
	do ix=1,ncx
	sumx=sumx+dx(ix)/2.
	enddo
	xg(ig)=2.*sumx + (dx(ncx)/2.)*(real(x(ilx))-1.)
	sumz=0.
	do iz=1,ncz
	sumz=sumz+dz(iz)/2.
	enddo
	zg(ig)=2.*(sumz-sumzn) + (dz(ncz)/2.)*(real(x(ilz))-1.)
	endif
		write(6,*)'----------------'
c
	enddo
	enddo
	enddo
	enddo	

	return
	end

	subroutine init
	character*80 b80,sfile,tfile
	parameter (ncmax=151)
	parameter (ncmaz=74)
c	Next line is the number of interior GLL points in 1D
	parameter (lmax=4)
c	N is the total number of GLL points in 1D including the points +-1
	parameter (N=lmax+2)
c	Next line is the maximum number of GLL points on the x-side of
c	the global grid times the maximum number on the z-side.
	parameter (nptmax=(N*ncmax-(ncmax-1))*(N*ncmaz-(ncmaz-1)))
c***
	common/grdpar/NX,NZ
	common/struc3/dx(ncmax),dz(ncmaz)
	common/glbgrd/igrd(ncmax,ncmaz,N*N)
	real*8 elat,elon,plat,plon,phiref,pi,rad
	real*8 dlat,dlon,cdelt,sdelt,sphi,cphi,delta,phi
	common/sphpar/dlat,plat,plon,phiref
	common/minfo/minc,mmax
cseis
	common/sinfo/nleno,dt,corper
c--
	parameter (m0max=100)
	parameter (bigr=6371.)
c***
	pi=3.14159265358979d0
	rad=180.d0/pi
	write(6,*)'INIT: Reading in simulation info'
cseis
c	open(2,file='simulation-spherical.info')
	open(2,file='seis2pt5d-spherical.info')
c--
	rewind(2)
	read(2,5) b80
	read(2,*) elat,elon
	read(2,5) b80
	read(2,*) dlat,dlon
c*	Determine (colat,lon) of pole of the spherical coordinate system.
	elat=pi/2.d0-elat/rad
	elon=elon/rad
	dlat=dlat/rad
	dlon=dlon/rad
        cdelt=dcos(elat)*dcos(dlat)+dsin(elat)*dsin(dlat)*dcos(dlon)
		write(6,*)'elat=',elat,'dlon=',dlon
		write(6,*)'cdelt=',cdelt
        delta=dacos(cdelt)
		write(6,*)'delta=',delta
	sdelt=dsin(delta)
        sphi=dsin(dlat)*dsin(dlon)/sdelt
        cphi=(dcos(dlat)-dcos(elat)*cdelt)/(dsin(elat)*sdelt)
	phi=datan2(sphi,cphi)
	plat=delta
	plon=elon+phi
		write(6,*)'INIT: elon,phi=',elon,phi
		write(6,*)'INIT: geographic plat,plon=',90.-rad*plat,rad*plon
cTE
c	Determine azimuth from the pole of the spherical coordinate system
c	to the reference point (theta=dlat,z=0) of the 2D grid in rad CCW from due N.
        sphi=dsin(elat)*dsin(dlon)/sdelt
        cphi=(dcos(elat)-dcos(dlat)*cdelt)/(dsin(dlat)*sdelt)
	phiref=datan2(sphi,cphi)
		write(6,*)'INIT: phi_ref=',rad*phiref
c*	
	read(2,5) b80
	read(2,*) NX
	read(2,5) b80
	read(2,*) (dx(i), i=1,NX)
	read(2,5) b80
	read(2,*) NZ
	read(2,5) b80
	read(2,*) (dz(i), i=1,NZ)
	read(2,5) b80
	read(2,*) ydist,wavmin
c	Structure and deformation is periodic with spatial periodicity ydist km, hence
	ydist=ydist/real(dsin(dlat))
	wavmin=wavmin/real(dsin(dlat))
	minc=int((2.*3.1415926535*bigr)/ydist)
	mmax=int((2.*3.1415926535*bigr)/wavmin)
	read(2,5) b80
	read(2,5) sfile
	read(2,5) b80
	read(2,5) tfile
cseis
	read(2,5) b90
	read(2,*) dt
	read(2,5) b90
	read(2,*) corper
	read(2,5) b90
	read(2,*) nleno
	nleno=nleno*2
c--
	  write(6,*)'Number of cells in theta-direction=',NX
	  write(6,*)'Number of cells in z-direction=',NZ
	  write(6,*)'Physical length of cell along theta-direction=',dx
	  write(6,*)'Physical length of cell along z-direction=',dz
	  write(6,*)'Structural parameter file is'
	  write(6,5) sfile
	  write(6,*)'Topography file is'
	  write(6,5) tfile
	close(2)
	write(6,*)'INIT: Finished reading in simulation info'
	write(6,*)'INIT: calling basis'
	call basis
	write(6,*)'INIT: out of basis'
	write(6,*)'INIT: calling globgrid'
	call globgrd(dlat)
	write(6,*)'INIT: out of globgrid'
	write(6,*)'INIT: calling gridgeom'
	call gridgeom(sfile)
	write(6,*)'INIT: out of gridgeom'
	write(6,*)'INIT: calling topogeom'
	call topogeom(tfile)
	write(6,*)'INIT: out of topogeom'
5	format(a80)
	return
	end

	subroutine matrel(w0,ky,nmat,Ap,Ai,jAi,AA0,AA1,AA2)
cccc	subroutine matrel(ome,ky,KU,KL)
c	Generate elements of G (and d) for given 
c	Laplace-transform parameter s
cccc	angular frequency ome
c	and y-wavenumber ky. KU and KL are the number of diagonals
c	above and below the central diagonal, respectively, with non-zero 
c	matrix elements.
c		real*8 xtmp
	character*80 dum80
	complex*16 w0,ui
cc	complex*16 w0p
	real*8 w0r,w0i,wmag,rterm,sterm
	real*8 ky,kyone
cccc	real*8 ome,ky
	parameter (ncmax=151)
	parameter (ncmaz=74)
c	Next line is the number of interior GLL points in 1D
	parameter (lmax=4)
c	N is the total number of GLL points in 1D including the points +-1
	parameter (N=lmax+2)
	parameter (N2=N*N)
c	Next line is the maximum number of GLL points on the x-side of
c	the global grid times the maximum number on the z-side.
	parameter (nptmax=(N*ncmax-(ncmax-1))*(N*ncmaz-(ncmaz-1)))
c	There are a maximum of nptmax GLL points in the 2D global grid,
c	each with three displacement components.
	common/glbxy/xg(nptmax),zg(nptmax)
c	igc stored in common will have the grid # of the last point
c	determined on the source plane -- representative of the source location.
	common/isrc/igc
c****
c	Every column of G will be expressed in a sparse matrix format
c	3*(N*OLDnptmax+N*N) is a slightly high estimate of KL and/or KU,
c	so icompr below is slightly higher than 2*KL+KU+1,
c	the intended leading dimension of G.
cOLD	parameter (icompr=3*(N*nptmax+N*N) * 3 + 1)
	complex*16 g
cOLD	dimension g(icompr,3*nptmax)
	complex*16 d
	common/datavec/d(3*nptmax,4)
cOLD	integer KL,KU,M,NRHS,LDAB,INFO
cOLD	integer IPIV(3*nptmax)
c****
cTOPO	The Jacobian J accounts for the distortion of elements due to topography
	real*8 jacm1
	common/struc5/jacm1(ncmax,ncmaz,N*N)
c****
	real*8 fac,fac1,fac2,cott
	real*8 wfac
	complex*16 dfac,dfacz,odfac,dprime
	parameter (bigr=6371.)
c***
	common/grdpar/NX,NZ
	real*8 vpval,vsval,rhoval,qbval,qkval
	real*8 vp,vs,rho,qb,qk
	real*8 drat,dratz
	complex*16 dxeff,dzeff
	common/struc3/dx(ncmax),dz(ncmaz)
	common/struc4/vp(ncmax,ncmaz,N*N),vs(ncmax,ncmaz,N*N),
     &	rho(ncmax,ncmaz,N*N),qb(ncmax,ncmaz,N*N),qk(ncmax,ncmaz,N*N)
	common/glbgrd/igrd(ncmax,ncmaz,N*N)
c*****
	parameter (kmax=ncmax*ncmaz*N2*N2*9)
	integer*8 nzmax,nmax
cUSED        parameter (nzmax = 32800000, nmax = 215000)
        parameter (nzmax = kmax, nmax = 3*(N*ncmax-(ncmax-1))*(N*ncmaz-(ncmaz-1))+1)
cUSED	integer*8 Ap,Ai,jAi
	integer*8 Ap
cUSED	dimension Ap(nmax),Ai(nzmax),jAi(kmax)
	dimension Ap(nmax)
	integer*8 Ai(nmat),jAi(nmat)
	complex*16 AA0(nmat),AA1(nmat),AA2(nmat)
cUSED	complex*16 AA (nzmax) , XX (nmax), BB (nmax), r (nmax)
	complex*16, ALLOCATABLE :: AA(:)
	double precision, ALLOCATABLE :: Ax(:)
	double precision, ALLOCATABLE :: Az(:)
	complex*16 XX (nmax), BB (nmax), r (nmax)
cUSED	double precision Ax (nzmax), xumf (nmax), b (nmax)
	double precision xumf (nmax), b (nmax)
	double precision control(20),info(90)
cUSED	double precision Az (nzmax), xz (nmax), bz (nmax)
	double precision xz (nmax), bz (nmax)
	integer*8 numf , numfz , numeric, symbolic, status, sys, filenum
c*****
	real*8 kappa0,mu0
	complex*16 lam2mu,lam,mu,kap
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
	allocate ( AA(nmat) )
	allocate ( Ax(nmat) )
	allocate ( Az(nmat) )

	write(6,*)'BB size of Ai=',size(Ai)
	write(6,*)'BB size of jAi=',size(jAi)

	write(6,*)'BB Ai(20)=',Ai(20)
	write(6,*)'BB jAi(20)=',jAi(20)
	pi=3.14159265
	twopi=2.*pi
	ui=dcmplx(0.d0,1.d0)
c	Assume a reference frequency of twopi rad/sec (i.e. 1 Hz).
	w0r=dreal(w0)/dble(twopi)
	w0i=dimag(w0)/dble(twopi)
	wmag=dsqrt(w0r*w0r+w0i*w0i)
	rterm=dlog(wmag)
	sterm=datan2(w0i,w0r)
		write(6,*)'MATREL: w0,w0r,w0i=',w0,w0r,w0i
		write(6,*)'MATREL: wmag=',wmag,'rterm,sterm=',rterm,sterm

c	Zero out the elements of AA
	ncol=3*(N*NX-(NX-1))*(N*NZ-(NZ-1))
		write(6,*)'MATREL: ncol=',ncol
	if(ky.lt.0.5d0) then
c	ky is 0 if arriving at this line
	write(6,*)'Zero-ing out elements of AA0, AA1, and AA2'
	write(6,*)'Ap(',ncol+1,')=',Ap(ncol+1)
	write(6,*)'nmat=',nmat
	do j=1,Ap(ncol+1)
	AA0(j)=dcmplx(0.,0.)
	AA1(j)=dcmplx(0.,0.)
	AA2(j)=dcmplx(0.,0.)
	enddo
	write(6,*)'Done zero-ing out elements of AA0, AA1, and AA2'
c	kyone is a placeholder for unit ky soo we can keep track of the
c	ky-terms that go into the AA1 and AA2 arrays
	kyone=1.d0
	else
	go to 35
	endif

	write(6,*)'Constructing G'
	k0=0
	do ncz=1,NZ
	do ncx=1,NX
	do ilz=1,N
	do ilx=1,N
	il=N*(ilz-1)+ilx
	i1=3*igrd(ncx,ncz,il)-2
	i2=i1+1
	i3=i2+1
c		write(6,*)'Working on rows',i1,i2,i3
c		write(6,*)'out of',3*(N*NX-(NX-1))*(N*NZ-(NZ-1))
c		write(6,*)'      '
c	i1,i2,i3 are the row numbers corresponding to the normal equation
c	resulting from multiplication of the momentum eqn with phi_[il]
c	in cell (ncx,ncz) and integration over that cell.
c	i1,i2,i3 refer to the x,y,z components.
c	[il] is equivalent to index [j] of notes.
	jl=0
10	jl=jl+1
	if(jl.gt.N2) go to 20
	jlz=(jl-1)/N+1
	jlx=jl-N*(jlz-1)
c	jl is the index for the expansion coefficients of displacement
c	[jl] is equivalent to index [k] of notes.
	j1=3*igrd(ncx,ncz,jl)-2
	j2=j1+1
	j3=j2+1
c	j1,j2,j3 are the column numbers, corresponding to the
c	a_jl, b-jl, c_jl expansion coeff of displacement in cell (ncx,ncz)
cNOTE	If the test function corresponds to
c	the bottom of the grid, then enforce zero-displacement BC
c	This assumes that the i1,i2,i3 elements of the data vector d are zero,
c	i.e. the source is not in one of the bottom cells.
c	if(ilz.eq.1.and.ncz.eq.1) then
c	write(6,*)'ilz=',ilz,'ncz=',ncz
c	write(6,*)'ncx,ncz=',ncx,ncz
c	write(6,*)'will go to 25'
c	endif
	if(ilz.eq.1.and.ncz.eq.1) go to 25
cNEW
	if(ilx.eq.1.and.ncx.eq.1) go to 30
	if(ilx.eq.N.and.ncx.eq.NX) go to 30
c--
	l=0
15	l=l+1
	if(l.gt.N2) then
	k0=k
c		write(6,*)'Check k0: go to 10 w/k0=',k
	go to 10
	endif
	lz=(l-1)/N+1
	lx=l-N*(lz-1)
		if(ilx.ne.lx.and.ilz.ne.lz) go to 15
c**
c	write(6,*)'ilx,lx=(',ilx,lx,')','ilz,lz=(',ilz,lz,')'
c	write(6,*)'dphi2x(',il,l,')=',dphi2x(il,l)
c	write(6,*)'dphi2z(',il,l,')=',dphi2z(il,l)
c	write(6,*)'-----------------------------'

	ig=igrd(ncx,ncz,l)
	rhoval=rho(ncx,ncz,l)
	vsval=vs(ncx,ncz,l)
	vpval=vp(ncx,ncz,l)
	qbval=qb(ncx,ncz,l)
	qkval=qk(ncx,ncz,l)

	mu0=0.1*rhoval*vsval**2
	kappa0=0.1*rhoval*vpval**2 - (4.d0/3.d0)*mu0

c	Assume a reference frequency of twopi rad/sec (i.e. 1 Hz).
c	This is done when computing the real and imaginary parts of log(w0/twopi)
c	in the main program, which are in rterm and sterm, respectively.
	mu=mu0 * 
     &	dcmplx(1.d0+(2.d0/pi)*dble(qbval)*rterm,dble(qbval)*(1.d0+(2.d0/pi)*sterm))
	kap=kappa0 * 
     &	dcmplx(1.d0+(2.d0/pi)*dble(qkval)*rterm,dble(qkval)*(1.d0+(2.d0/pi)*sterm))
	lam=kap-(2.d0/3.d0)*mu
	lam2mu=lam+2.d0*mu

cc	Apply correspondence principle for Maxwell solid.
c	mu=mu0*s/(s+mu0/etaval)
c	lam=kappa0-(2.d0/3.d0)*mu
c	lam2mu=lam+2.d0*mu

c	Eqn 7 of notes
c	darea=dble((dx/2.))*dble((dz/2.)) * (bigr+zg(ig))
	fac=bigr+zg(ig)
	fac1=1.d0/(bigr+zg(ig))
	fac2=fac1 * dble(1.0/sin(xg(ig)))
c		write(6,*)'fac1,fac2=',fac1,fac2,'xg(',ig,')=',xg(ig)
	cott=dble(cos(xg(ig))/sin(xg(ig)))
c		write(6,*)'cott=',cott
c*-*-*-*Absorbing layer
c	Define `effective' dx-value, which stands for (twice) d(theta)/dx, under PML coordinate stretching 
c	dxeff=dble(dx(ncx))
c	dzeff=dble(dz(ncz))
	drat=0.d0
	dratz=0.d0
	dfac=1.d0
	dfacz=1.d0
cOLD	odfac=1.d0
cTOPO
	odfac=jacm1(ncx,ncz,l)
c--
cc	dprime=0.d0
cc	w0p=w0
c	Absorbing layer at right hand and left hand sides of domain.
c	if(ncx.eq.(NX-4)) drat=0.2d0*(x(lx)-x(1))/(x(N)-x(1)) 
c	if(ncx.eq.(NX-3)) drat=0.2d0+0.2d0*(x(lx)-x(1))/(x(N)-x(1)) 
c	if(ncx.eq.(NX-2)) drat=0.4d0+0.2d0*(x(lx)-x(1))/(x(N)-x(1)) 
c	if(ncx.eq.(NX-1)) drat=0.6d0+0.2d0*(x(lx)-x(1))/(x(N)-x(1)) 
c	if(ncx.eq.NX) drat=0.8d0+0.2d0*(x(lx)-x(1))/(x(N)-x(1)) 
	if(ncx.eq.(NX-2)) drat=0.3333333d0*(x(lx)-x(1))/(x(N)-x(1)) 
	if(ncx.eq.(NX-1)) drat=0.333333d0+0.333333d0*(x(lx)-x(1))/(x(N)-x(1)) 
	if(ncx.eq.NX) drat=0.666667d0+0.333333d0*(x(lx)-x(1))/(x(N)-x(1)) 
c	if(ncx.eq.5) drat=0.2d0*(x(lx)-x(N))/(x(N)-x(1)) 
c	if(ncx.eq.4) drat=-0.2d0+0.2d0*(x(lx)-x(N))/(x(N)-x(1)) 
c	if(ncx.eq.3) drat=-0.4d0+0.2d0*(x(lx)-x(N))/(x(N)-x(1)) 
c	if(ncx.eq.2) drat=-0.6d0+0.2d0*(x(lx)-x(N))/(x(N)-x(1)) 
c	if(ncx.eq.1) drat=-0.8d0+0.2d0*(x(lx)-x(N))/(x(N)-x(1)) 
	if(ncx.eq.3) drat=0.333333d0*(x(lx)-x(N))/(x(N)-x(1)) 
	if(ncx.eq.2) drat=-0.333333d0+0.333333d0*(x(lx)-x(N))/(x(N)-x(1)) 
	if(ncx.eq.1) drat=-0.666667d0+0.333333d0*(x(lx)-x(N))/(x(N)-x(1)) 
c	Absorbing layer at bottom of domain.
	if(ncz.eq.3) dratz=0.333333d0*(x(lz)-x(N))/(x(N)-x(1)) 
	if(ncz.eq.2) dratz=-0.333333d0+0.333333d0*(x(lz)-x(N))/(x(N)-x(1)) 
	if(ncz.eq.1) dratz=-0.666667d0+0.333333d0*(x(lz)-x(N))/(x(N)-x(1)) 
	dfac=1.d0 - (ui/w0)*1.5d0*drat**2*(vpval/(3.d0*dx(ncx)*fac))*5.00
	dfacz=1.d0 - (ui/w0)*1.5d0*dratz**2*(vpval/(3.d0*dz(ncz)))*5.00
cNOTE	dprime == (d/dn)d(n)
	dxeff=dble(dx(ncx))*dfac
	dzeff=dble(dz(ncz))*dfacz
c		if(dabs(dimag(dzeff)).gt.1.d-09) then
c		write(6,*)'w0=',w0
c		write(6,*)'ncz=',ncz,'dfacz=',dfacz
c		write(6,*)'dzeff=',dzeff
c		endif
cc	dprime= 3.d0*drat*(0.2d0/(x(N)-x(1)))*(vpval/(3.d0*dx(ncx)*fac))*5.00
cc	odfac=1.d0/dfac
c		write(6,*)'ncx=',ncx,'lx=',lx,'drat=',drat
c		write(6,*)'dprime=',dprime,'odfac=',odfac
c		write(6,*)'- - - - - - - - - - - '
c	As (d/dx)d(n) 3 lines above is with respect to basis-function coordinates, it should be divided
c	by another factor of (dx(ncx)/2 * fac).  Another factor of (1/dfac) is omitted
c	from the dprime term of wfac below, so we are omitting a net factor of
c	((dx(ncx)/2)*fac*dfac)^(-1) in the dprime term of wfac.  But this term is implicitly
c	included in the AA(j) terms where wfac and other factors are multiplied out.   
c
cc	mu=mu/dfac**2
c	(1.d0 - (1.5d0/w0)*ui*max(drat,dratz)**2*(vpval/(3.d0*dx(ncx)*fac))*3.00)**2
cc	kap=kap/dfac**2
c	(1.d0 - (1.5d0/w0)*ui*max(drat,dratz)**2*(vpval/(3.d0*dx(ncx)*fac))*3.00)**2
cc	lam=kap-(2.d0/3.d0)*mu
cc	lam2mu=lam+2.d0*mu
cc	w0p=w0 - ui*max(drat**2,dratz**2)*(vpval/(3.d0*dx(ncx)*fac))*4.00
c	write(6,*)'dx(',ncx,')=',dx(ncx),'dxeff=',dxeff
c	write(6,*)'x(',lx,')=',x(lx),'quadratic factor=',((x(lx)-x(1))/(x(N)-x(1)))**2
c	write(6,*)'velocity=',vpval,'km s^-1.  Ratio=',(vpval/(dx(ncx)*fac))
c	write(6,*)'w0=',w0
c		write(6,*)'wfac 1=',dphi2x(il,l),'wfac 2=',phi2(il,l)*(ui/w0)*dprime
c		write(6,*)'phi2(',il,l,')=',phi2(il,l),'dprime=',dprime
c	write(6,*)'- - - - - - - - - - - - -'
cc	endif
c*-*-*-*
c	Terms arising from the theta- and r-parts of the gradient of the stress tensor
c	after integration by parts in theta and r.
	wfac = dphi2x(il,l)
cc	wfac = dphi2x(il,l) + phi2(il,l)*(ui/w0)*odfac*dprime
c		write(6,*)'Term1 dphi2x(',il,l,')=',dphi2x(il,l)
c		write(6,*)'Term 2=',phi2(il,l)*(ui/w0)*odfac*dprime
c		write(6,*)'ilx,lx=',ilx,lx
c		write(6,*)'ilz,lz=',ilz,lz
c		write(6,*)'- - - - - - - - - -'
c c	wfac = dphi2x(il,l) + phi2(il,l)*(ui/w0)*dprime
c		write(6,*)'wfac=',wfac
		if(dabs(wfac).gt.1.d-09) then
c	i1,j1 terms
	k=k0+1
	j=jAi(k)
	AA0(j)=AA0(j)
     &	- wt2(l)*wfac*lam2mu*dphi2x(jl,l) * dzeff/dxeff * fac1 * odfac
     &					    - wt2(l)*wfac*lam*phi2(jl,l)*cott * dzeff/2.d0 * fac1 * odfac

c	i1,j2 terms
	k=k+1
	j=jAi(k)
	AA1(j)=AA1(j)
     &	- wt2(l)*wfac*lam*ui*kyone*phi2(jl,l) * dzeff/2.d0 * fac2 * odfac

c	i1,j3 terms
	k=k+1
	j=jAi(k)
	AA0(j)=AA0(j)
     &	- wt2(l)*wfac*lam2mu*phi2(jl,l) * dzeff/2.d0 * fac1 * odfac
     &					    - wt2(l)*wfac*lam*dphi2z(jl,l) * 1.d0 * odfac
     &					    - wt2(l)*wfac*lam*phi2(jl,l) * dzeff/2.d0 * fac1 * odfac

c	i2,j1 terms
	k=k+1
	j=jAi(k)
	AA1(j)=AA1(j)
     &	- wt2(l)*wfac*mu*ui*kyone*phi2(jl,l) * dzeff/2.d0 * fac2 * odfac

c	i2,j2 terms
	k=k+1
	j=jAi(k)
	AA0(j)=AA0(j)
     &	- wt2(l)*wfac*mu*dphi2x(jl,l) * dzeff/dxeff * fac1 * odfac
     &					    + wt2(l)*wfac*mu*phi2(jl,l)*cott * dzeff/2.d0 * fac1 * odfac

c	i2,j3 terms
	k=k+1
	j=jAi(k)
c	AA(j)=AA(j)

c	i3,j1 terms
	k=k+1
	j=jAi(k)
	AA0(j)=AA0(j)
     &	- wt2(l)*wfac*mu*dphi2z(jl,l) * 1.d0 * odfac
     &					    + wt2(l)*wfac*mu*phi2(jl,l) * dzeff/2.d0 * fac1 * odfac

c	i3,j2 terms
	k=k+1
	j=jAi(k)
c	AA(j)=AA(j)

c	i3,j3 terms
	k=k+1
	j=jAi(k)
	AA0(j)=AA0(j)
     &	- wt2(l)*wfac*mu*dphi2x(jl,l) * dzeff/dxeff * fac1 * odfac
		endif

		if(dabs(dphi2z(il,l)).gt.1.d-09) then
c	i1,j1 terms
	k=k0+1
	j=jAi(k)
	AA0(j)=AA0(j)
     &		                            - wt2(l)*dphi2z(il,l)*mu*dphi2z(jl,l) * dxeff/dzeff * fac * odfac
     &		                            + wt2(l)*dphi2z(il,l)*mu*phi2(jl,l) * dxeff/2.0 * odfac

c	i1,j2 terms
	k=k+1
	j=jAi(k)
c	AA(j)=AA(j)

c	i1,j3 terms
	k=k+1
	j=jAi(k)
	AA0(j)=AA0(j)
     &		                            - wt2(l)*dphi2z(il,l)*mu*dphi2x(jl,l) * 1.d0 * odfac

c	i2,j1 terms
	k=k+1
	j=jAi(k)
c	AA(j)=AA(j)

c	i2,j2 terms
	k=k+1
	j=jAi(k)
	AA0(j)=AA0(j)
     &		                            - wt2(l)*dphi2z(il,l)*mu*dphi2z(jl,l) * dxeff/dzeff * fac * odfac
     &					    + wt2(l)*dphi2z(il,l)*mu*phi2(jl,l) * dxeff/2.d0 * odfac

c	i2,j3 terms
	k=k+1
	j=jAi(k)
	AA1(j)=AA1(j)
     &	- wt2(l)*dphi2z(il,l)*mu*ui*kyone*phi2(jl,l) * dxeff/2.d0 * fac2 * fac * odfac

c	i3,j1 terms
	k=k+1
	j=jAi(k)
	AA0(j)=AA0(j)
     &		                            - wt2(l)*dphi2z(il,l)*lam*dphi2x(jl,l) * 1.d0 * odfac
     &					    - wt2(l)*dphi2z(il,l)*lam*phi2(jl,l)*cott * dxeff/2.d0 * odfac

c	i3,j2 terms
	k=k+1
	j=jAi(k)
	AA1(j)=AA1(j)
     &	- wt2(l)*dphi2z(il,l)*lam*ui*kyone*phi2(jl,l) * dxeff/2.d0 * fac2 * fac * odfac

c	i3,j3 terms
	k=k+1
	j=jAi(k)
	AA0(j)=AA0(j)
     &		                            - wt2(l)*dphi2z(il,l)*lam2mu*dphi2z(jl,l) * dxeff/dzeff * fac * odfac
     &		                            - wt2(l)*dphi2z(il,l)*2.d0*lam*phi2(jl,l) * dxeff/2.d0 * odfac
		endif

c	Terms arising from the phi-part of the gradient of the stress tensor.
c	i1,j1 terms
		if(il.eq.l) then
	k=k0+1
	j=jAi(k)
	AA2(j)=AA2(j)
     &	+ wt2(l)*phi2(il,l)*fac2*ui*kyone *mu*ui*kyone*phi2(jl,l) 
     &	* dzeff*dxeff/4.d0 * fac2 * fac * odfac
	AA0(j)=AA0(j)
     &	- wt2(l)*phi2(il,l)*fac1*cott *lam2mu*phi2(jl,l)*cott * dzeff*dxeff/4.d0 * odfac
     &	- wt2(l)*phi2(il,l)*fac1*cott *lam*dphi2x(jl,l) * dzeff/2.d0 * odfac
     &	+ wt2(l)*phi2(il,l)*mu*dphi2z(jl,l) * dxeff/2.d0 * odfac
     &	- wt2(l)*phi2(il,l)*mu*phi2(jl,l) * dzeff*dxeff/4.d0 * fac1 * odfac
     &	+ 0.1d0*rhoval*w0**2*wt2(l)*phi2(il,l)*phi2(jl,l) * dzeff*dxeff/4.d0 * fac * odfac
c	Last term is inertial term
c	Above previous 2 lines: Additional terms arising from r- and theta-derivatives of unit vectors.

c	i1,j2 terms
	k=k+1
	j=jAi(k)
	AA1(j)=AA1(j)
     &	+ wt2(l)*phi2(il,l)*fac2*ui*kyone *mu*dphi2x(jl,l) * dzeff/2.d0 * odfac
     &	- wt2(l)*phi2(il,l)*fac2*ui*kyone *mu*phi2(jl,l)*cott * dzeff*dxeff/4.d0 * odfac
     &	- wt2(l)*phi2(il,l)*fac1*cott *lam2mu*ui*kyone*phi2(jl,l) 
     &	* dzeff*dxeff/4.d0 * fac2 * fac * odfac

c	i1,j3 terms
	k=k+1
	j=jAi(k)
	AA0(j)=AA0(j)
     &	- wt2(l)*phi2(il,l)*fac1*cott *lam2mu*phi2(jl,l) * dzeff*dxeff/4.d0 * odfac
     &	- wt2(l)*phi2(il,l)*fac1*cott *lam*phi2(jl,l) * dzeff*dxeff/4.d0 * odfac
     &	- wt2(l)*phi2(il,l)*fac1*cott *lam*dphi2z(jl,l) * dxeff/2.d0 * fac  * odfac
     &	+ wt2(l)*phi2(il,l)*mu*dphi2x(jl,l) * dzeff/2.d0 * fac1 * odfac
c	Above 1 line: Additional terms arising from r- and theta-derivatives of unit vectors.

c	i2,j1 terms
	k=k+1
	j=jAi(k)
	AA1(j)=AA1(j)
     &	+ wt2(l)*phi2(il,l)*fac2*ui*kyone *lam2mu*phi2(jl,l)*cott * dzeff*dxeff/4.d0 * odfac
     &	+ wt2(l)*phi2(il,l)*fac2*ui*kyone *lam*dphi2x(jl,l) * dzeff/2.d0 * odfac
     &	+ wt2(l)*phi2(il,l)*fac1*cott *mu*ui*kyone*phi2(jl,l) 
     &	* dzeff*dxeff/4.d0 * fac2 * fac * odfac

c	i2,j2 terms
	k=k+1
	j=jAi(k)
	AA2(j)=AA2(j)
     &	+ wt2(l)*phi2(il,l)*fac2*ui*kyone *lam2mu*ui*kyone*phi2(jl,l) 
     &	* dzeff*dxeff/4.d0 * fac2 * fac * odfac
	AA0(j)=AA0(j)
     &	+ wt2(l)*phi2(il,l)*fac1*cott *mu*dphi2x(jl,l) * dzeff/2.d0 * odfac
     &	- wt2(l)*phi2(il,l)*fac1*cott *mu*phi2(jl,l)*cott * dzeff*dxeff/4.d0 * odfac
     &	+ wt2(l)*phi2(il,l)*fac1 *mu*dphi2z(jl,l) * dxeff/2.d0 * fac * odfac
     &	- wt2(l)*phi2(il,l)*fac1 *mu*phi2(jl,l) * dzeff*dxeff/4.d0 * odfac
     &	+ 0.1d0*rhoval*w0**2*wt2(l)*phi2(il,l)*phi2(jl,l) * dzeff*dxeff/4.d0 * fac * odfac
c	Last term is inertial term

c	i2,j3 terms
	k=k+1
	j=jAi(k)
	AA1(j)=AA1(j)
     &	+ wt2(l)*phi2(il,l)*fac2*ui*kyone *lam2mu*phi2(jl,l) * dzeff*dxeff/4.d0 * odfac
     &	+ wt2(l)*phi2(il,l)*fac2*ui*kyone *lam*phi2(jl,l) * dzeff*dxeff/4.d0 * odfac
     &	+ wt2(l)*phi2(il,l)*fac2*ui*kyone *lam*dphi2z(jl,l) * dxeff/2.d0 * fac * odfac
     &	+ wt2(l)*phi2(il,l)*fac1 *mu*ui*kyone*phi2(jl,l) 
     &	* dzeff*dxeff/4.d0 * fac2 * fac * odfac

c	i3,j1 terms
	k=k+1
	j=jAi(k)
	AA0(j)=AA0(j)
     &	- wt2(l)*phi2(il,l)*fac1 *lam2mu*phi2(jl,l)*cott * dzeff*dxeff/4.d0 * odfac
     &	- wt2(l)*phi2(il,l)*fac1 *lam*dphi2x(jl,l) * dzeff/2.d0 * odfac
     &	- wt2(l)*phi2(il,l)*lam2mu*dphi2x(jl,l) * dzeff/2.d0 * fac1 * odfac
     &	- wt2(l)*phi2(il,l)*lam*phi2(jl,l)*cott * dzeff*dxeff/4.d0 * fac1 * odfac
c	Above 2 lines: Additional terms arising from r- and theta-derivatives of unit vectors.

c	i3,j2 terms
	k=k+1
	j=jAi(k)
	AA1(j)=AA1(j)
     &	+ wt2(l)*phi2(il,l)*fac2*ui*kyone *mu*dphi2z(jl,l) * dxeff/2.d0 * fac  * odfac
     &	- wt2(l)*phi2(il,l)*fac2*ui*kyone *mu*phi2(jl,l) * dzeff*dxeff/4.d0 * odfac
     &	- wt2(l)*phi2(il,l)*fac1 *lam2mu*ui*kyone*phi2(jl,l) 
     &	* dzeff*dxeff/4.d0 * fac2 * fac * odfac
     &	- wt2(l)*phi2(il,l)*lam*ui*kyone*phi2(jl,l) * dzeff*dxeff/4.d0 * fac2 * odfac
c	Above 1 line: Additional terms arising from r- and theta-derivatives of unit vectors.

c	i3,j3 terms
	k=k+1
	j=jAi(k)
	AA2(j)=AA2(j)
     &	+ wt2(l)*phi2(il,l)*fac2*ui*kyone *mu*ui*kyone*phi2(jl,l) 
     &	* dzeff*dxeff/4.d0 * fac2 * fac * odfac
	AA0(j)=AA0(j)
     &	- wt2(l)*phi2(il,l)*fac1 *lam2mu*phi2(jl,l) * dzeff*dxeff/4.d0 * odfac
     &	- wt2(l)*phi2(il,l)*fac1 *lam*phi2(jl,l) * dzeff*dxeff/4.d0 * odfac
     &	- wt2(l)*phi2(il,l)*fac1 *lam*dphi2z(jl,l) * dxeff/2.d0 * fac * odfac
     &	- wt2(l)*phi2(il,l)*lam2mu*phi2(jl,l) * dzeff*dxeff/4.d0 * fac1 * odfac
     &	- wt2(l)*phi2(il,l)*lam*dphi2z(jl,l) * dxeff/2.d0 * odfac
     &	- wt2(l)*phi2(il,l)*lam*phi2(jl,l) * dzeff*dxeff/4.d0 * fac1 * odfac
     &	+ 0.1d0*rhoval*w0**2*wt2(l)*phi2(il,l)*phi2(jl,l) * dzeff*dxeff/4.d0 * fac * odfac
c	Last term is inertial term
c	Above previous 3 lines: Additional terms arising from r- and theta-derivatives of unit vectors.
		endif

	go to 15
cNOTE	Enforce zero-displacement BC at the bottom of the grid.
25	continue
c	NOTE: if i1=j1, then it must also be true that i2=j2 and i3=j3
	if(i1.eq.j1) then
	k=k0+1
	j=jAi(k)
	AA0(j)=1.d0
	k=k+1
	j=jAi(k)
	AA0(j)=1.d0
	k=k+1
	j=jAi(k)
	AA0(j)=1.d0
	endif
c	if(i1.eq.j1) g(KU+KL+1+i1-j1,j1)=1.d0
c	if(i2.eq.j2) g(KU+KL+1+i2-j2,j2)=1.d0
c	if(i3.eq.j3) g(KU+KL+1+i3-j3,j3)=1.d0
c		write(6,*)'g(',KU+KL+1+i1-j1,j1,')=',g(KU+KL+1+i1-j1,j1)
c		write(6,*)'g(',KU+KL+1+i2-j2,j2,')=',g(KU+KL+1+i2-j2,j2)
c		write(6,*)'g(',KU+KL+1+i3-j3,j3,')=',g(KU+KL+1+i3-j3,j3)
	k0=k
	go to 10
30	continue
	igb=igrd(ncx,ncz,il)
cTE
c	roff=1.d0/(xg(igb)-xg(igc))
	roff=1.d0
c	Note the 0's below as well
c--
	k=k0+1
	j=jAi(k)
	AA0(j)=phi2(jl,il)*roff+dphi2x(jl,il)*(2.d0/dx(ncx))*0.
	k=k+1
	j=jAi(k)
	AA0(j)=phi2(jl,il)*roff+dphi2x(jl,il)*(2.d0/dx(ncx))*0.
	k=k+1
	j=jAi(k)
	AA0(j)=phi2(jl,il)*roff+dphi2x(jl,il)*(2.d0/dx(ncx))*0.
	k0=k
	go to 10
20	continue

	enddo
	enddo
	enddo
	enddo
35	continue
	write(6,*)'Assigning elements of AA'
	write(6,*)'Ap(',ncol+1,')=',Ap(ncol+1)
	do j=1,Ap(ncol+1)
	AA(j)=AA0(j) + ky*AA1(j) + ky*ky*AA2(j)
	enddo
	write(6,*)'Done assiging elements of AA'

40	continue

	write(6,*)'Done constructing G'
		write(6,*)'MATREL: end value of k=',k0
	numfz=Ap(ncol+1)
	do j=1,numfz
	Ax (j) = dreal(AA(j))
c		write(6,*)'Ax(',j,')=',Ax(j)
        Az (j) = dimag(AA(j))
c		write(6,*)'Az(',j,')=',Az(j)
	enddo

	numf=3*(N*NX-(NX-1))*(N*NZ-(NZ-1))

c       ----------------------------------------------------------------
c       factor the matrix and save to a file
c       ----------------------------------------------------------------

c       set default parameters
	write(6,*)'entering umf4zdef'
        call umf4zdef (control)

c       print control parameters.  set control (1) to 1 to print
c       error messages only
        control (1) = 2
	write(6,*)'entering umf4zpcon'
        call umf4zpcon (control)

c       pre-order and symbolic analysis
	write(6,*)'entering umf4zsym'
        call umf4zsym (numf, numf, Ap, Ai, Ax, Az, symbolic, control, info)

c       print statistics computed so far
c       call umf4zpinf (control, info) could also be done.
        print 80, info (1), info (16),
     $      (info (21) * info (4)) / 2**20,
     $      (info (22) * info (4)) / 2**20,
     $      info (23), info (24), info (25)
80      format ('symbolic analysis:',/,
     $      '   status:  ', f5.0, /,
     $      '   time:    ', e10.2, ' (sec)'/,
     $      '   estimates (upper bound) for numeric LU:', /,
     $      '   size of LU:    ', f10.2, ' (MB)', /,
     $      '   memory needed: ', f10.2, ' (MB)', /,
     $      '   flop count:    ', e10.2, /
     $      '   nnz (L):       ', f10.0, /
     $      '   nnz (U):       ', f10.0)

c       check umf4zsym error condition
        if (info (1) .lt. 0) then
            print *, 'Error occurred in umf4zsym: ', info (1)
            stop
        endif

c       numeric factorization
        call umf4znum (Ap, Ai, Ax, Az, symbolic, numeric, control, info)

c       print statistics for the numeric factorization
c       call umf4zpinf (control, info) could also be done.
        print 90, info (1), info (66),
     $      (info (41) * info (4)) / 2**20,
     $      (info (42) * info (4)) / 2**20,
     $      info (43), info (44), info (45)
90      format ('numeric factorization:',/,
     $      '   status:  ', f5.0, /,
     $      '   time:    ', e10.2, /,
     $      '   actual numeric LU statistics:', /,
     $      '   size of LU:    ', f10.2, ' (MB)', /,
     $      '   memory needed: ', f10.2, ' (MB)', /,
     $      '   flop count:    ', e10.2, /
     $      '   nnz (L):       ', f10.0, /
     $      '   nnz (U):       ', f10.0)

c       check umf4znum error condition
        if (info (1) .lt. 0) then
            print *, 'Error occurred in umf4znum: ', info (1)
            stop
        endif

c       save the symbolic analysis to the file s42.umf
c       note that this is not needed until another matrix is
c       factorized, below.
cU	filenum = 42
cU        call umf4zssym (symbolic, filenum, status)
cU        if (status .lt. 0) then
cU            print *, 'Error occurred in umf4zssym: ', status
cU            stop
cU        endif

c       save the LU factors to the file n0.umf
cU        call umf4zsnum (numeric, filenum, status)
cU        if (status .lt. 0) then
cU            print *, 'Error occurred in umf4zsnum: ', status
cU            stop
cU        endif

c       free the symbolic analysis
        call umf4zfsym (symbolic)

c       free the numeric factorization
cU        call umf4zfnum (numeric)

c       No LU factors (symbolic or numeric) are in memory at this point.

c	***** FIRST RHS *****
c       ----------------------------------------------------------------
c       load the LU factors back in, and solve the system
c       ----------------------------------------------------------------

c       At this point the program could terminate and load the LU
C       factors (numeric) from the n0.umf file, and solve the
c       system (see below).  Note that the symbolic object is not
c       required.

c       load the numeric factorization back in (filename: n0.umf)
cU        call umf4zlnum (numeric, filenum, status)
cU        if (status .lt. 0) then
cU            print *, 'Error occurred in umf4zlnum: ', status
cU            stop
cU        endif

c	Construct 1st rhs
	do i=1,numf
	BB(i)=d(i,1)
c		write(6,*)'BB(',i,')=',BB(i)
            b  (i) = dble (BB (i))
c		write(6,*)'b(',i,')=',b(i)
            bz (i) = imag (BB (i))
c		write(6,*)'bz(',i,')=',bz(i)
	enddo

c       solve Ax=b, without iterative refinement
        sys = 0
        call umf4zsol (sys, xumf, xz, b, bz, numeric, control, info)
        if (info (1) .lt. 0) then
            print *, 'Error occurred in umf4zsol: ', info (1)
            stop
        endif
        do i = 1,numf
            XX (i) = dcmplx (xumf (i), xz (i))
c	Save back into d(i,_) array
	    d(i,1)=XX(i)
	    write(6,*)'First rhs: XX(',i,')=',XX(i)
	enddo

c       free the numeric factorization
cU        call umf4zfnum (numeric)

c       No LU factors (symbolic or numeric) are in memory at this point.

c       print final statistics
c        call umf4zpinf (control, info)

c       print the residual.  [x (i) should be 1 + i/n]
c        call resid (numf, numfz, Ap, Ai, AA, XX, BB, r)

c	***** SECOND RHS *****
c       ----------------------------------------------------------------
c       load the LU factors back in, and solve the system
c       ----------------------------------------------------------------

c       At this point the program could terminate and load the LU
C       factors (numeric) from the n0.umf file, and solve the
c       system (see below).  Note that the symbolic object is not
c       required.

c       load the numeric factorization back in (filename: n0.umf)
cU        call umf4zlnum (numeric, filenum, status)
cU        if (status .lt. 0) then
cU            print *, 'Error occurred in umf4zlnum: ', status
cU            stop
cU        endif

c	Construct 2nd rhs
	do i=1,numf
	BB(i)=d(i,2)
c		write(6,*)'BB(',i,')=',BB(i)
            b  (i) = dble (BB (i))
c		write(6,*)'b(',i,')=',b(i)
            bz (i) = imag (BB (i))
c		write(6,*)'bz(',i,')=',bz(i)
	enddo

c       solve Ax=b, without iterative refinement
        sys = 0
        call umf4zsol (sys, xumf, xz, b, bz, numeric, control, info)
        if (info (1) .lt. 0) then
            print *, 'Error occurred in umf4zsol: ', info (1)
            stop
        endif
        do i = 1,numf
            XX (i) = dcmplx (xumf (i), xz (i))
c	Save back into d(i,_) array
	    d(i,2)=XX(i)
	    write(6,*)'Second rhs: XX(',i,')=',XX(i)
	enddo

c       free the numeric factorization
cU        call umf4zfnum (numeric)

c       No LU factors (symbolic or numeric) are in memory at this point.

c       print final statistics
c        call umf4zpinf (control, info)

c       print the residual.  [x (i) should be 1 + i/n]
c        call resid (numf, numfz, Ap, Ai, AA, XX, BB, r)

c	***** THIRD RHS *****
c       ----------------------------------------------------------------
c       load the LU factors back in, and solve the system
c       ----------------------------------------------------------------

c       At this point the program could terminate and load the LU
C       factors (numeric) from the n0.umf file, and solve the
c       system (see below).  Note that the symbolic object is not
c       required.

c       load the numeric factorization back in (filename: n0.umf)
cU        call umf4zlnum (numeric, filenum, status)
cU        if (status .lt. 0) then
cU            print *, 'Error occurred in umf4zlnum: ', status
cU            stop
cU        endif

c	Construct 3rd rhs
	do i=1,numf
	BB(i)=d(i,3)
c		write(6,*)'BB(',i,')=',BB(i)
            b  (i) = dble (BB (i))
c		write(6,*)'b(',i,')=',b(i)
            bz (i) = imag (BB (i))
c		write(6,*)'bz(',i,')=',bz(i)
	enddo

c       solve Ax=b, without iterative refinement
        sys = 0
        call umf4zsol (sys, xumf, xz, b, bz, numeric, control, info)
        if (info (1) .lt. 0) then
            print *, 'Error occurred in umf4zsol: ', info (1)
            stop
        endif
        do i = 1,numf
            XX (i) = dcmplx (xumf (i), xz (i))
c	Save back into d(i,_) array
	    d(i,3)=XX(i)
c	    write(6,*)'Third rhs: XX(',i,')=',XX(i)
	enddo

c       free the numeric factorization
cU        call umf4zfnum (numeric)

c       No LU factors (symbolic or numeric) are in memory at this point.

c       print final statistics
c        call umf4zpinf (control, info)

c       print the residual.  [x (i) should be 1 + i/n]
c        call resid (numf, numfz, Ap, Ai, AA, XX, BB, r)

c	***** FOURTH RHS *****
c       ----------------------------------------------------------------
c       load the LU factors back in, and solve the system
c       ----------------------------------------------------------------

c       At this point the program could terminate and load the LU
C       factors (numeric) from the n0.umf file, and solve the
c       system (see below).  Note that the symbolic object is not
c       required.

c       load the numeric factorization back in (filename: n0.umf)
cU        call umf4zlnum (numeric, filenum, status)
cU        if (status .lt. 0) then
cU            print *, 'Error occurred in umf4zlnum: ', status
cU            stop
cU        endif

c	Construct 4th rhs
	do i=1,numf
	BB(i)=d(i,4)
c		write(6,*)'BB(',i,')=',BB(i)
            b  (i) = dble (BB (i))
c		write(6,*)'b(',i,')=',b(i)
            bz (i) = imag (BB (i))
c		write(6,*)'bz(',i,')=',bz(i)
	enddo

c       solve Ax=b, without iterative refinement
        sys = 0
        call umf4zsol (sys, xumf, xz, b, bz, numeric, control, info)
        if (info (1) .lt. 0) then
            print *, 'Error occurred in umf4zsol: ', info (1)
            stop
        endif
        do i = 1,numf
            XX (i) = dcmplx (xumf (i), xz (i))
c	Save back into d(i,_) array
	    d(i,4)=XX(i)
c	    write(6,*)'Fourth rhs: XX(',i,')=',XX(i)
	enddo

c       free the numeric factorization
        call umf4zfnum (numeric)

c       No LU factors (symbolic or numeric) are in memory at this point.

c       print final statistics
c        call umf4zpinf (control, info)

c       print the residual.  [x (i) should be 1 + i/n]
c        call resid (numf, numfz, Ap, Ai, AA, XX, BB, r)

c		stop

c	open(2,file='matrel.solution1')
c	write(6,*)'solution vector 1='
c	do i=1,numf
c	write(6,*) i,d(i,1)
c	write(2,*) i,real(d(i,1))
c	enddo
c	close(2)

c	open(2,file='matrel.solution2')
c	write(6,*)'solution vector 2='
c	do i=1,numf
c	write(6,*) i,d(i,2)
c	write(2,*) i,real(d(i,2))
c	enddo
c	close(2)

c	open(2,file='matrel.solution3')
c	write(6,*)'solution vector 3='
c	do i=1,numf
c	write(6,*) i,d(i,3)
c	write(2,*) i,real(d(i,3))
c	enddo
c	close(2)

c	open(2,file='matrel.solution4')
c	write(6,*)'solution vector 4='
c	do i=1,numf
c	write(6,*) i,d(i,4)
c	write(2,*) i,real(d(i,4))
c	enddo
c	close(2)

	deallocate (AA)
	deallocate (Ax)
	deallocate (Az)

	return
	end

c=======================================================================
c== resid ==============================================================
c=======================================================================

c Compute the residual, r = Ax-b, its max-norm, and print the max-norm
C Note that A is zero-based.

        subroutine resid (n, nz, Ap, Ai, A, x, b, r)
        integer*8
     $      n, nz, Ap (n+1), Ai (n), j, i, p
        complex*16 A (nz), x (n), b (n), r (n), aij
	double precision rmax

        do 10 i = 1, n
            r (i) = -b (i)
10      continue

        do 30 j = 1,n
            do 20 p = Ap (j) + 1, Ap (j+1)
                i = Ai (p) + 1
                aij = A (p)
                r (i) = r (i) + aij * x (j)
20          continue
30      continue

        rmax = 0
        do 40 i = 1, n
            rmax = max (rmax, abs (r (i)))
40      continue

        print *, 'norm (A*x-b): ', rmax
        return
        end

	subroutine gridgeom(sfile)
	character*80 sfile
	parameter (ncmax=151)
	parameter (ncmaz=74)
c	Next line is the number of interior GLL points in 1D
	parameter (lmax=4)
c	N is the total number of GLL points in 1D including the points +-1
	parameter (N=lmax+2)
c	Next line is the maximum number of GLL points on the x-side of
c	the global grid times the maximum number on the z-side.
	parameter (nptmax=(N*ncmax-(ncmax-1))*(N*ncmaz-(ncmaz-1)))
c	xg and zg have the (x,z) coordinates at each of the global gridpoints
c	assuming dimensions of 2x2 for each cell (x and z each run from -1 to +1).
	common/glbxy/xg(nptmax),zg(nptmax)
	real*8 vp0(nptmax),vs0(nptmax),rho0(nptmax),
     &  qb0(nptmax),qk0(nptmax)
	dimension x0(nptmax),z0(nptmax)
	character*80 filin
	common/grdpar/NX,NZ
	real*8 vpval,vsval,rhoval,qbval,qkval
	real*8 vp,vs,rho,qb,qk
	common/struc3/dx(ncmax),dz(ncmaz)
	common/struc4/vp(ncmax,ncmaz,N*N),vs(ncmax,ncmaz,N*N),
     &	rho(ncmax,ncmaz,N*N),qb(ncmax,ncmaz,N*N),qk(ncmax,ncmaz,N*N)
	common/glbgrd/igrd(ncmax,ncmaz,N*N)
c
	pi=3.1415926535
	rad=180./pi
	open(2,file=sfile)
	rewind(2)
	write(6,*)'Reading in structural parameters'
	i=0
5	read(2,*,end=10) xin,zin,vpin,vsin,rhoin,qbin,qkin
	i=i+1
	vp0(i)=dble(vpin)
	vs0(i)=dble(vsin)
	rho0(i)=dble(rhoin)
c	Replace qbet and qkap with their inverses.
	qb0(i)=1.d0/qbin
	qk0(i)=1.d0/qkin
c	Convert (xin[deg.],zin[km]) to (xin[km],zin[km])
c	and reference to the initial x0-value
	x0(i)=(xin/rad)*6371.
	x0(i)=x0(i)-x0(1)
	z0(i)=zin
c		write(6,*)'x0 and z0(',i,')=',x0(i),z0(i)
	go to 5
10	continue
	close(2)
	imax=i
	write(6,*)'Finished reading structural parameters at',imax,'points'
c	Fill up model grid with values based on closest distance to just read-in values
c	igmax=(N*NX-(NX-1))*(N*NZ-(NZ-1))
c	do ig=1,igmax

	  do ncz=1,NZ
	  do ncx=1,NX
	  do ilz=1,N
	  do ilx=1,N
	il=N*(ilz-1)+ilx
	ig=igrd(ncx,ncz,il)

	xval=(xg(ig)-xg(1))*6371.
	if(ilx.eq.1) xval=xval+0.01*dx(ncx)*6371.
	if(ilx.eq.N) xval=xval-0.01*dx(ncx)*6371.
	zval=zg(ig)
	if(ilz.eq.1) zval=zval+0.01*dz(ncz)
	if(ilz.eq.N) zval=zval-0.01*dz(ncz)
c	find closest read-in gridpoint to (xg(ig),yg(ig))
	do i=1,imax
	dist=(xval-x0(i))**2 + (zval-z0(i))**2
c		if(ig.eq.94) write(6,*)'A latest dist=',dist
c		if(ig.eq.94) write(6,*)'A xg,zg=',xg(ig),zg(ig)
c		if(ig.eq.94) write(6,*)'A x0,z0=',x0(i),z0(i)
	if(i.eq.1) then
	distm=dist
	im=i
	endif
	if(i.gt.1.and.dist.lt.distm) then
	distm=dist
	im=i
c		if(ig.eq.94) write(6,*)'latest dist=',dist
c		if(ig.eq.94) write(6,*)'xg,zg=',xg(ig),zg(ig)
c		if(ig.eq.94) write(6,*)'x0,z0=',x0(i),z0(i)
c		if(ig.eq.94) write (6,*)' ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^ ^'
	endif
	enddo
c		if(ig.eq.94) write(6,*)'closest i=',im
		write(6,*)'NX=',NX,'NZ=',NZ
		write(6,*)'ig=',ig,'im=',im,'nptmax=',nptmax
	vp(ncx,ncz,il)=vp0(im)
	vs(ncx,ncz,il)=vs0(im)
	rho(ncx,ncz,il)=rho0(im)
	qb(ncx,ncz,il)=qb0(im)
	qk(ncx,ncz,il)=qk0(im)
c-----
c	Simulate absorbing boundaries at the lateral edges and bottom of the model
c	by prescribing very low Q layers.
c	if(ncx.le.2.or.ncx.ge.(NX-1)) then
c	qb(ig)=1.d0/3.d0
c	qk(ig)=1.d0/3.d0
c	endif
c	if(ncz.le.2) then
c	qb(ig)=1.d0/3.d0
c	qk(ig)=1.d0/3.d0
c	endif
c-----
		write(6,*)'GRIDGEOM: x,z=',xg(ig),'km',zg(ig),'km'
		write(6,*)'vp and vs(',ncx,ncz,il,')=',vp(ncx,ncz,il),vs(ncx,ncz,il)
		write(6,*)'qb(',ncx,ncz,il,')=',qb(ncx,ncz,il),'qk(',ncx,ncz,il,')=',qk(ncx,ncz,il)
		write(6,*)'------------------------'
	  enddo
	  enddo
	  enddo
	  enddo

c	enddo

	return
	end

	subroutine topogeom(tfile)
	character*80 tfile
	parameter (ncmax=151)
	parameter (ncmaz=74)
c	Next line is the number of interior GLL points in 1D
	parameter (lmax=4)
c	N is the total number of GLL points in 1D including the points +-1
	parameter (N=lmax+2)
c	Next line is the maximum number of GLL points on the x-side of
c	the global grid times the maximum number on the z-side.
	parameter (nptmax=(N*ncmax-(ncmax-1))*(N*ncmaz-(ncmaz-1)))
c	xg and zg have the (x,z) coordinates at each of the global gridpoints
c	assuming dimensions of 2x2 for each cell (x and z each run from -1 to +1).
	common/glbxy/xg(nptmax),zg(nptmax)
	common/struc3/dx(ncmax),dz(ncmaz)
	dimension indz(ncmaz)
	common/topo/dztopo(ncmaz,ncmax+1)
	real*8 jacm1
	common/struc5/jacm1(ncmax,ncmaz,N*N)
	real*8 wt,x
	common/glpts/wt(lmax+2),x(lmax+2)
	common/grdpar/NX,NZ
c	common/glbgrd/igrd(ncmax,ncmaz,N*N)
c	Initialize dztopo
	do j=1,NZ
	do i=1,NX+1
	dztopo(j,i)=0.
	enddo
	enddo
c
	open(2,file=tfile)
	rewind(2)
	read(2,5) b80
	read(2,*) nztopo
		if(nztopo.gt.NZ) then
		write(6,*)'# topographic layers',nztopo,'exceeds #z-layers=',NZ
		stop
		endif
	read(2,5) b80
	read(2,*) (indz(k), k=1,nztopo)
	read(2,5) b80
	do k=1,nztopo
	j=indz(k)
		if(j.gt.NZ) then
		write(6,*)'topographic index ind(',k,')=',j,'exceeds largest z-layer index=',NZ
		stop
		endif
	read(2,*) (dztopo(j,i), i=1,NX+1) 
	enddo
	close(2)
c	First set all topographic stretching Jacobian values to unity
c	dr * d(theta) = J dx * dz as employed in subroutine MATREL 
	do ncz=1,NZ
	do ncx=1,NX
	do ilz=1,N
	do ilx=1,N
	il=N*(ilz-1)+ilx
	jacm1(ncx,ncz,il)=1.d0
	enddo
	enddo
	enddo
	enddo

	do ncz=1,NZ
	do ncx=1,NX
	do ilz=1,N
	do ilx=1,N
	il=N*(ilz-1)+ilx
c	ig=igrd(ncx,ncz,il)
c	thfac1=x' of notes; x goes from -1 to 1, x' goes from 0 to 1 across the cell.
	thfac1=(x(ilx)-x(1))/(x(N)-x(1))
c	thfac2=1-x' of notes.
	thfac2=(x(N)-x(ilx))/(x(N)-x(1))
cOLD	if(ncz.gt.1) rfac=(dz(ncz)+dztopo(ncz,ncx+1)-dztopo(ncz-1,ncx+1))/(dz(ncz)+dztopo(ncz,ncx)-dztopo(ncz-1,ncx))
cOLD	if(ncz.eq.1) rfac=(dz(ncz)+dztopo(ncz,ncx+1))/(dz(ncz)+dztopo(ncz,ncx))
cOLD	jacm1(ncx,ncz,il)=dble(rfac*thfac1+thfac2)
	if(ncz.gt.1) then
	rfac1=1.+(dztopo(ncz,ncx+1)-dztopo(ncz-1,ncx+1))/dz(ncz)
	rfac2=1.+(dztopo(ncz,ncx)-dztopo(ncz-1,ncx))/dz(ncz)
	endif
	if(ncz.eq.1) then
	rfac1=1.+dztopo(ncz,ncx+1)/dz(ncz)
	rfac2=1.+dztopo(ncz,ncx)/dz(ncz)
	endif
	jacm1(ncx,ncz,il)=dble(thfac1*rfac1+thfac2*rfac2)
cNOTE	Additional factor of (dz(ncz)/2.)*(dx(ncx/2.) is already accounted for in MATREL.
	enddo
	enddo
	enddo
	enddo

5	format(a80)
	return
	end

	subroutine source(w0,ky)
	parameter (ncmax=151)
	parameter (ncmaz=74)
c	Next line is the number of interior GLL points in 1D
	parameter (lmax=4)
c	N is the total number of GLL points in 1D including the points +-1
	parameter (N=lmax+2)
c	Next line is the maximum number of GLL points on the x-side of
c	the global grid times the maximum number on the z-side.
	parameter (nptmax=(N*ncmax-(ncmax-1))*(N*ncmaz-(ncmaz-1)))
c	xg and zg have the (x,z) coordinates at each of the global gridpoints
c	assuming dimensions of 2x2 for each cell (x and z each run from -1 to +1).
	common/glbxy/xg(nptmax),zg(nptmax)
c	igc stored in common will have the grid # of the last point
c	determined on the source plane -- representative of the source location.
	common/isrc/igc
c
	character*80 b80,filin
c	Assume a maximum of 200 sources.
	dimension istyp(200)
c	istyp=0 for moment tensor source; istyp=1 for force source.
	common/grdpar/NX,NZ
	real*8 vp,vs,rho,qb,qk
	common/struc3/dx(ncmax),dz(ncmaz)
	common/struc4/vp(ncmax,ncmaz,N*N),vs(ncmax,ncmaz,N*N),
     &	rho(ncmax,ncmaz,N*N),qb(ncmax,ncmaz,N*N),qk(ncmax,ncmaz,N*N)
	common/topo/dztopo(ncmaz,ncmax+1)
	common/glbgrd/igrd(ncmax,ncmaz,N*N)
	real*8 dmax,dmin,dip,cdip,sdip,c2dip
	real*8 sstr,cstr,s2str,c2str,frak,srak,crak
	real*8 deltl,phil,deltd,delt0
	real*8 mu0
	real*8 shrm1,p1,p2,p3,p4,p5
	real*8 flat(200),flon(200),fbigl(200),fstr(200),frake(200),fwt(200)
	real*8 mrr0(200),mtt0(200),mpp0(200),mrt0(200),mrp0(200),mtp0(200)

	real*8, ALLOCATABLE :: fr0(: , :)
	real*8, ALLOCATABLE :: ft0(: , :)
	real*8, ALLOCATABLE :: fp0(: , :)
cOLD	real*8 fr0(200,2000),ft0(200,2000),fp0(200,2000)

	real*8 plat,plon,phiref,pi,rad
	real*8 dlat,cdelt,sdelt,delta,sphi,cphi,phisrc,phi1
	common/sphpar/dlat,plat,plon,phiref
	real*8 wt,x
	common/glpts/wt(lmax+2),x(lmax+2)
c***
	real*8 mxx,myy,mzz,mxy,mxz,myz
	real*8 mtt,mpp,mtp,mrt,mrp,mrr
	complex*16 fx,fy,fz
	complex*16 fr,ft,fp
	real*8 argr,argi,dts
	complex*16 efac0,efac,defac
	complex*16 d
	common/datavec/d(3*nptmax,4)
c***
	real*8 slat,slon
	real*8 phis,dphisx,dphisz
	real*8 phiv,dphivx,dphivz
	real*8 ky,fac1,fac2,cott
	complex*16 s,w0,ui
	complex*16 facc,facs
	parameter (bigr=6371.)
	pi=3.14159265358979d0
	rad=180.d0/pi
	ui=dcmplx(0.d0,1.d0)
	s=w0*ui

c	Read in source parameters.
	open(2,file='seis2pt5dsource-spherical.param')
	rewind(2)
	read(2,5) b80
	write(6,5) b80
	read(2,*) dmax,dmin,dip
	write(6,*) dmax,dmin,dip
	cdip=dcos(dip/rad)
	sdip=dsin(dip/rad)
	c2dip=cdip*cdip-sdip*sdip
	read(2,5) b80
	write(6,5) b80
	read(2,*) iseg
	write(6,*)'iseg=',iseg
	read(2,5) b80
	i=0
16	i=i+1
	read(2,5) b80
	read(b80,*,err=35) flat(i),flon(i),fbigl(i),fstr(i),frake(i),fwt(i)
	write(6,*) flat(i),flon(i),fbigl(i),fstr(i),frake(i),fwt(i)
	flat(i)=(pi/2.d0-flat(i)/rad)
	flon(i)=flon(i)/rad
	fstr(i)=fstr(i)/rad 
	frake(i)=frake(i)/rad 
	  istyp(i)=0
	  if(fwt(i).ge.0.0) go to 25
	  read(2,5) b80
          write(6,*)'Read in moment tensor components'
          write(6,*)'mrr,mtt,mpp,mrt,mrp,mtp (10^20 N m)'
	  write(6,*)'or time-reversed seismgrams implemented as point forces'
c	  write(6,*)'or force components fr,ft,fp (10^17 N)'
c         write(6,*)'where r=r hat=Up, t=theta hat=South, p=phi hat=East'
          read(b80,*) mrr0(i),mtt0(i),mpp0(i),mrt0(i),mrp0(i),mtp0(i)
25	continue
	if(i.lt.iseg) go to 16
	go to 40

35	continue
	  write(6,*)'Implement 3-component time-reversed seismograms at',iseg,' points'
	  write(6,*) iseg,' lat,lons on sites contributing seismograms?'
	  open(4,file=b80)
	  rewind(4)
	  do i=1,iseg
	  read(4,*)  flat(i),flon(i)
	  write(6,*)'flat(',i,')=',flat(i),'flon(',i,')=',flon(i)
	  flat(i)=(pi/2.d0-flat(i)/rad)
	  flon(i)=flon(i)/rad
	  fbigl(i)=1.e-6
	  fstr(i)=0.
	  frake(i)=0.
	  fwt(i)=-1.0
	  enddo
	  close(4)
	  write(6,*)'length of seismograms, time increment (s)?'
	  write(6,*)'Read in East-component, North-component, Up-component seismograms'
	  write(6,*)'Assume all start at time [dts] after t=0, and all seismograms have'
	  write(6,*)'the same time length and increment dts (determined by the program'
	  write(6,*)'from the seismograms)'
	  write(6,*)'length of seismograms?'
	  read(2,*) lens
	  write(6,*)'lens=',lens

	allocate ( fr0(iseg,lens) )
	allocate ( ft0(iseg,lens) )
	allocate ( fp0(iseg,lens) )

	  write(6,*)'File containing E, N, and Z-component seismograms?'
	  read(2,5) filin
	  open(4,file=filin)
	  rewind(4)	  
	  do i=1,iseg
	  do j=1,lens
	  read(4,*) tin,fp0(i,j),ft0(i,j),fr0(i,j)
c	  write(6,*) tin,fp0(i,j),ft0(i,j),fr0(i,j)
c	Cnange sign of ft0 to convert it to theta hat=South direction
	  ft0(i,j)=-ft0(i,j)
	  if(j.gt.1) dts=tin-t0
	  t0=tin
	  enddo
	  istyp(i)=1
	  enddo
	  close(4)
c		write(6,*)'dts=',dts
c		if(i.ne.9999) stop
40	continue
	read(2,5) b80
	read(2,*) idstr
	read(2,5) b80
	read(2,*) iddip
	close(2)

c	read(2,*) mtt,mpp,mzz,mtp,mtz,mpz
c	read(2,5) b80
c	read(2,*) slat,slon,zs
c	slat=pi/2.d0-slat/rad
c	slon=slon/rad	

c	Determine data vector of inverse problem
	write(6,*)'SOURCE: Determine data vector of inverse problem'
	write(6,*)'arising from the source term'

	write(6,*)'Max size of d=',3*nptmax
	write(6,*)'N=',N,'NX,NZ=',NX,NZ
	write(6,*)'Will zero out elements up to index',3*(N*NX-(NX-1))*(N*NZ-(NZ-1))
	do i=1,3*(N*NX-(NX-1))*(N*NZ-(NZ-1))
	d(i,1)=0.d0
	d(i,2)=0.d0
	d(i,3)=0.d0
	d(i,4)=0.d0
	enddo

	i=0
10	i=i+1
	if(i.gt.iseg) go to 30
	sstr=dsin(fstr(i))
	cstr=dcos(fstr(i))
	s2str=2.d0*sstr*cstr
	c2str=cstr*cstr-sstr*sstr 
	frak=frake(i)
	srak=dsin(frake(i))
	crak=dcos(frake(i))
c	Divide fault plane into [idstr] x [iddip] patches
	ilen=0
15	ilen=ilen+1
	if(ilen.gt.idstr) go to 10
c	Determine position (delta,phil) of running point on lower edge
c	This involves starting at (flat(i),flon(i)) and moving along
c	a direction striking fstr(i)-pi
	deltl=(fbigl(i)/dble(bigr))*(dble(ilen)-0.5d0)/dble(idstr)
        cdelt=dcos(flat(i))*dcos(deltl)+dsin(flat(i))*dsin(deltl)*(-cstr)
        delta=dacos(cdelt)
	sdelt=dsin(delta)
        sphi=dsin(deltl)*(-sstr)/sdelt
        cphi=(dcos(deltl)-dcos(flat(i))*cdelt)/(dsin(flat(i))*sdelt)
	phil=datan2(sphi,cphi)
c		write(6,*)'SOURCE: ilen=',ilen,'delta,phi=',90.-delta*rad,(phil+flon(i))*rad
	delt0=delta
	idip=0
20	idip=idip+1
	if(idip.gt.iddip) go to 15
c	Determine position (slat,slon) of running point along updip direction
c	This involves starting at (delt0,phil+flon(i)) and moving along
c	a direction striking fstr(i)-pi/2
	deltd=((dmax-dmin)*(cdip/sdip)/dble(bigr))*(dble(idip)-0.5d0)/dble(iddip)
	cdelt=dcos(delt0)*dcos(deltd)+dsin(delt0)*dsin(deltd)*(sstr)
        slat=dacos(cdelt)
	sdelt=sin(slat)
        sphi=dsin(deltd)*(-cstr)/sdelt
        cphi=(dcos(deltd)-dcos(delt0)*cdelt)/(dsin(delt0)*sdelt)
	slon=datan2(sphi,cphi)+phil+flon(i)
	zs=dmax+(dmin-dmax)*(dble(idip)-0.5d0)/dble(iddip)
	zs=-zs
c		write(6,*)'SOURCE: ilen,idip=',ilen,idip
c		write(6,*)'SOURCE: slat,slon,zs=',90.-slat*rad,slon*rad,zs
c		write(6,*)'------------'

c*	Determine angular distance (rad) and azimuth (rad CCW from due N)
c	of the source from the pole of the spherical coordinate system.
c	the azimuth is referenced to the azimuth of the origin (theta=dlat,z=0)
c	of the 2D grid.
        cdelt=dcos(slat)*dcos(plat)+dsin(slat)*dsin(plat)*dcos(plon-slon)
        delta=dacos(cdelt)
	sdelt=dsin(delta)
        sphi=dsin(slat)*dsin(plon-slon)/sdelt
        cphi=(dcos(slat)-dcos(plat)*cdelt)/(dsin(plat)*sdelt)
	phisrc=datan2(sphi,cphi)-phiref
cTE
c	phisrc=0.
c--
c		write(6,*)'SOURCE: phi_src=',rad*phisrc
c*	The theta-coord of the source is the angular distance from (plat,plon)
	xs=real(delta)
c		write(6,*)'SOURCE: xs=',xs
c	Determine which cell number the source is located in
	ncxs=0
	nczs=0
	do ncz=1,NZ
	do ncx=1,NX
c	lower left corner
	ilx=1
	ilz=1
	il=N*(ilz-1)+ilx
	ig=igrd(ncx,ncz,il)
	x1=xg(ig)
	z1=zg(ig)
c		write(6,*)'SOURCE: cell #',ncx,ncz,'dimensionless lower left corner=',x1,z1
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
	  ncxs=ncx
	  nczs=ncz
	  xsource=(2./dx(ncx))*(xs-x1)-1.
	  zsource=(2./dz(ncz))*(zs-z1)-1.
c		write(6,*)'xs,x1=',xs,x1
c		write(6,*)'zs,z1=',zs,z1
		write(6,*)'xsource,zsource=',xsource,zsource
	  ilc=N*((N/2)-1)+(N/2)
	  igc=igrd(ncx,ncz,ilc)
		write(6,*)'SOURCE: gridpoint of cell center=',igc
	  endif
	enddo
	enddo
	write(6,*)'SOURCE: Cell #s of source =(',ncxs,nczs,'). ilen,idip=',ilen,idip
	if(ncxs.eq.0.or.nczs.eq.0) then
	write(6,*)'SOURCE: Source location lies outside model domain'
	stop
	endif
c	write(6,*)'Local dimensionless source coords=(',xsource,zsource,')'
c*
c	Three cases for assigning source parameters:
c	1) moment tensor was read in directly
	  if(fwt(i).lt.0.0.and.istyp(i).eq.0) then
          mrr=mrr0(i)/dble(idstr*iddip)
          mtt=mtt0(i)/dble(idstr*iddip)
          mpp=mpp0(i)/dble(idstr*iddip)
          mrt=mrt0(i)/dble(idstr*iddip)
          mrp=mrp0(i)/dble(idstr*iddip)
          mtp=mtp0(i)/dble(idstr*iddip)
          endif
c	2) moment tensor needs to be computed from shear dislocation dimensions and slip
	  if(fwt(i).ge.0.0) then
c	Moment tensor from Ben Menahem and Singh, eqn. 4.115b for
c	shear dislocation, and derived from eqn. (4.101), (4.110), and
c	(4.113) for a tensile dislocation.
c		write(6,*)'strainA-momten: abs(frak)=',abs(frak)
cNOTE this line had previously mistakenly read:	mu0=0.1*rho(ncx,ncz,ilc)*vs(ncx,ncz,ilc)**2
	mu0=0.1*rho(ncxs,nczs,ilc)*vs(ncxs,nczs,ilc)**2
c	Next line is shear moment.
	shrm1=fwt(i)*fbigl(i)*((dmax-dmin)/sdip)*mu0*1.d-6 / dble(idstr*iddip)
c		write(6,*)'fwt=',fwt(i),'mu0=',mu0,'shrm1=',shrm1
	p1=srak*sdip*cdip*s2str+crak*sdip*c2str
	p2=crak*sdip*s2str-srak*sdip*cdip*c2str
	p3=-crak*cdip*sstr+srak*c2dip*cstr
	p4=srak*c2dip*sstr+crak*cdip*cstr
	p5=srak*sdip*cdip
c		write(6,*)'p1-5=',p1,p2,p3,p4,p5
	mrr=shrm1*2.d0*p5
	mtt=-shrm1*(p2+p5)
	mpp=shrm1*(p2-p5)
	mrt=-shrm1*p4
	mrp=-shrm1*p3
	mtp=-shrm1*p1
	endif
c	3) force vector was read in directly
	  if(istyp(i).eq.1) then
	  fr=0.d0
	  ft=0.d0
	  fp=0.d0
c	Take the complex conjugate of the Fourier transform of each time series
c	to yield the FT of the time-reversed seismograms. We have a exp(-ui*w0*t) convention
c	fore the FT, so the sum below uses efac ~ exp(+ui*w0*t)=exp(s*(j*dts))
	  argr=dreal(s*dts)
	  argi=dimag(s*dts)
	  defac=dexp(argr)*(dcos(argi)+ui*dsin(argi))
c	Assign a `origin' time of the adjoint seismograms to be one-half of the
c	way into the first input seismogram.
	  if(i.eq.1) then
	  argr=dreal(s*dts*dble(lens)/2.d0)
	  argi=dimag(s*dts*dble(lens)/2.d0)
	  efac0=dexp(-argr)*(dcos(argi)-ui*dsin(argi))
	  endif
c	  efac=1.d0
	  efac=efac0
	  do j=1,lens
	  efac=efac*defac
	  fr=fr+fr0(i,j)*efac*dts
	  ft=ft+ft0(i,j)*efac*dts
	  fp=fp+fp0(i,j)*efac*dts
	  enddo
          fr=fr/dble(idstr*iddip)
          ft=ft/dble(idstr*iddip)
          fp=fp/dble(idstr*iddip)
	  endif
cNOTE	At this point, the moment tensor or force vector, which should be read in wrt local
c	geographic r,theta,phi coordinates, should be rotated by an angle phi1
c	about the vertical axis at the source point. 
        sphi=dsin(plat)*dsin(plon-slon)/sdelt
        cphi=(dcos(plat)-dcos(slat)*cdelt)/(dsin(slat)*sdelt)
	phi1=datan2(sphi,cphi)
c		write(6,*)'SOURCE: phi1=',rad*phi1
c		if(i.ne.9999) stop
c	Rotate moment tensor components an amount phi1 about a vertical axis
c	at the source.
	if(istyp(i).eq.0) then
	mxx=mtt*cphi*cphi + mpp*sphi*sphi - 2.d0*mtp*cphi*sphi
	myy=mpp*cphi*cphi + mtt*sphi*sphi + 2.d0*mtp*cphi*sphi
	mxy=(mtt-mpp)*cphi*sphi + mtp*(cphi*cphi-sphi*sphi)
	mxz=mrt*cphi - mrp*sphi
	myz=mrp*cphi + mrt*sphi
	mzz=mrr
	else
c	Rotate force components an amount phi1 about a vertical axis
c       at the source.
	fx=ft*cphi - fp*sphi
	fy=fp*cphi + ft*sphi
	fz=fr
	endif
cTE
c		mxx=0.
c		myy=0.
c		mxy=0.
c		mxz=2.d-2
c		myz=0.
c		mzz=0.
c--
c		write(6,*)'SOURCE: original mtt,mpp,mrr,mtp,mrt,mrp=',mtt,mpp,mrr,mtp,mrt,mrp
c		write(6,*)'SOURCE: rotated  mxx,myy,mzz,mxy,mxz,myz=',mxx,myy,mzz,mxy,mxz,myz 
c*
cTE
c		if(i.ne.9999) stop
cTOPO
c	thfac1=x' of notes; x goes from -1 to 1, x' goes from 0 to 1 across the cell.
	thfac1=(xsource-x(1))/(x(N)-x(1))
c	thfac2=1-x' of notes.
	thfac2=(x(N)-xsource)/(x(N)-x(1))
	if(nczs.gt.1) then
	rfac1=1.+(dztopo(nczs,ncxs+1)-dztopo(nczs-1,ncxs+1))/dz(nczs)
	rfac2=1.+(dztopo(nczs,ncxs)-dztopo(nczs-1,ncxs))/dz(nczs)
	endif
	if(nczs.eq.1) then
	rfac1=1.+dztopo(nczs,ncxs+1)/dz(nczs)
	rfac2=1.+dztopo(nczs,ncxs)/dz(nczs)
	endif
	sjacm1=thfac1*rfac1+thfac2*rfac2
c--
	do ilz=1,N
	do ilx=1,N
	il=N*(ilz-1)+ilx
	i1=3*igrd(ncxs,nczs,il)-2
	i2=i1+1
	i3=i2+1
c	i1,i2,i3 are the row numbers corresponding to the normal equation
c	resulting from multiplication of the momentum eqn with phi_[il]
c	in cell (ncxs,nczs) and integration over that cell.
c	i1,i2,i3 refer to the x,y,z components.
c	Next interpolate the basis function phi2_[il] at x=xsource,z=zsource
	phis=phiv(ilx,ilz,xsource,zsource)
c	Next interpolate the basis function d(phi2_[il])/dx at x=xsource,z=zsource
	dphisx=dphivx(ilx,ilz,xsource,zsource) * dble(2./dx(ncxs))
c	Next interpolate the basis function d(phi2_[il])/dz at x=xsource,z=zsource
	dphisz=dphivz(ilx,ilz,xsource,zsource) * dble(2./dz(nczs))
c	Eqn 10 of notes
c		write(6,*)'phis=',phis
c		write(6,*)'dphisx=',dphisx
c		write(6,*)'dphisz=',dphisz
c		write(6,*)'i1,i2,i3=',i1,i2,i3
c		write(6,*)'3*nptmax=',3*nptmax
c		write(6,*)'term1=',(-mxx*dphisx+ui*ky*mxy*phis-mxz*dphisz)
c		write(6,*)'term2=',(-mxy*dphisx+ui*ky*myy*phis-myz*dphisz)
c		write(6,*)'term3=',(-mxz*dphisx+ui*ky*myz*phis-mzz*dphisz)
	facc=dcos(ky*phisrc)
	facs=-ui*dsin(ky*phisrc)
	fac1=1.d0/(bigr+zs)
	fac2=1.d0/((bigr+zs)*sdelt)
	cott=cdelt/sdelt
c	d(i1,1)=(-mxx*(dphisx-phis*cott)*fac1+ui*ky*mxy*phis*fac2-mxz*(dphisz-phis*fac1))/s * fac * fac2 
c	d(i2,1)=(-mxy*(dphisx-phis*cott)*fac1+ui*ky*myy*phis*fac2-myz*(dphisz-phis*fac1))/s * fac * fac2
c	d(i3,1)=(-mxz*(dphisx-phis*cott)*fac1+ui*ky*myz*phis*fac2-mzz*(dphisz-phis*fac1))/s * fac * fac2
c	Matrix elements representing the source vector.  Note that the minc/twopi factor
c	that appears in eqn B-6 of Pollitz (2014 GJI) is accounted for elsewhere (facp,facs
c	in subroutine addm, etc.).
c	Moment tensor source with H(t) (step function) time dependence.
	if(istyp(i).eq.0) then
	d(i1,1)=d(i1,1) + (-mxx*(dphisx)*fac1-mxz*(dphisz))/s * facc * fac2 * sjacm1
	d(i2,1)=d(i2,1) + (ui*ky*myy*phis*fac2)/s * facc * fac2 * sjacm1
	d(i3,1)=d(i3,1) + (-mxz*(dphisx)*fac1-mzz*(dphisz))/s * facc * fac2 * sjacm1
	d(i1,2)=d(i1,2) + (ui*ky*mxy*phis*fac2)/s * facc * fac2 * sjacm1
	d(i2,2)=d(i2,2) + (-mxy*(dphisx)*fac1-myz*(dphisz))/s * facc * fac2 * sjacm1
	d(i3,2)=d(i3,2) + (ui*ky*myz*phis*fac2)/s * facc * fac2 * sjacm1
	d(i1,3)=d(i1,3) + (-mxx*(dphisx)*fac1-mxz*(dphisz))/s * facs * fac2 * sjacm1
	d(i2,3)=d(i2,3) + (ui*ky*myy*phis*fac2)/s * facs * fac2 * sjacm1
	d(i3,3)=d(i3,3) + (-mxz*(dphisx)*fac1-mzz*(dphisz))/s * facs * fac2 * sjacm1
	d(i1,4)=d(i1,4) + (ui*ky*mxy*phis*fac2)/s * facs * fac2 * sjacm1
	d(i2,4)=d(i2,4) + (-mxy*(dphisx)*fac1-myz*(dphisz))/s * facs * fac2 * sjacm1
	d(i3,4)=d(i3,4) + (ui*ky*myz*phis*fac2)/s * facs * fac2 * sjacm1
	else
c	Force vector source with delta(t) (Dirac delta function) time dependence.
	d(i1,1)=d(i1,1) - (fx*phis) * facc * fac2 * sjacm1
	d(i2,1)=0.d0
	d(i3,1)=d(i3,1) - (fz*phis) * facc * fac2 * sjacm1
	d(i1,2)=0.
	d(i2,2)=d(i2,2) - (fy*phis) * facc * fac2 * sjacm1
	d(i3,2)=0.
	d(i1,3)=d(i1,3) - (fx*phis) * facs * fac2 * sjacm1
	d(i2,3)=0.d0
	d(i3,3)=d(i3,3) - (fz*phis) * facs * fac2 * sjacm1
	d(i1,4)=0.
	d(i2,4)=d(i2,4) - (fy*phis) * facs * fac2 * sjacm1
	d(i3,4)=0.
	endif
c		write(6,*)'d(',i1,1,')=',d(i1,1)
c		write(6,*)'d(',i2,1,')=',d(i2,1)
c		write(6,*)'d(',i3,1,')=',d(i3,1)
c		write(6,*)'d(',i1,2,')=',d(i1,2)
c		write(6,*)'d(',i2,2,')=',d(i2,2)
c		write(6,*)'d(',i3,2,')=',d(i3,2)
	enddo
	enddo
cTE
c		if(ilen.eq.idstr.and.idip.eq.iddip) stop
	go to 20
30	continue
	if(istyp(1).eq.1) then
	deallocate (fr0)
	deallocate (ft0)
	deallocate (fp0)
	endif
	
5	format(a80)
	return
	end

	subroutine invft(nleno,dt,obeta,corper)
c	Obtain inverse Fourier transform of spectrum in xo-array with
c	time increment dt.  The
c	Filter the time series xo of length nleno and time step dt
c	using corner period corper.
c	Result output into array xt.
c
c	The input time series is specified in the frequency domain
c	at freq. samples that are a distance obeta below the real omega axis.
c	The input spectrum is assumed to be sampled in intervals
c	of angular frequency equal to ommin, as specified below.
c
	parameter (maxom=681)
	parameter (iommax=3*maxom/2)
	real*8 s(2048)
	complex*8 xo
	complex*8 xval,ui
	real*4 xt,ommin,cost0,cosdt0,sint0,sindt0
	real*4 facf(iommax)
	real*8 xor(2048),xoi(2048)
	common/spec/xo(iommax),xt(1024)
	common/invtyp/ityp
	ui=cmplx(0.,1.)
	pi=3.14159265
c**	take inverse FT with corner freq. at period [corper] sec.
c	(low pass filter).
	ommin=2.*pi/(dt*real(nleno))
	jmax1=int((2.*pi*2./corper)/real(ommin))
c		write(6,*)'invft: jmax1=',jmax1
	jmax0=jmax1/2
	do j=1,jmax1
	facf(j)=1.0
	if(j.eq.1) facf(j)=0.5
c	Above line: we pick up zero frequency only once, so
c	want it multiplied by (1/(2*pi)) instead of (1/pi) in j-loop below.
	if(j.gt.jmax0) facf(j)=0.5-0.5*cos((real(j-jmax1)/real(
     &	jmax1-jmax0))*(pi))
	xo(j)=xo(j)*(1./pi)*ommin*facf(j)
c	Revise xo-array for velocity seismograms -- mult by ui*(real(j-1)*ommin - ui*obeta)
	if(ityp.eq.1) xo(j)=xo(j)*(ui*(real(j-1)*ommin - ui*obeta))
	xor(j)=dble(real(xo(j)))
	xoi(j)=dble(imag(xo(j)))
	enddo
c		write(6,*)'xor(',100 and 101,')=',xor(100),xor(101)
c	Pad time domain xor and xoi arrays with zeroes up to 2048
	do j=jmax1+1,2048
	xor(j)=0.d0
	xoi(j)=0.d0
	enddo
	call fft(xor,2048)
	call fft(xoi,2048)
	delt=dt*real(nleno)/2048.
	de0=exp(obeta*delt)
	e0=1./de0
	do j=1,1024
        k1=2*j-1
        k2=k1+1
c	Restore the exponential factor arising from integrating
c	below the real omega axis.
	e0=e0*de0
        xt(j)=(xor(k1)+xoi(k2))*e0
        enddo
c* * *
	return
	end

	subroutine fft(s,n)
c  FFT calculates the fast fourier transform of series s and places
c  the result in xr(cos transform) and xi(sine transform).  This is
c  a modification of the routine found in P. Bloomfield's book.
	implicit real*8 (a-h,o-z)
	dimension s(n),xr(2048),xi(2048),ur(15),ui(15)
	do 10 i=1,n
	xi(i) = 0.d0
  10    xr(i) = s(i)
        ur(1) = 0.d0
        ui(1) = 1.d0
        do 110 i=2,15
	ur(i) = dsqrt((1.d0+ur(i-1)) / 2.d0)
  110   ui(i) = ui(i-1) / (2.d0*ur(i))
        n0 = 1
        ii = 0
  140   n0 = n0+n0
        ii = ii+1
        if(n0.lt.n) go to 140
        i1 = n0/2
        i3 = 1
        i0 = ii
        do 260 i4 = 1,ii
        do 250 k = 1,i1
        wr = 1.d0
        wi = 0.d0
        kk = k-1
        do 230 i=1,i0
        if(kk.eq.0) go to 240
        if(mod(kk,2).eq.0) go to 230
        j0 = i0-i
        ws = wr*ur(j0)-wi*ui(j0)
        wi = wr*ui(j0)+wi*ur(j0)
        wr = ws
  230   kk = kk/2
  240   wi = -wi
        l = k
	do 250 j=1,i3
        l1 = l+i1
        zr = xr(l)+xr(l1)
        zi = xi(l)+xi(l1)
        z = wr*(xr(l)-xr(l1))-wi*(xi(l)-xi(l1))
        xi(l1) = wr*(xi(l)-xi(l1))+wi*(xr(l)-xr(l1))
        xr(l1) = z
        xr(l) = zr
        xi(l) = zi
  250   l = l1+i1
        i0 = i0-1
        i3 = i3+i3
  260   i1 = i1/2
        um = 1.d0/dfloat(n0)
	do 310 j=1,n0
        k = 0
        j1 = j-1
	do 320 i=1,ii
        k = 2*k+mod(j1,2)
  320   j1 = j1/2
        k = k+1
        if(k.lt.j) go to 310
        zr = xr(j)
        zi = xi(j)
        xr(j) = xr(k)
        xi(j) = xi(k)
        xr(k) = zr
        xi(k) = zi
  310   continue
	n2=n/2
	do 330 j=1,n2
	k=2*j-1
	s(k)=xr(j)
	s(k+1)=xi(j)
330	continue 
	return
	end

	subroutine addmr(nleno,nrec,phir)
c	Do the summation of transformed displacement
c	over azimuthal order number using
c	the response functions which are stored in 'rec-j-m'
	parameter (ncmax=151)
	parameter (ncmaz=74)
c	Next line is the number of interior GLL points in 1D
	parameter (lmax=4)
c	N is the total number of GLL points in 1D including the points +-1
	parameter (N=lmax+2)
c	Next line is the maximum number of GLL points on the x-side of
c	the global grid times the maximum number on the z-side.
	parameter (nptmax=(N*ncmax-(ncmax-1))*(N*ncmaz-(ncmaz-1)))
	parameter (maxom=681)
	common/grdpar/NX,NZ
	real*8 ang
	complex*16 dsrec,dsrecp
	complex*16 ds1,ds2,ds3
	complex*16 ui,facp,facm
	common/minfo/minc,mmax
	parameter (maxrec=300)
	dimension phir(maxrec)
	common/msumr/dsrec(3*N*N*maxrec,maxom),dsrecp(3*N*N*maxrec,maxom)
	real*8 twopi
c	m0tot=1+2*(mmax/minc)
c
	twopi=2.d0*3.14159265358979d0
	ui=dcmplx(0.,1.)
	open(4,file='rec-j-m',form='unformatted')
	rewind(4)
	do j=1,nleno/3
		write(6,*)'addmr: Doing j=',j,'out of',nleno/3
c	Zero out the dsrec-array and dsrecp-array
	do ic=1,N*N*nrec
	i1=3*ic-2
	i2=i1+1
	i3=i2+1
	dsrec(i1,j)=0.d0
	dsrec(i2,j)=0.d0
	dsrec(i3,j)=0.d0
	dsrecp(i1,j)=0.d0
	dsrecp(i2,j)=0.d0
	dsrecp(i3,j)=0.d0
	enddo
	  do m=minc,mmax+minc,minc
	ic=0
	do nr=1,nrec
	phival=phir(nr)
	  ang=dble(phival)*dble(m-minc)
	  facp=(dcos(ang)+ui*dsin(ang))*dble(minc)/twopi
	  facm=(dcos(ang)-ui*dsin(ang))*dble(minc)/twopi
	do idum=1,N*N
	ic=ic+1
	i1=3*ic-2
	i2=i1+1
	i3=i2+1
	read(4) ds1,ds2,ds3
c		write(6,*)'addmr: ic=',ic,'ds1,ds2,ds3=',ds1,ds2,ds3
	dsrec(i1,j)=dsrec(i1,j)+ds1*facp
	dsrec(i2,j)=dsrec(i2,j)+ds2*facp
	dsrec(i3,j)=dsrec(i3,j)+ds3*facp
	dsrecp(i1,j)=dsrecp(i1,j)+ds1*facp*ui*dble(m-minc)
	dsrecp(i2,j)=dsrecp(i2,j)+ds2*facp*ui*dble(m-minc)
	dsrecp(i3,j)=dsrecp(i3,j)+ds3*facp*ui*dble(m-minc)
	if(m.gt.minc) then
	read(4) ds1,ds2,ds3
c		write(6,*)'addmr again: ic=',ic,'ds1,ds2,ds3=',ds1,ds2,ds3
	dsrec(i1,j)=dsrec(i1,j)+ds1*facm
	dsrec(i2,j)=dsrec(i2,j)+ds2*facm
	dsrec(i3,j)=dsrec(i3,j)+ds3*facm
	dsrecp(i1,j)=dsrecp(i1,j)-ds1*facm*ui*dble(m-minc)
	dsrecp(i2,j)=dsrecp(i2,j)-ds2*facm*ui*dble(m-minc)
	dsrecp(i3,j)=dsrecp(i3,j)-ds3*facm*ui*dble(m-minc)
	endif
	enddo
	enddo
	  enddo
	enddo
	close(4)
	return
	end

	subroutine addm(nleno,pscal,nphi)
c	Do the summation of transformed displacement
c	over azimuthal order number using
c	the response functions which are stored in 'displ-j-m'
	parameter (ncmax=151)
	parameter (ncmaz=74)
c	Next line is the number of interior GLL points in 1D
	parameter (lmax=4)
c	N is the total number of GLL points in 1D including the points +-1
	parameter (N=lmax+2)
c	Next line is the maximum number of GLL points on the x-side of
c	the global grid times the maximum number on the z-side.
	parameter (nptmax=(N*ncmax-(ncmax-1))*(N*ncmaz-(ncmaz-1)))
	parameter (maxom=681)
	  parameter (maxphi=261)
	common/grdpar/NX,NZ
	real*8 ang
	complex*16 ds
	complex*16 ds1,ds2,ds3
	complex*16 ui,facp,facm
	real*8 dlat,plat,plon,phiref
	common/sphpar/dlat,plat,plon,phiref
	common/minfo/minc,mmax
	common/msum/ds(maxphi,3*ncmax*N,maxom)
	real*8 twopi
c	m0tot=1+2*(mmax/minc)
c
	twopi=2.d0*3.14159265358979d0
	ui=dcmplx(0.,1.)
	  nphio2=nphi/2
	open(4,file='displ-j-m',form='unformatted')
	rewind(4)
	do j=1,nleno/3
c	Zero out the ds-array
	do ic=1,NX*N
	i1=3*ic-2
	i2=i1+1
	i3=i2+1
	do iphi=1,maxphi
	ds(iphi,i1,j)=0.d0
	ds(iphi,i2,j)=0.d0
	ds(iphi,i3,j)=0.d0
	enddo
	enddo
	  do m=minc,mmax+minc,minc
		write(6,*)'j=',j,'out of',nleno/3,'m=',m,'out of',mmax+minc
	do ic=1,NX*N
	i1=3*ic-2
	i2=i1+1
	i3=i2+1
	read(4) ds1,ds2,ds3
	  do iphi=1,nphi
	  phival=pscal*real(twopi)/(real(minc))*real(iphi-(nphio2+1))/real(nphio2)
c		write(6,*)'iphi=',iphi,'phival=',phival
	  ang=dble(phival)*dble(m-minc)
	  facp=(dcos(ang)+ui*dsin(ang))*dble(minc)/twopi
	  facm=(dcos(ang)-ui*dsin(ang))*dble(minc)/twopi
	ds(iphi,i1,j)=ds(iphi,i1,j)+ds1*facp
	ds(iphi,i2,j)=ds(iphi,i2,j)+ds2*facp
	ds(iphi,i3,j)=ds(iphi,i3,j)+ds3*facp
	  enddo
	if(m.gt.minc) then
	read(4) ds1,ds2,ds3
	  do iphi=1,nphi
	  phival=pscal*real(twopi)/(real(minc))*real(iphi-(nphio2+1))/real(nphio2)
	  ang=dble(phival)*dble(m-minc)
	  facp=(dcos(ang)+ui*dsin(ang))*dble(minc)/twopi
	  facm=(dcos(ang)-ui*dsin(ang))*dble(minc)/twopi
	ds(iphi,i1,j)=ds(iphi,i1,j)+ds1*facm
	ds(iphi,i2,j)=ds(iphi,i2,j)+ds2*facm
	ds(iphi,i3,j)=ds(iphi,i3,j)+ds3*facm
	  enddo
	endif
	enddo
	  enddo
	enddo
	close(4)
	return
	end

	subroutine readp(nleno)
c	Read in the response functions which are stored in 'vertp-j'
	parameter (ncmax=151)
	parameter (ncmaz=74)
c	Next line is the number of interior GLL points in 1D
	parameter (lmax=4)
c	N is the total number of GLL points in 1D including the points +-1
	parameter (N=lmax+2)
c	Next line is the maximum number of GLL points on the x-side of
c	the global grid times the maximum number on the z-side.
	parameter (nptmax=(N*ncmax-(ncmax-1))*(N*ncmaz-(ncmaz-1)))
	parameter (maxom=681)
	common/grdpar/NX,NZ
	complex*16 ds1,ds2,ds3
	complex*16 dsvert,dsvrtp,dsvrts
	common/vertpr/dsvert(3*nptmax,maxom),dsvrtp(3*nptmax,maxom),dsvrts(3*nptmax,maxom)
	open(4,file='vertp-j',form='unformatted')
	rewind(4)
	igmax=(N*NX-(NX-1))*(N*NZ-(NZ-1))
	do j=1,nleno/3
	do ig=1,igmax
	i1=3*ig-2
	i2=i1+1
	i3=i2+1
	read(4) ds1,ds2,ds3
	dsvert(i1,j)=ds1
	dsvert(i2,j)=ds2
	dsvert(i3,j)=ds3
	read(4) ds1,ds2,ds3
	dsvrtp(i1,j)=ds1
	dsvrtp(i2,j)=ds2
	dsvrtp(i3,j)=ds3
	read(4) ds1,ds2,ds3
	dsvrts(i1,j)=ds1
	dsvrts(i2,j)=ds2
	dsvrts(i3,j)=ds3
c		write(6,*)'readp: dsvert(',i1,j,')=',dsvert(i1,j)
c		write(6,*)'readp: dsvert(',i2,j,')=',dsvert(i2,j)
	enddo
	enddo
	close(4)
	return
	end

      SUBROUTINE SORT(X,Mx,N,Y)
C
C     PURPOSE--THIS SUBROUTINE SORTS (IN ASCENDING ORDER)
C              THE N ELEMENTS OF THE SINGLE PRECISION VECTOR X
C              AND PUTS THE RESULTING N SORTED VALUES INTO THE
C              SINGLE PRECISION VECTOR Y.
C     INPUT  ARGUMENTS--X      = THE SINGLE PRECISION VECTOR OF
C                                OBSERVATIONS TO BE SORTED. 
C                     --Mx     = THE SIZE OF X AND Y
C                     --N      = THE INTEGER NUMBER OF OBSERVATIONS
C                                IN THE VECTOR X. 
C     OUTPUT ARGUMENTS--Y      = THE SINGLE PRECISION VECTOR
C                                INTO WHICH THE SORTED DATA VALUES
C                                FROM X WILL BE PLACED.
C     OUTPUT--THE SINGLE PRECISION VECTOR Y
C             CONTAINING THE SORTED
C             (IN ASCENDING ORDER) VALUES
C             OF THE SINGLE PRECISION VECTOR X.
C     PRINTING--NONE UNLESS AN INPUT ARGUMENT ERROR CONDITION EXISTS. 
C     RESTRICTIONS--THE DIMENSIONS OF THE VECTORS IL AND IU 
C                   (DEFINED AND USED INTERNALLY WITHIN
C                   THIS SUBROUTINE) DICTATE THE MAXIMUM
C                   ALLOWABLE VALUE OF N FOR THIS SUBROUTINE.
C                   IF IL AND IU EACH HAVE DIMENSION K,
C                   THEN N MAY NOT EXCEED 2**(K+1) - 1.
C                   FOR THIS SUBROUTINE AS WRITTEN, THE DIMENSIONS
C                   OF IL AND IU HAVE BEEN SET TO 36,
C                   THUS THE MAXIMUM ALLOWABLE VALUE OF N IS
C                   APPROXIMATELY 137 BILLION.
C                   SINCE THIS EXCEEDS THE MAXIMUM ALLOWABLE
C                   VALUE FOR AN INTEGER VARIABLE IN MANY COMPUTERS,
C                   AND SINCE A SORT OF 137 BILLION ELEMENTS
C                   IS PRESENTLY IMPRACTICAL AND UNLIKELY,
C                   THEN THERE IS NO PRACTICAL RESTRICTION
C                   ON THE MAXIMUM VALUE OF N FOR THIS SUBROUTINE.
C                   (IN LIGHT OF THE ABOVE, NO CHECK OF THE 
C                   UPPER LIMIT OF N HAS BEEN INCORPORATED
C                   INTO THIS SUBROUTINE.)
C     OTHER DATAPAC   SUBROUTINES NEEDED--NONE.
C     FORTRAN LIBRARY SUBROUTINES NEEDED--NONE.
C     MODE OF INTERNAL OPERATIONS--SINGLE PRECISION.
C     LANGUAGE--ANSI FORTRAN. 
C     COMMENT--THE SMALLEST ELEMENT OF THE VECTOR X
C              WILL BE PLACED IN THE FIRST POSITION
C              OF THE VECTOR Y,
C              THE SECOND SMALLEST ELEMENT IN THE VECTOR X
C              WILL BE PLACED IN THE SECOND POSITION
C              OF THE VECTOR Y, ETC.
C     COMMENT--THE INPUT VECTOR X REMAINS UNALTERED.
C     COMMENT--IF THE ANALYST DESIRES A SORT 'IN PLACE',
C              THIS IS DONE BY HAVING THE SAME
C              OUTPUT VECTOR AS INPUT VECTOR IN THE CALLING SEQUENCE. 
C              THUS, FOR EXAMPLE, THE CALLING SEQUENCE
C              CALL SORT(X,N,X)
C              IS ALLOWABLE AND WILL RESULT IN
C              THE DESIRED 'IN-PLACE' SORT.
C     COMMENT--THE SORTING ALGORTHM USED HEREIN
C              IS THE BINARY SORT.
C              THIS ALGORTHIM IS EXTREMELY FAST AS THE
C              FOLLOWING TIME TRIALS INDICATE.
C              THESE TIME TRIALS WERE CARRIED OUT ON THE
C              UNIVAC 1108 EXEC 8 SYSTEM AT NBS
C              IN AUGUST OF 1974.
C              BY WAY OF COMPARISON, THE TIME TRIAL VALUES
C              FOR THE EASY-TO-PROGRAM BUT EXTREMELY
C              INEFFICIENT BUBBLE SORT ALGORITHM HAVE
C              ALSO BEEN INCLUDED--
C              NUMBER OF RANDOM        BINARY SORT       BUBBLE SORT
C               NUMBERS SORTED
C                N = 10                 .002 SEC          .002 SEC
C                N = 100                .011 SEC          .045 SEC
C                N = 1000               .141 SEC         4.332 SEC
C                N = 3000               .476 SEC        37.683 SEC
C                N = 10000             1.887 SEC      NOT COMPUTED
C     REFERENCES--CACM MARCH 1969, PAGE 186 (BINARY SORT ALGORITHM
C                 BY RICHARD C. SINGLETON).
C               --CACM JANUARY 1970, PAGE 54.
C               --CACM OCTOBER 1970, PAGE 624.
C               --JACM JANUARY 1961, PAGE 41.
C     WRITTEN BY--JAMES J. FILLIBEN
C                 STATISTICAL ENGINEERING LABORATORY (205.03)
C                 NATIONAL BUREAU OF STANDARDS
C                 WASHINGTON, D. C. 20234
C                 PHONE--301-921-2315
C     ORIGINAL VERSION--JUNE      1972. 
C     UPDATED         --NOVEMBER  1975. 
C
C---------------------------------------------------------------------
C
      INTEGER X(Mx),Y(Mx)
	INTEGER HOLD,AMED
      DIMENSION IU(36),IL(36) 
C
      IPR=6
C
C     CHECK THE INPUT ARGUMENTS FOR ERRORS
C
      IF(N.LT.1)GOTO50
      IF(N.EQ.1)GOTO55
      HOLD=X(1)
      DO60I=2,N
      IF(X(I).NE.HOLD)GOTO90
60 	CONTINUE
cc      WRITE(IPR, 9)HOLD
      DO61I=1,N
      Y(I)=X(I)
61 	CONTINUE
      RETURN
50 	CONTINUE
cc	WRITE(IPR,15) 
cc      WRITE(IPR,47)N
      RETURN
55 	CONTINUE
cc	WRITE(IPR,18) 
      Y(1)=X(1)
      RETURN
90 	CONTINUE
cc9	FORMAT(1H ,108H***** NON-FATAL DIAGNOSTIC--THE FIRST  INPUT ARGUME)
cc     1NT (A VECTOR) TO THE SORT   SUBROUTINE HAS ALL ELEMENTS = ,E15.8,6
cc     1H *****)
cc15 	FORMAT(1H , 91H***** FATAL ERROR--THE SECOND INPUT ARGUMENT TO THE
cc     1 SORT   SUBROUTINE IS NON-POSITIVE *****)
cc18 	FORMAT(1H ,100H***** NON-FATAL DIAGNOSTIC--THE SECOND INPUT ARGUME
cc     1NT TO THE SORT   SUBROUTINE HAS THE VALUE 1 *****)
cc47 	FORMAT(1H , 35H***** THE VALUE OF THE ARGUMENT IS ,I8   ,6H *****)
C
C-----START POINT-----------------------------------------------------
C
C     COPY THE VECTOR X INTO THE VECTOR Y
      DO100I=1,N
      Y(I)=X(I)
100 	CONTINUE
C
C     CHECK TO SEE IF THE INPUT VECTOR IS ALREADY SORTED
C
      NM1=N-1
      DO200I=1,NM1
      IP1=I+1
      IF(Y(I).LE.Y(IP1))GOTO200
      GOTO250
200 	CONTINUE
      RETURN
250 	M=1 
      I=1 
      J=N 
305 	IF(I.GE.J)GOTO370
310 	K=I 
      MID=(I+J)/2
      AMED=Y(MID)
      IF(Y(I).LE.AMED)GOTO320 
      Y(MID)=Y(I)
      Y(I)=AMED
      AMED=Y(MID)
320 	L=J 
      IF(Y(J).GE.AMED)GOTO340 
      Y(MID)=Y(J)
      Y(J)=AMED
      AMED=Y(MID)
      IF(Y(I).LE.AMED)GOTO340 
      Y(MID)=Y(I)
      Y(I)=AMED
      AMED=Y(MID)
      GOTO340
330 	Y(L)=Y(K)
      Y(K)=TT
340 	L=L-1
      IF(Y(L).GT.AMED)GOTO340 
      TT=Y(L)
350 	K=K+1
      IF(Y(K).LT.AMED)GOTO350 
      IF(K.LE.L)GOTO330
      LMI=L-I
      JMK=J-K
      IF(LMI.LE.JMK)GOTO360
      IL(M)=I
      IU(M)=L
      I=K 
      M=M+1
      GOTO380
360 	IL(M)=K
      IU(M)=J
      J=L 
      M=M+1
      GOTO380
370 	M=M-1
      IF(M.EQ.0)RETURN
      I=IL(M)
      J=IU(M)
380 	JMI=J-I
      IF(JMI.GE.11)GOTO310
      IF(I.EQ.1)GOTO305
      I=I-1
390 	I=I+1
      IF(I.EQ.J)GOTO370
      AMED=Y(I+1)
      IF(Y(I).LE.AMED)GOTO390 
      K=I 
395 	Y(K+1)=Y(K)
      K=K-1
      IF(AMED.LT.Y(K))GOTO395 
      Y(K+1)=AMED
      GOTO390
      END 
