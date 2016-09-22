      program netsdf2cmp
	
!	include '/work/irudeva/include/netcdf.inc'
	include '/usr/local/include/netcdf.inc'
        
        character*120 ncfile,cmpfile,command,head*80
	character hs*2,uvar*4,vvar*4
	character*4 ouvar,ovvar,wvar
	real miss

	integer ncid,errorcode
	character dimname*10,time(4)*4
	integer londimid,latdimid,lonvarid,latvarid,uvarid,vvarid
	integer lonsize,latsize,slat,elat
	integer startcoo(3),countcoo(3)
	real*4, allocatable:: longitudes(:),latitudes(:)
	real,allocatable:: v(:,:),avv(:,:,:),u(:,:),avu(:,:,:),w(:,:)
     *                     ,w12(:,:,:,:)
        integer*2,allocatable:: sv(:,:),su(:,:)
	real n(12)
	real*8 uscale,uoffset, vscale,voffset
	integer scaleid, iffsetid, atype, alen


	integer y,m,sy,sm,ey,em,t,m1,m2
	character*2 chm,chd
	integer system,status

	integer i,j,nw

	read (*,'(a)') ncfile
	read (*,'(a)') hs
	read (*,'(i4,i2)') sy,sm;write(*,'("start: ",i4,i2)') sy,sm
	read (*,'(i4,i2)') ey,em;write(*,'("end:   ",i4,i2)')ey,em
	read (*,'(a)') uvar
	read (*,'(a)') vvar
	read (*,'(a)') cmpfile
	read (*,'(a)') ouvar
	read (*,'(a)') ovvar
	write(*,'(''CMP file: '',a)') cmpfile	
	open (20,file=cmpfile,form='unformatted')	



	ncid=ncopn(ncfile,ncnowrit,rcode)
	
	LONDIMID=NCDID(ncid,'longitude',ERRORCODE)
	LATDIMID=NCDID(ncid,'latitude',ERRORCODE)

	call ncdinq(ncid,londimid,dimname,lonsize, errorcode)
	call ncdinq(ncid,latdimid,dimname,latsize, errorcode)
	if(latsize==121)then
	if(hs=='0N')then
	 slat=1
	 elat=(latsize-1)/2+1  !!!!for NH!!!!
	elseif(hs=='0S')then
	 slat=(latsize-1)/2+1 
	 elat=latsize 
	else
	 write(*,'("ERROR: Check hs")')
	endif
	else
	 slat=1
	 elat=latsize
	endif

	
	lonvarid=ncvid(ncid,'longitude',errorcode)
	latvarid=ncvid(ncid,'latitude',errorcode)
	uvarid = ncvid(ncid,uvar,errorcode)
	vvarid = ncvid(ncid,vvar,errorcode)

	call ncainq(ncid,uvarid,"scale_factor",atype,
     &             alen,errorcode)
	if(atype.ne.6)then
	 write(*,'("!!!ERROR: check type of scale_factor!!!" )')
	 stop
	endif
	call ncainq(ncid,uvarid,"add_offset",atype,
     &             alen,errorcode)
	if(atype.ne.6)then
	 write(*,'("!!!ERROR: check type of add_offset!!!" )')
	 stop
	endif

	call ncagt(ncid,uvarid,"scale_factor",uscale,errorcode)
	call ncagt(ncid,uvarid,"add_offset",uoffset,errorcode)
	print*, "uscale=",uscale
	print*, "uoffset=",uoffset
	call ncagt(ncid,vvarid,"scale_factor",vscale,errorcode)
	call ncagt(ncid,vvarid,"add_offset",voffset,errorcode)
	print*, "vscale=",vscale
	print*, "voffset=",voffset
	
	allocate(longitudes(lonsize),latitudes(latsize))
	startcoo(1)=1
	countcoo(1)=lonsize
        call ncvgt(ncid,lonvarid,startcoo(1),countcoo(1),longitudes,
     &		   errorcode)
	startcoo(1)=1
	countcoo(1)=latsize
        call ncvgt(ncid,latvarid,startcoo(1),countcoo(1),latitudes,
     &		   errorcode)

	print*, "start out"
!	if(all(latitudes<=0.))then
	 write(20)elat-slat+1
	 write(21,'(i10)')elat-slat+1

	if(hs=="0N")then
	 j1=slat
	 j2=elat
	 j0=1
	elseif(hs=="0S")then
	 j1=elat
	 j2=slat
	 j0=-1
	endif
	

	  write(20)(latitudes(i),i=j1,j2,j0)
	  write(21,'(100f10.2)')(latitudes(i),i=j1,j2,j0)

	if(all(longitudes>=0.).and.all(longitudes<360.))then
	 write(20)lonsize
	 write(20)(longitudes(i),i=1,lonsize)
	 write(21,'(i10)')lonsize
	 write(21,'(240f10.2)')(longitudes(i),i=1,lonsize)
	elseif(all(longitudes>=-180.).and.all(longitudes<180.))then
	 write(*,'("ERROR: Check longitudes")')
	 stop
	else
	 write(20)lonsize
	 write(20)(longitudes(i),i=1,lonsize)
	 write(21,'(i10)')lonsize
	 write(21,'(240f10.2)')(longitudes(i),i=1,lonsize)
	endif 
	

	t=0;n=0.
         do 1 y=sy,ey	

	m1=1;m2=12
	if(y==sy)m1=sm
	if(y==ey)m2=em
	
	do 1 m=m1,m2
	t=t+1;print*, y,m,t


	if(y==sy.and.m==sm)then  !first field
	
	allocate(u(lonsize,latsize),v(lonsize,latsize)
     &          ,su(lonsize,latsize),sv(lonsize,latsize)
     &          ,avu(lonsize,latsize,12),avv(lonsize,latsize,12)
     &          ,w(lonsize,latsize),w12(lonsize,latsize,2,12))

        startcoo(1)=1
	countcoo(1)=lonsize
        startcoo(2)=1
	countcoo(2)=latsize
	countcoo(3)=1

	print*,"end 1"	
	endif  !for the first time step
        
	startcoo(3)=t

!******end first step *************************************************
	

	call ncvgt(ncid,uvarid,startcoo,countcoo,su
     &              ,errorcode)
	call ncvgt(ncid,vvarid,startcoo,countcoo,sv
     &              ,errorcode)
	u=float(su)*uscale+uoffset
	v=float(sv)*vscale+voffset
	
	avu(:,:,m)=avu(:,:,m)+u
	avv(:,:,m)=avv(:,:,m)+v
	n(m)=n(m)+1.

1	continue
	call ncclos(ncid,errorcode)

	do 2 m=1,12

	write(chm,'(i2)') m
	if(chm(1:1)==' ')chm(1:1)="0"


	do 2 nw=1,2
	if(nw==1)then
	 w(:,:)=avu(:,:,m)/n(m)
	 wvar=ouvar
	else
	 w(:,:)=avv(:,:,m)/n(m)
	 wvar=ovvar
	endif
!*******head **********************************************************
	write(head,'(a5,"clim",21x,a5,5x,4x,a,7x,a7,10x,a10)')wvar
     &         ,"ERAIN",chm,"m/s","1.5x1.5DEG"
	write(20) head
	write(21,'(a80)') head
	write(*,*) head
!*******end head ******************************************************
!	if(m==1.and.nw==1)then
!	write(22,'(3f10.2)') ((longitudes(i),latitudes(j),
!     &           w(i,j),i=1,lonsize),j=1,latsize)
!	endif


	 write(20)((w(i,j),i=1,lonsize),j=j1,j2,j0)	
	 write(21,'(240f10.2)')((w(i,j),i=1,lonsize),j=j1,j2,j0)
	w12(:,:,nw,m)=w
		
2	continue

	!close(21,status='delete')
	close(21)
	close(20)

!!!!climatology check!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	open(22,file="u500.erain.clim."//hs//".dat",action='write')
	open(23,file="v500.erain.clim."//hs//".dat",action='write')
	

	 do j=j1,j2,j0;print*, j	
	 do i=1,lonsize
	
	write(22,'(14f10.2)')longitudes(i),latitudes(j)
     &                        ,(w12(i,j,1,m),m=1,12)
	write(23,'(14f10.2)')longitudes(i),latitudes(j)
     &                        ,(w12(i,j,2,m),m=1,12)
			
	 enddo
	 enddo

	close(22)
	close(23)
	end
