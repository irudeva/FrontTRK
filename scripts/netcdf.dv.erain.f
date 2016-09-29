      program netsdf2cmp

	include '/opt/local/include/netcdf.inc'

        character*120 Dnc,ncfile2,ncfile,cmpfile,command,head*80
	character hs*2,uvar*4,vvar*4,dvvarcmp*5,lmon*12
	real miss

	integer ncid,errorcode,crit
	character dimname*10,time(4)*4,time0*4,lev*3
	integer londimid,latdimid,lonvarid,latvarid,uvarid,vvarid
	integer lonsize,latsize,slat,elat
	integer startcoo(3),countcoo(3)
	real, allocatable:: longitudes(:),latitudes(:)
	real,allocatable:: v(:,:),v1(:,:),u(:,:),u1(:,:),dv(:,:)
	integer*2, allocatable:: su(:,:),sv(:,:)
 	real*8 uscale, uoffset,vscale,voffset
	integer*2 missv

	integer y,m,d,sy,sm,sd,ey,em,ed,y0
	character*2 chm,chd,chm0,chd0
	integer system,status

	integer t,th,i,j,t1,t2,j1,j2,j0

        data time/"0000","0600","1200","1800"/
	crit=0


	read (*,'(a)') hs
        read (*,'(a)') lev
	read (*,'(i4,2i2)') sy,sm,sd
	read (*,'(i4,2i2)') ey,em,ed
	read (*,'(i4)') ly
	read (*,'(a)') lmon
	read (*,'(a)') uvar
	read (*,'(a)') vvar
	read (*,'(a)') Dnc
	read (*,'(a)') cmpfile
	read (*,'(a)') dvvarcmp
	read (*,*) miss
	write(*,'(''CMP file: '',a)') cmpfile
	open (20,file=cmpfile,form='unformatted')

	do 1 y=sy,ey;print*, y
	 if (mod(y,4)==0)then
	 open(10,file='dates366.dat',action='read')
	 else
	 open(10,file='dates365.dat',action='read')
	 endif


        write(ncfile,'(a,i4,".nc")')Dnc(:len_trim(Dnc)),y
        if(y.gt.ly)then
         write(ncfile,'(2a,i4,".nc")')Dnc(:len_trim(Dnc)),trim(lmon),y
	end if


	write(command,'("ls ",a,".gz")')ncfile(:len_trim(ncfile))
	status=system(command)
	print*, command
	if(status==0)then
	 write(command,'("gunzip ",a,".gz")')ncfile(:len_trim(ncfile))
	STATUS=SYSTEM(command)
	 if(STATUS==1)write(*,'("Gunzip failed")')
	endif

	ncid=ncopn(ncfile,ncnowrit,rcode)
        print*, "NCOPN: file open ", ncfile

!*******first step ****************************************************
	if(y==sy)then
	LONDIMID=NCDID(ncid,'lon',ERRORCODE)
	LATDIMID=NCDID(ncid,'lat',ERRORCODE)

	call ncdinq(ncid,londimid,dimname,lonsize, errorcode)
	call ncdinq(ncid,latdimid,dimname,latsize, errorcode)
	print*, "ok"
	allocate (longitudes(lonsize),latitudes(latsize))
	if(latsize==121)then
	if(hs=='0N')then
	 slat=1
	 elat=(latsize-1)/2+1  !!!!for SH!!!!
	 j0=1
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


	if(hs=="0N")then
	 j1=slat
	 j2=elat
	 j0=1
	elseif(hs=="0S")then
	 j1=elat
	 j2=slat
	 j0=-1
	endif


	lonvarid=ncvid(ncid,'lon',errorcode)
	latvarid=ncvid(ncid,'lat',errorcode)
	uvarid = ncvid(ncid,uvar,errorcode)
	vvarid = ncvid(ncid,vvar,errorcode)

	print*, "ok"
	startcoo(1)=1
	countcoo(1)=lonsize
        call ncvgt(ncid,lonvarid,startcoo(1),countcoo(1),longitudes,
     &		   errorcode)
	startcoo(1)=1
	countcoo(1)=latsize
        call ncvgt(ncid,latvarid,startcoo(1),countcoo(1),latitudes,
     &		   errorcode)


	 write(20)elat-slat+1
	 write(21,'(i10)')elat-slat+1
	 write(20)(latitudes(i),i=j1,j2,j0)
	 write(21,'(61f10.2)')(latitudes(i),i=j1,j2,j0)

	 if(all(longitudes>=0.).and.all(longitudes<360.))then
	 write(20)lonsize
	 write(20)(longitudes(i),i=1,lonsize)
	 write(21,'(i10)')lonsize
	 write(21,'(540f10.2)')(longitudes(i),i=1,lonsize)
	elseif(all(longitudes>=-180.).and.all(longitudes<180.))then
	 write(20)lonsize
	 write(20)(longitudes(i),i=lonsize/2+1,lonsize),
     &            (longitudes(i)+360.,i=1,lonsize/2)
	 write(21,'(i10)')lonsize
	 write(21,'(240f10.2)')(longitudes(i),i=lonsize/2+1,lonsize),
     &		  (longitudes(i)+360.,i=1,lonsize/2)
	else
	 write(*,'("ERROR: Check longitudes")')
	 stop
	endif

	print*, "ok"
	allocate(u(lonsize,latsize),v(lonsize,latsize)
     &          ,su(lonsize,latsize),sv(lonsize,latsize)
     &          ,u1(lonsize,latsize),v1(lonsize,latsize)
     &          ,dv(lonsize,latsize))

        startcoo(1)=1
	countcoo(1)=lonsize
        startcoo(2)=1
	countcoo(2)=latsize
	countcoo(3)=1
	endif   !y==sy
!**************************************************************************

	call ncagt(ncid,uvarid,"scale_factor",uscale,errorcode)
	call ncagt(ncid,uvarid,"add_offset",uoffset,errorcode)
	print*, "uscale=",uscale
	print*, "uoffset=",uoffset
	call ncagt(ncid,vvarid,"scale_factor",vscale,errorcode)
	call ncagt(ncid,vvarid,"add_offset",voffset,errorcode)
	print*, "vscale=",vscale
	print*, "voffset=",voffset

	call ncagt(ncid,uvarid,'_FillValue',missv,
     &              errorcode)

	t1=1
	t2=365
	 if (mod(y,4)==0)t2=366
	if(y==sy)then
	t1=0
	 loop_a: do
	  read(10,'(2i2)',end=9)m,d
	  t1=t1+1
	  if(m==sm.and.d==sd)then
	   backspace(10)
	   exit loop_a
	  endif
	 enddo loop_a
	elseif(y==ey)then
	t2=0
	 loop_b: do
	  read(10,'(2i2)',end=9)m,d
	  t2=t2+1
	  if(m==em.and.d==ed)then
	   rewind(10)
	   exit loop_b
	  endif
	 enddo loop_b
	endif

	print*, y,t1,t2

	do 11 t=t1,t2
	 read(10,'(2a)')chm,chd
	 read(chm,'(i2)')m
	 read(chd,'(i2)')d
	do 13 th=1,4

!        startcoo(3)=(t-1)*4+th
!	call ncvgt(ncid,uvarid,startcoo,countcoo,su,
!     &              errorcode)
!	call ncvgt(ncid,vvarid,startcoo,countcoo,sv,
!     &              errorcode)
!	u=float(su)*uscale+uoffset
!	v=float(sv)*vscale+voffset
!	do j=j1,j2,j0
!	do i=1,lonsize
!	if(su(i,j)==missv)u(i,j)=miss
!	if(sv(i,j)==missv)v(i,j)=miss
!	enddo
!	enddo


!	!write(25,'(4f10.2)') ((longitudes(i),latitudes(j),u(i,j),v(i,j),
!	write(25,'(2f10.2,2e14.6)')
!     & ((longitudes(i),latitudes(j),u(i,j),v(i,j),
!     & i=1,lonsize),j=slat,elat)



!tmpppppp
!        write(ncfile2,'("../../merra/MERRA",i3
!     *  ".prod.assim.inst6_3d_ana_Np.",i4,2a,".SUB.nc")')a,y,chm,chd

!	ncid2=ncopn(ncfile2,ncnowrit,rcode)
!	call ncvgt(ncid2,uvarid,startcoo,countcoo,u1,
!     &              errorcode)


!	do j=slat,elat
!	do i=1,lonsize
!	if(longitudes(i)<0.)then
!	write(22,'(3f10.2)')longitudes(i)+360.,latitudes(j),u1(i,j)/100.
!	else
!	write(22,'(3f10.2)')longitudes(i),latitudes(j),u1(i,j)/100.
!	endif
!	enddo
!	enddo
!endtmp

!******end first step *************************************************


        startcoo(3)=(t-1)*4+th;print*, startcoo(3)
	call ncvgt(ncid,uvarid,startcoo,countcoo,su
     &              ,errorcode)
	call ncvgt(ncid,vvarid,startcoo,countcoo,sv
     &              ,errorcode)

	u1=float(su)*uscale+uoffset
	v1=float(sv)*vscale+voffset
	do j=j1,j2,j0
	do i=1,lonsize
	if(su(i,j)==missv)u1(i,j)=miss
	if(sv(i,j)==missv)v1(i,j)=miss
	enddo
	enddo

	print*, u1(1,j1)

	if(crit/=0)then
      dv=miss
	do 15 j=slat,elat
	do 15 i=1,lonsize
	 if (u(i,j)>=0..and.u1(i,j)>=0.)then
	  if(hs=='0S'.and.v(i,j)<=0..and.v1(i,j)>0.)then
           if(v1(i,j)/=miss.and.v(i,j).ne.miss)
     &	 		    dv(i,j)=v1(i,j)-v(i,j)
	  endif
	  if(hs=='0N'.and.v(i,j)>=0..and.v1(i,j)<0.)then
           if(v1(i,j)/=miss.and.v(i,j).ne.miss)
     &                      dv(i,j)=v(i,j)-v1(i,j)
	  endif
	 endif
15	continue

!*******head **********************************************************
	write(head,'(a5,25x,a5,5x,i4,2a,x,a4,a7,10x,a10)')dvvarcmp
     &         ,"ERAIN",y0,chm0,chd0,time0,"m/s","1.5x1.5DEG"
	write(20) head
	write(21,'(a80)') head
	write(*,*) head
	print*, y0,chm0,chd0,time0
!*******end head ******************************************************

!tmp
!	if(t==1.and.y==2009)then
!	print*, 'start writing 23', hs
!	if(hs=='0S')print*,'hs is OK'
!	do j=slat,elat
!	do i=1,lonsize
!	if(longitudes(i)<0.)then
!	write(23,'(6f10.2)')longitudes(i)+360.,latitudes(j),
!     &    v(i,j),v1(i,j),v(i,j)-v1(i,j),dv(i,j)
!	else
!	write(23,'(6f10.2)')longitudes(i),latitudes(j),
!     &    v(i,j),v1(i,j),v(i,j)-v1(i,j),dv(i,j)
!	endif
!	enddo
!	enddo
!	stop
!
!endtmp
!	write(23,'(7f10.2)')((longitudes(i),latitudes(j),u(i,j),u1(i,j),
!     & v(i,j),v1(i,j),dv(i,j),i=1,lonsize),j=j1,j2,j0)
!	stop
!	endif

	 write(20)((dv(i,j),i=1,lonsize),j=j1,j2,j0)
	 write(21,'(240f10.2)')((dv(i,j),i=1,lonsize),j=j1,j2,j0)
!	 write(22,'(240f10.2)')((v(i,j),i=lonsize/2+1,lonsize),j=j1,j2,j0)
!	 write(23,'(240f10.2)')((v1(i,j),i=lonsize/2+1,lonsize),j=j1,j2,j0)
	endif !crit/=0

	 u=u1
	 v=v1
	 y0=y
	 chm0=chm
	 chd0=chd
	 time0=time(th)
	 crit=1


13	enddo !hr

	if(y==ey.and.m==em.and.d==ed)then
	print*, "current  ", y,m,d
	print*, "end  ",ey,em,ed
	print*,"stopped"
	goto 88
	endif

11	enddo !time

88	close(10)
	call ncclos(ncid,errorcode)

!	 write(command,'("gzip ",a)')ncfile(:len_trim(ncfile))
!	print*, command
!	 STATUS=SYSTEM(command)
!	 if(STATUS==1)write(*,'("Gzip failed")')

1	enddo !year

	close(20)
	!close(21,status='delete')
	close(21)


8	stop

9      write(*,'("ERROR: Check dates!!!!")')

	end
