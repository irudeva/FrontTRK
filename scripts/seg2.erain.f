      module mcoor
      parameter(mlon= 300, mlat= 100) ! Grid size
! mlat=180./dlat + 1 (73 for NCEP1)
! mlon=360./dlon + 1 (145 for NCEP1)
      end module mcoor
      
      program seg2
c
c * Notes
c
c (1) 1/4/2008: Added -1 option: If there are just two points in an object and these are on the same latitude
c     then replace by the most eastward (single) point
c (2) 8/4/2008: Added dv info in seg1
c (3) 1/12/2009: Corrected determination of most eastward pixel for the case of a frontal object straddling 
c     the Greenwich meridian. For a sequence of longitude indices e.g. 239 240 1 2 where 241 <--> 360 deg.
c     the most eastward pixel was set to 240; now it is 2.
c (4)  Compilation:
c      g77 -ffortran-bounds-check -o pc/seg2.erain seg2.erain.f ./velleman_and_hoaglin/libeda.a
c
      
      use mcoor
      real,parameter:: pi=3.14159265
      parameter (maxt=10000)
      integer t
      parameter (maxi=3000)
      integer ilon(maxi),ilat(maxi)
      integer ilon2(maxi),ilat2(maxi)
      real dv(maxi)
      real dv2(maxi)
      character*2 chm,chd,hs,cdh*4
c
      integer islon(maxi)
      real slon(maxi)
      real xrlons(maxi),xrlats(maxi),cglat,cglon,xrlonsn(maxi)
      real sx(maxi),sy(maxi),x,y
      real rx(maxi),ry(maxi)
c
      integer jslon(maxi), jslat(maxi)
      real sdv(maxi)
c
      real xlons(mlon),xlats(mlat)
c
      character*80 infile,outfile
      character*80 head,ch80,chrec
      character*80 optarg
      integer err
c
       if(iargc().eq.0)then
        write(*,*)'Usage: seg2 -i infile -o outfile'
	write(*,*)'Output: fort.4 (fort.3: diagnostic)'
        stop
      else
        narg= iargc()
	i=1
 	do while (i.le.narg)
         call getarg(i,optarg)
	 if(optarg.eq.'-i')then
          i=i+1
          call getarg(i,optarg)
          infile=optarg
	 elseif(optarg.eq.'-o')then
          i=i+1
          call getarg (i,outfile)
	 elseif(optarg.eq.'-h')then
          i=i+1
          call getarg (i,hs)
         endif
         i=i+1
	enddo
	write(*,*)'Input file: ',infile(:ilen(infile))
	write(*,*)'Output file: ',outfile(:ilen(outfile))
      endif
      open (1,file=infile)
      open (4,file=outfile)
c
c * *** For NCEP/NCEP2/ERA40/JRA25 reanalysis data
c
ckk      nlats= 73
ckk      do i=1,nlats
ckk        xlats(i)= -90 + (i-1)*2.5
ckk      enddo
ckk      nlons= 145
ckk      do i=1,nlons
ckk        xlons(i)= (i-1)*2.5
ckk      enddo
c
c * ERA-Interim global grid
c
	read(1,'(7x,i5)') nlons
	read(1,'(7x,i5)') nlats
	write(*,'("Nlats: ",i5,"Nlons: ", i5)') nlats,nlons

	if(nlats>mlat.or.nlons>mlon)then
	write(*,'("ERROR: Check coordinates!!!!")')
	stop
	endif
ci      nlats= mlat
ci      nlons= mlon
ci      dlat= 1.5
ci      dlon= 1.5
	dlat=90./float(nlats-1)
	dlon=360./float(nlons)
      do i=1,nlats
	if(hs=='0N')then
        xlats(i)= 90. - (i-1)*dlat
	elseif(hs=='0S')then
        xlats(i)= -90. + (i-1)*dlat
	endif
      enddo
      do i=1,nlons
        xlons(i)= (i-1)*dlon
      enddo
	write(*,*)nlats, xlats(1:nlats)
	write(*,*)nlons,xlons(1:nlons)
c*******************************************************
c
      do t=1,maxt
        write(*,*)'t=',t
        read(1,'(A)',end=9)chrec
	read(chrec,'(5x,I4,1x,I4,a2,a2,1x,a4)')irec,iy,chm,chd,chh
	write(*,'(''Rec: '',I4,1x,I4,a2,a2,1x,a4)')t,iy,chm,chd,chh
ci	write(3,'(''Rec: '',I4,1x,I4,I2,I2,1x,I4)')t,iy,im,id,ih
	write(4,'(''Rec: '',I4,1x,I4,a2,a2,1x,a4)')t,iy,chm,chd,chh
	read(1,'(6x,I4)')nf
	write(*,'(''No.F: '',I4)')nf
ci	write(3,'(''No.F: '',I4)')nf
	write(4,'(''No.F: '',I4)')nf
       do k=1,nf
         read(1,'(7x,I4,7x,I4)')kf,nkf
ci	 write(3,'(''Front: '',I4,1x,''No.P: '',I4)')kf,nkf
	 write(*,'(''Front: '',I4,1x,''No.P: '',I4)')kf,nkf
	 do j=1,nkf
	   read(1,'(3(I4,1x),18x,1x,F8.3)')jd,ilat(j),ilon(j),dv(j)
ci	   write(3,'(3(I4,1x),'':'',F8.3)')j,ilat(j),ilon(j),dv(j)
	 enddo
c
c * Get most eastward point at each latitude
c
         if(nkf.eq.1)then
	   ilat2(1)= ilat(1)
	   if(ilon(1).eq.nlons)then
	     ilon2(1)= 2
	   else
	     ilon2(1)= ilon(1) +1
	   endif
	   dv2(1)= dv(1)
	   nkf2= 1
           icross0= 0
c
	 else  ! Fronts where nkf >= 2
           icross0= 0
c
c * Determine most eastern gridpoint for each latitude
c
c * Check for longitude 360 -> 0
c
c * Add an extra (dummy) row to list to assist with loop processing
            ilat(nkf+1)= -999
	    ilon(nkf+1)= -999
	    dv(nkf+1)= -999.
c
	     klat0= ilat(1)
	     klon0= ilon(1)
	     dv0= dv(1)
	     kk= 0
	     jj= 1
	     jslat(jj)= klat0
	     jslon(jj)= klon0
	     sdv(jj)= dv0
             do j=2,nkf+1
	       klat= ilat(j)
	       klon= ilon(j)
	       dvj= dv(j)
	       if(klat.eq.klat0)then
	         jj= jj +1
	         jslat(jj)= klat
	         jslon(jj)= klon
	         sdv(jj)= dvj
	       else
	         kk= kk +1
	         njj= jj
ci                 write(14,*)'njj: ',njj,' ',(jslon(jj),jj=1,njj)
                 jlonr= abs(jslon(njj)-jslon(1))
	         if(jslon(njj).eq.nlons-1.or.jlonr.gt.0.5*nlons)then ! nlons corresponds to 360
ci	           write(12,'(''Rec: '',I4,1x,I4,I2,I2,1x,I4)')
ci     &                   t,iy,im,id,ih
ci                   write(12,*)'njj: ',njj,' ',(jslon(jj),jj=1,njj)
	           do jj=1,njj
		     if(jslon(jj).lt.0.5*nlons)then
		       jslon(jj)= jslon(jj) +nlons -1
		     endif
	           enddo
ci                   write(12,*)'Adj njj: ',njj,' ',(jslon(jj),jj=1,njj)
		   if(njj.gt.1)then
                     call indexi (njj,jslon,islon)
	             jn= islon(njj)
		   else
	             jn= 1
		   endif
c * jn is index of max. (eastward) longitude
	           ilat2(kk)= jslat(jn)
		   ilon2(kk)= jslon(jn)
		   if(ilon2(kk).gt.nlons) then
		     ilon2(kk)= ilon2(kk) -nlons +1
		   endif
		   dv2(kk)= sdv(jn)
	         else
		   ilat2(kk)= jslat(njj) 
		   ilon2(kk)= jslon(njj)
		   dv2(kk)= sdv(njj)
	         endif
	         jj= 1
		 klat0= klat
	         jslat(jj)= klat
	         jslon(jj)= klon
	         sdv(jj)= dvj
	       endif
	     enddo
	     nkf2= kk
	 endif

ci	 do j=1,nkf2
ci	   write(3,'(3(I4,1x))')j,ilat2(j),ilon2(j)
ci	 enddo
c
c * For case of front at t+1, front should be one pixel to east
c
         if(nkf.gt.1)then  ! If original no. lat. points >1
	 do j=1,nkf2
	   ilon2a= ilon2(j)
	   if (ilon2a.eq.nlons)then
	     ilon2a= 2  ! Assumes 1 == nlons so next lon is 2
	   else
	     ilon2a= ilon2a +1
	   endif
	   ilon2(j)= ilon2a
ci	   write(3,'(3(I4,1x))')j,ilat2(j),ilon2(j)
	 enddo
	 endif
c

c####
c
         if(nkf2.ge.7)then  ! Need N >= 7 for rsm
      do j=1,nkf2
        slon(j)= real(ilon2(j))
      enddo
      call indexx (nkf2,slon,islon)
      lon_min= slon(islon(1))
      lon_max= slon(islon(nkf2))
      if(mod(nkf2,2).eq.0)then
        kmed= nkf2/2
      else
        kmed= (nkf2+1)/2
      endif
      lon_med= slon(islon(kmed))
ci      write(3,*)'lon_min: ',lon_min,' lon_max: ',lon_max,
ci     * ' lon_med: ',lon_med
      icross0= 0
      if(lon_min.lt.0.2*nlons.and.
     *      lon_max.gt.0.8*nlons)then
        icross0= 1
      endif
ci      write(3,*)'icross0: ',icross0
c
	 do j=1,nkf2
	   if(icross0.eq.1)then
	     xrlats(j)= xlats(ilat2(j))
	     if(xlons(ilon2(j)).gt.180.and.xlons(ilon2(j)).le.360)then
	       xrlons(j)= xlons(ilon2(j)) -360.
		print*, j,ilon2(j),xlons(ilon2(j)),xrlons(j)
	     else
	       xrlons(j)= xlons(ilon2(j))
	     endif
	   else
	     xrlons(j)= xlons(ilon2(j))
	     xrlats(j)= xlats(ilat2(j))
	   endif
ckk	   dv2(j)= dv(j)
ci	   write(3,'(3(I4,1x),2(F7.2,1x))')j,ilat2(j),ilon2(j),
ci     * xrlats(j),xrlons(j)
	 enddo
	 do j=1,nkf2
	   sx(j)= 0.
	   rx(j)= 0.
	   sy(j)= 0.
	   ry(j)= 0.
	 enddo
	 call rsm (xrlons,nkf2,sx,rx,2,err)
	 if(err.gt.0)write(*,*)'ERROR in rsm lons: ',err
	 call rsm (xrlats,nkf2,sy,ry,2,err)
	 if(err.gt.0)write(*,*)'ERROR in rsm lats: ',err
ci	 open(2,file='j.grz')
ci	 write(2,*)'Raw'
ci	 write(2,*)nkf2
ci	 do i=1,nkf2
ci	   write(2,*)xrlons(i),xrlats(i)
ci	 enddo
ci	 write(2,*)'Smooth'
ci	 write(2,*)nkf2
ci	 do i=1,nkf2
ci	   write(2,*)sx(i),sy(i)
ci	 enddo
ci	 close(2)
ci	 write(3,*)'After rsm'
	 do j=1,nkf2
	   if(icross0.eq.1)then
	     xrlats(j)= sy(j)
	     if(sx(j).lt.0)then
	       xrlons(j)= sx(j) +360.
	     else
	       xrlons(j)= sx(j)
	     endif
	   else
	     xrlons(j)= sx(j)
	     xrlats(j)= sy(j)
	   endif
ci	   write(3,'(3(I4,1x),2(F7.2,1x))')j,ilat2(j),ilon2(j),
ci     * xrlats(j),xrlons(j)
	 enddo
	 else
c * N < 7 : Can't do rsm
	 do j=1,nkf2
	   xrlons(j)= xlons(ilon2(j))
	   xrlats(j)= xlats(ilat2(j))
	 enddo
c
	 endif
c
c * Check for crossing lon 0 == 360
c
         if(nkf2.gt.1)then
           call indexx (nkf2,xrlons,islon)
           xlon_min= xrlons(islon(1))
           xlon_max= xrlons(islon(nkf2))
	   icross0= 0
	   if((xlon_max-xlon_min).ge.180.)then
	     icross0= 1
	   endif
	 endif
c
c * Check for lons >360 or lons <0 due to smoothing
c
         if(nkf2.ge.7)then
           call indexx (nkf2,xrlons,islon)
           xlon_min= xrlons(islon(1))
           xlon_max= xrlons(islon(nkf2))
	   if(xlon_max.gt.360.)then
	     icross0= 1
	     do j=1,nkf2
	       if(xrlons(j).gt.360.)then
	         xrlons(j)= xrlons(j) -360.
	       endif
	     enddo
	   endif
	   if(xlon_min.lt.0.)then
	     icross0= 1
	     do j=1,nkf2
	       if(xrlons(j).lt.0.)then
	         xrlons(j)= xrlons(j) +360.
	       endif
	     enddo
	   endif
	 endif

c  cglat/cglon

        xrlonsn=xrlons
	if(icross0.eq.1)then
        do j=1,nkf2
         if(xrlonsn(j)>180.)then
          xrlonsn(j)=xrlonsn(j)-360.
         end if
        end do
	end if

        cglat=sum(xrlats(1:nkf2))/real(nkf2)

        y=sum(cos(pi*xrlats(1:nkf2)/180.))
        x=sum(cos(pi*xrlats(1:nkf2)/180.)*xrlonsn(1:nkf2))
        cglon=x/y
        if(cglon<0.)cglon=cglon+360.
c  end cglon/cglat
c
	 if(t==1.and.icross0/=1)write(9,'(i5,'',1," "'')')nkf2
	 write(4,'(''Front: '',I4,1x,''No.P: '',I4,'' Cross0: '',I1,
     *''  CG: '',2f7.2)')
     * kf,nkf2,icross0,cglat,cglon
	 do j=1,nkf2
	 if(t==1.and.icross0/=1)write(9,'(3f10.2)') 
     & xrlons(j),xrlats(j),dv2(j)
	   write(4,'(3(I4,1x),2(F7.2,1x),F8.3)')j,ilat2(j),ilon2(j),
     * xrlats(j),xrlons(j),dv2(j)
	 enddo
c####
       enddo
      enddo
      stop
c
9     continue
c
c
      nt= t -1
      write(*,*)'No. of time records: ',nt
c
      end

c ))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))
c ((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((

      SUBROUTINE INDEXX(N,ARRIN,INDX)
      DIMENSION ARRIN(N),INDX(N)
      DO 11 J=1,N
        INDX(J)=J
11    CONTINUE
      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          INDXT=INDX(L)
          Q=ARRIN(INDXT)
        ELSE
          INDXT=INDX(IR)
          Q=ARRIN(INDXT)
          INDX(IR)=INDX(1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            INDX(1)=INDXT
            RETURN
          ENDIF
        ENDIF
        I=L
        J=L+L
20      IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF(ARRIN(INDX(J)).LT.ARRIN(INDX(J+1)))J=J+1
          ENDIF
          IF(Q.LT.ARRIN(INDX(J)))THEN
            INDX(I)=INDX(J)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        INDX(I)=INDXT
      GO TO 10
      END
       
c ))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))
c ((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((

      integer function ilen(string)
c
      character*1 c
      character*(*) string
c
      do i=len(string),1,-1
        c= string(i:i)
        if (c.ne.' ') goto 10
      enddo
c * String is wholly blank
      ilen= -1
      return
c
   10 ilen= i
c
      return
      end

c ))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))
c ((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((

      SUBROUTINE INDEXI(N,ARRIN,INDX)
      integer ARRIN(N),INDX(N)
      DO 11 J=1,N
        INDX(J)=J
11    CONTINUE
      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          INDXT=INDX(L)
          Q=ARRIN(INDXT)
        ELSE
          INDXT=INDX(IR)
          Q=ARRIN(INDXT)
          INDX(IR)=INDX(1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            INDX(1)=INDXT
            RETURN
          ENDIF
        ENDIF
        I=L
        J=L+L
20      IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF(ARRIN(INDX(J)).LT.ARRIN(INDX(J+1)))J=J+1
          ENDIF
          IF(Q.LT.ARRIN(INDX(J)))THEN
            INDX(I)=INDX(J)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        INDX(I)=INDXT
      GO TO 10
      END
       
c ))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))
c ((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((

