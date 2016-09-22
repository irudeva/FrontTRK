      module mcoor
       parameter (mlon=241,mlat=61)
      end module mcoor


      program front
	
	use mcoor
	implicit none
c	First program to identify fronts using cmp-files of v-component of wind

       character *80 fdv,farea,optarg
	
       character *80 head,ch42*42
       integer nlat,nlon
       real xlat(mlat),xlon(mlon)
       real dv(mlon,mlat),dat(mlon,mlat)


	real miss,cdv


	integer t,i,j,narg,nk,kk,ik,klab2


	narg=iargc()
	 if(narg.eq.0)then
	  write(*,*) 'Usage: front -C cdv -dv filedv -a areafile'
	 else
	  fdv=''
	  farea=''
	  cdv=2.  ! Critical dv
 	  i=1
	  do while(i.le.narg)
	  call getarg(i,optarg)
	  if(optarg.eq.'-C')then
	   i=i+1
	   call getarg(i,optarg)
	   read(optarg,*) cdv
	   write(*,*) "Note: cdv = ", cdv
	  elseif(optarg.eq.'-dv')then
	   i=i+1
	   call getarg(i,optarg)
	   fdv=optarg
	  elseif(optarg.eq.'-a')then
	   i=i+1
	   call getarg(i,optarg)
	   farea=optarg
	  endif
	  i=i+1
	  enddo
	 endif

         open(13,file=farea)
                                                     
         open(1,file=fdv,form='unformatted')
         write(*,*)'dv file: ',fdv(:len_trim(fdv))
   	                                                            
         read(1)nlat;print*, nlat
         read(1)(xlat(i),i=1,nlat)
         read(1)nlon
         read(1)(xlon(i),i=1,nlon)

	 write(13,'(''Nlon : '',i5)') nlon
	 write(13,'(''Nlat : '',i5)') nlat

	 write(*,'(''Nlat : '',i5)') nlat
	 write(*,'(''LatS '',2f10.2,'' ... '',2f10.2)') 
     &              xlat(1:2),xlat((nlat-1):nlat)
	 write(*,'(''Nlon : '',i5)') nlon
	 write(*,'(''LonS '',2f10.2,'' ... '',2f10.2)') 
     &              xlon(1:2),xlon((nlon-1):nlon)

	 t=0
	 do	 
         read(1,end=9)head
         read(1)((dv(i,j),i=1,nlon),j=1,nlat) 
         t=t+1

	 dv(nlon+1,:)=dv(1,:)     	

        dat=0.
	do i=1,nlon+1
	  do j=1,nlat
	    if(dv(i,j)>=cdv) dat(i,j)=1.
	  enddo
	enddo
	write(*,'(''Start connect'',i5)') t  
	call connect8c (nlon+1,nlat,dat,klab2)

	write(13,'(''Rec: '',i4,x,a)')t,head(41:53)
        write(13,'(''No.F:'',i4)')klab2
      do nk=1,klab2 
        kk= 0
	open(131,file='dummy')
        do j=1,nlat
	  do i=1,nlon
	    if(dat(i,j).eq.real(nk))then
	      kk= kk +1
	      write(131,'(3(i4,x),2(f8.2,x),'':'',f8.3)')kk,j,i,
     &                    xlat(j), xlon(i), dv(i,j)
	    endif
	  enddo
	enddo
	rewind(131)
        write(13,'(''Front: '',i4,x,''No.P: '',i4)') nk,kk
	do ik=1,kk
	 read(131,'(a)') ch42
	 write(13,'(a)') ch42
        enddo
        close(131,status='delete')
      enddo
	
	
	
	 enddo
9	 close(1)
         close(13)

	end

c ((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((

      subroutine connect8c (nlons,nlats,dat2,klab2)
      use mcoor
c
c * Author: Kevin Keay  Date: Jan 10 2008
c
c   (1) 24/2/2010; Routine: connect8c: Corrected the code that detects
c       the equivalence of objects that straddle lon 0 == 360 (wrapped
c       e.g. if 8 == 10, 10 == 12 and 12 == 17 then (8,10,12) -> 17
c   (2) 4/3/2010: The correction in (1) required further tweaking
c
      integer nlons,nlats,nr,nc,i,j
c
      real dat2(mlon,mlat)
c
      parameter (maxr=200,maxc=600)
      real b(maxr,maxc)
      real be(0:maxr+1,0:maxc+1)
      real cp(maxr,maxc)
      real cw(maxr,0:1)  ! For wrapped case
c
      parameter (maxlab=24000)  ! Assume 20% of grid could be individual objects
      real e(maxlab)
      real ef(maxlab)
      real ef2(maxlab)
      integer ie(maxlab)
      integer ko(maxlab)
c
      real lab,l
      integer klab2
c
ci      character*600 ch400  ! To hold (300I2) where maxc=300
c
c
c * Turn on column wrapping
c
      iwrap= 1  ! For wrapping (j= 1 == j= nc)
c
c * Convert dat1(mapk,i,j) -> b(i,j)
c
      nr= nlats
      nc= nlons
       do i=1,nr
        do j=1,nc
	  b(i,j)= dat2(j,i)
	  if(b(i,j).gt.0.)b(i,j)= 1.
	enddo
      enddo
c
c * Add a boundary of thickness one zero pixel to image
c
      do i=0,nr+1
        do j=0,nc+1
	  if (i.gt.0.and.i.lt.nr+1.and.j.gt.0.and.j.lt.nc+1)then
	    be(i,j)= b(i,j)
	  else
	    be(i,j)= 0.
	  endif
	enddo
      enddo
c
c
c * Initialise labels
c
      do i=1,nr
        do j=1,nc
	  cp(i,j)= 0.
	enddo
      enddo
c
c * Equivalence label array
c
      do k=1,maxlab
        e(k)= real(k)
      enddo
c
c * (1) Find components
c
      i= 1
      lab= 0.
      do while (i.le.nr)
        j= 1
        do while (j.le.nc)
	  p= be(i,j)
	  if(p.eq.0.)then  ! p is a background pixel
	    j= j +1  ! Move to next pixel
	  else  ! p is a foreground pixel
c
c d u f
c l p 
c
	    d= be(i-1,j-1)
	    u= be(i-1,j)
	    f= be(i-1,j+1)  ! 'e' in some texts e.g. McAndrew
	    l= be(i,j-1)
	    psum= u +l +d +f  ! Sum of pixels
c
c * All 4 pixels are zero - new component
c
            if(psum.eq.0)then
	      lab= lab +1 ! New label for new component
	      cp(i,j)= lab ! Assign label to pixel
c
c * There is 1 foreground pixel
c   The 4 combinations are: d u f l
c
	    elseif(psum.eq.1)then
	      if(d.eq.1)then
	        cp(i,j)= cp(i-1,j-1) ! Assign label of d to p
	      elseif(u.eq.1)then
	        cp(i,j)= cp(i-1,j)   ! Assign label of u to p
	      elseif(f.eq.1)then
	        cp(i,j)= cp(i-1,j+1) ! Assign label of f to p
	      elseif(l.eq.1)then
	        cp(i,j)= cp(i,j-1)   ! Assign label of l to p
	      endif
c
c * --------------
c
c * There are 2 foreground pixels
c   The 6 combinations are: du df dl uf ul fl
c
	    elseif(psum.eq.2)then
              if(d.eq.1.and.u.eq.1)then
	        call etab2 (i,j,i-1,j-1,i-1,j,maxr,maxc,maxlab,cp,e)   ! du
              elseif(d.eq.1.and.f.eq.1)then
	        call etab2 (i,j,i-1,j-1,i-1,j+1,maxr,maxc,maxlab,cp,e) ! df
              elseif(d.eq.1.and.l.eq.1)then
	        call etab2 (i,j,i-1,j-1,i,j-1,maxr,maxc,maxlab,cp,e)   ! dl
              elseif(u.eq.1.and.f.eq.1)then
	        call etab2 (i,j,i-1,j,i-1,j+1,maxr,maxc,maxlab,cp,e)   ! uf
              elseif(u.eq.1.and.l.eq.1)then
	        call etab2 (i,j,i-1,j,i,j-1,maxr,maxc,maxlab,cp,e)     ! ul
              elseif(f.eq.1.and.l.eq.1)then
	        call etab2 (i,j,i-1,j+1,i,j-1,maxr,maxc,maxlab,cp,e)   ! fl
	      endif
c
c * There are 3 foreground pixels
c   The 4 combinations are: duf dul dfl ufl
c
	    elseif(psum.eq.3)then
              if(d.eq.1.and.u.eq.1.and.f.eq.1)then
                call etab3 (i,j,i-1,j-1,i-1,j,i-1,j+1,
     *		  maxr,maxc,maxlab,cp,e)   ! duf
              elseif(d.eq.1.and.u.eq.1.and.l.eq.1)then
                call etab3 (i,j,i-1,j-1,i-1,j,i,j-1,
     *		  maxr,maxc,maxlab,cp,e)   ! dul
              elseif(d.eq.1.and.f.eq.1.and.l.eq.1)then
                call etab3 (i,j,i-1,j-1,i-1,j+1,i,j-1,
     *		  maxr,maxc,maxlab,cp,e)   ! dfl
              elseif(u.eq.1.and.f.eq.1.and.l.eq.1)then
                call etab3 (i,j,i-1,j,i-1,j+1,i,j-1,
     *		  maxr,maxc,maxlab,cp,e)   ! ufl
	      endif
c
c * There are 4 foreground pixels i.e. a singe combinaton, dufl
c
	    elseif(psum.eq.4)then
	      call etab4 (i,j,i-1,j-1,i-1,j,i-1,j+1,i,j-1,
     *	        maxr,maxc,maxlab,cp,e) ! dufl
c
	    endif
c
c * --------------
c	    
	    j= j +1
	  endif
	enddo
	i= i +1
      enddo
      nlab= lab  ! No. of labels
      write(*,*)'nlab: ',nlab
      if(nlab.eq.0)then
        write(*,*)'WARNING: No components detected'
	stop
      endif
c
ci      write(8,*)'First components:'
ci      do i=1,nr
ci        write(ch400,'(300I2)')(int(cp(i,j)),j=1,nc)
ci	write(8,'(A)')ch400(:ilen(ch400))
ci      enddo
c
c * (2) Use equivalence information to determine final component labels
c
c * Firstly, may need to sort array e e.g. 1 == 3, 2 == 4 
c
      do k=1,nlab  ! Sequence of labels
        ko(k)= k
      enddo
      call indexx (nlab,e,ie)
ci      write(8,*)'Equivalence array: raw and sorted'
ci      do k=1,nlab
ci        write(8,*)ko(ie(k)),e(ie(k))
ci      enddo
c
      elab= e(ie(1))
      klab= 1
      ef(ie(1))= klab
      do k=2,nlab
        if(e(ie(k)).ne.elab)then
	  klab= klab +1  ! New unique label
	  elab= e(ie(k))
	endif
	ef(ie(k))= klab
      enddo
c
ci      write(8,*)'Equivalence table:'
ci      write(8,*)'nlab: ',nlab
ci      write(8,*)'klab: ',klab
ci      write(*,*)'klab: ',klab
ci      write(8,*)'k,ie(k),ko(),ef():'
ci      do k=1,nlab
ci        write(8,*)k,ie(k),ko(ie(k)),ef(ie(k))
ci      enddo
c
      do i=1,nr
        do j=1,nc
	  c= cp(i,j) ! Current component label
	  if(c.gt.0)then
	    ic= int(c)
	    cnew= ef(ic)  ! New component label
	    cp(i,j)= cnew
	  endif
	enddo
      enddo
c
ci      write(8,*)'Final components (no wrapping):'
ci      do i=1,nr
ci        write(ch400,'(400I2)')(int(cp(i,j)),j=1,nc)
ci	write(8,'(A)')ch400(:ilen(ch400))
ci      enddo
c
c * Check for wrapping of data i.e. j= 1 == j= nc
c
      if(iwrap.eq.1)then
        write(*,*)'iwrap: ',iwrap
ci        write(8,*)'iwrap: ',iwrap
c
c * Reset labels to 1 ... klab
c   These represent the current set of objects
c
      do k=1,klab
        ef2(k)= real(k)
      enddo
c
c * Copy columns 1 and nc to a work array since these will be updated
c   in the wrapping equivalence process
c
      do i=1,nr
        cw(i,0)= cp(i,1)
        cw(i,1)= cp(i,nc)
      enddo
c
      do i=1,nr
        c1= cw(i,0)  ! cp(i,1)
        cnc= cw(i,1) ! cp(i,nc)
	if(c1.gt.0.and.cnc.gt.0)then
	  if(c1.lt.cnc)then
	    c= c1
	    c0= cnc
	    ic0= int(c0)
	  else  ! cnc < c1
	    c= cnc
	    c0= c1
	    ic0= int(c0)
	  endif
ckk          write(8,*)'---',i,c,c0
c * Replace larger colour by smaller colour
	  do k=1,klab
	    if(ef2(k).eq.c0)then
	      ef2(k)= c
	    endif
	  enddo
	  do ii=1,nr
	    do jj=0,1
	      if(cw(ii,jj).eq.c0)cw(ii,jj)= c 
	    enddo
	  enddo
	endif
      enddo

c * Check for smallest colour
ckk      do k=1,klab
ckk	kw= int(ef2(k))
ckk        if(kw.ne.k)then
ckk	  if(ef2(kw).lt.ef2(k))then
ckk	    ef2(k)= ef2(kw)
ckk	  endif
ckk	endif
ckk      enddo
c
      call indexx (klab,ef2,ie)
c
      elab= ef2(ie(1))
      klab2= 1
      ef(ie(1))= klab2
      do k=2,klab
        if(ef2(ie(k)).ne.elab)then
	  klab2= klab2 +1  ! New unique label
	  elab= ef2(ie(k))
	endif
	ef(ie(k))= klab2
      enddo
c
      do i=1,nr
        do j=1,nc
	  c= cp(i,j) ! Current component label
	  if(c.gt.0)then
	    ic= int(c)
	    cnew= ef(ic)  ! New component label
	    cp(i,j)= cnew
	  endif
	enddo
      enddo
c
      write(*,*)'Unique labels after wrapping:'
      write(*,*)'klab2: ',klab2
c
c
      endif
c
c
c
      do i=1,nlons
        do j=1,nlats
	  dat2(i,j)= cp(j,i)
	enddo
      enddo
      
      return
      end

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

      subroutine etab2 (i,j,i1,j1,i2,j2,maxr,maxc,maxlab,cp,e)
c
c * Author: Kevin Keay  Date: Jan 9 2008
c
      real cp(maxr,maxc)
      real e(maxlab)
c
      real lab1,lab2
c
      lab1= cp(i1,j1)
      lab2= cp(i2,j2)
      ilab1= int(lab1)
      ilab2= int(lab2)
c * Note equivalence of labels; assign smaller label to pixel
      if(lab1.lt.lab2)then
        e(ilab2)= e(ilab1)
        cp(i,j)= lab1
      else
        e(ilab1)= e(ilab2)
        cp(i,j)= lab2
      endif
c
      return
      end
       
c ))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))
c ((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((

      subroutine etab3 (i,j,i1,j1,i2,j2,i3,j3,maxr,maxc,maxlab,cp,e)
c
c * Author: Kevin Keay  Date: Jan 9 2008
c
      real cp(maxr,maxc)
      real e(maxlab)
c
      parameter (nx=3)
      integer ix(nx)
      real x(nx)
      real lab1,lab2,lab3
      equivalence (lab1,x(1))
      equivalence (lab2,x(2))
      equivalence (lab3,x(3))
c
      real lab
c
      lab1= cp(i1,j1)
      lab2= cp(i2,j2)
      lab3= cp(i3,j3)
      ilab1= int(lab1)
      ilab2= int(lab2)
      ilab3= int(lab3)
 
c * Sort labels into ascending order
      call indexx (nx,x,ix)

c * Note equivalence of labels; assign smallest label to pixel

      lab= x(ix(1))  ! Smallest label
      ilab= int(lab)
c
      ilab2= int(x(ix(2)))
      e(ilab2)= e(ilab)
      ilab3= int(x(ix(3)))
      e(ilab3)= e(ilab)
c
      cp(i,j)= lab
      if(lab.eq.0)then
        write(*,*)'i,j,i1,j1,i2,j2,i3,j3:'
        write(*,*)i,j,i1,j1,i2,j2,i3,j3
	write(*,*)lab,lab1,lab2,lab3
      endif
c
      return
      end
       
c ))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))
c ((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((

      subroutine etab4 (i,j,i1,j1,i2,j2,i3,j3,i4,j4,
     * maxr,maxc,maxlab,cp,e)
c
c * Author: Kevin Keay  Date: Jan 9 2008
c
      real cp(maxr,maxc)
      real e(maxlab)
c
      parameter (nx=4)
      integer ix(nx)
      real x(nx)
      real lab1,lab2,lab3,lab4
      equivalence (lab1,x(1))
      equivalence (lab2,x(2))
      equivalence (lab3,x(3))
      equivalence (lab4,x(4))
c
      real lab
c
      lab1= cp(i1,j1)
      lab2= cp(i2,j2)
      lab3= cp(i3,j3)
      lab4= cp(i4,j4)
      ilab1= int(lab1)
      ilab2= int(lab2)
      ilab3= int(lab3)
      ilab4= int(lab4)
 
c * Sort labels into ascending order
      call indexx (nx,x,ix)

c * Note equivalence of labels; assign smallest label to pixel

      lab= x(ix(1))  ! Smallest label
      ilab= int(lab)
c
      ilab2= int(x(ix(2)))
      e(ilab2)= e(ilab)
      ilab3= int(x(ix(3)))
      e(ilab3)= e(ilab)
      ilab4= int(x(ix(4)))
      e(ilab4)= e(ilab)
c
      cp(i,j)= lab
      if(lab.eq.0)then
        write(*,*)'x:',x
        write(*,*)'ix:',ix
        write(*,*)'i,j,i1,j1,i2,j2,i3,j3,i4,j4:'
        write(*,*)i,j,i1,j1,i2,j2,i3,j3,i4,j4
	write(*,*)lab,lab1,lab2,lab3,lab4
      endif
c
      return
      end
       
c ))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))
