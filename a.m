	subroutine elmt40(d,ul,xl,ix,tl,s,p,ndf,ndm,nst,isw)

c       Alexandros z Richtung
	
	!Bearbeitet 2D ELEM --> 3D ELEM
	
!	
c	Bearbeitet : 3 bilineare 8-Knoten Elemente
c	Erweiterung zum dG-Interface-Element
c	alpha & additional alpha
c
c		ebener Verzerrungszustand
c		geometrisch linear
c		isotropes Materialverhalten
c	Materialparameter:
c		d(1)   : E-Modul E
c		d(2)   : Querkontraktionszahl nu
c		d(3)   : Theta
c		d(4)   : Switch
c		d(5)   : ElementSize
c		d(6)   : additional alpha
c	Voigt-Notation:
c		1.Index: Komponente 11
c		2.Index: Komponente 22
c		3.Index: Komponente 12
c		3.Index: Komponente 33
c	DOF:
c		1.Index: Verschiebung in x-Richtung u
c		2.Index: Verschiebung in y-Richtung v

c	---------------------
c	         DEKLARATION
c	--------------------- 

	implicit none
	include 'iofile.h'
	include 'eldata.h'
	include 'cdata.h'
	include 'comblk.h'
	include 'prstrs.h'
	include 'pdata5.h'
	include 'pdata6.h'

c	FEAP-Variablen
	integer ix(*),ndf,ndm,lint,nst,isw,inf
! 	real*8 d(6),ul(3,16),xl(3,16),tl(*),s(nst,nst),p(nst) 
	real*8 d(4),ul(48),xl(ndm,*),tl(*),s(nst,*),p(*)
c	Schleifenzaehler
	integer ig,i,j,k,l,ii,jj,kk,ll,i1,j1,r,f,t,a,m
c	Materialmatrix
	real*8 cmat(6,6)
c	Spannungen/Dehnungen
	real*8 sigmaGP(4,4),sig(6),eps(6) !sigmadg nicht angepasst
	real*8 sigmapos(6),sigmaneg(6),epspos(6),epsneg(6)
c	Shape-Funktionen
	real*8 shp(16),shpX(16),shpY(16),shpZ(16)
	real*8 Bmatpos(6,3,8),Bmatneg(6,3,8)
	real*8 Nmatpos(3,3,8),Nmatneg(3,3,8)
	real*8 Naddpos(3,3,8),Naddneg(3,3,8)
	real*8 Bmatpos2(6,24),Bmatneg2(6,24)
	real*8 Nmatpos2(3,24),Nmatneg2(3,24)
	real*8 w,detJ,detJ2D

c	dG zusatz
	integer theta
	real*8 ds,beta
	real*8 normal(3,6),rint,elmtsize,alpha
	real*8 Kcon11(24,24),Kcon12(24,24),Kcon21(24,24),Kcon22(24,24)
	real*8 Kcon(48,48),KconT(48,48)
	real*8 Kalp11(24,24),Kalp12(24,24),Kalp21(24,24),Kalp22(24,24)
	real*8 Kalp(48,48)
	real*8 Kalpadd1(24,24),Kalpadd2(24,24),Kalpadd3(24,24)
	real*8 Kalpadd4(24,24)
	real*8 Kalpadd(48,48)

c	Ausgabe
	logical isda
	character*20 outfile
	integer outf
	character*20 Kconmitmatlab,temp(48)
	real*8 zahl


c	---------------------
c	         FEM
c	---------------------

c---------------Materialparameter einlesen
	if(isw.eq. 1) then

	call pinput(d(1), 4)
	write(*,*)'# # # # # # # # # # # # # # # # # # # # # # #'
	write(*,*)'#         I N P U T P A R A M E T E R       #'
	write(*,*)'# # # # # # # # # # # # # # # # # # # # # # #'
	write(*,*)'plain strain 8-node dG-element'
	write(*,*)'  E-Modul E			: ',d(1)
	write(*,*)'  Querkontraktion nu		: ',d(2)
	write(*,*)'  Theta			: ',d(3)
	write(*,*)'  additional Penalty		: ',d(4)

	write(iow,*)'plain strain 8-node dG-element'
	write(iow,*)'  E-Modul E	  : ',d(1)
	write(iow,*)'  Querkontraktion nu : ',d(2)
	write(iow,*)'  Theta		  : ',d(3)
	write(iow,*)'  additional Penalty : ',d(4)

	endif
	
	  call pzero(Kcon11,24*24)
	  call pzero(Kcon12,24*24)
	  call pzero(Kcon21,24*24)
	  call pzero(Kcon22,24*24)
	  call pzero(Kcon,48*48)
	  call pzero(KconT,48*48)
	  call pzero(Kalp11,24*24)
	  call pzero(Kalp12,24*24)
	  call pzero(Kalp21,24*24)
	  call pzero(Kalp22,24*24)
	  call pzero(Kalp,48*48)
	  call pzero(Kalpadd1,24*24)
	  call pzero(Kalpadd2,24*24)
	  call pzero(Kalpadd3,24*24)
	  call pzero(Kalpadd4,24*24)
	  call pzero(Kalpadd,48*48)
	  call pzero(Bmatpos,6*3*8)
	  call pzero(Bmatneg,6*3*8)
          call pzero(Nmatpos,3*3*8)
          call pzero(Nmatneg,3*3*8)
	  
c---------------Kontrolle
	if(isw .eq. 2) then
	
	endif

c-----------------------------------------------------------------------
c                                     Interface Formulation with 8 nodes
c-----------------------------------------------------------------------
c---------------STEMA/ Residuum
	if((isw.eq.3).or.(isw.eq.6)) then

	call pzero(s,48*48)
	call pzero(p,48)
	
C---------------------------------------------------------
c	Materialmatrix C-Matrix 3D

	call mate3D40(d,cmat)

	alpha = d(3)
	
C---------------------------------------------------------
c       Normal matrix

c	get normal matrix of face 2 on positive element between nodes 2 and 3
	call normalmat40(xl,normal)
	
C---------------------------------------------------------
c	loop over gauss points
        do ig=1,4
         
C---------------------------------------------------------
c	shape functions and derivatives wrt. xsi and eta
	    call shape40(ig,xl,shp,Bmatpos,Bmatneg,Nmatpos,Nmatneg,detJ
     + ,rint)
        
! !         rint=1.0d0
	 
c-------------------------------------------------------
c        !Original Kcon 
        
	do i = 1,8 !
	 do j = 1,8 !
	  do l = 1,6 !
	   do ii = 1,3 !
	    do jj = 1,3 !
	     do kk = 1,6 !
	      do ll = 1,3 !
	Kcon11(3*(i-1)+ii,3*(j-1)+ll) = Kcon11(3*(i-1)+ii,3*(j-1)+ll) +
     +	Nmatpos(jj,ii,i) * normal(jj,kk) * cmat(kk,l) * 
     +  Bmatpos(l,ll,j) * rint * (-0.5d0) 
	
	Kcon12(3*(i-1)+ii,3*(j-1)+ll) = Kcon12(3*(i-1)+ii,3*(j-1)+ll) +
     +	Nmatpos(jj,ii,i) * normal(jj,kk) * cmat(kk,l) * 
     +  Bmatneg(l,ll,j) * rint * (-0.5d0)

	Kcon21(3*(i-1)+ii,3*(j-1)+ll) = Kcon21(3*(i-1)+ii,3*(j-1)+ll) +
     +	Nmatneg(jj,ii,i) * normal(jj,kk) * cmat(kk,l) *
     +  Bmatpos(l,ll,j) * rint * (0.5d0)

	Kcon22(3*(i-1)+ii,3*(j-1)+ll) = Kcon22(3*(i-1)+ii,3*(j-1)+ll) +
     +	Nmatneg(jj,ii,i) * normal(jj,kk) * cmat(kk,l) * 
     +  Bmatneg(l,ll,j) * rint * (0.5d0)
 
	      enddo ! ll
	     enddo ! kk
	    enddo ! jj
	   enddo ! ii
	  enddo ! l
	 enddo ! j
	enddo ! i
c----------------------------------------------------------

c	stiffness loop Kalpha 11 : Kalpha11(24,24)=Kalp11(24,24)+ Nmatpos(3,3,8) * Nmatpos(3,3,8) * alpha * rint
	do i = 1,8
	 do ii = 1,3
	  do j = 1,8
	   do jj = 1,3
	    do kk = 1,3
	Kalp11(3*(i-1)+ii,3*(j-1)+kk) = Kalp11(3*(i-1)+ii,3*(j-1)+kk) +
     +  Nmatpos(jj,ii,i) * Nmatpos(jj,kk,j) * alpha * rint
     
	Kalp12(3*(i-1)+ii,3*(j-1)+kk) = Kalp12(3*(i-1)+ii,3*(j-1)+kk) +  
     +  Nmatpos(jj,ii,i) * Nmatneg(jj,kk,j) * alpha * rint * (-1.d0)
   
	Kalp21(3*(i-1)+ii,3*(j-1)+kk) = Kalp21(3*(i-1)+ii,3*(j-1)+kk) +
     +  Nmatneg(jj,ii,i) * Nmatpos(jj,kk,j) * alpha * rint * (-1.d0)
  
	Kalp22(3*(i-1)+ii,3*(j-1)+kk) = Kalp22(3*(i-1)+ii,3*(j-1)+kk) +
     +  Nmatneg(jj,ii,i) * Nmatneg(jj,kk,j) * alpha * rint
          
	    enddo
	   enddo
	  enddo
	 enddo
	enddo          

	enddo !4gp/face  
 	
 	ig = 1
 	
!  c	Zusätzlichen Anteil zur Elementsteifigkeitsmatrix mit alpha*/beta bilden
  	call shp1gp40(ig,xl,shp,Naddpos,Naddneg,w,rint)
 
!  c	Integration & Vorfaktor alpha* bzw beta
        !rint=1.0d0
  	ds = rint * w
  	beta = d(4)
c	write(*,*)'rint', rint
 
c	stiffness loop Kalphaadd 11
 	do i = 1,8
 	 do ii = 1,3
 	  do j = 1,8
 	   do jj = 1,3
 	    do kk = 1,3
 	Kalpadd1(3*(i-1)+ii,3*(j-1)+kk) = 
     +	Kalpadd1(3*(i-1)+ii,3*(j-1)+kk) + Naddpos(jj,ii,i) *
     +  Naddpos(jj,kk,j) * beta * ds
 
 	Kalpadd2(3*(i-1)+ii,3*(j-1)+kk) = 
     +	Kalpadd2(3*(i-1)+ii,3*(j-1)+kk) + Naddpos(jj,ii,i) *
     +  Naddneg(jj,kk,j) * beta * ds * (-1.d0)
 
 	Kalpadd3(3*(i-1)+ii,3*(j-1)+kk) = 
     +	Kalpadd3(3*(i-1)+ii,3*(j-1)+kk) + Naddneg(jj,ii,i) *
     +  Naddpos(jj,kk,j) * beta * ds * (-1.d0)
 
 	Kalpadd4(3*(i-1)+ii,3*(j-1)+kk) = 
     +	Kalpadd4(3*(i-1)+ii,3*(j-1)+kk) + Naddneg(jj,ii,i) *
     +  Naddneg(jj,kk,j) * beta * ds
 	    enddo
 	   enddo
 	  enddo
 	 enddo
 	enddo

c	Untermatrizen in Elementsteifigkeitsmatrix einsortieren element stiffness Kcon  

c       Kcon(48,48), Kconii(24,24), Kalp ...

	do i = 1,24 
	 do j = 1,24
	 
	 Kcon(i,j) = Kcon11(i,j)
	 Kcon(i,j+24) = Kcon12(i,j)
	 Kcon(i+24,j) = Kcon21(i,j)
	 Kcon(i+24,j+24) = Kcon22(i,j)

	 Kalp(i,j) = Kalp11(i,j)
	 Kalp(i,j+24) = Kalp12(i,j)
	 Kalp(i+24,j) = Kalp21(i,j)
	 Kalp(i+24,j+24) = Kalp22(i,j)
 
 	 Kalpadd(i,j) = Kalpadd1(i,j)
 	 Kalpadd(i,j+24) = Kalpadd2(i,j)
 	 Kalpadd(i+24,j) = Kalpadd3(i,j)
 	 Kalpadd(i+24,j+24) = Kalpadd4(i,j)
 	 
	 enddo
	enddo
! 	
! 	open(unit=777,file='Kalp_3Dz_alex.txt',status='unknown')
! 	 do ii = 1,48
! 	  write(777,'(48(E15.7,A1))')(Kalp(ii,jj),'  ',jj = 1,48)
! 	 enddo
! 	close(777)

c-------------------------------------------------------------------------------
c      Lösung mit Matlab
! 	
!        open(unit=333,file='/home/hiwits01/feap_alexandros/Kconmitmatlab'
!      +  ,status='old',action='read')
!         
!         do i=1,48
!            read(333,*) (temp(j),j=1,48)
!            do j=1,48
!            read(temp(j),*) zahl 
!            Kcon(i,j)=zahl
!          
!           enddo
!         enddo
!         
c-----------------------------------------------------------------------------------
        
!   call mprint(Kcon(1,1),48,48,48,'Kcon')
!   call mprint(Kcon11(1,1),24,24,24,'Kcon11')
      call mprint(Kalp(1,1),48,48,48,'Kalp_z')
     
c	element stiffness total
	do i = 1,48
	 do j = 1,48
	
          s(i,j) = Kcon(i,j) + Kalp(i,j) + Kalpadd(i,j)
!             s(i,j) = Kalp(i,j) + Kcon(i,j) 
!           s(i,j) = Kalp(i,j)

	 enddo
	enddo
	
 	do i = 1,48
	 do j = 1,48
	  p(i) = p(i) - s(i,j) * ul(j)
	 enddo
	enddo
	
C---------------------------------------------------------

c	neue gausspunkt schleife zur bestimmung des residuums
! 	do ig = 1,4                                          
! C---------------------------------------------------------
! c	Dehnungen Epsilon berechnen
! 
! 	call eps40(ig,xl,ul,epspos,epsneg)  
! 	!! mesa sthn routine eps40 kalloume thn shape40 oste na ipologisoume ta epspos kai epsneg
! 	!! otan omos pame na upologisoume to residuum den xanakaloume 
! 	!! thn shape40 me apotelesma na ginontai oloi oi ipologismoi 
! 	!! me ton teleutaio Nmatpos(ig=4) apo thn schleife opou ipologisthike o Kcon
! 	call shape40(ig,xl,shp,Bmatpos,Bmatneg,Nmatpos,Nmatneg,detJ,
!      +  rint)
!         rint=0.125d0
! 
! 	call pzero(sigmapos,6)
! 	call pzero(sigmaneg,6)
! C---------------------------------------------------------
! c       Spannungen Sigma berechnen
! 	do i = 1,6
! 	 do j = 1,6
! 	  sigmapos(i) = sigmapos(i) + cmat(i,j) * epspos(j)
! 	  sigmaneg(i) = sigmaneg(i) + cmat(i,j) * epsneg(j)
! 	 enddo
! 	enddo
! C---------------------------------------------------------
! c	Residuum berechnen (Anteil aus Kcon) Nmatpos(3,3,8), normal(3,6), sigmapos(6)
! 	do i = 1,3
! 	 do j = 1,3
! 	  do k = 1,6
! 	   do ii = 1,8
! 	   p(3*(ii-1)+i) = p(3*(ii-1)+i) - (
!      +		  Nmatpos(j,i,ii) * normal(j,k) * sigmapos(k) * (-0.5d0)
!      +		+ Nmatpos(j,i,ii) * normal(j,k) * sigmaneg(k) * (-0.5d0)
!      +					    ) * rint
! 	   p(24+(3*(ii-1)+i)) = p(24+(3*(ii-1)+i)) - (
!      +            Nmatneg(j,i,ii) * normal(j,k) * sigmapos(k) * (0.5d0)
!      +		+ Nmatneg(j,i,ii) * normal(j,k) * sigmaneg(k) * (0.5d0)
!      +					          ) * rint
! 	   enddo
! 	  enddo
! 	 enddo
! 	enddo
! 	
! c	Anteil aus Kalpha und Kalphaadd
! ! 	do i = 1,48
! ! 	 do j = 1,48
! !  	  p(i) = p(i) - (Kalp(i,j) + Kalpadd(i,j)) * ul(j)
! !           p(i) = p(i) - Kalp(i,j)*ul(j)
! ! 	 enddo
! ! 	enddo
! 	
! 	enddo
C---------------------------------------------------------
      endif !isw 3     
C---------------------------------------------------------
 
! c---------------Postprozess (Text)
! 	if(isw .eq. 4) then
! 
! ! 	do ig = 1,4
! ! 	  call norm08(ig,nel,xl,ndm,ul,d,n)
! ! 	enddo
! 	
!      	endif
! 
! c---------------Postprozess (Graphisch)
! 	if(isw .eq. 8) then
! 
! ! 	call stcn08(ix,xl,ul,d,hr(nph),hr(nph+numnp),lint,nel,numnp)
! 
! 	endif
! 
 	end

c=======================================================================
c---------   S U B R O U T I N E S   -----------------------------------
c=======================================================================    

c-----------------------------------------------------------------------
      subroutine mate3D40(d,cmat)
c-----------------------------------------------------------------------
c
c      Purpose: calculate linear isotropic elasticity matrix
c
c      Inputs:
c         v         - Poissons ratio
c         E         - Youngs modulus
c
c      Outputs:
c         cmat(3,3) - Material Matrix
c
c-----------------------------------------------------------------------
      implicit none
      real*8 d(*),cmat(6,6),E,v,fact
      
      call pzero(cmat,6*6)
      
c.... cmat plain strain
	fact=d(1)*(1.0d0-d(2))/((1.0d0+d(2))*(1.0d0-2.0d0*d(2)))

	cmat(1,1)=fact
	cmat(1,2)=fact*d(2)/(1.0d0-d(2))
	cmat(1,3)=cmat(1,2)
	
	cmat(2,1)=cmat(1,2)
	cmat(2,2)=fact
	cmat(2,3)=cmat(1,2)
	
	cmat(3,1)=cmat(1,2)	
	cmat(3,2)=cmat(1,2)
	cmat(3,3)=fact
	
	cmat(4,4)=fact*(1-2*d(2))/(2*(1-d(2)))
	cmat(5,5)=cmat(4,4)
	cmat(6,6)=cmat(4,4)

      return
      end  
      
c-----------------------------------------------------------------------
!       subroutine eps40(ig,xl,ul,epspos,epsneg)
! c-----------------------------------------------------------------------
! c
! c	VERZERRUNGEN BESTIMMEN
! c	Uebergabeparameter:
! c		ig    : Gausspunkt
! c		xl    : Knotenkoordinaten
! c		ul    : Knotenverschiebungen
! c		eps   : Verzerrungen		->Rueckgabe (Voigt-Notation)
! c
! c-----------------------------------------------------------------------
!       implicit none
!       integer ig,i,j,jj
!       real*8 xl(3,16),ul(48),epspos(6),epsneg(6)
!       real*8 shp(16),Bmatpos(6,3,8),Bmatneg(6,3,8),Nmatpos(3,3,8)
!       real*8 Nmatneg(3,3,8),detJ,rint
!       
! 	  call shape40(ig,xl,shp,Bmatpos,Bmatneg,Nmatpos,Nmatneg,detJ,
!      +    rint)
! 	  call pzero(epspos,6)
! 	  call pzero(epsneg,6)
! 	  
! 	  do i = 1,6
! 	   do jj = 1,3
! 	    do j = 1,8
! 	    epspos(i) = epspos(i) + Bmatpos(i,jj,j) * ul(3*(j-1)+jj)
! 	    epsneg(i) = epsneg(i) + Bmatneg(i,jj,j) *ul(24+(3*(j-1)+jj))
! 	    enddo
! 	   enddo
! 	  enddo
! 	  
! 	  return
! 	  end                                                            

c-----------------------------------------------------------------------
      subroutine normalmat40(xl,normal)                                      
c-----------------------------------------------------------------------
c
c      Purpose: compute normal vector of every face of a 4-noded element
c
c      Inputs:
c         xl(ndm,*)	- Nodal coordinates for element
c
c      
c      
c
c-----------------------------------------------------------------------
      implicit none
      real*8 xl(3,16),normvec(3,1),va(3,1),vb(3,1),vc(3,1),crvec(3,1)                                
      real*8 normal(3,6),magnitudea, magnitudeb, magnitudec 
      real*8 fi, h, surface, betragnormvec 
      integer i  !, betragnormvec  
      

      
      call pzero(normal,3*6)

c.... create list for the distance between the nodal coordinates of node 2-3 and 2-6 
	  
	  va(1,1) = xl(1,6)-xl(1,5)
	  va(2,1) = xl(2,6)-xl(2,5)
 	  va(3,1) = xl(3,6)-xl(3,5)
	
	  vb(1,1) = xl(1,8)-xl(1,5)
          vb(2,1) = xl(2,8)-xl(2,5)
          vb(3,1) = xl(3,8)-xl(3,5)

          
! c.... create list of the normalvector n = a x b  , a-Vektor from node 2->3, b-Vektor from node 2->6
! 
        normvec(1,1) = (va(2,1)*vb(3,1)-vb(2,1)*va(3,1))
        normvec(2,1) = (vb(1,1)*va(3,1)-va(1,1)*vb(3,1))
        normvec(3,1) = (va(1,1)*vb(2,1)-vb(1,1)*va(2,1))
        
        betragnormvec = sqrt(normvec(1,1)*normvec(1,1)+normvec(2,1)*
     +  normvec(2,1)+normvec(3,1)*normvec(3,1))
        
        
!          write(*,*)'normvec(1,1)',normvec(1,1)
!          write(*,*)'normvec(2,1)',normvec(2,1)
!          write(*,*)'normvec(3,1)',normvec(3,1)
!          write(*,*)'betragnormvec', betragnormvec
         
	
c.... normalmatrix on positive element (face 2)
! 	
! 	  normal(1,1) = normvec(1,1)/betragnormvec
! 	  normal(2,2) = normvec(2,1)/betragnormvec
! 	  normal(3,3) = normvec(3,1)/betragnormvec  
! 	  normal(1,4) = normvec(2,1)/betragnormvec
! 	  normal(2,4) = normvec(1,1)/betragnormvec  
! 	  normal(2,5) = normvec(3,1)/betragnormvec       
!           normal(3,5) = normvec(2,1)/betragnormvec  
!           normal(1,6) = normvec(3,1)/betragnormvec
!           normal(3,6) = normvec(1,1)/betragnormvec
          
          normal(1,1) = normvec(1,1)/betragnormvec
	  normal(2,2) = normvec(2,1)/betragnormvec
	  normal(3,3) = normvec(3,1)/betragnormvec 
	  normal(1,4) = normvec(2,1)/betragnormvec
	  normal(2,4) = normvec(1,1)/betragnormvec 
	  normal(2,5) = normvec(3,1)/betragnormvec      
          normal(3,5) = normvec(2,1)/betragnormvec 
          normal(1,6) = normvec(3,1)/betragnormvec
          normal(3,6) = normvec(1,1)/betragnormvec
          
!          
!         write(*,*)'normal(1,1)',normal(1,1)
!           write(*,*)'normal(2,2)',normal(2,2)
!           write(*,*)'normal(3,3)',normal(3,3)
!           write(*,*)'normal(1,4)',normal(1,4)
!           write(*,*)'normal(2,4)',normal(2,4)
!           write(*,*)'normal(2,5)',normal(2,5)
!           write(*,*)'normal(3,5)',normal(3,5)
!           write(*,*)'normal(1,6)',normal(1,6)
!           write(*,*)'normal(3,6)',normal(3,6)
          
! c--------------------------------
! c        Matlab beispiel: 
        
!         normal(1,1) = 0.0d0
!         normal(2,2) = 0.0d0
!         normal(3,3) = 1.0d0
!         normal(1,4) = 0.0d0
!         normal(2,4) = 0.0d0
!         normal(2,5) = 1.0d0
!         normal(3,5) = 0.0d0
!         normal(1,6) = 1.0d0
!         normal(3,6) = 0.0d0
c------------------------------------
        
!         call mprint(normal(1,1),3,6,3,'normal')
        
! c.... <<<<< Calculating the Surface of a Parallelogramm >>>>
! 
!         vc(1,1) = xl(1,7)-xl(1,5)
!         vc(2,1) = xl(2,7)-xl(2,5)
!         vc(3,1) = xl(3,7)-xl(3,5)
!         
!         magnitudea = dsqrt(va(1,1)**2+va(2,1)**2+va(3,1)**2)
!         magnitudeb = dsqrt(vb(1,1)**2+vb(2,1)**2+vb(3,1)**2)
!         magnitudec = dsqrt(vc(1,1)**2+vc(2,1)**2+vc(3,1)**2)
!         
!         fi = acos((magnitudea**2+magnitudeb**2-magnitudec**2)
!      +   /(2*magnitudea*magnitudeb))
!         
!         h = sin(fi)*magnitudeb
!         
!         surface = magnitudea*h
! 	
! 	rint = surface  

      return              
      end

c-----------------------------------------------------------------------
      subroutine shape40(ig,xl,shp,Bmatpos,Bmatneg,Nmatpos,Nmatneg,detJ,
     + rint)
c-----------------------------------------------------------------------
c
c      Purpose: 
c
c         ig         - gauss point
c         xl         - nodal coordinates
c	  shp	     - shape functionss Ziel verfolgt, einfach zu handhabende Standardprogramme zur Verfügung zu stellen, auf dass der Benutzer weniger häufig darauf angewiesen ist, das Terminal zu benutzen.
c	  shpX	     - derivative of shape function in X
c	  shpY	     - derivative of shape function in Y
c	  Bmat	     - B-matrix
c	  w	     - weighting
c	  detJ	     - determinant of Jacobian
c
c-----------------------------------------------------------------------

      implicit none
      integer i,ig
      real*8 xl(3,16),shp(16),shpX(16),shpY(16),shpZ(16),invJ(3,3)
      real*8 w,detJ,CofJ(3,3),c,d,e,rint,Jzeta,Jzeta1,Jzeta2,detJ2D
      real*8 gp(8,3),shpXi(16),shpEta(16),shpZeta(16),J(3,3)
      real*8 Bmatpos(6,3,8),Bmatneg(6,3,8)
      real*8 Nmatpos(3,3,8),Nmatneg(3,3,8)

      w = 1.0d0  

      gp(1,1)= -1.0d0/dsqrt(3.0d0)           ! GAUSS PUNKTE Xi ELEM 1
      gp(2,1)=  1.0d0/dsqrt(3.0d0)
      gp(3,1)=  1.0d0/dsqrt(3.0d0)
      gp(4,1)= -1.0d0/dsqrt(3.0d0)
      
      gp(5,1)= -1.0d0/dsqrt(3.0d0)              ! GAUSS PUNKTE Xi ELEM 2
      gp(6,1)=  1.0d0/dsqrt(3.0d0)
      gp(7,1)=  1.0d0/dsqrt(3.0d0)
      gp(8,1)= -1.0d0/dsqrt(3.0d0)
      
      gp(1,2)= -1.0d0/dsqrt(3.0d0)        ! GAUSS PUNKTE ETA ELEM 1
      gp(2,2)= -1.0d0/dsqrt(3.0d0)
      gp(3,2)=  1.0d0/dsqrt(3.0d0)
      gp(4,2)=  1.0d0/dsqrt(3.0d0)
      
      gp(5,2)= -1.0d0/dsqrt(3.0d0)        ! GAUSS PUNKTE ETA ELEM 2
      gp(6,2)= -1.0d0/dsqrt(3.0d0)
      gp(7,2)=  1.0d0/dsqrt(3.0d0)
      gp(8,2)=  1.0d0/dsqrt(3.0d0)
      
      gp(1,3)= 1.0d0                       ! GAUSS PUNKTE ZETA ELEM 1
      gp(2,3)= 1.0d0 
      gp(3,3)= 1.0d0  
      gp(4,3)= 1.0d0  
      
      gp(5,3)= -1.0d0                      ! GAUSS PUNKTE ZETA ELEM 2
      gp(6,3)= -1.0d0
      gp(7,3)= -1.0d0
      gp(8,3)= -1.0d0
                                                                                              ! SHAPE FUKTIONS 8 KNOTEN ELEMENT FÜR ELEM 1
	  shp(1) = 0.125d0 * (1.0d0 - gp(ig,1)) * (1.0d0 - gp(ig,2))
     +  * (1.0d0 - gp(ig,3))       
	  shp(2) = 0.125d0 * (1.0d0 + gp(ig,1)) * (1.0d0 - gp(ig,2)) 
     +  * (1.0d0 - gp(ig,3))
          shp(3) = 0.125d0 * (1.0d0 + gp(ig,1)) * (1.0d0 + gp(ig,2))
     +  * (1.0d0 - gp(ig,3))
	  shp(4) = 0.125d0 * (1.0d0 - gp(ig,1)) * (1.0d0 + gp(ig,2))
     +  * (1.0d0 - gp(ig,3))
          shp(5) = 0.125d0 * (1.0d0 - gp(ig,1)) * (1.0d0 - gp(ig,2))
     +  * (1.0d0 + gp(ig,3))
	  shp(6) = 0.125d0 * (1.0d0 + gp(ig,1)) * (1.0d0 - gp(ig,2))
     +  * (1.0d0 + gp(ig,3))
	  shp(7) = 0.125d0 * (1.0d0 + gp(ig,1)) * (1.0d0 + gp(ig,2))
     +  * (1.0d0 + gp(ig,3))
	  shp(8) = 0.125d0 * (1.0d0 - gp(ig,1)) * (1.0d0 + gp(ig,2))
     +  * (1.0d0 + gp(ig,3))

	shpXi(1) = -0.125d0 *  (1.0d0 - gp(ig,2)) * (1.0d0 - gp(ig,3))                ! ABLEITUNG NACH XI
	shpXi(2) = +0.125d0 *  (1.0d0 - gp(ig,2)) * (1.0d0 - gp(ig,3))
	shpXi(3) = +0.125d0 *  (1.0d0 + gp(ig,2)) * (1.0d0 - gp(ig,3))
	shpXi(4) = -0.125d0 *  (1.0d0 + gp(ig,2)) * (1.0d0 - gp(ig,3))
	shpXi(5) = -0.125d0 *  (1.0d0 - gp(ig,2)) * (1.0d0 + gp(ig,3))
	shpXi(6) = +0.125d0 *  (1.0d0 - gp(ig,2)) * (1.0d0 + gp(ig,3))
	shpXi(7) = +0.125d0 *  (1.0d0 + gp(ig,2)) * (1.0d0 + gp(ig,3))
	shpXi(8) = -0.125d0 *  (1.0d0 + gp(ig,2)) * (1.0d0 + gp(ig,3))

	shpEta(1) =  -0.125d0 * (1.0d0 - gp(ig,1)) * (1.0d0 - gp(ig,3))                ! ABLEITUNG NACH ETA 
	shpEta(2) =  -0.125d0 * (1.0d0 + gp(ig,1)) * (1.0d0 - gp(ig,3))
	shpEta(3) =  +0.125d0 * (1.0d0 + gp(ig,1)) * (1.0d0 - gp(ig,3))
	shpEta(4) =  +0.125d0 * (1.0d0 - gp(ig,1)) * (1.0d0 - gp(ig,3))
	shpEta(5) =  -0.125d0 * (1.0d0 - gp(ig,1)) * (1.0d0 + gp(ig,3))
	shpEta(6) =  -0.125d0 * (1.0d0 + gp(ig,1)) * (1.0d0 + gp(ig,3))
	shpEta(7) =  +0.125d0 * (1.0d0 + gp(ig,1)) * (1.0d0 + gp(ig,3))
	shpEta(8) =  +0.125d0 * (1.0d0 - gp(ig,1)) * (1.0d0 + gp(ig,3))
	
	shpZeta(1) = -0.125d0 * (1.0d0 - gp(ig,1)) * (1.0d0 - gp(ig,2))                 ! ABLEUTUNG NACH ZETA
	shpZeta(2) = -0.125d0 * (1.0d0 + gp(ig,1)) * (1.0d0 - gp(ig,2))
	shpZeta(3) = -0.125d0 * (1.0d0 + gp(ig,1)) * (1.0d0 + gp(ig,2))
	shpZeta(4) = -0.125d0 * (1.0d0 - gp(ig,1)) * (1.0d0 + gp(ig,2))
	shpZeta(5) = +0.125d0 * (1.0d0 - gp(ig,1)) * (1.0d0 - gp(ig,2))
	shpZeta(6) = +0.125d0 * (1.0d0 + gp(ig,1)) * (1.0d0 - gp(ig,2))
	shpZeta(7) = +0.125d0 * (1.0d0 + gp(ig,1)) * (1.0d0 + gp(ig,2))
	shpZeta(8) = +0.125d0 * (1.0d0 - gp(ig,1)) * (1.0d0 + gp(ig,2))
	
	  call pzero(J,3*3)                                                          ! JACOBI LAUT ELEM 1 
	  do i = 1,8
	    J(1,1) = J(1,1) + shpXi(i) * xl(1,i)                                 !!!!!! Richtig
	    J(1,2) = J(1,2) + shpXi(i) * xl(2,i)
	    J(1,3) = J(1,3) + shpXi(i) * xl(3,i)
	    
	    J(2,1) = J(2,1) + shpEta(i) * xl(1,i)
	    J(2,2) = J(2,2) + shpEta(i) * xl(2,i)
	    J(2,3) = J(2,3) + shpEta(i) * xl(3,i)
	    
	    J(3,1) = J(3,1) + shpZeta(i) * xl(1,i)
	    J(3,2) = J(3,2) + shpZeta(i) * xl(2,i)
	    J(3,3) = J(3,3) + shpZeta(i) * xl(3,i)
	  enddo

	detJ =J(1,1)*(J(2,2)*J(3,3)-J(2,3)*J(3,2))-J(1,2)*(J(2,1)*J(3,3)      ! determinant von jacobian matrix
     +          -J(2,3)*J(3,1))+J(1,3)*(J(2,1)*J(3,2)-J(2,2)*J(3,1))
     
      detJ2D = J(1,1)*J(2,2)-J(1,2)*J(2,1)
    
      rint = detJ2D * w
      
      !write(*,*)

	  CofJ(1,1)=J(2,2)*J(3,3)-J(2,3)*J(3,2)        
	  CofJ(1,2)=J(1,3)*J(3,2)-J(1,2)*J(3,3)
	  CofJ(1,3)=J(1,2)*J(2,3)-J(1,3)*J(2,2)
	  
	  CofJ(2,1)=J(2,3)*J(3,1)-J(2,1)*J(3,3)
          CofJ(2,2)=J(1,1)*J(3,3)-J(1,3)*J(3,1)
          CofJ(2,3)=J(1,3)*J(2,1)-J(1,1)*J(2,3)
          
          CofJ(3,1)=J(2,1)*J(3,2)-J(2,2)*J(3,1)
          CofJ(3,2)=J(1,2)*J(3,1)-J(1,1)*J(3,2)
          CofJ(3,3)=J(1,1)*J(2,2)-J(1,2)*J(2,1)
          
          call pzero(invJ,3*3)
          
          invJ(1,1)=CofJ(1,1)/detJ                         ! Inverse jacobian matrix
          invJ(1,2)=CofJ(1,2)/detJ
          invJ(1,3)=CofJ(1,3)/detJ
          
          invJ(2,1)=CofJ(2,1)/detJ
          invJ(2,2)=CofJ(2,2)/detJ
          invJ(2,3)=CofJ(2,3)/detJ
          
          invJ(3,1)=CofJ(3,1)/detJ
          invJ(3,2)=CofJ(3,2)/detJ
          invJ(3,3)=CofJ(3,3)/detJ
                    
	    call pzero(Bmatpos,6*3*8)                   ! Bmat(6,3,8) und  Nmat(3,3,8) 
	    call pzero(Nmatpos,3*3*8)
	    do i = 1,8

	      Nmatpos(1,1,i) = shp(i)                         
	      Nmatpos(2,2,i) = shp(i)
	      Nmatpos(3,3,i) = shp(i) 
	      
	      shpX(i) = invJ(1,1)*shpXi(i)+invJ(1,2)*shpEta(i)+invJ(1,3)         !!!!!! Richtig
     +        * shpZeta(i)
	      shpY(i) = invJ(2,1)*shpXi(i)+invJ(2,2)*shpEta(i)+invJ(2,3)
     +        * shpZeta(i)
	      shpZ(i) = invJ(3,1)*shpXi(i)+invJ(3,2)*shpEta(i)+invJ(3,3)
     +        * shpZeta(i)

	      Bmatpos(1,1,i) = shpX(i)        !Bmatpos                           !!!!!! Richtig
	      Bmatpos(2,2,i) = shpY(i)
	      Bmatpos(3,3,i) = shpZ(i)
	      
	      Bmatpos(4,1,i) = shpY(i)
	      Bmatpos(4,2,i) = shpX(i)
	      
	      Bmatpos(5,2,i) = shpZ(i)
	      Bmatpos(5,3,i) = shpY(i)
	      
	      Bmatpos(6,1,i) = shpZ(i)
	      Bmatpos(6,3,i) = shpX(i)
	      
	    enddo
	    
	    ! SHAPE FUKTIONS 8 KNOTEN ELEMENT FÜR ELEM 2
	 shp(9)  = 0.125d0 * (1.0d0 - gp(ig+4,1)) * (1.0d0 - gp(ig+4,2))
     +  *(1.0d0 - gp(ig+4,3))        
	 shp(10) = 0.125d0 * (1.0d0 + gp(ig+4,1)) * (1.0d0 - gp(ig+4,2))
     +  *(1.0d0 - gp(ig+4,3))
	 shp(11) = 0.125d0 * (1.0d0 + gp(ig+4,1)) * (1.0d0 + gp(ig+4,2))
     +  *(1.0d0 - gp(ig+4,3))
	 shp(12) = 0.125d0 * (1.0d0 - gp(ig+4,1)) * (1.0d0 + gp(ig+4,2))
     +  *(1.0d0 - gp(ig+4,3))
         shp(13) = 0.125d0 * (1.0d0 - gp(ig+4,1)) * (1.0d0 - gp(ig+4,2))
     +  *(1.0d0 + gp(ig+4,3))
	 shp(14) = 0.125d0 * (1.0d0 + gp(ig+4,1)) * (1.0d0 - gp(ig+4,2))
     +  *(1.0d0 + gp(ig+4,3))
	 shp(15) = 0.125d0 * (1.0d0 + gp(ig+4,1)) * (1.0d0 + gp(ig+4,2))
     +  *(1.0d0 + gp(ig+4,3))
	 shp(16) = 0.125d0 * (1.0d0 - gp(ig+4,1)) * (1.0d0 + gp(ig+4,2))
     +  *(1.0d0 + gp(ig+4,3))
      
	shpXi(9)  = -0.125d0 * (1.0d0 - gp(ig+4,2)) *(1.0d0 -gp(ig+4,3))                              ! ABLEITUNG NACH XI    ELEM 2
	shpXi(10) = +0.125d0 * (1.0d0 - gp(ig+4,2)) *(1.0d0 -gp(ig+4,3))
	shpXi(11) = +0.125d0 * (1.0d0 + gp(ig+4,2)) *(1.0d0 -gp(ig+4,3))
	shpXi(12) = -0.125d0 * (1.0d0 + gp(ig+4,2)) *(1.0d0 -gp(ig+4,3))
	shpXi(13) = -0.125d0 * (1.0d0 - gp(ig+4,2)) *(1.0d0 +gp(ig+4,3))
	shpXi(14) = +0.125d0 * (1.0d0 - gp(ig+4,2)) *(1.0d0 +gp(ig+4,3))
	shpXi(15) = +0.125d0 * (1.0d0 + gp(ig+4,2)) *(1.0d0 +gp(ig+4,3))
	shpXi(16) = -0.125d0 * (1.0d0 + gp(ig+4,2)) *(1.0d0 +gp(ig+4,3))
	
	shpEta(9)  =  -0.125d0 * (1.0d0 - gp(ig+4,1))*(1.0d0-gp(ig+4,3))                    ! ABLEITUNG NACH ETA    ELEM 2
	shpEta(10) =  -0.125d0 * (1.0d0 + gp(ig+4,1))*(1.0d0-gp(ig+4,3))                                          
	shpEta(11) =  +0.125d0 * (1.0d0 + gp(ig+4,1))*(1.0d0-gp(ig+4,3))
	shpEta(12) =  +0.125d0 * (1.0d0 - gp(ig+4,1))*(1.0d0-gp(ig+4,3))
	shpEta(13) =  -0.125d0 * (1.0d0 - gp(ig+4,1))*(1.0d0+gp(ig+4,3))
        shpEta(14) =  -0.125d0 * (1.0d0 + gp(ig+4,1))*(1.0d0+gp(ig+4,3))
	shpEta(15) =  +0.125d0 * (1.0d0 + gp(ig+4,1))*(1.0d0+gp(ig+4,3))
	shpEta(16) =  +0.125d0 * (1.0d0 - gp(ig+4,1))*(1.0d0+gp(ig+4,3))
	
	shpZeta(9)  = -0.125d0 * (1.0d0 - gp(ig+4,1))*(1.0d0-gp(ig+4,2))                  ! ABLEUTUNG NACH ZETA  ELEM 2
	shpZeta(10) = -0.125d0 * (1.0d0 + gp(ig+4,1))*(1.0d0-gp(ig+4,2)) 
	shpZeta(11) = -0.125d0 * (1.0d0 + gp(ig+4,1))*(1.0d0+gp(ig+4,2)) 
	shpZeta(12) = -0.125d0 * (1.0d0 - gp(ig+4,1))*(1.0d0+gp(ig+4,2)) 
	shpZeta(13) = +0.125d0 * (1.0d0 - gp(ig+4,1))*(1.0d0-gp(ig+4,2)) 
	shpZeta(14) = +0.125d0 * (1.0d0 + gp(ig+4,1))*(1.0d0-gp(ig+4,2)) 
	shpZeta(15) = +0.125d0 * (1.0d0 + gp(ig+4,1))*(1.0d0+gp(ig+4,2)) 
	shpZeta(16) = +0.125d0 * (1.0d0 - gp(ig+4,1))*(1.0d0+gp(ig+4,2)) 
	                                                                     
	  call pzero(J,3*3)
	  call pzero(CofJ,3*3)
	  call pzero(invJ,3*3)                                                                ! JACOBI LAUT ELEM 2
	  
	  do i = 9,16
	    J(1,1) = J(1,1) + shpXi(i) * xl(1,i)                           !!!!!! Richtig
	    J(1,2) = J(1,2) + shpXi(i) * xl(2,i)
	    J(1,3) = J(1,3) + shpXi(i) * xl(3,i)
	    
	    J(2,1) = J(2,1) + shpEta(i) * xl(1,i)
	    J(2,2) = J(2,2) + shpEta(i) * xl(2,i)
	    J(2,3) = J(2,3) + shpEta(i) * xl(3,i)
	    
	    J(3,1) = J(3,1) + shpZeta(i) * xl(1,i)
	    J(3,2) = J(3,2) + shpZeta(i) * xl(2,i)
	    J(3,3) = J(3,3) + shpZeta(i) * xl(3,i)
	  enddo
	  	  
          detJ=J(1,1)*(J(2,2)*J(3,3)-J(2,3)*J(3,2))-J(1,2)*  ! determinant von jacobian matrix
     +	  (J(2,1)*J(3,3)-J(2,3)*J(3,1))+
     +    J(1,3)*(J(2,1)*J(3,2)-J(2,2)*J(3,1)) 
          
	  CofJ(1,1)=J(2,2)*J(3,3)-J(2,3)*J(3,2)        
	  CofJ(1,2)=J(1,3)*J(3,2)-J(1,2)*J(3,3)
	  CofJ(1,3)=J(1,2)*J(2,3)-J(1,3)*J(2,2)
	  
      CofJ(2,1)=J(2,3)*J(3,1)-J(2,1)*J(3,3)
	  CofJ(2,2)=J(1,1)*J(3,3)-J(1,3)*J(3,1)
	  CofJ(2,3)=J(1,3)*J(2,1)-J(1,1)*J(2,3)
	  
	  CofJ(3,1)=J(2,1)*J(3,2)-J(2,2)*J(3,1)
	  CofJ(3,2)=J(1,2)*J(3,1)-J(1,1)*J(3,2)
	  CofJ(3,3)=J(1,1)*J(2,2)-J(1,2)*J(2,1)
	  
	  invJ(1,1)=CofJ(1,1)/detJ                       ! Inverse jacobian matrix
	  invJ(1,2)=CofJ(1,2)/detJ
	  invJ(1,3)=CofJ(1,3)/detJ
	  
	  invJ(2,1)=CofJ(2,1)/detJ
	  invJ(2,2)=CofJ(2,2)/detJ
	  invJ(2,3)=CofJ(2,3)/detJ
	  
	  invJ(3,1)=CofJ(3,1)/detJ
	  invJ(3,2)=CofJ(3,2)/detJ
	  invJ(3,3)=CofJ(3,3)/detJ
	  
	    call pzero(Bmatneg,6*3*8)                    ! Bmatneg und Nmatneg 
	    call pzero(Nmatneg,3*3*8)
	    do i = 9,16                                                                   

	      Nmatneg(1,1,i-8) = shp(i)
	      Nmatneg(2,2,i-8) = shp(i)
	      Nmatneg(3,3,i-8) = shp(i)
                                                                                        !!!!!! Richtig            ! Shape FUKTIONS nach X Y Z
	      shpX(i) = invJ(1,1)*shpXi(i)+invJ(1,2)*shpEta(i)+invJ(1,3)                   
     +        * shpZeta(i)       
	      shpY(i) = invJ(2,1)*shpXi(i)+invJ(2,2)*shpEta(i)+invJ(2,3)
     +        * shpZeta(i)
              shpZ(i) = invJ(3,1)*shpXi(i)+invJ(3,2)*shpEta(i)+invJ(3,3)
     +        * shpZeta(i)
              
                  
	     Bmatneg(1,1,i-8) = shpX(i)            ! Bmatneg                !!!!!! Richtig
	     Bmatneg(2,2,i-8) = shpY(i)
	     Bmatneg(3,3,i-8) = shpZ(i)
	                 
	     Bmatneg(4,1,i-8) = shpY(i) 
	     Bmatneg(4,2,i-8) = shpX(i) 
	                 
	     Bmatneg(5,2,i-8) = shpZ(i) 
	     Bmatneg(5,3,i-8) = shpY(i) 
	                 
	     Bmatneg(6,1,i-8) = shpZ(i) 
	     Bmatneg(6,3,i-8) = shpX(i) 	      	      
	      
	    enddo
	    	    	   
      return
      end

c-----------------------------------------------------------------------
      subroutine shp1gp40(ig,xl,shp,Naddpos,Naddneg,w,rint)
c-----------------------------------------------------------------------
c
c      Purpose: 
c
c         ig         - gauss point
c         xl         - nodal coordinates
c	  shp	     - shape functions
c	  shpX	     - derivative of shape function in X
c	  shpY	     - derivative of shape function in Y
c	  Bmat	     - B-matrix
c	  w	     - weighting
c	  detJ	     - determinant of Jacobian
c
c-----------------------------------------------------------------------

        implicit none                                                                     
        integer i,ig
        real*8 xl(3,16),shp(16),w,CofJ(3,3)
        real*8 gp(2,3),rint,shpXi(16),shpEta(16),shpZeta(16),J(3,3)
        real*8 Naddpos(3,3,8),Naddneg(3,3,8),Jzeta,c,d,e,detJ,invJ(3,3)

        w = 2.0d0  ! (w = 4.0d0)

        gp(1,1) =  0.0d0
        gp(2,1) =  0.0d0
      
        gp(1,2) =  0.0d0
        gp(2,2) =  0.0d0
      
        gp(1,3) =  1.0d0
        gp(2,3) = -1.0d0

	shp(1) = 0.125d0 * (1.0d0 - gp(ig,1)) * (1.0d0 - gp(ig,2))
     +  * (1.0d0 - gp(ig,3))
	shp(2) = 0.125d0 * (1.0d0 + gp(ig,1)) * (1.0d0 - gp(ig,2))
     +  * (1.0d0 - gp(ig,3))
	shp(3) = 0.125d0 * (1.0d0 + gp(ig,1)) * (1.0d0 + gp(ig,2))
     +  * (1.0d0 - gp(ig,3))
	shp(4) = 0.125d0 * (1.0d0 - gp(ig,1)) * (1.0d0 + gp(ig,2))
     +  * (1.0d0 - gp(ig,3))
        shp(5) = 0.125d0 * (1.0d0 - gp(ig,1)) * (1.0d0 - gp(ig,2))
     +  * (1.0d0 + gp(ig,3))
	shp(6) = 0.125d0 * (1.0d0 + gp(ig,1)) * (1.0d0 - gp(ig,2)) 
     +  * (1.0d0 + gp(ig,3))
	shp(7) = 0.125d0 * (1.0d0 + gp(ig,1)) * (1.0d0 + gp(ig,2))
     +  * (1.0d0 + gp(ig,3))
	shp(8) = 0.125d0 * (1.0d0 - gp(ig,1)) * (1.0d0 + gp(ig,2))
     +  * (1.0d0 + gp(ig,3))
	
	    call pzero(Naddpos,3*3*8)
	    do i = 1,8

	      Naddpos(1,1,i) = shp(i)
	      Naddpos(2,2,i) = shp(i)
	      Naddpos(3,3,i) = shp(i)
	      
	    enddo
!         
       shpXi(1) = -0.125d0 *  (1.0d0 - gp(ig,2)) * (1.0d0 - gp(ig,3))                   ! ABLEITUNG NACH XI
 	shpXi(2) = +0.125d0 *  (1.0d0 - gp(ig,2)) * (1.0d0 - gp(ig,3))
 	shpXi(3) = +0.125d0 *  (1.0d0 + gp(ig,2)) * (1.0d0 - gp(ig,3))
 	shpXi(4) = -0.125d0 *  (1.0d0 + gp(ig,2)) * (1.0d0 - gp(ig,3))
 	shpXi(5) = -0.125d0 *  (1.0d0 - gp(ig,2)) * (1.0d0 + gp(ig,3))
 	shpXi(6) = +0.125d0 *  (1.0d0 - gp(ig,2)) * (1.0d0 + gp(ig,3))
 	shpXi(7) = +0.125d0 *  (1.0d0 + gp(ig,2)) * (1.0d0 + gp(ig,3))
 	shpXi(8) = -0.125d0 *  (1.0d0 + gp(ig,2)) * (1.0d0 + gp(ig,3))
 
 	shpEta(1) =  -0.125d0 * (1.0d0 - gp(ig,1)) * (1.0d0 - gp(ig,3))                ! ABLEITUNG NACH ETA 
 	shpEta(2) =  -0.125d0 * (1.0d0 + gp(ig,1)) * (1.0d0 - gp(ig,3))
 	shpEta(3) =  +0.125d0 * (1.0d0 + gp(ig,1)) * (1.0d0 - gp(ig,3))
 	shpEta(4) =  +0.125d0 * (1.0d0 - gp(ig,1)) * (1.0d0 - gp(ig,3))
 	shpEta(5) =  -0.125d0 * (1.0d0 - gp(ig,1)) * (1.0d0 + gp(ig,3))
 	shpEta(6) =  -0.125d0 * (1.0d0 + gp(ig,1)) * (1.0d0 + gp(ig,3))
 	shpEta(7) =  +0.125d0 * (1.0d0 + gp(ig,1)) * (1.0d0 + gp(ig,3))
 	shpEta(8) =  +0.125d0 * (1.0d0 - gp(ig,1)) * (1.0d0 + gp(ig,3))
 	
 	shpZeta(1) = -0.125d0 * (1.0d0 - gp(ig,1)) * (1.0d0 - gp(ig,2))                 ! ABLEUTUNG NACH ZETA
 	shpZeta(2) = -0.125d0 * (1.0d0 + gp(ig,1)) * (1.0d0 - gp(ig,2))
 	shpZeta(3) = -0.125d0 * (1.0d0 + gp(ig,1)) * (1.0d0 + gp(ig,2))
 	shpZeta(4) = -0.125d0 * (1.0d0 - gp(ig,1)) * (1.0d0 + gp(ig,2))
 	shpZeta(5) = +0.125d0 * (1.0d0 - gp(ig,1)) * (1.0d0 - gp(ig,2))
 	shpZeta(6) = +0.125d0 * (1.0d0 + gp(ig,1)) * (1.0d0 - gp(ig,2))
 	shpZeta(7) = +0.125d0 * (1.0d0 + gp(ig,1)) * (1.0d0 + gp(ig,2))
 	shpZeta(8) = +0.125d0 * (1.0d0 - gp(ig,1)) * (1.0d0 + gp(ig,2))
 	
 	  call pzero(J,3*3)                                                          ! JACOBI LAUT ELEM 1 
 	  do i = 1,8
 	    J(1,1) = J(1,1) + shpXi(i) * xl(1,i)
 	    J(1,2) = J(1,2) + shpXi(i) * xl(2,i)
 	    J(1,3) = J(1,3) + shpXi(i) * xl(3,i)
 	    
 	    J(2,1) = J(2,1) + shpEta(i) * xl(1,i)
 	    J(2,2) = J(2,2) + shpEta(i) * xl(2,i)
 	    J(2,3) = J(2,3) + shpEta(i) * xl(3,i)
 	    
 	    J(3,1) = J(3,1) + shpZeta(i) * xl(1,i)
 	    J(3,2) = J(3,2) + shpZeta(i) * xl(2,i)
 	    J(3,3) = J(3,3) + shpZeta(i) * xl(3,i)
 	  enddo
 
 	detJ =J(1,1)*(J(2,2)*J(3,3)-J(2,3)*J(3,2))-J(1,2)*(J(2,1)*J(3,3)      ! determinant von jacobian matrix
     +          -J(2,3)*J(3,1))+J(1,3)*(J(2,1)*J(3,2)-J(2,2)*J(3,1))
     
        rint = detJ
  !      write(*,*)'rint shp1gp40',rint
     
 	  CofJ(1,1)=J(2,2)*J(3,3)-J(2,3)*J(3,2)        
 	  CofJ(1,2)=J(1,3)*J(3,2)-J(1,2)*J(3,3)
 	  CofJ(1,3)=J(1,2)*J(2,3)-J(1,3)*J(2,2)
 	  
 	  CofJ(2,1)=J(2,3)*J(3,1)-J(2,1)*J(3,3)
           CofJ(2,2)=J(1,1)*J(3,3)-J(1,3)*J(3,1)
           CofJ(2,3)=J(1,3)*J(2,1)-J(1,1)*J(2,3)
           
           CofJ(3,1)=J(2,1)*J(3,2)-J(2,2)*J(3,1)
           CofJ(3,2)=J(1,2)*J(3,1)-J(1,1)*J(3,2)
           CofJ(3,3)=J(1,1)*J(2,2)-J(1,2)*J(2,1)
           
           
           invJ(1,1)=CofJ(1,1)/detJ                         ! Inverse jacobian matrix
           invJ(1,2)=CofJ(1,2)/detJ
           invJ(1,3)=CofJ(1,3)/detJ
           
           invJ(2,1)=CofJ(2,1)/detJ
           invJ(2,2)=CofJ(2,2)/detJ
           invJ(2,3)=CofJ(2,3)/detJ
           
           invJ(3,1)=CofJ(3,1)/detJ
           invJ(3,2)=CofJ(3,2)/detJ
           invJ(3,3)=CofJ(3,3)/detJ

        shp(9)  = 0.125d0 * (1.0d0 - gp(ig+1,1)) * (1.0d0 - gp(ig+1,2))
     +  * (1.0d0 - gp(ig+1,3))
	shp(10) = 0.125d0 * (1.0d0 + gp(ig+1,1)) * (1.0d0 - gp(ig+1,2))
     +  * (1.0d0 - gp(ig+1,3))
	shp(11) = 0.125d0 * (1.0d0 + gp(ig+1,1)) * (1.0d0 + gp(ig+1,2))
     +  * (1.0d0 - gp(ig+1,3))
	shp(12) = 0.125d0 * (1.0d0 - gp(ig+1,1)) * (1.0d0 + gp(ig+1,2))
     +  * (1.0d0 - gp(ig+1,3))
	shp(13) = 0.125d0 * (1.0d0 - gp(ig+1,1)) * (1.0d0 - gp(ig+1,2))
     +  * (1.0d0 + gp(ig+1,3))
        shp(14) = 0.125d0 * (1.0d0 + gp(ig+1,1)) * (1.0d0 - gp(ig+1,2))
     +  * (1.0d0 + gp(ig+1,3))
	shp(15) = 0.125d0 * (1.0d0 + gp(ig+1,1)) * (1.0d0 + gp(ig+1,2)) 
     +  * (1.0d0 + gp(ig+1,3))
	shp(16) = 0.125d0 * (1.0d0 - gp(ig+1,1)) * (1.0d0 + gp(ig+1,2)) 
     +  * (1.0d0 + gp(ig+1,3))
	
	    call pzero(Naddneg,3*3*8)
	    do i = 9,16

	      Naddneg(1,1,i-8) = shp(i)
	      Naddneg(2,2,i-8) = shp(i)
	      Naddneg(3,3,i-8) = shp(i)
	    enddo
	!         
      shpXi(1) = -0.125d0 *  (1.0d0 - gp(ig+1,2)) * (1.0d0 - gp(ig+1,3))                   ! ABLEITUNG NACH XI
      shpXi(2) = +0.125d0 *  (1.0d0 - gp(ig+1,2)) * (1.0d0 - gp(ig+1,3))
      shpXi(3) = +0.125d0 *  (1.0d0 + gp(ig+1,2)) * (1.0d0 - gp(ig+1,3))
      shpXi(4) = -0.125d0 *  (1.0d0 + gp(ig+1,2)) * (1.0d0 - gp(ig+1,3))
      shpXi(5) = -0.125d0 *  (1.0d0 - gp(ig+1,2)) * (1.0d0 + gp(ig+1,3))
      shpXi(6) = +0.125d0 *  (1.0d0 - gp(ig+1,2)) * (1.0d0 + gp(ig+1,3))
      shpXi(7) = +0.125d0 *  (1.0d0 + gp(ig+1,2)) * (1.0d0 + gp(ig+1,3))
      shpXi(8) = -0.125d0 *  (1.0d0 + gp(ig+1,2)) * (1.0d0 + gp(ig+1,3))
 
      shpEta(1) = -0.125d0 * (1.0d0 - gp(ig+1,1)) * (1.0d0 - gp(ig+1,3))                ! ABLEITUNG NACH ETA 
      shpEta(2) = -0.125d0 * (1.0d0 + gp(ig+1,1)) * (1.0d0 - gp(ig+1,3))
      shpEta(3) = +0.125d0 * (1.0d0 + gp(ig+1,1)) * (1.0d0 - gp(ig+1,3))
      shpEta(4) = +0.125d0 * (1.0d0 - gp(ig+1,1)) * (1.0d0 - gp(ig+1,3))
      shpEta(5) = -0.125d0 * (1.0d0 - gp(ig+1,1)) * (1.0d0 + gp(ig+1,3))
      shpEta(6) = -0.125d0 * (1.0d0 + gp(ig+1,1)) * (1.0d0 + gp(ig+1,3))
      shpEta(7) = +0.125d0 * (1.0d0 + gp(ig+1,1)) * (1.0d0 + gp(ig+1,3))
      shpEta(8) = +0.125d0 * (1.0d0 - gp(ig+1,1)) * (1.0d0 + gp(ig+1,3))
     
      shpZeta(1) =-0.125d0 * (1.0d0 - gp(ig+1,1)) * (1.0d0 - gp(ig+1,2))                 ! ABLEUTUNG NACH ZETA
      shpZeta(2) =-0.125d0 * (1.0d0 + gp(ig+1,1)) * (1.0d0 - gp(ig+1,2))
      shpZeta(3) =-0.125d0 * (1.0d0 + gp(ig+1,1)) * (1.0d0 + gp(ig+1,2))
      shpZeta(4) =-0.125d0 * (1.0d0 - gp(ig+1,1)) * (1.0d0 + gp(ig+1,2))
      shpZeta(5) =+0.125d0 * (1.0d0 - gp(ig+1,1)) * (1.0d0 - gp(ig+1,2))
      shpZeta(6) =+0.125d0 * (1.0d0 + gp(ig+1,1)) * (1.0d0 - gp(ig+1,2))
      shpZeta(7) =+0.125d0 * (1.0d0 + gp(ig+1,1)) * (1.0d0 + gp(ig+1,2))
      shpZeta(8) =+0.125d0 * (1.0d0 - gp(ig+1,1)) * (1.0d0 + gp(ig+1,2))
 	
 	  call pzero(J,3*3)                                                          ! JACOBI LAUT ELEM 1 
 	  do i = 1,8
 	    J(1,1) = J(1,1) + shpXi(i) * xl(1,i)
 	    J(1,2) = J(1,2) + shpXi(i) * xl(2,i)
 	    J(1,3) = J(1,3) + shpXi(i) * xl(3,i)
 	    
 	    J(2,1) = J(2,1) + shpEta(i) * xl(1,i)
 	    J(2,2) = J(2,2) + shpEta(i) * xl(2,i)
 	    J(2,3) = J(2,3) + shpEta(i) * xl(3,i)
 	    
 	    J(3,1) = J(3,1) + shpZeta(i) * xl(1,i)
 	    J(3,2) = J(3,2) + shpZeta(i) * xl(2,i)
 	    J(3,3) = J(3,3) + shpZeta(i) * xl(3,i)
 	  enddo
 
 	detJ =J(1,1)*(J(2,2)*J(3,3)-J(2,3)*J(3,2))-J(1,2)*(J(2,1)*J(3,3)      ! determinant von jacobian matrix
     +          -J(2,3)*J(3,1))+J(1,3)*(J(2,1)*J(3,2)-J(2,2)*J(3,1))
      
 	   CofJ(1,1)=J(2,2)*J(3,3)-J(2,3)*J(3,2)        
 	   CofJ(1,2)=J(1,3)*J(3,2)-J(1,2)*J(3,3)
 	   CofJ(1,3)=J(1,2)*J(2,3)-J(1,3)*J(2,2)
 	   
 	   CofJ(2,1)=J(2,3)*J(3,1)-J(2,1)*J(3,3)
           CofJ(2,2)=J(1,1)*J(3,3)-J(1,3)*J(3,1)
           CofJ(2,3)=J(1,3)*J(2,1)-J(1,1)*J(2,3)
           
           CofJ(3,1)=J(2,1)*J(3,2)-J(2,2)*J(3,1)
           CofJ(3,2)=J(1,2)*J(3,1)-J(1,1)*J(3,2)
           CofJ(3,3)=J(1,1)*J(2,2)-J(1,2)*J(2,1)
           
           
           invJ(1,1)=CofJ(1,1)/detJ                         ! Inverse jacobian matrix
           invJ(1,2)=CofJ(1,2)/detJ
           invJ(1,3)=CofJ(1,3)/detJ
           
           invJ(2,1)=CofJ(2,1)/detJ
           invJ(2,2)=CofJ(2,2)/detJ
           invJ(2,3)=CofJ(2,3)/detJ
           
           invJ(3,1)=CofJ(3,1)/detJ
           invJ(3,2)=CofJ(3,2)/detJ
           invJ(3,3)=CofJ(3,3)/detJ
 
        return
        end
