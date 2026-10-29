c program for backpropagation method
      implicit none
c
      double precision x,xx,y,snet,tpot_t,tpot0,delpot,tpotmin,w
      double precision ttneti,y0,u,z,ll,yy,azi,arnd,ran3,OO,s
      double precision netj,netjj,error,aaexp,ttt,Out,omega
      double precision temp0,tempt,err_0,errort,err_sq,wnew,error1
c
      integer nr,nc,nn,npatrn,kk,i,j,niter,n,k
      integer epoc,epocmax,ik,aM,as,idum,jdum,itmm
      integer jj,ii,kdum,mdum,ndum,itemp,patt,kkdum,imk
c
      parameter(nr=3,nc=3,npatrn=4)
      parameter(n=9,nn=9,epocmax=10)
c
c
      dimension w(52),xx(nn,npatrn)
      dimension azi(52),err_0(npatrn)
      dimension netj(4),netjj(4),ttt(4,4),Out(4),OO(4)
      dimension wnew(52)
c
       open(2,file="tpot0.dat")
       open(3,file="wij.dat")
       open(71,file="converged-patt.dat")
       open(73,file="all-graph.dat")
c
c===========assigning the values of w=============================
c
        jdum=-54
        do i=1,52
        azi(i)=ran3(jdum)
        azi(i)=azi(i)-0.5d0
        w(i)=azi(i)
c        
       write(74,*)i,azi(i),jdum,w(i)
       enddo  
c============== end of assigning the values of w====================        
c
       call pattern(xx,ttt)
c
       epoc=1
c
 120   continue
c
       temp0=100.0d0
       s=0.01d0
       itmm=20
       do  80 patt=1,4
c
       call init_err(err_0)
       errort=err_0(patt)
      write(*,*)"errort"
      write(*,*) errort
c
  90   continue
c
       do 81 itemp=1,itmm
       tempt=temp0-5.0d0*(i-1)
       omega=tempt
c 
c
       kdum=-55
       mdum=-110
       ndum=-180
       kkdum=-250
c
  70   continue
c
        do 82 kk=1,10     
c
       imk=1+51*ran3(kkdum)
       w(imk)=w(imk)+(1.0d0-2.0d0*ran3(kdum))*s*ran3(mdum)
      write(*,*)"imk,w(imk),ran3(kkdum),ran3(mdum)"
      write(*,*) imk,w(imk),ran3(kkdum),ran3(mdum)
c
c       kdum=kdum+1
c       mdum=mdum+1
c       ndum=ndum+1
c       kkdum=kkdum+1
c
       do jj=1,4 
       netj(jj)=0.0d0
       if (jj.eq.1) then
       do ii=1,9       
       netj(jj)=netj(jj)+w(ii)*xx(ii,patt)    
       enddo
       endif
       if (jj.eq.2) then
       do ii=10,18
       netj(jj)=netj(jj)+w(ii)*xx(ii-9,patt)
       enddo
       endif
c
       if (jj.eq.3) then
       do ii=19,27
       netj(jj)=netj(jj)+w(ii)*xx(ii-18,patt)
       enddo
       endif
c
       if (jj.eq.4) then
       do ii=28,36
       netj(jj)=netj(jj)+w(ii)*xx(ii-27,patt)
       enddo
       endif
c
       enddo
c
       do jj=1,4
       OO(jj)=1/(1+exp(-netj(jj)))
       enddo
       do jj=1,4
       netjj(jj)=0.0d0
c
       if(jj.eq.1)then
       do ii=1,4
       netjj(jj)=netjj(jj)+w(36+ii)*OO(ii)
       enddo
       endif
c
       if(jj.eq.2)then
       do ii=1,4
       netjj(jj)=netjj(jj)+w(40+ii)*OO(ii)
      enddo
      endif
c
       if(jj.eq.3)then
       do ii=1,4
       netjj(jj)=netjj(jj)+w(44+ii)*OO(ii)
      enddo
      endif
c
       if(jj.eq.4)then
       do ii=1,4
       netjj(jj)=netjj(jj)+w(48+ii)*OO(ii)
       enddo
       endif
c
       enddo
c
      error=0.0d0
      do jj=1,4
      Out(jj)=1/(1+exp(-netjj(jj)))
c
      error=error+0.5d0*(ttt(jj,patt)-Out(jj))**2  
      enddo
      if(error.lt.errort)then 
      wnew(imk)=w(imk)
      write(4,*)"error,errort"
      write(4,*) wnew(imk),error,errort
      error1=errort
      errort=error
      go to  70
      else
      if(error.gt.errort)then
c      errort=error
      err_sq=sqrt(error)
      aaexp=exp(-omega/err_sq) 
      else
      endif
      endif
      if(aaexp.gt.ran3(ndum)) then 
      wnew(imk)=w(imk)
c
      errort=error
      write(4,*)"wnew(imk),omega,error,errort"
      write(4,*) wnew(imk),error,errort
c
      go to 70
c      
      endif
c 
 82   continue
c
c
 81   continue
c
      if(abs(error).lt.0.00001d0)then
      go to 60
      else
      temp0=temp0/2
      itmm=itmm/2
      itmm=int(itmm)
      go to 90
      endif
c
 60   continue
 80   continue
c
c
      epoc=epoc+1
c
      if(epoc.gt.10)then
c
      go to 110
      endif
      write(*,*)"epoc"
      write(*,*)epoc
c
      go to 120     
c      
 110  continue
c 
      stop
      end
c end of main program
c================================================================
c subroutine to define the pattern
      subroutine pattern(xx,ttt)
      implicit none
c
      double precision x,xx,z,ll,ttt
      integer nr,nc,nn,npatrn,kk,i,j,k,ik
c
      parameter(nr=3,nc=3,nn=9,npatrn=4)
c
      dimension xx(nn,npatrn),z(nr,nc,npatrn),ll(nr,nc)
      dimension ttt(4,4)
c
      do k=1,npatrn
c k=1 is for matrix-1, k=2 is matrix-2, k=3 is matrix-3, k=4 is unknown matrix-4
      if (k.eq.1)then

      xx(1,k)=1.0d0
      xx(2,k)=2.0d0
      xx(3,k)=3.0d0
      xx(4,k)=2.0d0
      xx(5,k)=5.0d0
      xx(6,k)=6.0d0
      xx(7,k)=3.0d0
      xx(8,k)=6.0d0
      xx(9,k)=15.0d0
      ttt(1,k)=0.1392252d0
      ttt(2,k)=2.3249987d0
      ttt(3,k)=18.535776d0
      ttt(4,k)=3.2d0
   
      endif 

c
      if (k.eq.2)then

      xx(1,k)=11.0d0
      xx(2,k)=19.0d0
      xx(3,k)=36.0d0
      xx(4,k)=19.0d0
      xx(5,k)=25.0d0
      xx(6,k)=9.0d0
      xx(7,k)=36.0d0
      xx(8,k)=9.0d0
      xx(9,k)=45.0d0
      ttt(1,k)=-15.271218d0
      ttt(2,k)=21.870458d0
      ttt(3,k)=74.400759d0
      ttt(4,k)=4.5d0

      endif
      if (k.eq.3)then

      xx(1,k)=10.0d0
      xx(2,k)=17.0d0
      xx(3,k)=15.0d0
      xx(4,k)=17.0d0
      xx(5,k)=6.0d0
      xx(6,k)=49.0d0
      xx(7,k)=15.0d0
      xx(8,k)=49.0d0
      xx(9,k)=91.0d0
      ttt(1,k)=-19.598669d0
      ttt(2,k)=9.21939d0
      ttt(3,k)=117.37928d0
      ttt(4,k)=1.7d0
      endif
      if (k.eq.4) then
      xx(1,k)=51.0d0
      xx(2,k)=85.0d0
      xx(3,k)=66.0d0
      xx(4,k)=85.0d0
      xx(5,k)=7.0d0
      xx(6,k)=19.0d0
      xx(7,k)=66.0d0
      xx(8,k)=19.0d0
      xx(9,k)=100.0d0
c
      ttt(1,k)=-63.375719d0
      ttt(2,k)=47.387331d0
      ttt(3,k)=173.98839d0
      ttt(4,k)=10.2d0
c
      endif
c
      enddo
c
       ik =0
       do i =1,nr
       do j =1,nc
       ik=ik+1
       ll(i,j)=ik
       enddo
       enddo
c
      do k=1,npatrn
      do i=1,nr
      do j=1,nc
      kk=ll(i,j)
      z(i,j,k)=xx(kk,k)
      enddo
      enddo
      enddo
c 
      return
      end
c
      subroutine pack(i,j,ij)
      integer i,j,ij,II,JJ
      II=10*i
      JJ=j
      ij=II+JJ
      return
      end
c*********end of subroutine*********************
c=========subroutine-calculate-initial error=================
c
      subroutine init_err(err_0)
c
      implicit none
c
      double precision x,xx,y,snet,tpot_t,tpot0,delpot,tpotmin,w
      double precision ttneti,y0,u,z,ll,yy,azi,arnd,ran3,OO,s
      double precision netj,netjj,error,aaexp,ttt,Out,omega
      double precision temp0,tempt,err_0,errort,err_sq,wnew,error1
c
      integer nr,nc,nn,npatrn,kk,i,j,niter,n,k
      integer epoc,epocmax,ik,aM,as,idum,jdum,itmm
      integer jj,ii,kdum,mdum,ndum,itemp,patt,kkdum,imk
c
      parameter(nr=3,nc=3,npatrn=4)
      parameter(n=9,nn=9,epocmax=10)
c
c
      dimension w(52),xx(nn,npatrn)
      dimension azi(52),err_0(npatrn)
      dimension netj(4),netjj(4),ttt(4,4),Out(4),OO(4)
      dimension wnew(52)
c
c===========assigning the values of w=============================
c
        do i=1,52
        jdum=-54+i
        azi(i)=ran3(jdum)
        azi(i)=azi(i)-0.5d0
        w(i)=azi(i)
c
       enddo
       call pattern(xx,ttt)       
       do patt=1,npatrn
c
       do jj=1,4
       netj(jj)=0.0d0
       if (jj.eq.1) then
       do ii=1,9
       netj(jj)=netj(jj)+w(ii)*xx(ii,patt)
       enddo
       endif
       if (jj.eq.2) then
       do ii=10,18
       netj(jj)=netj(jj)+w(ii)*xx(ii-9,patt)
       enddo
       endif
c
       if (jj.eq.3) then
       do ii=19,27
       netj(jj)=netj(jj)+w(ii)*xx(ii-18,patt)
       enddo
       endif
c
       if (jj.eq.4) then
       do ii=28,36
       netj(jj)=netj(jj)+w(ii)*xx(ii-27,patt)
       enddo
       endif
c
       enddo
c
       do jj=1,4
       OO(jj)=1/(1+exp(-netj(jj)))
       enddo
       do jj=1,4
       netjj(jj)=0.0d0
c
       if(jj.eq.1)then
       do ii=1,4
       netjj(jj)=netjj(jj)+w(36+ii)*OO(ii)
       enddo
       endif
c
       if(jj.eq.2)then
       do ii=1,4
       netjj(jj)=netjj(jj)+w(40+ii)*OO(ii)
      enddo
      endif
c
       if(jj.eq.3)then
       do ii=1,4
       netjj(jj)=netjj(jj)+w(44+ii)*OO(ii)
      enddo
      endif
c
       if(jj.eq.4)then
       do ii=1,4
       netjj(jj)=netjj(jj)+w(48+ii)*OO(ii)
       enddo
       endif
c
       enddo
c
      err_0(patt)=0.0d0
      do jj=1,4
      Out(jj)=1/(1+exp(-netjj(jj)))
c
      err_0(patt)=err_0(patt)+0.5d0*(ttt(jj,patt)-Out(jj))**2
      enddo
c
      enddo
c
      return
      end
c=========end of subroutine==================================
************************************************************************
* Filename: ran3.f
************************************************************************
*
* Random number generator from Numerical Recipes (no modifications)
* "ran3" takes a -ve number to initialize a set of random
* numbers (only needed once). After that, the argument that gets
* passed in is "1" and does not initialize the function again.
*
************************************************************************
       double precision FUNCTION ran3(idum)
        IMPLICIT none
C Global Variable Declarations
         INTEGER idum
C Local Variable Declarations
       INTEGER MBIG,MSEED,MZ
       DOUBLE PRECISION FAC
       PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1.0d0/MBIG)
       INTEGER i,iff,ii,inext,inextp,k
       INTEGER mj,mk,ma(55)
       SAVE iff,inext,inextp,ma
       DATA iff /0/
*************************************************************************
       if(idum.lt.0.or.iff.eq.0)then
       iff=1
       mj=MSEED-iabs(idum)
       mj=mod(mj,MBIG)
       ma(55)=mj
       mk=1
       do 11 i=1,54
       ii=mod(21*i,55)
       ma(ii)=mk
       mk=mj-mk
       if(mk.lt.MZ)mk=mk+MBIG
       mj=ma(ii)
 11    continue
       do 13 k=1,4
      do 12 i=1,55
      ma(i)=ma(i)-ma(1+mod(i+30,55))
      if(ma(i).lt.MZ)ma(i)=ma(i)+MBIG
 12   continue
 13   continue
      inext=0
      inextp=31
c idum=1
      endif
      inext=inext+1
      if(inext.eq.56)inext=1
       inextp=inextp+1
       if(inextp.eq.56)inextp=1
       mj=ma(inext)-ma(inextp)
       if(mj.lt.MZ)mj=mj+MBIG
       ma(inext)=mj
       ran3 = mj * FAC
       IF (ran3 .LT. 0. .OR. ran3 .GT. 1.) THEN
       WRITE (*,*) '*** ran3.f: something went wrong!'
       WRITE (*,*) '*** ran3.f: ran3 =', ran3
      ENDIF
      RETURN
      END
