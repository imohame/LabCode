      subroutine CG(a,d,x,xo,n)
      integer n
      real*8 a(n,n),d(n,1),x(n,1),xo(n,1),xstore(n,1)
      real*8 p(n,1),c,r(n,1),rT(1,n),pT(1,n),top(1,1),bottom(1,1)
      real*8 roT(1,n),test,ro(n,1),beta(1,1),rstore,test1(n,1)
      
   10 r=d-matmul(A,xo)
      rstore=maxval(r)
      p=r
      if (abs(maxval(r)).lt.1E-16) then
      x=xo
      return
      endif
      do i=1,1000
c      write(*,*) maxval(r)
      call matrixtrans(p,pT,n,1)
      call matrixtrans(r,rT,n,1)
      
c      write(*,*) 'before top'
      top=matmul(rT,r)
c      write(*,*) 'before bottom'
      bottom=matmul(pT,matmul(a,p))
c      write(*,*) bottom
      
      c=top(1,1)/bottom(1,1)
c      write(*,*) c

      
      x=xo+c*p
      ro=r
      roT=rT
      xo=x
c      write(*,*) 'before r'
      r=ro-c*matmul(A,p)
      if (c.eq.0) then
      test1=x-xstore
      test=abs(maxval(test1))
c      write(*,*) test
      if (abs(maxval(d-matmul(a,x))).ge.abs(rstore)) then

c      write(*,*) maxval(abs(d-matmul(a,x)))
      return
      endif
      if (test.le.1.E-8) then
      return
      endif
      xstore=x
c      write(*,*) maxval(abs(d-matmul(a,x))),rstore
      
      goto 10

      endif

      
      if (maxval(r).lt.1E-16)then
      return
      endif
      
      
      call matrixtrans(r,rT,n,1)
c      write(*,*) 'before beta'
      beta=matmul(rT,r)/matmul(roT,ro)
      p=r+beta*p
      enddo
      

      end
      
      
      
      
      
      
      