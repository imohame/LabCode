      subroutine convck(d,dp,dtol,conv,nv,iter)
c     implicit double precision (a-h,o-z)                                    dp
      logical conv
      character*4 mess
      common/bk14/lfna(15),lfnt(6)
      common/bk29/numfrq,clengt
      dimension d(*),dp(*),dtol(*)
      data tol/1.e-9/,its/50/
      do 10 n=1,nv
      dtol(n)=abs((d(n)-dp(n))/d(n))
   10 dp(n)=d(n)
      write(lfnt(2),30) iter,(d(n),n=1,nv)
c
      call intrup  (mess,nwr,1)
      if (mess.eq.'sw2.') write (*,30) iter,(d(n),n=1,nv)
c
      write(lfnt(2),40) (dtol(n),n=1,nv)
c     call empty (17)                                                   cray1
      do 20 n=1,numfrq
      if (dtol(n).gt.tol) return
   20 continue
      conv=.true.
      return
c
   30 format(//'iteration',i5,' current eigenvalues'/(4e20.8))
   40 format(/5x,'current residuals'/(4e20.8))
      end
