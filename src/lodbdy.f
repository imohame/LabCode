      subroutine lodbdy(fval,xm)
c     implicit double precision (a-h,o-z)                                    dp
c
c     calulate body force loads
c
       use mod_parameters
      common/range/mft,mlt,lft,llt,nftm1
      
      common/vect1/r1(nelemg),r2(nelemg),r3(nelemg),
     > r4(nelemg),r5(nelemg),r6(nelemg),
     > r7(nelemg),r8(nelemg),mtype(nelemg),mte(nelemg)
      
      common/vect2/
     1 ed1(nelemg),ed2(nelemg),ed3(nelemg),ed4(nelemg),
     2 ed5(nelemg),ed6(nelemg),ed7(nelemg),ed8(nelemg),
     3 fd1(nelemg),fd2(nelemg),fd3(nelemg),fd4(nelemg),
     4 fd5(nelemg),fd6(nelemg),fd7(nelemg),fd8(nelemg)
      common/vect9/scale(8,nelemg),yz(8,nelemg)
      
      common/bk27/nlcur,nptst,nthpy,nthpz,nthps,xmy,xmz,xms,nload,nptm
      dimension fval(*),xm(4,*)
      
      if (nthpy.ne.0) then
      facy=fval(nthpy)*xmy
      do 20 i=lft,llt
      r1(i)=r1(i)+xm(1,i)*facy
      r3(i)=r3(i)+xm(2,i)*facy
      r5(i)=r5(i)+xm(3,i)*facy
      r7(i)=r7(i)+xm(4,i)*facy
   20 continue
      endif
      if (nthpz.ne.0) then
      facz=fval(nthpz)*xmz
      do 30 i=lft,llt
      r2(i)=r2(i)+xm(1,i)*facz
      r4(i)=r4(i)+xm(2,i)*facz
      r6(i)=r6(i)+xm(3,i)*facz
      r8(i)=r8(i)+xm(4,i)*facz
   30 continue
      endif
      if (nthps.eq.0) return
      ispny=0
      if (nthps.lt.0) ispny=1
      nthps=abs(nthps)
      facs=(fval(nthps)*xms)**2
      do 60 i=lft,llt
      r1(i)=r1(i)-xm(1,i)*(yz(1,i)+ed1(i))*facs
      r3(i)=r3(i)-xm(2,i)*(yz(3,i)+ed3(i))*facs
      r5(i)=r5(i)-xm(3,i)*(yz(5,i)+ed5(i))*facs
      r7(i)=r7(i)-xm(4,i)*(yz(7,i)+ed7(i))*facs
   60 continue
      if (ispny.ne.0) then
      do 70 i=lft,llt
      r2(i)=r2(i)-xm(1,i)*(yz(2,i)+ed2(i))*facs
      r4(i)=r4(i)-xm(2,i)*(yz(4,i)+ed4(i))*facs
      r6(i)=r6(i)-xm(3,i)*(yz(6,i)+ed6(i))*facs
      r8(i)=r8(i)-xm(4,i)*(yz(8,i)+ed8(i))*facs
   70 continue
      endif
      if (ispny.eq.1) nthps=-nthps
      return
      end
