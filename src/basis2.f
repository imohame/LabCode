      subroutine basis2 (h,p,xj,det,xx)
c     implicit double precision (a-h,o-z)                                    dp
c
      common/bk01/h4(4,5),p14(4,5),p24(4,5)
      common/bk02/ioofc,iphase,imass,lpar(9)
      common/bk14/lfna(15),lfnt(6)
      common/bk48/stress(4),strain(4),d(4,4),ipt,nel,nstate
      dimension h(*),p(2,*),xj(2,2),xx(2,*)
c
     
      do 10 i=1,4
      h(i)=h4(i,ipt)
      p(1,i)=p14(i,ipt)
   10 p(2,i)=p24(i,ipt)
      do 30 i=1,2
      do 30 j=1,2
      xj(i,j)=0.0
      do 20 k=1,4
   20 xj(i,j)=xj(i,j)+p(i,k)*xx(j,k)
   30 continue
      det=xj(1,1)*xj(2,2)-xj(2,1)*xj(1,2)
      if (det.gt.0.0) go to 40
      write(lfnt(2),50) nel
      write (*,50) nel
   
      call bye (2)
   40 continue
      return
c
   50 format(/'zero or negative jacobian for element',i5)
      end
