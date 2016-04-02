      subroutine lodmbf(fval,freep,xm,den,matp)
c     implicit double precision (a-h,o-z)                                    dp
c
c     calulate magnetic field body force loads
c
       use mod_parameters
      real*8 hed                                                        vax750
      common/bk07/mbfc,nelpg,hed(12)
      common/range/mft,mlt,lft,llt,nftm1
      common/vect1/ r(nelemg,10)
      common/vect9/scl(8,nelemg),yz(8,nelemg)
      dimension fval(*),freep(5,*),xm(4,*),den(*),matp(1)
      if (mbfc.eq.0) return
      do 10 i=lft,llt
      bfmgy=freep(1,i)
      bfmgz=freep(2,i)
      if (bfmgy.eq.0.0.and.bfmgz.eq.0.0) go to 10
      rho=den(matp(i))
      facy=fval(mbfc)*bfmgy/rho
      facz=fval(mbfc)*bfmgz/rho
      r(i,1)=r(i,1)-xm(1,i)*facy
      r(i,2)=r(i,2)-xm(1,i)*facz
      r(i,3)=r(i,3)-xm(2,i)*facy
      r(i,4)=r(i,4)-xm(2,i)*facz
      r(i,5)=r(i,5)-xm(3,i)*facy
      r(i,6)=r(i,6)-xm(3,i)*facz
      r(i,7)=r(i,7)-xm(4,i)*facy
      r(i,8)=r(i,8)-xm(4,i)*facz
   10 continue
      return
      end
