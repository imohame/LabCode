      subroutine promul(a,b,c,jdiag,neq)
c     implicit double precision (a-h,o-z)                                    dp
      dimension a(*),b(*),c(*),jdiag(*)
      js=1
      do 30 j=1,neq
      jd=jdiag(j)
      if (js.gt.jd) go to 30
      bj=b(j)
      ab=a(jd)*bj
      if (js.eq.jd) go to 20
      jb=j-jd
      je=jd-1
      do 10 jj=js,je
   10 c(jj+jb)=c(jj+jb)+a(jj)*bj
      ab=ab+fdot(a(js),b(js+jb),jd-js)
   20 c(j)=c(j)+ab
   30 js=jd+1
      return
      end
